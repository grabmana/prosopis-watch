import ee
import geemap

ee.Authenticate()
ee.Initialize(project='your-gcp-project-id')  # e.g. 'my-gcp-project'


CONFIG = {
    'prosopis_asset': 'path-to-your-prosopis-presence-points',  # e.g. 'users/yourusername/prosopis_presence'
    'absence_asset':  'path-to-your-absence-points',   # e.g. 'users/yourusername/prosopis_absence'
    'bbox': [38.10164648805493, 6.898560072409943, 43.06746680055493, 14.694609699887478],
    'dry_months': (11, 2),
    'wet_months': (6, 9),
    'max_cloud': 20,
    'input_bands': ['NDVI_dry', 'EVI_dry', 'elevation', 'dist_river'],
    'n_trees': 50,
    'train_split': 0.7,
    'random_seed': 42,
    'pseudo_absence_n': 200,
    'exclude_lc_classes': [1, 6, 7, 8],
    'export_folder': 'GEE_Exports',
    'export_filename': 'prosopis_afar',
    'export_scale': 30,
    'export_crs': 'EPSG:4326',
}

VIS = {
    'ndvi':       {'min': -0.1, 'max': 0.6, 'palette': ['brown', 'white', 'darkgreen']},
    'amplitude':  {'min': 0, 'max': 0.5, 'palette': ['darkgreen', 'white', 'brown']},
    'classified': {'min': 0, 'max': 1, 'palette': ['lightgrey', 'darkgreen']},
    'landcover':  {'min': 1, 'max': 8, 'palette': ['#c4281b', '#e49634', '#397e48',
                   '#DFC35A', '#88af52', '#429ae4', '#7c85c9', '#a6a6a6']},
}


def make_study_area(bbox):
    """
    Build a rectangular study area from a bounding box and return the Afar admin boundary.

    Args:
        bbox (list): [lon_min, lat_min, lon_max, lat_max]

    Returns:
        tuple: (ee.Geometry.Polygon, ee.FeatureCollection)
    """
    lon_w, lat_s, lon_e, lat_n = bbox
    geom = ee.Geometry.Polygon(
        [[[lon_w, lat_n], [lon_w, lat_s], [lon_e, lat_s], [lon_e, lat_n]]], None, False)
    admin = (ee.FeatureCollection('FAO/GAUL/2015/level1')
             .filter(ee.Filter.eq('ADM0_NAME', 'Ethiopia'))
             .filter(ee.Filter.eq('ADM1_NAME', 'Afar')))
    return geom, admin


def load_s2_median(study_area, month_start, month_end, max_cloud=20, bands=None):
    """
    Load a cloud-filtered Sentinel-2 SR median composite for a given calendar month range.

    Args:
        study_area (ee.Geometry): Region to filter and clip the collection.
        month_start (int): Start month (1-12).
        month_end (int): End month (1-12). Can wrap across year-end (e.g. 11 to 2).
        max_cloud (int): Maximum cloud cover percentage. Defaults to 20.
        bands (list): List of band names to select. Defaults to ['B2','B4','B5','B8','B8A','B11'].

    Returns:
        ee.Image: Median composite clipped to study_area.
    """
    if bands is None:
        bands = ['B2', 'B4', 'B5', 'B8', 'B8A', 'B11']
    return (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
            .filterBounds(study_area)
            .filter(ee.Filter.calendarRange(month_start, month_end, 'month'))
            .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', max_cloud))
            .select(bands)
            .median()
            .clip(study_area))


def compute_ndvi(image, nir='B8', red='B4', name='NDVI'):
    """
    Compute Normalised Difference Vegetation Index.

    Args:
        image (ee.Image): Input image with NIR and Red bands.
        nir (str): NIR band name. Defaults to 'B8'.
        red (str): Red band name. Defaults to 'B4'.
        name (str): Output band name. Defaults to 'NDVI'.

    Returns:
        ee.Image: Single-band NDVI image.
    """
    return image.normalizedDifference([nir, red]).rename(name)


def compute_evi(image, name='EVI'):
    """
    Compute Enhanced Vegetation Index. Less prone to saturation than NDVI in dense canopy.

    Args:
        image (ee.Image): Input image with B2, B4, B8 bands.
        name (str): Output band name. Defaults to 'EVI'.

    Returns:
        ee.Image: Single-band EVI image.
    """
    return image.expression(
        '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))',
        {'NIR': image.select('B8'), 'RED': image.select('B4'), 'BLUE': image.select('B2')}
    ).rename(name)


def compute_ndre(image, name='NDRE'):
    """
    Compute Red-Edge NDVI using Sentinel-2 B8A and B5. More sensitive to
    chlorophyll content than standard NDVI.

    Args:
        image (ee.Image): Input image with B8A and B5 bands.
        name (str): Output band name. Defaults to 'NDRE'.

    Returns:
        ee.Image: Single-band NDRE image.
    """
    return image.normalizedDifference(['B8A', 'B5']).rename(name)


def compute_ndvi_amplitude(ndvi_wet, ndvi_dry, name='NDVI_amplitude'):
    """
    Compute seasonal NDVI amplitude (wet minus dry). Prosopis juliflora remains
    green year-round, so low amplitude values are indicative of its presence.

    Args:
        ndvi_wet (ee.Image): NDVI from the wet season composite.
        ndvi_dry (ee.Image): NDVI from the dry season composite.
        name (str): Output band name. Defaults to 'NDVI_amplitude'.

    Returns:
        ee.Image: Single-band amplitude image.
    """
    return ndvi_wet.subtract(ndvi_dry).rename(name)


def load_elevation(study_area):
    """
    Load SRTM 30 m elevation data clipped to the study area.

    Args:
        study_area (ee.Geometry): Region to clip the DEM.

    Returns:
        ee.Image: Single-band elevation image in metres.
    """
    return ee.Image('USGS/SRTMGL1_003').clip(study_area).rename('elevation')


def load_dist_to_river(study_area, threshold=30, kernel=256):
    """
    Derive approximate distance-to-river (metres) from JRC Global Surface Water.
    Pixels with water occurrence above the threshold are treated as permanent water.

    Args:
        study_area (ee.Geometry): Region to clip the output.
        threshold (int): Minimum water occurrence percentage to classify as water. Defaults to 30.
        kernel (int): Neighbourhood size for distance transform. Defaults to 256.

    Returns:
        ee.Image: Single-band distance-to-river image in metres.
    """
    water = (ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
             .select('occurrence')
             .clip(study_area)
             .gt(threshold))
    return (water.fastDistanceTransform(kernel)
            .sqrt()
            .multiply(ee.Image.pixelArea().sqrt())
            .rename('dist_river')
            .clip(study_area))


def build_feature_stack(study_area, cfg):
    """
    Build a multi-band feature image combining spectral indices and environmental layers.
    Bands included: NDVI_dry, NDVI_wet, NDVI_amplitude, NDRE_dry, EVI_dry, elevation, dist_river.

    Args:
        study_area (ee.Geometry): Study area geometry.
        cfg (dict): Configuration dictionary with season windows and cloud threshold.

    Returns:
        ee.Image: Stacked feature image.
    """
    dry_s, dry_e = cfg['dry_months']
    wet_s, wet_e = cfg['wet_months']

    s2_dry = load_s2_median(study_area, dry_s, dry_e, cfg['max_cloud'])
    # s2_wet = load_s2_median(study_area, wet_s, wet_e, cfg['max_cloud'], bands=['B4', 'B8'])

    ndvi_dry = compute_ndvi(s2_dry, name='NDVI_dry')
    # ndvi_wet = compute_ndvi(s2_wet, name='NDVI_wet')

    return (ndvi_dry
            # .addBands(ndvi_wet)
            # .addBands(compute_ndvi_amplitude(ndvi_wet, ndvi_dry))
            # .addBands(compute_ndre(s2_dry, name='NDRE_dry'))
            .addBands(compute_evi(s2_dry, name='EVI_dry'))
            .addBands(load_elevation(study_area))
            .addBands(load_dist_to_river(study_area)))


def apply_landcover_mask(stack, study_area, exclude_classes):
    """
    Mask out unwanted land-cover classes from the feature stack using ESA WorldCover 2021.
    Original WorldCover values are remapped to: 1=Urban, 2=Cropland, 3=Forest,
    4=Shrub, 5=Grass, 6=Water, 7=Wetland, 8=Bare.

    Args:
        stack (ee.Image): Feature stack to mask.
        study_area (ee.Geometry): Used to clip the land-cover image for display.
        exclude_classes (list): Remapped class integers to exclude (e.g. [1, 6, 7, 8]).

    Returns:
        tuple: (masked ee.Image, land_cover ee.Image)
    """
    old_lc = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]
    new_lc = [ 3,  4,  5,  2,  1,  8,  9,  6,  7,  7,   5]
    land_cover = ee.Image('ESA/WorldCover/v200/2021').remap(old_lc, new_lc).rename('constant')
    mask = land_cover.remap(exclude_classes, ee.List.repeat(0, len(exclude_classes)), 1)
    return stack.updateMask(mask), land_cover


def load_ground_truth(study_area, land_cover, cfg):
    """
    Load and merge presence and absence points for classifier training.
    Absence points include both field-collected data and pseudo-absences
    sampled from cropland, forest, shrub, and grassland cover classes.

    Args:
        study_area (ee.Geometry): Region for pseudo-absence sampling.
        land_cover (ee.Image): Remapped WorldCover image.
        cfg (dict): Config with asset paths, pseudo_absence_n, and random_seed.

    Returns:
        tuple: (merged FeatureCollection, presence FeatureCollection, absence FeatureCollection)
    """
    presence = (ee.FeatureCollection(cfg['prosopis_asset'])
                .map(lambda f: f.set('label', 1)))
    absence = (ee.FeatureCollection(cfg['absence_asset'])
               .map(lambda f: f.set('label', 0)))
    pseudo = (land_cover
              .updateMask(land_cover.eq(2).Or(land_cover.eq(3))
                          .Or(land_cover.eq(4)).Or(land_cover.eq(5)))
              .sample(region=study_area, scale=30,
                      numPixels=cfg['pseudo_absence_n'],
                      seed=cfg['random_seed'], geometries=True)
              .map(lambda f: f.set('label', 0)))
    return presence.merge(absence).merge(pseudo), presence, absence


def sample_features(masked_stack, ground_truth, input_bands):
    """
    Sample the feature stack at ground-truth point locations.
    Drops any points that fall outside the masked area (null values),
    matching the implicit behaviour of GEE JavaScript sampleRegions.

    Args:
        masked_stack (ee.Image): Masked feature stack.
        ground_truth (ee.FeatureCollection): Points with a 'label' property (1=presence, 0=absence).
        input_bands (list): Band names used for classification — used to filter nulls.

    Returns:
        ee.FeatureCollection: Sampled training data with band values and labels.
    """
    sampled = masked_stack.sampleRegions(
        collection=ground_truth, properties=['label'], scale=30, tileScale=4)
    return sampled.filter(ee.Filter.notNull(input_bands))


def train_classifier(training, input_bands, n_trees=50):
    """
    Train a Random Forest classifier.

    Args:
        training (ee.FeatureCollection): Sampled training data.
        input_bands (list): Band names to use as input features.
        n_trees (int): Number of trees. Defaults to 50.

    Returns:
        ee.Classifier: Trained classifier.
    """
    return (ee.Classifier.smileRandomForest(numberOfTrees=n_trees)
            .train(features=training, classProperty='label', inputProperties=input_bands))


def classify_image(masked_stack, clf, study_area):
    """
    Apply a trained classifier to the feature stack.

    Args:
        masked_stack (ee.Image): Masked feature stack.
        clf (ee.Classifier): Trained classifier.
        study_area (ee.Geometry): Region to clip the output.

    Returns:
        ee.Image: Binary classified image (1=Prosopis, 0=other).
    """
    return masked_stack.classify(clf).clip(study_area)


def accuracy_assessment(training, input_bands, n_trees=50, train_split=0.7, random_seed=42):
    """
    Evaluate classifier performance using a random hold-out split.

    Args:
        training (ee.FeatureCollection): Full training dataset.
        input_bands (list): Band names used for classification.
        n_trees (int): Number of trees for the validation classifier. Defaults to 50.
        train_split (float): Fraction of data used for training. Defaults to 0.7.
        random_seed (int): Seed for reproducible splitting. Defaults to 42.

    Returns:
        dict: Keys are confusion_matrix, overall_accuracy, kappa,
              producers_accuracy, consumers_accuracy.
    """
    split     = training.randomColumn('random', random_seed)
    train_set = split.filter(ee.Filter.lt('random', train_split))
    test_set  = split.filter(ee.Filter.gte('random', train_split))
    clf       = (ee.Classifier.smileRandomForest(n_trees)
                 .train(features=train_set, classProperty='label', inputProperties=input_bands))
    conf_mat  = test_set.classify(clf).errorMatrix('label', 'classification')
    return {
        'confusion_matrix':   conf_mat,
        'overall_accuracy':   conf_mat.accuracy(),
        'kappa':              conf_mat.kappa(),
        'producers_accuracy': conf_mat.producersAccuracy(),
        'consumers_accuracy': conf_mat.consumersAccuracy(),
    }


def print_accuracy(metrics):
    """
    Print accuracy metrics returned by accuracy_assessment().

    Args:
        metrics (dict): Output from accuracy_assessment().
    """
    print('Confusion Matrix:')
    print(metrics['confusion_matrix'].getInfo())
    print(f"Overall Accuracy : {metrics['overall_accuracy'].getInfo():.4f}")
    print(f"Kappa            : {metrics['kappa'].getInfo():.4f}")
    print('Producers Accuracy:', metrics['producers_accuracy'].getInfo())
    print('Consumers Accuracy:', metrics['consumers_accuracy'].getInfo())


def export_map(classified, study_area, cfg):
    """
    Start a Google Drive export task for the classified image.

    Args:
        classified (ee.Image): Binary Prosopis classification output.
        study_area (ee.Geometry): Export region.
        cfg (dict): Config with export_folder, export_filename, export_scale, export_crs.

    Returns:
        ee.batch.Task: Started export task.
    """
    task = ee.batch.Export.image.toDrive(
        image=classified,
        description='Prosopis_Afar_Map',
        folder=cfg['export_folder'],
        fileNamePrefix=cfg['export_filename'],
        region=study_area,
        scale=cfg['export_scale'],
        maxPixels=int(1e9),
        crs=cfg['export_crs']
    )
    task.start()
    return task


def add_layers(m, study_area, admin, land_cover, ndvi_dry,
               ndvi_wet, ndvi_amp, classified, presence, absence):
    """
    Add all map layers to a geemap.Map instance.

    Args:
        m (geemap.Map): Target map object.
        study_area (ee.Geometry): Study area geometry.
        admin (ee.FeatureCollection): Afar admin boundary.
        land_cover (ee.Image): Remapped WorldCover image.
        ndvi_dry (ee.Image): Dry-season NDVI.
        ndvi_wet (ee.Image): Wet-season NDVI.
        ndvi_amp (ee.Image): NDVI amplitude.
        classified (ee.Image): Binary Prosopis classification.
        presence (ee.FeatureCollection): Prosopis presence points.
        absence (ee.FeatureCollection): Absence points.
    """
    m.centerObject(study_area, 9)
    m.addLayer(admin,                       {'color': '888888'}, 'Afar Region')
    m.addLayer(land_cover.clip(study_area), VIS['landcover'],    'WorldCover 2021', False)
    m.addLayer(ndvi_dry,                    VIS['ndvi'],         'NDVI dry',        False)
    m.addLayer(ndvi_wet,                    VIS['ndvi'],         'NDVI wet',        False)
    m.addLayer(ndvi_amp,                    VIS['amplitude'],    'NDVI Amplitude',  False)
    m.addLayer(classified,                  VIS['classified'],   'Prosopis Map')
    m.addLayer(presence,                    {'color': 'lime'},   'Prosopis presence')
    m.addLayer(absence,                     {'color': 'red'},    'Absence points')


def run_pipeline(cfg):
    """
    Run the full Prosopis mapping pipeline end-to-end.

    Steps: study area → feature stack → land-cover mask → ground truth →
           classifier training → classification → accuracy assessment → visualisation.
    Export is commented out by default — uncomment to activate.

    Args:
        cfg (dict): Configuration dictionary (see CONFIG at top of file).

    Returns:
        tuple: (classified ee.Image, interactive geemap.Map)
    """
    study_area, admin = make_study_area(cfg['bbox'])

    raw_stack = build_feature_stack(study_area, cfg)
    masked_stack, land_cover = apply_landcover_mask(raw_stack, study_area, cfg['exclude_lc_classes'])
    print('Bands:', raw_stack.bandNames().getInfo())

    ground_truth, presence, absence = load_ground_truth(study_area, land_cover, cfg)
    training = sample_features(masked_stack, ground_truth, cfg['input_bands'])

    print(f"Prosopis : {presence.size().getInfo()}")
    print(f"Absence  : {absence.size().getInfo()}")
    print(f"Training : {training.size().getInfo()}")

    clf = train_classifier(training, cfg['input_bands'], cfg['n_trees'])
    importance = clf.explain().get('importance').getInfo()
    print('\nFeature importance:')
    for band, score in sorted(importance.items(), key=lambda x: -x[1]):
        print(f'  {band:<22} {score:.4f}')

    classified = classify_image(masked_stack, clf, study_area)

    print('\nAccuracy (70/30 hold-out):')
    print_accuracy(accuracy_assessment(
        training, cfg['input_bands'],
        n_trees=cfg['n_trees'],
        train_split=cfg['train_split'],
        random_seed=cfg['random_seed']
    ))

    # export_map(classified, study_area, cfg)

    dry_s, dry_e = cfg['dry_months']
    wet_s, wet_e = cfg['wet_months']
    s2_dry   = load_s2_median(study_area, dry_s, dry_e, cfg['max_cloud'])
    s2_wet   = load_s2_median(study_area, wet_s, wet_e, cfg['max_cloud'], bands=['B4', 'B8'])
    ndvi_dry = compute_ndvi(s2_dry, name='NDVI_dry')
    ndvi_wet = compute_ndvi(s2_wet, name='NDVI_wet')
    ndvi_amp = compute_ndvi_amplitude(ndvi_wet, ndvi_dry)

    m = geemap.Map()
    add_layers(m, study_area, admin, land_cover,
               ndvi_dry, ndvi_wet, ndvi_amp, classified, presence, absence)

    return classified, m


if __name__ == '__main__':
    classified_map, interactive_map = run_pipeline(CONFIG)
    interactive_map