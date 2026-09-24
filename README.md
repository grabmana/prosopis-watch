# Prosopis Watch

Open [`index.html`](index.html) in a modern browser with internet access. It is a
single-file satellite screening and decision-support app: **no build step,
backend, API keys, or installed packages**. Styles and application logic are
inline; Leaflet 1.9.4 loads from unpkg. The existing research datasets and
Earth Engine/Python model versions below are preserved and are **not loaded or
executed by this app**.

## Five capabilities

1. **Map & change:** place search, dated Sentinel-2 scene selection with actual
   catalog cloud percentages, NDVI + NDMI screening, clustered candidate pins,
   and new/stable/loss comparison. Two independent, synchronized Leaflet maps
   let you swipe between baseline and comparison pins with mouse, touch, or
   keyboard. When every in-window scene is at least 5% cloudy, the app searches
   ±45 days and recommends up to three clearer outside-period observations.
2. **Regrowth vs recovery:** twelve separate monthly searches, one least-cloudy
   scene per month, with local SCL cloud/shadow rejection. Woody/grass/inconclusive
   phenology is assessed only with at least five usable months per location
   spanning at least 240 days. Requested and observed date ranges are displayed;
   absent months are not interpolated.
3. **Re-invasion risk:** 16/24/32-cell weighted composite of distance to mapped
   candidate seed sources, NDMI, Copernicus DEM slope, and observed MODIS burns
   during the last 24 completed calendar months. This is **not a calibrated
   probability or time-to-invasion prediction**. Missing observations and zero
   fire coverage are reported. With no MODIS scenes, the fire term is zero as an
   explicit fallback; a fire-only weighting remains unscored because there is
   no measured input contributing to that composite.
4. **Targeting:** relative ranking using nearby seasonal grass signals, NDWI
   water proximity, slope, and sampled candidate density, with CSV export.
   Requires matching screening and seasonal analysis extents.
5. **MRV quality control:** paste `lat,lon,YYYY-MM-DD,activity` CSV records
   (optional header; quoted fields supported; up to 30 records). Clearing/removal
   claims are screened using the least-cloudy scene before and after the claim
   date, a default 60-day window on each side, and a default 0.10 NDVI drop.
   Results and CSV audit include scene IDs, dates, cloud cover, measurements,
   thresholds, reasons, and check time. Other activities are inconclusive.

### Scientific interpretation

**This app is a spectral indicator, not a trained classifier.** It does not use
the repository's trained research models or claim their accuracy. Defaults are
NDVI ≥ 0.35 and NDMI ≥ 0.05; both are adjustable. Persistent dry-season greenness
and moisture can distinguish woody cover from senescent grass, but other moist
woody species, native trees, and irrigated crops can produce the same signal.
**All findings require field validation.**

- Screening samples every eighth pixel across at most 16 lossless PNG tiles;
  zoom is reduced to cover the viewport. Grouping within ~0.0035° produces
  candidate locations, not stand polygons or measured invasion area. Displayed
  confidence is heuristic spectral strength, not species probability.
- SCL retains vegetation, bare ground and water classes only. Residual clouds,
  mixed pixels, 20 m SWIR resolution, coarse sampling, and ~0.0078 index
  quantization limit interpretation. Bands are used as delivered; there is no
  bespoke radiometric/processing-baseline harmonization.
- Compare matching seasons. New means newly meeting both thresholds with a
  sufficient NDVI rise; loss means a baseline candidate with sufficient decline;
  stable means a remaining candidate in both. Unresolved transitions and
  missing paired observations are not assigned a confident change class.
- Phenology: woody if minimum NDVI ≥ 0.28 and amplitude ≤ 0.28; grass if minimum
  < 0.22 and amplitude > 0.25; otherwise inconclusive. Sparse monthly sampling
  can miss seasonal extrema. These are locally uncalibrated rules.
- DEM PNGs use a coarse first pass followed by locally bounded rescaling to
  reduce elevation quantization. Slopes span approximately 300 m at Afar;
  the actual maximum quantization step is reported. This is not DEM accuracy
  or a fine-scale erosion measurement. MODIS is approximately 500 m and is
  sampled at grid centers; incomplete fire observations can miss burns.
- Risk transforms and weights are assumptions, not fitted parameters. Land use
  and climate are omitted because no reliable free browser source is integrated
  at this resolution. ESA WorldCover is a legitimate future land-use addition.
- Target scores are **relative within the current extent, not absolute** benefits.
  Grass share is not forage yield, candidate sample count is not canopy density,
  and a spectral water signal is not guaranteed accessible water. Clearing on
  slopes can worsen erosion without restoration. Costs, tenure and biodiversity
  tradeoffs are not modeled.
- **“Not verified” means the satellite record does not support the claim — it
  is not proof the work was not done.** Cloud, scene gaps and small plots cause
  false negatives. This is a screening flag for follow-up. “Verified” means a
  spectral drop only: drought, seasonal change, fire and cropping can also
  explain it. No usable before/after pair means inconclusive.

The baseline captures an analysis extent. Comparison, seasonal, risk and
targeting results retain that extent even if you pan elsewhere; re-screen A to
change it. No results are generated from placeholder imagery or random data.

## Sources, attribution and service limits

| Purpose | Anonymous endpoint / attribution |
| --- | --- |
| STAC searches | `https://planetarycomputer.microsoft.com/api/stac/v1/search` — Microsoft Planetary Computer |
| PNG data and dated imagery | `https://planetarycomputer.microsoft.com/api/data/v1/item/tiles/WebMercatorQuad/{z}/{x}/{y}@1x` |
| Spectral bands, SCL, visual imagery | Copernicus Sentinel-2 L2A, collection `sentinel-2-l2a` |
| Terrain | Copernicus DEM GLO-30, collection `cop-dem-glo-30`, asset `data` |
| Fire | NASA MODIS MCD64A1 v6.1, collection `modis-64A1-061`, asset `Burn_Date` |
| Context basemap | `https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}` — Esri, Maxar, Earthstar Geographics and GIS User Community |
| Place search | `https://nominatim.openstreetmap.org/search` — © OpenStreetMap contributors |
| Maps | Leaflet 1.9.4 (BSD-2-Clause), JavaScript/CSS from `https://unpkg.com/leaflet@1.9.4/` with integrity hashes |

Planetary Computer anonymous access is **rate-limited and intended for
research/non-commercial use**; consult its current service terms and each data
provider's license before deployment or commercial use. A displayed Esri basemap
is context imagery, not an observation from either selected date.

STAC calls use GET, including pagination; datetime slashes and expression plus
signs are URL-encoded. Data tiles request PNG, band expressions set
`asset_as_band=True`, and transparent pixels are skipped. Out-of-footprint tile
404s are normal and ignored; other failures are surfaced, not replaced with
synthetic measurements. Catalog results are never silently truncated. Requests
have timeouts, transient rate-limit retries, cancellation, and in-memory caches.

Nominatim search is submit-only, with a queued **hard 1100 ms minimum** between
requests in each page instance, stale submissions dropped, and visible refusal
messages. This is not a multi-user, cross-tab or hosting-wide rate limiter.
Public deployments must comply with the
[Nominatim usage policy](https://operations.osmfoundation.org/policies/nominatim/);
a heavily used public service may need a different geocoder, which is outside
this no-backend app's scope. Never use this as a bulk geocoding client.

Queries, map extents and record coordinates are sent only to the named services.
Field activity text stays in the page and exported CSV; nothing is persistently
stored by the app. Internet access and service CORS support are required even
when opening the file locally. Service availability, coverage, and quotas can
change.

## Deploy to GitHub Pages

The app needs only `/index.html` and this `/README.md`; existing research files
can remain in the repository.

1. Commit the files to `main`.
2. Open **Settings → Pages → Source: Deploy from a branch**.
3. Select **`main` / `(root)` → Save**.
4. Visit **https://soilwatch.github.io/Prosopis_mapping/** after Pages finishes
   publishing (usually a few minutes).

No Actions workflow, secrets, backend or build is needed. This change does not
itself enable Pages in repository settings. **Do not share a OneDrive/SharePoint
preview link:** its sandbox can block scripts and fetch requests and render a
blank app. Use GitHub Pages or open the HTML directly in a normal browser.

## Validation and real-data acceptance checklist

The repository has no existing web build, lint or test suite. Static JavaScript
and local browser checks are possible without adding dependencies. **Live Afar
acceptance remains unverified in the implementation environment:** the listed
public services failed DNS resolution before HTTP responses could be obtained.
This is an environment limitation, not evidence that those services are down.
No real scene counts, detected stands, geographic patterns, or successful MRV
claims are asserted here.

Local checks performed: inline JavaScript syntax, CSV/calendar/encoding/masking/
clustering invariants, isolated synthetic change and MRV outcomes, and twelve
monthly queries with refusal below five usable months. Headless Chromium checks
covered direct-file and local-HTTP startup, mobile width, the unavailable-CDN
fallback, and absence of uncaught exceptions on those paths. Synthetic fixtures
were validation inputs only, never app data or evidence of scientific accuracy.
Full map interaction and real-data checks below still require network access.

Before operational use, check Afar at **11.78 N, 40.96 E, zoom 13**:

- Place search returns results and moves both maps.
- Scene cards show actual dates, IDs and cloud percentages; verify the <5%
  outside-period recommendation rule when triggered.
- Screen a 2024 baseline: inspect thousands of sampled pixels and candidate
  clusters against watercourses/settlements and field evidence. Do not tune
  thresholds to force a predetermined count.
- Compare a later, seasonally matched observation and inspect new/stable/loss
  pins; a class may legitimately be absent. No observations are invented to
  force all three.
- Request January–December of a completed year; check all twelve month rows,
  actual observation range, missing months and per-location minimum coverage.
- Run the risk grid, inspect raw proxy values and coverage in cell popups, then
  run targeting and inspect the ranking/CSV. Fire can legitimately be zero.
- Use documented clearing records and documented non-clearing controls to
  exercise MRV outcomes and inspect the audit CSV. Do not manufacture claims
  solely to obtain a “verified” label.
- Drag the swipe divider with mouse/touch/keyboard; pan and zoom either side,
  resize the viewport, and confirm identical centers/zooms and independent pins.
- Inspect the console for uncaught exceptions; expected footprint tile 404s are
  not analysis errors. Check cancellation, service refusal and missing data.

---

# Research archive updates

# OUR OWN FIELD DATASET
- Find prosopis/non-prosopis polygons and centroids of polygons in different formats (.csv, .xlsx, .geojson) in the *"Ground_truth_data"* folder. This cames from an ODK form. The date of the file download can be found in its name (year_month_day)
- Find also a set of prosopis/non-prosopis gps points collected by SoilWatch team in different trips ("trips_gps_points"). 

# Datasets for ground truth: 
- prosopis juliflora (world registers download from https://www.gbif.org/)  - presence
- prosopis_juliflora_ethiopia (only Ethiopia registers download from https://www.gbif.org/)  - presence
- acacia-tortilis_balanites-aegyptiaca (Acacia tortilis and Balanites aegyptiaca only Ethiopia registers download from https://www.gbif.org/)  - used for abscense of prosopis

- prosopis_juliflora_combined_from2016_africa-occAsia: is the prosopis_juliflora dataset from GBIF combined with Biendata.org dataset. Keeped all Africa and occidental Asia points.
  
In this version of dataset GBIF data was filtered by:
1. valid coordinates, 
2. human observations or occurrence, 
3. equal or more than 4 individuals, 
4. removed points in the ocean, 
5. keep only 2016 to 2026 points.

Biendata.org dataset was filtered by:
1. no GBIF data, for not repeating, 
2. observation type= checklist, occurrence and trait occurence. 
  
# Models (run in Google Earth Engine) - ordered from new to old:
- version_20b_01-09-26: uses prosopis_juliflora_combined_from2016_africa-occAsia and acacia-tortilis_balanites-aegyptiaca (filtered dataset).
  In this version:
  - using new field data --> 13 hard negatives added in calibration and the rest added to validation
  - changed the land covers from where pseudoabscenses are generated
  - changed the use of FeatureMasked
  - added S1 and dynamic world variables
  - rainfall variable was fixed

- version_19_19-08-26: uses prosopis_juliflora_combined_from2016_africa-occAsia and acacia-tortilis_balanites-aegyptiaca (filtered dataset).
  In this version:
  - study area: whole Ethiopia and Kenya countries
  - EMIT was removed as it hasnt got complete cover in our new study area
  - pseudo-abscense data is taken at similar distance to roads as our own data (to limit the dist_road overfitting)
  - uses ONLINE data to calibrate and OWN GROUND TRUTH data to validate.
  - Dry and wet season by country: Ethiopia: Nov-Feb dry/Jun-Sep wet; Kenya: Jan-Feb dry/Mar-May wet
  - New spectral indexes and variables were added (as VIRE, TCARI, dist_built, etc.)

- version_6_13-05-26: uses prosopis_juliflora_combined_from2016_africa-occAsia and acacia-tortilis_balanites-aegyptiaca (filtered dataset).
  In this version:
  - a bigger dataset was used,
  - a buffer of 20m radius was used to sample locations and uses mean values of circular polygon for all bands;
  - climate (rainfall), soil data (soc) and distance to roads were included;
  - S2 cloudmask were added;
  - S2 NDWI was added;
  - EMIT bands and indices for dry season were added;
  - areas calculations were added.
    
 - version_1_17-04-26: uses prosopis_juliflora_ethiopia and acacia-tortilis_balanites-aegyptiaca.

    
# Python version (can be run in VS Code, Jupyter, or Google Colab):
- version_1_17-04-26: prosopis_mapping.py — direct port of version_1, modular function-based pipeline
- entry point: run_pipeline(CONFIG)

  
