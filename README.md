# Updates

# OUR OWN FIELD DATASET
- Find prosopis/non-prosopis polygons and centroids of polygons in different formats (.csv, .xlsx, .geojson) in the *"Ground_truth_data"* folder. 
- The date of the file download can be found in its name (year_month_day)

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
  
# Models (run in Google Earth Engine):
- version_1_17-04-26: uses prosopis_juliflora_ethiopia and acacia-tortilis_balanites-aegyptiaca.
- version_6_13-05-26: uses prosopis_juliflora_combined_from2016_africa-occAsia and acacia-tortilis_balanites-aegyptiaca (filtered dataset).

  In this version:
  - a bigger dataset was used,
  - a buffer of 20m radius was used to sample locations and uses mean values of circular polygon for all bands;
  - climate (rainfall), soil data (soc) and distance to roads were included;
  - S2 cloudmask were added;
  - S2 NDWI was added;
  - EMIT bands and indices for dry season were added;
  - areas calculations were added.

# Python version (can be run in VS Code, Jupyter, or Google Colab):
- version_1_17-04-26: prosopis_mapping.py — direct port of version_1, modular function-based pipeline
- entry point: run_pipeline(CONFIG)

# Future changes to make to v6 model:
- adjust variables of model when no land cover mask is used to improve accuracy
- use our own ground truth data!
- try using different land cover layer (dynamic world instead of worl cover for example)
- try Sentinel1 (radar)
- try GEDI
  
