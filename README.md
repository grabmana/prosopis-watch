Updates

Datasets for ground truth: 
- prosopis juliflora (world registers download from https://www.gbif.org/)  - presence
- prosopis_juliflora_ethiopia (only Ethiopia registers download from https://www.gbif.org/)  - presence
- acacia-tortilis_balanites-aegyptiaca (Acacia tortilis and Balanites aegyptiaca only Ethiopia registers download from https://www.gbif.org/ for )  - used for abscense of prosopis

Models (run in Google Earth Engine):
- version_1_17-04-26: uses prosopis_juliflora_ethiopia and acacia-tortilis_balanites-aegyptiaca.


Future changes to make to model:
- use a bigger dataset (world or africa) instead only Ethiopia
- remove very old registers
- make buffer of 20m radius to sample locations and use mean values of circular polygons for all bands
- try new dataset: https://www.biendata.org/
- add climate or soil data?
- add cloudmask of s2
- add distance to roads
- try using ndvi of dryest month

  
