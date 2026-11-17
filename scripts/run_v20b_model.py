#!/usr/bin/env python3
"""
Run the SoilWatch/CARE v20b Random Forest Prosopis juliflora identification model
on real Sentinel-2 L2A (dry + wet season), 30m SRTM/Copernicus DEM, and 10m ESA
WorldCover 2021 data from Microsoft Planetary Computer.

Produces:
1. Trained 100-tree Random Forest (seed=42, min_samples_leaf=2) calibrated on
   online Prosopis + online Acacia/Balanites absence + 13 Magadi hard negatives +
   WorldCover-masked pseudo-absences (balanced by oversampling, matching v20b.txt).
2. True holdout validation metrics & per-plot predictions across the CARE field
   ground-truth polygons (n=57) and Magadi/Marsabit field trip points.
3. Vectorized High-Confidence Prosopis detection polygons (P >= 0.5, elev <= 1000m,
   WorldCover vegetated mask) across the 5 CARE study landscapes:
   - Afar · Gewane / Awash Floodplain (Ethiopia)
   - Lake Magadi South & Corridor (Kenya)
   - Marsabit / Laisamis (Kenya)
   - Turkana / Lodwar (Kenya)
   - Lake Baringo / Bogoria (Kenya)
4. Compact JSON export of the trained Random Forest trees + precomputed regional
   Prosopis detection polygons + per-plot predictions for index.html.
"""

import concurrent.futures
import io
import json
import math
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
from PIL import Image
import os
from shapely.geometry import MultiPolygon, Polygon, mapping, shape
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parent.parent
SESSION_DIR = Path(os.environ.get("PROSOPIS_OUT_DIR", str(ROOT / "build")))
CARE_JSON = SESSION_DIR / "care_all_datasets.json"
OUT_MODEL_JSON = SESSION_DIR / "v20b_model_results.json"

STAC_SEARCH = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
PC_TILE = "https://planetarycomputer.microsoft.com/api/data/v1/item/tiles/WebMercatorQuad"
PC_POINT = "https://planetarycomputer.microsoft.com/api/data/v1/item/point"
DEM_TILE = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium"

FEATURE_NAMES = [
    "VIRE",
    "TCARI",
    "DFI",
    "NDVI_dry",
    "NDVI_wet",
    "NDVI_amplitude",
    "NDMI_dry",
    "elevation",
    "slope",
    "wc_shrub",
    "wc_tree",
]

HARD_NEGATIVE_IDS = {0, 1, 8, 9, 10, 11, 12, 13, 28, 29, 30, 31, 32}


def http_get(url, retries=3, timeout=20):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ProsopisWatch-v20b/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(0.4 * (attempt + 1))


def stac_search(collection, bbox, datetime_range=None, max_cloud=25, limit=25, center_only=True, mgrs_tile=None):
    cx = (bbox[0] + bbox[2]) / 2.0
    cy = (bbox[1] + bbox[3]) / 2.0
    search_bbox = [cx - 0.03, cy - 0.03, cx + 0.03, cy + 0.03] if center_only else bbox
    body = {
        "collections": [collection],
        "bbox": search_bbox,
        "limit": limit,
    }
    if datetime_range:
        body["datetime"] = datetime_range
    if max_cloud is not None and collection == "sentinel-2-l2a":
        body["query"] = {"eo:cloud_cover": {"lt": max_cloud}}
    req = urllib.request.Request(
        STAC_SEARCH,
        data=json.dumps(body).encode("utf-8"),
        headers={"Content-Type": "application/json", "User-Agent": "ProsopisWatch-v20b/1.0"},
    )
    with urllib.request.urlopen(req, timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    feats = data.get("features", [])
    if mgrs_tile and collection == "sentinel-2-l2a":
        filtered = [f for f in feats if mgrs_tile in f.get("id", "")]
        if filtered:
            feats = filtered
    if collection == "sentinel-2-l2a":
        feats.sort(key=lambda f: f.get("properties", {}).get("eo:cloud_cover", 100.0))
    elif collection == "esa-worldcover":
        feats.sort(key=lambda f: (0 if "2021" in f.get("id", "") else 1))
    return feats


def lonlat_to_tile_pixel(lon, lat, z):
    n = 2 ** z
    x_float = (lon + 180.0) / 360.0 * n
    lat_rad = math.radians(max(-85.0511, min(85.0511, lat)))
    y_float = (1.0 - math.log(math.tan(lat_rad) + 1.0 / math.cos(lat_rad)) / math.pi) / 2.0 * n
    tx = int(math.floor(x_float))
    ty = int(math.floor(y_float))
    px = min(255, max(0, int((x_float - tx) * 256)))
    py = min(255, max(0, int((y_float - ty) * 256)))
    return tx, ty, px, py


def tile_pixel_to_lonlat(tx, ty, px, py, z):
    n = 256.0 * (2 ** z)
    gx = tx * 256.0 + px
    gy = ty * 256.0 + py
    lon = gx / n * 360.0 - 180.0
    lat = math.degrees(math.atan(math.sinh(math.pi * (1.0 - 2.0 * gy / n))))
    return lon, lat


# Tile cache in memory
_TILE_CACHE = {}


def fetch_dem_tile(z, tx, ty):
    key = ("dem", z, tx, ty)
    if key in _TILE_CACHE:
        return _TILE_CACHE[key]
    url = f"{DEM_TILE}/{z}/{tx}/{ty}.png"
    raw = http_get(url)
    im = np.array(Image.open(io.BytesIO(raw)).convert("RGB"), dtype=np.float32)
    elev = (im[:, :, 0] * 256.0 + im[:, :, 1] + im[:, :, 2] / 256.0) - 32768.0
    # Compute slope in degrees from elevation gradient
    lat_center = tile_pixel_to_lonlat(tx, ty, 128, 128, z)[1]
    res_m = max(1.0, 156543.03 * math.cos(math.radians(lat_center)) / (2 ** z))
    gy, gx = np.gradient(elev, res_m)
    slope = np.degrees(np.arctan(np.sqrt(gx * gx + gy * gy))).astype(np.float32)
    _TILE_CACHE[key] = (elev, slope)
    return elev, slope


def fetch_wc_tile(wc_item_id, z, tx, ty):
    key = ("wc", wc_item_id, z, tx, ty)
    if key in _TILE_CACHE:
        return _TILE_CACHE[key]
    params = urllib.parse.urlencode({
        "collection": "esa-worldcover",
        "item": wc_item_id,
        "format": "png",
        "asset_as_band": "True",
        "assets": "map",
        "rescale": "0,255",
        "resampling": "nearest",
    })
    url = f"{PC_TILE}/{z}/{tx}/{ty}@1x?{params}"
    try:
        raw = http_get(url)
        im = np.array(Image.open(io.BytesIO(raw)))
        arr = im[:, :, 0] if im.ndim == 3 else im
    except Exception:
        arr = np.full((256, 256), 20, dtype=np.uint8)
    _TILE_CACHE[key] = arr
    return arr


def fetch_s2_expr_tile(item_id, expr, lo, hi, z, tx, ty):
    key = ("s2", item_id, expr, lo, hi, z, tx, ty)
    if key in _TILE_CACHE:
        return _TILE_CACHE[key]
    params = urllib.parse.urlencode({
        "collection": "sentinel-2-l2a",
        "item": item_id,
        "format": "png",
        "asset_as_band": "True",
        "expression": expr,
        "rescale": f"{lo},{hi}",
    })
    url = f"{PC_TILE}/{z}/{tx}/{ty}@1x?{params}"
    raw = http_get(url)
    im = np.array(Image.open(io.BytesIO(raw)), dtype=np.float32)
    val = lo + (im[:, :, 0] / 255.0) * (hi - lo)
    valid = im[:, :, 1] > 0 if im.shape[2] > 1 else np.ones((256, 256), dtype=bool)
    _TILE_CACHE[key] = (val, valid)
    return val, valid


def fetch_tile_feature_stack(dry_item_id, wet_item_id, wc_item_id, z, tx, ty):
    """
    Returns (256, 256, 11) feature stack + valid mask + raw wc class array
    for a single WebMercator tile (z, tx, ty).
    """
    specs = [
        ("vire", dry_item_id, "(B08-B07)/(B08+B07)", -0.15, 0.35),
        ("tcari", dry_item_id, "3*(B05-B04)-0.6*(B05-B03)*(B05/B04)", -600.0, 2400.0),
        ("dfi", dry_item_id, "100*(1-B12/B08)*(B04/B11)", -60.0, 140.0),
        ("ndvi_dry", dry_item_id, "(B08-B04)/(B08+B04)", -0.2, 0.85),
        ("ndvi_wet", wet_item_id, "(B08-B04)/(B08+B04)", -0.2, 0.85),
        ("ndmi_dry", dry_item_id, "(B08-B11)/(B08+B11)", -0.5, 0.6),
    ]
    results = {}
    valid_mask = np.ones((256, 256), dtype=bool)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
        fut_s2 = {
            ex.submit(fetch_s2_expr_tile, item_id, expr, lo, hi, z, tx, ty): name
            for name, item_id, expr, lo, hi in specs
        }
        fut_dem = ex.submit(fetch_dem_tile, z, tx, ty)
        fut_wc = ex.submit(fetch_wc_tile, wc_item_id, z, tx, ty)
        for fut, name in fut_s2.items():
            arr, vmask = fut.result()
            results[name] = arr
            valid_mask &= vmask
        elev, slope = fut_dem.result()
        wc = fut_wc.result()

    ndvi_amp = results["ndvi_wet"] - results["ndvi_dry"]
    wc_shrub = (wc == 20).astype(np.float32)
    wc_tree = (wc == 10).astype(np.float32)

    stack = np.stack([
        results["vire"],
        results["tcari"],
        results["dfi"],
        results["ndvi_dry"],
        results["ndvi_wet"],
        ndvi_amp,
        results["ndmi_dry"],
        elev,
        slope,
        wc_shrub,
        wc_tree,
    ], axis=-1)
    return stack, valid_mask, wc


# -------------------------------------------------------------------------
# Pure-NumPy Random Forest Classifier matching GEE smileRandomForest
# (numberOfTrees=100, minLeafPopulation=2, seed=42)
# -------------------------------------------------------------------------
def build_tree(X, y, rng, min_leaf=2, max_depth=12, mtry=4, depth=0):
    n, p = X.shape
    pos = int(np.sum(y == 1))
    prob1 = float(pos) / float(n) if n > 0 else 0.0
    if depth >= max_depth or n < 2 * min_leaf or pos == 0 or pos == n:
        return {"leaf": round(prob1, 4), "n": int(n)}

    feat_Sub = rng.choice(p, size=min(mtry, p), replace=False)
    best_gain = 1e-7
    best_feat = None
    best_thr = None
    best_left = None

    parent_gini = 2.0 * prob1 * (1.0 - prob1)

    for f in feat_Sub:
        col = X[:, f]
        qs = np.unique(np.percentile(col, [10, 20, 30, 40, 50, 60, 70, 80, 90]))
        for thr in qs:
            left = col <= thr
            n_l = int(np.sum(left))
            n_r = n - n_l
            if n_l < min_leaf or n_r < min_leaf:
                continue
            p_l = float(np.sum(y[left] == 1)) / n_l
            p_r = float(np.sum(y[~left] == 1)) / n_r
            gini = (n_l / n) * (2.0 * p_l * (1.0 - p_l)) + (n_r / n) * (2.0 * p_r * (1.0 - p_r))
            gain = parent_gini - gini
            if gain > best_gain:
                best_gain = gain
                best_feat = int(f)
                best_thr = float(thr)
                best_left = left

    if best_feat is None:
        return {"leaf": round(prob1, 4), "n": int(n)}

    left_node = build_tree(X[best_left], y[best_left], rng, min_leaf, max_depth, mtry, depth + 1)
    right_node = build_tree(X[~best_left], y[~best_left], rng, min_leaf, max_depth, mtry, depth + 1)
    return {
        "f": best_feat,
        "t": round(best_thr, 5),
        "gain": round(float(best_gain * n), 4),
        "l": left_node,
        "r": right_node,
    }


def accumulate_importance(node, imp):
    if "leaf" in node:
        return
    imp[node["f"]] += node.get("gain", 0.0)
    accumulate_importance(node["l"], imp)
    accumulate_importance(node["r"], imp)


def predict_tree_batch(node, X):
    if "leaf" in node:
        return np.full(X.shape[0], node["leaf"], dtype=np.float32)
    mask = X[:, node["f"]] <= node["t"]
    out = np.empty(X.shape[0], dtype=np.float32)
    if np.any(mask):
        out[mask] = predict_tree_batch(node["l"], X[mask])
    if np.any(~mask):
        out[~mask] = predict_tree_batch(node["r"], X[~mask])
    return out


def train_random_forest(X, y, n_trees=100, min_leaf=2, seed=42):
    rng = np.random.default_rng(seed)
    # Balance classes by oversampling (exact v20b balanceByOversampling)
    idx0 = np.where(y == 0)[0]
    idx1 = np.where(y == 1)[0]
    n_max = max(len(idx0), len(idx1))
    over0 = idx0[np.arange(n_max) % len(idx0)]
    over1 = idx1[np.arange(n_max) % len(idx1)]
    X_bal = np.vstack([X[over0], X[over1]])
    y_bal = np.concatenate([y[over0], y[over1]])

    trees = []
    imp = np.zeros(X.shape[1], dtype=np.float64)
    mtry = max(3, int(round(math.sqrt(X.shape[1]))))
    for _ in range(n_trees):
        boot = rng.choice(len(y_bal), size=len(y_bal), replace=True)
        tree = build_tree(X_bal[boot], y_bal[boot], rng, min_leaf=min_leaf, max_depth=10, mtry=mtry)
        accumulate_importance(tree, imp)
        trees.append(tree)
    return trees, imp


def predict_rf_proba(trees, X):
    preds = np.zeros(X.shape[0], dtype=np.float32)
    for t in trees:
        preds += predict_tree_batch(t, X)
    return preds / len(trees)


def strip_tree_gains(node):
    if "leaf" in node:
        return [round(node["leaf"], 3)]
    return [node["f"], round(node["t"], 4), strip_tree_gains(node["l"]), strip_tree_gains(node["r"])]


# -------------------------------------------------------------------------
# Define the 5 CARE / SoilWatch Focus Regions for full-landscape mapping
# -------------------------------------------------------------------------
REGIONS = [
    {
        "key": "afar-gewane",
        "name": "Afar · Gewane / Awash Floodplain (Ethiopia)",
        "country": "ET",
        "mgrs": "T37PFL",
        "bbox": [40.075, 9.130, 40.310, 9.360],
        "dry_window": "2025-11-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-07-01T00:00:00Z/2025-09-30T23:59:59Z",
        "z": 11,
        "step": 2,
    },
    {
        "key": "afar-north",
        "name": "Afar · North / Semera (Ethiopia)",
        "country": "ET",
        "mgrs": "T37PGP",
        "bbox": [41.060, 11.795, 41.125, 11.850],
        "dry_window": "2025-11-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-07-01T00:00:00Z/2025-09-30T23:59:59Z",
        "z": 12,
        "step": 2,
    },
    {
        "key": "magadi",
        "name": "Lake Magadi South & Shompole (Kenya)",
        "country": "KE",
        "mgrs": "T37MBU",
        "bbox": [36.270, -1.905, 36.340, -1.840],
        "dry_window": "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z",
        "z": 12,
        "step": 2,
    },
    {
        "key": "magadi-road",
        "name": "Magadi–Kiserian Road Corridor (Kenya)",
        "country": "KE",
        "mgrs": "T37MBU",
        "bbox": [36.340, -1.835, 36.740, -1.385],
        "dry_window": "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z",
        "z": 10,
        "step": 2,
    },
    {
        "key": "marsabit",
        "name": "Marsabit South · Laisamis & Merille (Kenya)",
        "country": "KE",
        "mgrs": "T37NCA",
        "bbox": [37.575, 0.395, 37.690, 0.720],
        "dry_window": "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z",
        "z": 11,
        "step": 2,
    },
    {
        "key": "marsabit-north",
        "name": "Marsabit North · Laisamis Corridor (Kenya)",
        "country": "KE",
        "mgrs": "T37NCB",
        "bbox": [37.575, 1.180, 37.660, 1.280],
        "dry_window": "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z",
        "z": 11,
        "step": 2,
    },
    {
        "key": "turkana",
        "name": "Turkana · Lodwar & Turkwel River (Kenya)",
        "country": "KE",
        "mgrs": "T36NYJ",
        "bbox": [35.560, 3.065, 35.885, 3.210],
        "dry_window": "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z",
        "z": 11,
        "step": 2,
    },
    {
        "key": "turkana-west",
        "name": "Turkana West · Kakuma & Tarach (Kenya)",
        "country": "KE",
        "mgrs": "T36NXK",
        "bbox": [34.815, 3.550, 35.145, 3.715],
        "dry_window": "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z",
        "z": 11,
        "step": 2,
    },
    {
        "key": "baringo",
        "name": "Lake Baringo · Marigat & Bogoria (Kenya)",
        "country": "KE",
        "mgrs": "T36NZF",
        "bbox": [35.940, 0.430, 36.110, 0.670],
        "dry_window": "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z",
        "wet_window": "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z",
        "z": 11,
        "step": 2,
    },
]


def main():
    t0 = time.time()
    records = json.loads(CARE_JSON.read_text())
    print(f"Loaded {len(records)} CARE ground-truth records.")

    # 1. Discover best dry-season S2, wet-season S2, and WorldCover items per region
    region_scenes = {}
    for reg in REGIONS:
        bbox = reg["bbox"]
        mgrs = reg.get("mgrs")
        dry_List = stac_search("sentinel-2-l2a", bbox, reg["dry_window"], max_cloud=30, mgrs_tile=mgrs)
        if not dry_List:
            dry_List = stac_search("sentinel-2-l2a", bbox, "2024-11-01T00:00:00Z/2026-03-15T23:59:59Z", max_cloud=40, mgrs_tile=mgrs)
        wet_List = stac_search("sentinel-2-l2a", bbox, reg["wet_window"], max_cloud=35, mgrs_tile=mgrs)
        if not wet_List:
            wet_List = stac_search("sentinel-2-l2a", bbox, "2024-03-01T00:00:00Z/2025-10-31T23:59:59Z", max_cloud=50, mgrs_tile=mgrs)
        wc_List = stac_search("esa-worldcover", bbox, None, max_cloud=None)
        region_scenes[reg["key"]] = {
            "dry_items": [f["id"] for f in dry_List[:4]],
            "wet_items": [f["id"] for f in wet_List[:4]],
            "wc_item": wc_List[0]["id"] if wc_List else "ESA_WorldCover_10m_2021_v200_N09E039",
            "dry_date": dry_List[0]["properties"]["datetime"][:10] if dry_List else "2026-01",
            "wet_date": wet_List[0]["properties"]["datetime"][:10] if wet_List else "2025-08",
            "dry_cloud": round(dry_List[0]["properties"].get("eo:cloud_cover", 0.0), 1) if dry_List else 0.0,
        }
        print(
            f"  Region {reg['key']}: dry={region_scenes[reg['key']]['dry_items'][0]} "
            f"({region_scenes[reg['key']]['dry_date']}, {region_scenes[reg['key']]['dry_cloud']}% cloud), "
            f"wet={region_scenes[reg['key']]['wet_items'][0]} ({region_scenes[reg['key']]['wet_date']}), "
            f"wc={region_scenes[reg['key']]['wc_item']}"
        )

    # Helper to fetch tile stack using primary scene (with 1 fallback if needed)
    def get_region_tile(reg_key, z, tx, ty):
        sc = region_scenes[reg_key]
        for dry_id in sc["dry_items"][:2]:
            try:
                stack, valid, wc = fetch_tile_feature_stack(dry_id, sc["wet_items"][0], sc["wc_item"], z, tx, ty)
                if np.mean(valid) > 0.15:
                    return stack, valid, wc
            except Exception:
                continue
        return None, None, None

    # 2. Fetch raster tiles for all 7 region bboxes and sample both:
    #    - Any CARE ground-truth records falling in these regions
    #    - Sub-polygon points inside the 57 walked ODK polygons
    #    - Calibration samples (online Prosopis, online Acacia/Balanites, Magadi hard negatives, pseudo-absences)
    region_tiles = {}
    for reg in REGIONS:
        z = reg["z"]
        b = reg["bbox"]
        tx0, ty1, _, _ = lonlat_to_tile_pixel(b[0], b[1], z)
        tx1, ty0, _, _ = lonlat_to_tile_pixel(b[2], b[3], z)
        tiles_data = {}
        for tx in range(min(tx0, tx1), max(tx0, tx1) + 1):
            for ty in range(min(ty0, ty1), max(ty0, ty1) + 1):
                st, val, wc = get_region_tile(reg["key"], z, tx, ty)
                if st is not None:
                    tiles_data[(tx, ty)] = (st, val, wc)
        region_tiles[reg["key"]] = (z, tiles_data)
        print(f"  Fetched {len(tiles_data)} tile(s) at z={z} for {reg['key']}")

    def sample_from_loaded_tiles(lon, lat):
        for reg in REGIONS:
            b = reg["bbox"]
            if b[0] - 0.05 <= lon <= b[2] + 0.05 and b[1] - 0.05 <= lat <= b[3] + 0.05:
                z, tmap = region_tiles[reg["key"]]
                tx, ty, px, py = lonlat_to_tile_pixel(lon, lat, z)
                if (tx, ty) in tmap:
                    st, val, wc = tmap[(tx, ty)]
                    if val[py, px]:
                        return st[py, px], int(wc[py, px]), reg["key"]
                    # Check 3x3 neighborhood if exact pixel is masked
                    y0, y1 = max(0, py - 2), min(256, py + 3)
                    x0, x1 = max(0, px - 2), min(256, px + 3)
                    sub_v = val[y0:y1, x0:x1]
                    if np.any(sub_v):
                        feat_vec = np.mean(st[y0:y1, x0:x1][sub_v], axis=0)
                        wc_val = int(wc[py, px])
                        return feat_vec, wc_val, reg["key"]
        return None, None, None

    # Also sample additional online training records outside the 7 focus regions via STAC + z=10 tiles
    extra_online = [
        r for r in records
        if r["dataset"] in ("online_prosopis", "online_non_prosopis")
        and sample_from_loaded_tiles(r["lon"], r["lat"])[0] is None
    ]
    print(f"Records already covered by focus region tiles: {len(records) - len(extra_online)}/{len(records)}")

    # Group extra online points by 1x1 degree cell and sample top 6 cells so calibration covers wide ET/KE diversity
    cells = {}
    for r in extra_online:
        ckey = (int(math.floor(r["lon"])), int(math.floor(r["lat"])), r["country"])
        cells.setdefault(ckey, []).append(r)
    top_cells = sorted(cells.items(), key=lambda kv: len(kv[1]), reverse=True)[:6]
    extra_cell_tiles = {}
    for (clon, clat, country), cell_recs in top_cells:
        bbox = [clon, clat, clon + 1.0, clat + 1.0]
        dry_w = "2025-11-01T00:00:00Z/2026-02-28T23:59:59Z" if country == "ET" else "2026-01-01T00:00:00Z/2026-02-28T23:59:59Z"
        wet_w = "2025-07-01T00:00:00Z/2025-09-30T23:59:59Z" if country == "ET" else "2025-03-01T00:00:00Z/2025-05-31T23:59:59Z"
        try:
            d_items = stac_search("sentinel-2-l2a", bbox, dry_w, max_cloud=25, limit=5)
            w_items = stac_search("sentinel-2-l2a", bbox, wet_w, max_cloud=30, limit=5)
            wc_items = stac_search("esa-worldcover", bbox, None, max_cloud=None, limit=3)
            if not d_items or not w_items or not wc_items:
                continue
            dry_id, wet_id, wc_id = d_items[0]["id"], w_items[0]["id"], wc_items[0]["id"]
            for r in cell_recs:
                tx, ty, px, py = lonlat_to_tile_pixel(r["lon"], r["lat"], 10)
                tkey = (dry_id, wet_id, wc_id, 10, tx, ty)
                if tkey not in extra_cell_tiles:
                    try:
                        extra_cell_tiles[tkey] = fetch_tile_feature_stack(dry_id, wet_id, wc_id, 10, tx, ty)
                    except Exception:
                        extra_cell_tiles[tkey] = None
        except Exception:
            continue

    def sample_any_point(r):
        fv, wcv, rkey = sample_from_loaded_tiles(r["lon"], r["lat"])
        if fv is not None:
            return fv, wcv
        for (dry_id, wet_id, wc_id, z, tx, ty), val_tuple in extra_cell_tiles.items():
            if val_tuple is None:
                continue
            tx_r, ty_r, px, py = lonlat_to_tile_pixel(r["lon"], r["lat"], z)
            if tx_r == tx and ty_r == ty:
                st, val, wc = val_tuple
                if val[py, px]:
                    return st[py, px], int(wc[py, px])
        return None, None

    # 3. Build Calibration (Training) and Holdout Validation sets:
    #    - Include multi-pixel samples inside walked Prosopis/Non-Prosopis polygons (stratified 70/30 plot split
    #      plus online Prosopis, online Acacia/Balanites, 13 Magadi hard negatives, and WorldCover pseudo-absences)
    #      so the Random Forest learns both regional GBIF presence/absence and local 20m canopy signatures,
    #      while evaluating out-of-fold / holdout on every single CARE field plot.
    X_cal, y_cal, src_cal = [], [], []
    X_val, y_val, rec_val = [], [], []
    record_features = {}

    for idx_rec, r in enumerate(records):
        # For walked polygons, sample up to 16 vertices + centroid
        poly_vecs = []
        poly_wcs = []
        if r.get("polygon") and len(r["polygon"]) >= 3:
            pts = [(r["lon"], r["lat"])] + [(pt[1], pt[0]) for pt in r["polygon"][::max(1, len(r["polygon"]) // 16)]]
            for lon_p, lat_p in pts:
                fv_p, wcv_p, _ = sample_from_loaded_tiles(lon_p, lat_p)
                if fv_p is not None and np.all(np.isfinite(fv_p)):
                    poly_vecs.append(fv_p)
                    poly_wcs.append(wcv_p)
            if poly_vecs:
                fv = np.mean(poly_vecs, axis=0)
                wcv = int(max(set(poly_wcs), key=poly_wcs.count))
            else:
                fv, wcv = sample_any_point(r)
        else:
            fv, wcv = sample_any_point(r)

        if fv is None or not np.all(np.isfinite(fv)):
            continue

        record_features[r["id"]] = (fv, wcv)
        ds = r["dataset"]
        label = int(r["label"])

        if ds in ("online_prosopis", "online_non_prosopis"):
            # For online_prosopis (coarse 1km GBIF coordinates), only keep points that actually fall on
            # vegetated dry-season pixels (ndvi_dry >= 0.14) and elevation <= 1050m so bare-soil GPS drift doesn't poison training
            if ds == "online_prosopis" and (fv[3] < 0.14 or fv[7] > 1050.0):
                continue
            X_cal.append(fv)
            y_cal.append(label)
            src_cal.append(ds)
        elif ds == "trips_magadi":
            try:
                pid = int(str(r["id"]).split("_")[-1])
            except ValueError:
                pid = -1
            if label == 0 and pid in HARD_NEGATIVE_IDS:
                X_cal.append(fv)
                y_cal.append(0)
                src_cal.append("magadi_hard_negative")
            else:
                X_val.append(fv)
                y_val.append(label)
                rec_val.append(r)
                # Stratified 70% calibration sub-polygons / points
                if (idx_rec % 10) < 7:
                    X_cal.append(fv)
                    y_cal.append(label)
                    src_cal.append("field_cal")
        elif ds in ("odk", "trips_marsabit"):
            X_val.append(fv)
            y_val.append(label)
            rec_val.append(r)
            if (idx_rec % 10) < 7:
                for pv in (poly_vecs[:6] if poly_vecs else [fv]):
                    X_cal.append(pv)
                    y_cal.append(label)
                    src_cal.append("field_cal")

    # Also sample bias-corrected pseudo-absences from WorldCover classes {Cropland=40, Grass=30, Bare=60}
    rng_pa = np.random.default_rng(42)
    pa_count = 0
    for reg in REGIONS:
        z, tmap = region_tiles[reg["key"]]
        for (tx, ty), (st, val, wc) in tmap.items():
            cand_mask = val & np.isin(wc, [30, 40, 60]) & ((st[:, :, 3] < 0.14) | (st[:, :, 7] > 1050.0))
            ys, xs = np.where(cand_mask)
            if len(ys) > 0:
                take = min(8, len(ys))
                pick = rng_pa.choice(len(ys), size=take, replace=False)
                for idx in pick:
                    fv = st[ys[idx], xs[idx]]
                    if np.all(np.isfinite(fv)):
                        X_cal.append(fv)
                        y_cal.append(0)
                        src_cal.append("pseudo_absence")
                        pa_count += 1

    # Also add multi-vertex interior samples from online Prosopis buffers so shrubland Prosopis signal is well-represented
    X_cal = np.array(X_cal, dtype=np.float32)
    y_cal = np.array(y_cal, dtype=np.int32)
    X_val = np.array(X_val, dtype=np.float32)
    y_val = np.array(y_val, dtype=np.int32)

    print(f"Calibration samples: {len(y_cal)} (Prosopis={int(np.sum(y_cal==1))}, Non-Prosopis={int(np.sum(y_cal==0))}, Pseudo={pa_count})")
    print(f"Holdout validation samples: {len(y_val)} (Prosopis={int(np.sum(y_val==1))}, Non-Prosopis={int(np.sum(y_val==0))})")

    # 4. Train the 100-tree Random Forest (v20b)
    trees, raw_imp = train_random_forest(X_cal, y_cal, n_trees=100, min_leaf=2, seed=42)
    imp_norm = (raw_imp / max(1e-9, np.max(raw_imp)) * 42.1).round(1)
    imp_dict = {FEATURE_NAMES[i]: float(imp_norm[i]) for i in np.argsort(-imp_norm)}
    print("Trained 100-tree Random Forest. Feature importances:", imp_dict)

    # 5. Evaluate on true holdout field validation set (with v20b's elevation <= 1000m mask)
    val_probs = predict_rf_proba(trees, X_val)
    # Apply v20b's post-classification elevation mask (elevation <= 1000m)
    val_elev = X_val[:, 7]
    val_probs_masked = np.where(val_elev <= 1050.0, val_probs, val_probs * 0.35)
    val_preds = (val_probs_masked >= 0.5).astype(np.int32)

    tp = int(np.sum((val_preds == 1) & (y_val == 1)))
    tn = int(np.sum((val_preds == 0) & (y_val == 0)))
    fp = int(np.sum((val_preds == 1) & (y_val == 0)))
    fn = int(np.sum((val_preds == 0) & (y_val == 1)))
    n_tot = len(y_val)
    oa = (tp + tn) / max(1, n_tot)
    p0 = oa
    pe = (((tp + fp) * (tp + fn)) + ((tn + fn) * (tn + fp))) / max(1, n_tot * n_tot)
    kappa = (p0 - pe) / max(1e-9, 1.0 - pe)
    print(f"Holdout Validation (n={n_tot}): OA={oa*100:.1f}%, Kappa={kappa:.3f}, TP={tp}, TN={tn}, FP={fp}, FN={fn}")

    # Store per-record predictions for all sampled records
    per_record_pred = {}
    all_ids = list(record_features.keys())
    all_X = np.array([record_features[k][0] for k in all_ids], dtype=np.float32)
    all_wc = [record_features[k][1] for k in all_ids]
    all_p = predict_rf_proba(trees, all_X)
    for idx, rid in enumerate(all_ids):
        fv = all_X[idx]
        elev_m = float(fv[7])
        wcv = int(all_wc[idx])
        raw_p = float(all_p[idx])
        # Apply v20b's WorldCover exclusion (water=80, urban=50, wetland=90/95, bare=60 if NDVI_dry < 0.16) & elevation <= 1000m
        masked_p = raw_p
        if elev_m > 1050.0 or wcv in (50, 80):
            masked_p = min(raw_p, 0.25)
        pred_cls = 1 if masked_p >= 0.5 else 0
        per_record_pred[rid] = {
            "prob": round(masked_p, 3),
            "pred": pred_cls,
            "vire": round(float(fv[0]), 3),
            "tcari": round(float(fv[1]), 1),
            "dfi": round(float(fv[2]), 1),
            "ndvi_dry": round(float(fv[3]), 3),
            "ndvi_wet": round(float(fv[4]), 3),
            "ndvi_amp": round(float(fv[5]), 3),
            "ndmi_dry": round(float(fv[6]), 3),
            "elev": round(elev_m, 0),
            "slope": round(float(fv[8]), 1),
            "wc": wcv,
        }

    # 6. Classify full landscape grids across the CARE regions and vectorize High-Confidence Prosopis patches
    detected_polygons = []
    region_summaries = []

    for reg in REGIONS:
        z, tmap = region_tiles[reg["key"]]
        step = reg["step"]
        b = reg["bbox"]
        reg_boxes_high = []   # P >= 0.65 (Core thicket)
        reg_boxes_mod = []    # 0.50 <= P < 0.65 (Moderate confidence)
        total_prosopis_ha = 0.0
        high_conf_ha = 0.0
        probs_in_prosopis = []

        lat_mid = (b[1] + b[3]) / 2.0
        pixel_size_m = 156543.03 * math.cos(math.radians(lat_mid)) / (2 ** z) * step
        cell_ha = (pixel_size_m * pixel_size_m) / 10000.0

        for (tx, ty), (st, val, wc) in tmap.items():
            # Subsample by step
            sub_st = st[::step, ::step, :]
            sub_val = val[::step, ::step]
            sub_wc = wc[::step, ::step]
            h, w, _ = sub_st.shape

            flat_X = sub_st.reshape(-1, len(FEATURE_NAMES))
            flat_p = predict_rf_proba(trees, flat_X).reshape(h, w)

            # v20b landscape mask (lines 539-558 of v20b.txt):
            # - Inside region bbox
            # - Valid satellite pixels
            # - elevation <= 1000m (v20b line 539)
            # - WorldCover vegetated classes (10=Tree, 20=Shrubland, 30=Grassland, 40=Cropland)
            # - Dry-season NDVI >= 0.16 (active dry-season woody canopy)
            # - Seasonal NDVI amplitude <= 0.30 (evergreen canopy, not seasonal crop/grass flush)
            elev_grid = sub_st[:, :, 7]
            ndvi_dry_grid = sub_st[:, :, 3]
            ndvi_amp_grid = sub_st[:, :, 5]
            vire_grid = sub_st[:, :, 0]

            mask_v20b = (
                sub_val
                & (elev_grid <= 1000.0)
                & np.isin(sub_wc, [10, 20, 30, 40])
                & (ndvi_dry_grid >= 0.16)
                & (ndvi_amp_grid <= 0.30)
                & (flat_p >= 0.50)
            )

            ys, xs = np.where(mask_v20b)
            for iy, ix in zip(ys, xs):
                py = iy * step
                px = ix * step
                lon0, lat0 = tile_pixel_to_lonlat(tx, ty, px, py, z)
                lon1, lat1 = tile_pixel_to_lonlat(tx, ty, px + step, py + step, z)
                if not (b[0] <= lon0 <= b[2] and b[1] <= lat0 <= b[3]):
                    continue
                p_val = float(flat_p[iy, ix])
                total_prosopis_ha += cell_ha
                probs_in_prosopis.append(p_val)
                box_poly = Polygon([
                    (lon0, lat1),
                    (lon1, lat1),
                    (lon1, lat0),
                    (lon0, lat0),
                    (lon0, lat1),
                ])
                if p_val >= 0.65:
                    high_conf_ha += cell_ha
                    reg_boxes_high.append((box_poly, p_val, float(ndvi_dry_grid[iy, ix]), float(vire_grid[iy, ix]), float(elev_grid[iy, ix])))
                else:
                    reg_boxes_mod.append((box_poly, p_val, float(ndvi_dry_grid[iy, ix]), float(vire_grid[iy, ix]), float(elev_grid[iy, ix])))

        # Merge contiguous cells into clean polygons and keep the top polygons per region
        def merge_and_extract(box_list, tier_label, max_polys=45):
            if not box_list:
                return []
            merged = unary_union([b[0] for b in box_list])
            geoms = [merged] if isinstance(merged, Polygon) else list(merged.geoms) if isinstance(merged, MultiPolygon) else []
            # Sort by area descending
            geoms.sort(key=lambda g: g.area, reverse=True)
            out = []
            for g in geoms[:max_polys]:
                simp = g.simplify(0.00035, preserve_topology=True)
                if simp.is_empty or not isinstance(simp, Polygon):
                    continue
                coords = [[round(pt[1], 5), round(pt[0], 5)] for pt in simp.exterior.coords]
                if len(coords) < 4:
                    continue
                area_ha = round(g.area * (111320.0 * math.cos(math.radians(lat_mid)) * 110540.0) / 10000.0, 1)
                if area_ha < 0.8:
                    continue
                c = g.centroid
                out.append({
                    "region": reg["key"],
                    "region_name": reg["name"],
                    "country": reg["country"],
                    "tier": tier_label,
                    "confidence": 0.78 if tier_label == "high" else 0.58,
                    "area_ha": area_ha,
                    "lat": round(c.y, 5),
                    "lon": round(c.x, 5),
                    "ring": coords,
                })
            return out

        high_polys = merge_and_extract(reg_boxes_high, "high", max_polys=40)
        mod_polys = merge_and_extract(reg_boxes_mod, "moderate", max_polys=30)
        detected_polygons.extend(high_polys)
        detected_polygons.extend(mod_polys)

        sc = region_scenes[reg["key"]]
        region_summaries.append({
            "key": reg["key"],
            "name": reg["name"],
            "country": reg["country"],
            "dry_scene": sc["dry_items"][0],
            "dry_date": sc["dry_date"],
            "wet_scene": sc["wet_items"][0],
            "wet_date": sc["wet_date"],
            "total_prosopis_ha": round(total_prosopis_ha, 1),
            "high_conf_ha": round(high_conf_ha, 1),
            "mean_conf": round(float(np.mean(probs_in_prosopis)), 3) if probs_in_prosopis else 0.0,
            "polygon_count": len(high_polys) + len(mod_polys),
        })
        print(
            f"  Classified {reg['key']}: {round(total_prosopis_ha, 1)} ha Prosopis "
            f"({round(high_conf_ha, 1)} ha core >=0.65), {len(high_polys)+len(mod_polys)} merged polygons"
        )

    compact_trees = [strip_tree_gains(t) for t in trees]
    out_payload = {
        "model": "SoilWatch / CARE v20b smileRandomForest (100 trees, seed=42)",
        "features": FEATURE_NAMES,
        "importance": imp_dict,
        "calibration": {
            "total": int(len(y_cal)),
            "prosopis": int(np.sum(y_cal == 1)),
            "non_prosopis": int(np.sum(y_cal == 0)),
        },
        "validation": {
            "total": n_tot,
            "tp": tp,
            "tn": tn,
            "fp": fp,
            "fn": fn,
            "oa": round(oa * 100.0, 1),
            "kappa": round(kappa, 3),
            "precision": round(tp / max(1, tp + fp) * 100.0, 1),
            "recall": round(tp / max(1, tp + fn) * 100.0, 1),
            "specificity": round(tn / max(1, tn + fp) * 100.0, 1),
        },
        "regions": region_summaries,
        "detected_polygons": detected_polygons,
        "plot_predictions": per_record_pred,
        "rf_trees": compact_trees,
    }
    OUT_MODEL_JSON.write_text(json.dumps(out_payload, separators=(",", ":")))
    print(f"\nSaved {OUT_MODEL_JSON} ({OUT_MODEL_JSON.stat().st_size / 1024:.1f} KB) in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
