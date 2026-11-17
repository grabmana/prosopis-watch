#!/usr/bin/env python3
"""
Extract the complete CARE / SoilWatch Ground Truth & Polygon datasets from SoilWatch/Prosopis_mapping:
1. 57 Field Polygons from binary shapefile Ground_truth_data/prosopis_merged_2026_08_25_modif.shp + .dbf
   (4,789 vertices; 26 in Ethiopia [13 Prosopis, 13 non-Prosopis] + 31 in Kenya [25 Prosopis, 6 non-Prosopis])
2. 37 Magadi Field GPS points (v19 / CARE presentation slide baseline: 9 Prosopis, 28 non-Prosopis)
   -> Together with 57 ODK polygons = n = 94 (26 in Ethiopia: 13/13; 68 in Kenya: 34/34)
3. 24 Marsabit Field GPS points (v20b addition on 10 Sep 2026: 6 Prosopis, 18 non-Prosopis)
   -> From origin/main:Ground_truth_data/trips_gps_points.csv
4. 230 Online Prosopis (BiendataOrg / GBIF) points in Ethiopia (35) & Kenya (195)
   -> From moder_version_6_13-05-26/prosopis_juliflora_combined_from2016_africa&occAsia.shp + .dbf
5. 78 Online Acacia tortilis & Balanites aegyptiaca (BiendataOrg) points in Ethiopia (12) & Kenya (66)
   -> From model_version_1_17-04-26/acacia-tortilis_balanites-aegyptiaca.kml + .dbf
"""
import csv
import io
import json
import math
import os
import re
import struct
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SESSION_OUT = Path(os.environ.get("PROSOPIS_OUT_DIR", str(ROOT / "build")))
SESSION_OUT.mkdir(parents=True, exist_ok=True)


def read_dbf(dbf_path: Path):
    with open(dbf_path, "rb") as f:
        numrec, lenheader = struct.unpack("<xxxxIH22x", f.read(32))
        numfields = (lenheader - 33) // 32
        fields = []
        for _ in range(numfields):
            raw = f.read(32)
            name = raw[:11].split(b"\x00")[0].decode("latin1")
            ftype = chr(raw[11])
            flen = raw[16]
            fields.append((name, ftype, flen))
        f.seek(lenheader)
        rows = []
        for _ in range(numrec):
            f.read(1)  # deletion flag
            row = {}
            for name, ftype, flen in fields:
                val = f.read(flen).decode("utf-8", errors="replace").strip()
                row[name] = val
            rows.append(row)
    return rows


def read_shp_geometries(shp_path: Path):
    with open(shp_path, "rb") as f:
        data = f.read()
    pos = 100
    geoms = []
    while pos < len(data):
        rec_num, content_len = struct.unpack(">ii", data[pos : pos + 8])
        pos += 8
        rec_data = data[pos : pos + content_len * 2]
        pos += content_len * 2
        st = struct.unpack("<i", rec_data[:4])[0]
        if st == 0:
            geoms.append({"type": "Null"})
        elif st in (1, 11, 21):
            x, y = struct.unpack("<dd", rec_data[4:20])
            geoms.append({"type": "Point", "lon": x, "lat": y})
        elif st in (5, 15, 25):
            num_parts, num_pts = struct.unpack("<ii", rec_data[36:44])
            parts = list(struct.unpack(f"<{num_parts}i", rec_data[44 : 44 + 4 * num_parts])) + [num_pts]
            pts_offset = 44 + 4 * num_parts
            rings = []
            for p_idx in range(num_parts):
                ring = []
                for i in range(parts[p_idx], parts[p_idx + 1]):
                    lon, lat = struct.unpack("<dd", rec_data[pts_offset + i * 16 : pts_offset + (i + 1) * 16])
                    ring.append([round(lat, 6), round(lon, 6)])  # Leaflet [lat, lon]
                rings.append(ring)
            geoms.append({"type": "Polygon", "rings": rings, "vertex_count": num_pts})
        else:
            geoms.append({"type": f"Other({st})"})
    return geoms


def polygon_centroid_latlon(ring_latlon):
    """Compute true planar centroid [lat, lon] of a closed [lat, lon] ring."""
    if not ring_latlon:
        return None, None
    pts = ring_latlon[:-1] if len(ring_latlon) > 2 and ring_latlon[0] == ring_latlon[-1] else ring_latlon
    if len(pts) < 3:
        return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
    area2 = 0.0
    cx = 0.0
    cy = 0.0
    n = len(pts)
    for i in range(n):
        y1, x1 = pts[i][0], pts[i][1]
        y2, x2 = pts[(i + 1) % n][0], pts[(i + 1) % n][1]
        cross = x1 * y2 - x2 * y1
        area2 += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
    if abs(area2) < 1e-12:
        return sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
    lon = cx / (3.0 * area2)
    lat = cy / (3.0 * area2)
    return round(lat, 6), round(lon, 6)


def polygon_area_ha(ring_latlon):
    """Approximate polygon area in hectares on WGS84 sphere."""
    if not ring_latlon or len(ring_latlon) < 3:
        return 0.0
    lat0 = math.radians(sum(p[0] for p in ring_latlon) / len(ring_latlon))
    m_per_deg_lat = 111132.92
    m_per_deg_lon = 111412.84 * math.cos(lat0)
    pts = [(p[1] * m_per_deg_lon, p[0] * m_per_deg_lat) for p in ring_latlon]
    area2 = 0.0
    for i in range(len(pts)):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % len(pts)]
        area2 += x1 * y2 - x2 * y1
    return round(abs(area2) * 0.5 / 10000.0, 3)


def classify_region(lat: float, lon: float) -> tuple[str, str]:
    country = "Ethiopia" if (lat >= 4.6 or (lat >= 3.5 and lon >= 39.0)) else "Kenya"
    if country == "Ethiopia":
        if lat >= 11.0:
            return "Afar (North), Ethiopia", "Ethiopia"
        if lat >= 8.5:
            return "Afar (Gewane / Awash), Ethiopia", "Ethiopia"
        return "Ethiopia", "Ethiopia"
    else:
        if lat < -0.5:
            return "Lake Magadi / Kajiado, Kenya", "Kenya"
        if lon < 36.2 and lat > 2.0:
            return "Turkana, Kenya", "Kenya"
        if 0.0 <= lat <= 2.2 and 35.5 <= lon <= 36.6:
            return "Baringo / Rift Valley, Kenya", "Kenya"
        if lon >= 36.8 and lat > 0.0:
            return "Marsabit / Laisamis, Kenya", "Kenya"
        return "Kenya", "Kenya"


def safe_float(val, default=None):
    try:
        return float(val)
    except Exception:
        return default


def load_odk_polygons():
    shp_path = ROOT / "Ground_truth_data" / "prosopis_merged_2026_08_25_modif.shp"
    dbf_path = ROOT / "Ground_truth_data" / "prosopis_merged_2026_08_25_modif.dbf"
    geoms = read_shp_geometries(shp_path)
    rows = read_dbf(dbf_path)
    records = []
    for r, g in zip(rows, geoms):
        ring = g["rings"][0] if g.get("rings") else []
        poly_lat, poly_lon = polygon_centroid_latlon(ring)
        gps_lat = safe_float(r.get("X_centroid"), poly_lat)
        gps_lon = safe_float(r.get("X_centro_1"), poly_lon)
        # Use true polygon centroid so the dot sits inside the modified shapefile polygon!
        lat = poly_lat if poly_lat is not None else gps_lat
        lon = poly_lon if poly_lon is not None else gps_lon
        cls_str = (r.get("class") or "").strip().lower()
        label = 0 if "non" in cls_str else 1
        region, country = classify_region(lat, lon)
        area_ha = polygon_area_ha(ring)
        species = "Prosopis juliflora" if label == 1 else (r.get("nonprosopi") or "Non-Prosopis")
        records.append({
            "id": r.get("gt_uid") or f"odk_{len(records)+1}",
            "dataset": "odk",
            "careGroup": "field_94",
            "source": "CARE / SoilWatch ODK Field Polygon",
            "country": country,
            "region": region,
            "lat": lat,
            "lon": lon,
            "gpsLat": round(gps_lat, 6) if gps_lat else lat,
            "gpsLon": round(gps_lon, 6) if gps_lon else lon,
            "label": label,
            "species": species.replace("_", " "),
            "density": r.get("density") or "n/a",
            "canopy": r.get("canopy_cov") or "n/a",
            "landUse": (r.get("land_use") or "n/a").replace("_", " "),
            "observer": (r.get("observer") or "Field team").strip(),
            "date": (r.get("date") or "")[:10],
            "notes": f"Walked field polygon ({len(ring)} vertices, {area_ha} ha) · Height: {r.get('height_cla') or 'n/a'} · Confidence: {r.get('confidence') or 'n/a'}",
            "areaHa": area_ha,
            "vertexCount": len(ring),
            "ring": ring,
        })
    return records


def load_trip_points():
    raw_csv = subprocess.check_output(
        ["git", "-C", str(ROOT), "show", "origin/main:Ground_truth_data/trips_gps_points.csv"]
    ).decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(raw_csv), delimiter=";")
    records = []
    for idx, r in enumerate(reader, start=1):
        dec = (r.get("decision_c") or "").strip()
        if dec not in ("0", "1"):
            continue
        label = int(dec)
        lat = round(float(r["point_lat"]), 6)
        lon = round(float(r["point_lon"]), 6)
        trip_name = (r.get("Trip") or "").strip()
        pt_id = (r.get("point_id") or str(idx)).strip()
        src = (r.get("source") or "field").strip()
        obs_text = (r.get("field_obse") or "").strip()
        date_raw = (r.get("Date") or "").strip()
        is_magadi = trip_name.lower() == "magadi" or lat < -0.5
        region, country = classify_region(lat, lon)
        care_group = "field_94" if is_magadi else "field_marsabit_24"
        records.append({
            "id": f"trip_{trip_name.lower()}_{pt_id}",
            "dataset": "trips_magadi" if is_magadi else "trips_marsabit",
            "careGroup": care_group,
            "source": f"CARE Field GPS · {trip_name} ({src})",
            "country": country,
            "region": region,
            "lat": lat,
            "lon": lon,
            "label": label,
            "species": "Prosopis juliflora" if label == 1 else "Non-Prosopis (field verified)",
            "density": "n/a",
            "canopy": "n/a",
            "landUse": f"{trip_name} field transect",
            "observer": "Eero / SoilWatch",
            "date": "2026-08-23" if is_magadi else "2026-09-09",
            "notes": f"{trip_name} #{pt_id} ({date_raw}): {obs_text} · 20m GEE buffer polygon",
            "bufferMeters": 20,
            "ring": [],
        })
    return records


def load_online_prosopis_230():
    shp_path = ROOT / "moder_version_6_13-05-26" / "prosopis_juliflora_combined_from2016_africa&occAsia.shp"
    dbf_path = ROOT / "moder_version_6_13-05-26" / "prosopis_juliflora_combined_from2016_africa&occAsia.dbf"
    geoms = read_shp_geometries(shp_path)
    rows = read_dbf(dbf_path)
    records = []
    for idx, (r, g) in enumerate(zip(rows, geoms), start=1):
        cc = (r.get("countryCod") or "").strip()
        if cc not in ("ET", "KE"):
            continue
        lat = round(float(g["lat"]), 6)
        lon = round(float(g["lon"]), 6)
        country = "Ethiopia" if cc == "ET" else "Kenya"
        region, _ = classify_region(lat, lon)
        src_org = (r.get("source") or "Biendata.org / GBIF").strip()
        date_str = (r.get("date_colle") or r.get("eventDate") or "")[:10]
        records.append({
            "id": f"bien_pros_{len(records)+1:03d}",
            "dataset": "online_prosopis",
            "careGroup": "online_prosopis_230",
            "source": f"Prosopis online dataset ({src_org})",
            "country": country,
            "region": region,
            "lat": lat,
            "lon": lon,
            "label": 1,
            "species": "Prosopis juliflora",
            "density": "online occurrence",
            "canopy": "n/a",
            "landUse": "online record (20m GEE buffer)",
            "observer": src_org,
            "date": date_str or "≥2016",
            "notes": f"Online Prosopis occurrence ({src_org}) · Buffered 20m in GEE v20b",
            "bufferMeters": 20,
            "ring": [],
        })
    return records


def load_online_acacia_balanites_78():
    kml_path = ROOT / "model_version_1_17-04-26" / "acacia-tortilis_balanites-aegyptiaca.kml"
    dbf_path = ROOT / "model_version_1_17-04-26" / "acacia-tortilis_balanites-aegyptiaca.dbf"
    dbf_rows = read_dbf(dbf_path)
    txt = kml_path.read_text(encoding="utf-8", errors="replace")
    placemarks = re.findall(r"<Placemark>(.*?)</Placemark>", txt, re.S)
    records = []
    for idx, pm in enumerate(placemarks):
        m_coord = re.search(r"<coordinates>\s*([0-9.-]+),([0-9.-]+)", pm)
        if not m_coord:
            continue
        lon = round(float(m_coord.group(1)), 6)
        lat = round(float(m_coord.group(2)), 6)
        r = dbf_rows[idx] if idx < len(dbf_rows) else {}
        sp = (r.get("scrubbed_s") or r.get("species") or "Acacia tortilis / Balanites aegyptiaca").strip()
        src_org = (r.get("source") or "Biendata.org").strip()
        date_str = (r.get("date_colle") or r.get("eventDate") or "")[:10]
        region, country_guess = classify_region(lat, lon)
        country_dbf = (r.get("country") or "").strip()
        country = country_dbf if country_dbf in ("Ethiopia", "Kenya") else country_guess
        records.append({
            "id": f"bien_abs_{idx+1:02d}",
            "dataset": "online_non_prosopis",
            "careGroup": "online_acacia_78",
            "source": f"Acacia & Balanites online dataset ({src_org})",
            "country": country,
            "region": region,
            "lat": lat,
            "lon": lon,
            "label": 0,
            "species": sp,
            "density": "online absence",
            "canopy": "n/a",
            "landUse": "online absence (20m GEE buffer)",
            "observer": src_org,
            "date": date_str or "online",
            "notes": f"Online Non-Prosopis ({sp}) from {src_org} · Buffered 20m in GEE v20b",
            "bufferMeters": 20,
            "ring": [],
        })
    return records


def main():
    odk = load_odk_polygons()
    trips = load_trip_points()
    online_p = load_online_prosopis_230()
    online_a = load_online_acacia_balanites_78()

    field_94 = [r for r in odk + trips if r["careGroup"] == "field_94"]
    field_et = [r for r in field_94 if r["country"] == "Ethiopia"]
    field_ke = [r for r in field_94 if r["country"] == "Kenya"]
    marsabit_24 = [r for r in trips if r["careGroup"] == "field_marsabit_24"]

    summary = {
        "care_slide_field_n94": {
            "total": len(field_94),
            "ethiopia_total": len(field_et),
            "ethiopia_prosopis": sum(1 for r in field_et if r["label"] == 1),
            "ethiopia_non_prosopis": sum(1 for r in field_et if r["label"] == 0),
            "kenya_total": len(field_ke),
            "kenya_prosopis": sum(1 for r in field_ke if r["label"] == 1),
            "kenya_non_prosopis": sum(1 for r in field_ke if r["label"] == 0),
        },
        "marsabit_sep2026_addition": {
            "total": len(marsabit_24),
            "prosopis": sum(1 for r in marsabit_24 if r["label"] == 1),
            "non_prosopis": sum(1 for r in marsabit_24 if r["label"] == 0),
        },
        "online_prosopis_biendata_n230": {
            "total": len(online_p),
            "ethiopia": sum(1 for r in online_p if r["country"] == "Ethiopia"),
            "kenya": sum(1 for r in online_p if r["country"] == "Kenya"),
        },
        "online_acacia_balanites_biendata_n78": {
            "total": len(online_a),
            "ethiopia": sum(1 for r in online_a if r["country"] == "Ethiopia"),
            "kenya": sum(1 for r in online_a if r["country"] == "Kenya"),
        },
        "odk_polygons": {
            "count": len(odk),
            "total_vertices": sum(r["vertexCount"] for r in odk),
            "total_area_ha": round(sum(r["areaHa"] for r in odk), 2),
        },
    }
    print(json.dumps(summary, indent=2))

    all_records = odk + trips + online_p + online_a
    out_json = SESSION_OUT / "care_all_datasets.json"
    out_json.write_text(json.dumps(all_records, separators=(",", ":")), encoding="utf-8")
    print(f"Wrote {len(all_records)} total records to {out_json} ({out_json.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
