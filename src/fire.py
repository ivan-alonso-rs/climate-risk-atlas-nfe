"""Static wildfire susceptibility index for forest areas.

Framework (Chuvieco et al. 2010): combine human ignition probability with fuel
conditions (fuel type + fuel moisture). All inputs are resampled to one grid.
"""

import math

import ee
import numpy as np
import requests
import rasterio
from rasterio.features import rasterize
from rasterio.transform import Affine
from rasterio.warp import Resampling, reproject
from scipy.ndimage import distance_transform_edt
from shapely.geometry import LineString

S2_SR = "COPERNICUS/S2_SR_HARMONIZED"
WORLDCOVER = "ESA/WorldCover/v200/2021"
# DLR Tree Species Germany 2022 (10 m, EPSG:3035, CC BY 4.0), read remotely as COG
TREE_SPECIES = "/vsicurl/https://download.geoservice.dlr.de/TREE_SPECIES_DE/files/treespecies_de_2022.tif"
TREE_SPECIES_NODATA = 999
OVERPASS = "https://overpass-api.de/api/interpreter"

# Relative flammability per DLR class. Expert assumption: pine > other conifers > broadleaves;
# canopy-cover loss (clear-cuts with slash and dry grass) rated high.
FUEL_SCORE = {0: 1.0,                                  # pine
              1: 0.6, 2: 0.6, 3: 0.6, 4: 0.6,          # spruce, Douglas fir, larch, fir
              5: 0.3, 6: 0.3, 7: 0.3, 8: 0.3, 9: 0.3,  # beech, oak, birch, alder, other
              666: 0.8}                                # canopy cover loss

# OSM highway values treated as forest tracks/paths (the rest count as roads)
TRACKS = {"track", "path", "footway", "cycleway", "bridleway", "steps"}


# ---------- grid ----------

def make_grid(geom, crs, res=20, buffer_m=1500):
    """Snapped raster grid (transform, width, height) covering geom + buffer."""
    b = geom.transform(crs, 1).buffer(buffer_m, 1).bounds(1, crs).coordinates().get(0).getInfo()
    xs, ys = [p[0] for p in b], [p[1] for p in b]
    xmin, ymax = math.floor(min(xs) / res) * res, math.ceil(max(ys) / res) * res
    width = math.ceil((max(xs) - xmin) / res)
    height = math.ceil((ymax - min(ys)) / res)
    return Affine(res, 0, xmin, 0, -res, ymax), width, height


def ee_to_array(img, crs, transform, width, height):
    """Download a single-band ee.Image onto the grid as a float32 numpy array."""
    url = img.toFloat().getDownloadURL({
        "crs": crs, "crs_transform": list(transform)[:6],
        "dimensions": f"{width}x{height}", "format": "GEO_TIFF"})
    with rasterio.MemoryFile(requests.get(url, timeout=300).content) as mem, mem.open() as src:
        return src.read(1, masked=True).filled(np.nan).astype("float32")


# ---------- fuel moisture (Sentinel-2 NDMI anomaly) ----------

def _mask_scl(img):
    scl = img.select("SCL")
    return img.select(["B8", "B11"]).updateMask(scl.eq(4).Or(scl.eq(5)))


def summer_ndmi(geom, year):
    col = (ee.ImageCollection(S2_SR).filterBounds(geom)
           .filterDate(f"{year}-06-01", f"{year}-08-31")
           .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 40))
           .map(_mask_scl))
    return col.median().normalizedDifference(["B8", "B11"]).rename("NDMI")


def ndmi_zscore(geom, year, years):
    """Per-pixel z-score of summer NDMI in `year` relative to `years`.

    Each pixel is compared with its own history, so species-specific NDMI levels cancel out.
    """
    stack = ee.ImageCollection([summer_ndmi(geom, y) for y in years])
    mean, sd = stack.mean(), stack.reduce(ee.Reducer.stdDev())
    return summer_ndmi(geom, year).subtract(mean).divide(sd).rename("NDMI_z")


# ---------- fuel type (DLR tree species) ----------

def tree_species_on_grid(crs, transform, width, height):
    """Dominant tree species resampled (majority) to the grid; NODATA outside forest."""
    dst = np.full((height, width), TREE_SPECIES_NODATA, dtype="uint16")
    with rasterio.open(TREE_SPECIES) as src:
        reproject(rasterio.band(src, 1), dst, dst_transform=transform, dst_crs=crs,
                  dst_nodata=TREE_SPECIES_NODATA, resampling=Resampling.mode)
    return dst


def fuel_score(species):
    out = np.full(species.shape, np.nan, dtype="float32")
    for cls, score in FUEL_SCORE.items():
        out[species == cls] = score
    return out


# ---------- human ignition (OSM roads/tracks, built-up) ----------

def osm_highways(bbox_wgs84):
    """OSM highway lines in [west, south, east, north] as (geometry, highway) pairs."""
    w, s, e, n = bbox_wgs84
    q = f'[out:json][timeout:180];way["highway"]({s},{w},{n},{e});out geom;'
    r = requests.post(OVERPASS, data={"data": q}, timeout=300,
                      headers={"User-Agent": "climate-risk-atlas-nfe (research/portfolio)"})
    r.raise_for_status()
    return [(LineString([(p["lon"], p["lat"]) for p in el["geometry"]]), el["tags"]["highway"])
            for el in r.json()["elements"] if len(el.get("geometry", [])) > 1]


def distance_to(mask, res):
    """Euclidean distance (m) from every cell to the nearest True cell."""
    return distance_transform_edt(~mask) * res


def lines_to_mask(lines, crs, transform, width, height):
    import geopandas as gpd
    g = gpd.GeoSeries(lines, crs="EPSG:4326").to_crs(crs)
    return rasterize(((geom, 1) for geom in g), out_shape=(height, width),
                     transform=transform, all_touched=True, dtype="uint8").astype(bool)


def decay(d, d0):
    """Exponential distance decay: 1 at the source, 0.37 at d0."""
    return np.exp(-d / d0)


# ---------- index ----------

def minmax(a, valid):
    lo, hi = np.nanpercentile(a[valid], [2, 98])
    return np.clip((a - lo) / (hi - lo), 0, 1)


def susceptibility(factors, weights):
    """Weighted linear combination of factors (dict name -> 0..1 array)."""
    total = sum(weights.values())
    return sum(factors[k] * w for k, w in weights.items()) / total


def rank_correlation(a, b):
    """Spearman correlation computed element-wise (no BLAS needed)."""
    from scipy.stats import rankdata
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean(); rb -= rb.mean()
    return float((ra * rb).sum() / np.sqrt((ra * ra).sum() * (rb * rb).sum()))


def quantile_classes(index, valid, q=(0.5, 0.8)):
    """Relative classes: 1 = low (< median), 2 = moderate, 3 = high (top 20 %)."""
    t = np.nanquantile(index[valid], q)
    cls = np.digitize(index, t) + 1
    return np.where(valid, cls, 0).astype("uint8")
