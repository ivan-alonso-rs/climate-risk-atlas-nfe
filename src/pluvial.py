"""Topographic screening for pluvial (heavy-rain) flooding.

Depression analysis on a DEM: fill sinks, keep depressions deeper than
MIN_DEPTH, compute their depth/area/volume, and a D8 flow accumulation.
A screening of where surface water concentrates; it is not a hydrodynamic model.
"""
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from scipy import ndimage
from pysheds.grid import Grid
from pysheds.sview import Raster

if not hasattr(np, "in1d"):  # pysheds 0.5 still calls np.in1d, removed in NumPy 2
    np.in1d = np.isin

# Minimum depression depth (m). Berlin's Senkenanalyse also uses 0.20 m;
# shallower sinks are mostly DEM artefacts.
MIN_DEPTH = 0.20


def read_dem(path, bounds=None):
    """Read DEM as float32 with NaN as nodata; optionally clip to (xmin, ymin, xmax, ymax) in the DEM's CRS."""
    with rasterio.open(path) as src:
        win = from_bounds(*bounds, transform=src.transform) if bounds else None
        dem = src.read(1, window=win).astype("float32")
        tfm = src.window_transform(win) if win else src.transform
        nodata, crs = src.nodata, src.crs
    if nodata is not None:
        dem[dem == nodata] = np.nan
    return dem, tfm, crs


def depressions(dem, res, min_depth=MIN_DEPTH):
    """Fill sinks and return (depth, labels, table).

    depth: filled DEM minus DEM (m), only inside kept depressions, else 0.
    labels: integer id per depression (0 = none).
    table: dict of arrays (id, area_m2, max_depth_m, volume_m3).
    """
    from skimage.morphology import reconstruction
    # Priority-flood style fill by morphological reconstruction (erosion from the border)
    dem0 = np.where(np.isnan(dem), np.nanmax(dem), dem)
    seed = dem0.copy()
    seed[1:-1, 1:-1] = dem0.max()
    filled = reconstruction(seed, dem0, method="erosion")
    depth = (filled - dem0).astype("float32")
    depth[np.isnan(dem)] = 0

    labels, n = ndimage.label(depth > 0)
    idx = np.arange(1, n + 1)
    max_d = ndimage.maximum(depth, labels, idx)
    keep = idx[max_d >= min_depth]
    mask = np.isin(labels, keep)
    labels = np.where(mask, labels, 0)
    depth = np.where(mask, depth, 0).astype("float32")

    area = ndimage.sum(np.ones_like(depth), labels, keep) * res * res
    vol = ndimage.sum(depth, labels, keep) * res * res
    mx = ndimage.maximum(depth, labels, keep)
    return depth, labels, {"id": keep, "area_m2": area, "max_depth_m": mx, "volume_m3": vol}


def flow_accumulation(dem_path):
    """D8 flow accumulation (upstream cell count) from a DEM GeoTIFF.

    Sinks are filled and flats resolved first, so water routes through
    depressions to the outlet. Returns a float32 array aligned with the DEM.
    """
    grid = Grid.from_raster(dem_path)
    dem = grid.read_raster(dem_path)
    flooded = grid.fill_pits(dem)
    inflated = grid.resolve_flats(grid.fill_depressions(flooded))
    fdir = grid.flowdir(inflated)
    acc = grid.accumulation(fdir)
    return np.asarray(acc, dtype="float32")


def osm_buildings_mask(bounds, shape, transform, crs="EPSG:25832"):
    """Rasterise OpenStreetMap building footprints onto the DEM grid (True = building).

    bounds: (xmin, ymin, xmax, ymax) in `crs`. Multipolygon relations are skipped
    (simple closed ways only), so coverage is slightly incomplete.
    """
    import requests
    from pyproj import Transformer
    from rasterio.features import rasterize
    from shapely.geometry import Polygon
    to_ll = Transformer.from_crs(crs, 4326, always_xy=True)
    to_utm = Transformer.from_crs(4326, crs, always_xy=True)
    w, s = to_ll.transform(bounds[0], bounds[1])
    e, n = to_ll.transform(bounds[2], bounds[3])
    q = f'[out:json][timeout:180];way["building"]({s},{w},{n},{e});out geom;'
    r = requests.post("https://overpass-api.de/api/interpreter", data={"data": q}, timeout=240,
                      headers={"User-Agent": "climate-risk-atlas/0.1 (EO portfolio project)"})
    r.raise_for_status()
    polys = []
    for el in r.json()["elements"]:
        pts = [to_utm.transform(p["lon"], p["lat"]) for p in el.get("geometry", [])]
        if len(pts) >= 4:
            polys.append(Polygon(pts))
    mask = rasterize([(p, 1) for p in polys if p.is_valid], out_shape=shape,
                     transform=transform, fill=0, dtype="uint8").astype(bool)
    return mask, len(polys)


def rank_depressions(dem, res, acc, labels, table):
    """Add catchment area and fill ratio to the depression table.

    catchment_m2: largest upstream area (flow accumulation) inside the depression.
    fill_ratio: catchment_m2 / volume_m3. High = a lot of runoff for little storage.
    """
    ids = table["id"]
    acc_max = ndimage.maximum(acc, labels, ids)
    table = dict(table)
    table["catchment_m2"] = acc_max * res * res
    table["fill_ratio"] = table["catchment_m2"] / np.maximum(table["volume_m3"], 1.0)
    return table


def osm_water_mask(bounds, shape, transform, crs="EPSG:25832"):
    """Rasterise OSM permanent water (natural=water, riverbanks) onto the DEM grid.

    Relations are assembled from their outer member ways with shapely.polygonize.
    Returns a boolean mask (True = water). Used to drop rivers/lakes from the depression ranking.
    """
    import requests
    from pyproj import Transformer
    from rasterio.features import rasterize
    from shapely.geometry import LineString, Polygon
    from shapely.ops import polygonize, unary_union
    to_ll = Transformer.from_crs(crs, 4326, always_xy=True)
    to_utm = Transformer.from_crs(4326, crs, always_xy=True)
    w, s = to_ll.transform(bounds[0], bounds[1])
    e, n = to_ll.transform(bounds[2], bounds[3])
    bb = f"({s},{w},{n},{e})"
    q = (f'[out:json][timeout:180];(way["natural"="water"]{bb};relation["natural"="water"]{bb};'
         f'way["waterway"="riverbank"]{bb};relation["waterway"="riverbank"]{bb};);out geom;')
    r = requests.post("https://overpass-api.de/api/interpreter", data={"data": q}, timeout=240,
                      headers={"User-Agent": "climate-risk-atlas/0.1 (EO portfolio project)"})
    r.raise_for_status()

    def xy(geom):
        return [to_utm.transform(p["lon"], p["lat"]) for p in geom]

    polys = []
    for el in r.json()["elements"]:
        if el["type"] == "way":
            pts = xy(el.get("geometry", []))
            if len(pts) >= 4 and pts[0] == pts[-1]:
                polys.append(Polygon(pts))
        else:
            lines = [LineString(xy(m["geometry"])) for m in el.get("members", [])
                     if m.get("role") == "outer" and len(m.get("geometry", [])) >= 2]
            if lines:
                polys.extend(polygonize(unary_union(lines)))
    shapes = [(p, 1) for p in polys if p.is_valid and not p.is_empty]
    if not shapes:
        return np.zeros(shape, bool), 0
    return rasterize(shapes, out_shape=shape, transform=transform, fill=0, dtype="uint8").astype(bool), len(shapes)
