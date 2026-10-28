"""Reusable Sentinel-2 helpers for Google Earth Engine (Python API)."""

import ee

S2_SR = "COPERNICUS/S2_SR_HARMONIZED"
CLOUD_SCORE = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"


def s2_clean(aoi, start, end, cs_threshold=0.6):
    """Sentinel-2 surface reflectance, cloud-masked with Cloud Score+.

    cs_cdf goes from 0 (cloudy) to 1 (clear). Pixels below the threshold are masked.
    """
    cs = ee.ImageCollection(CLOUD_SCORE)
    return (
        ee.ImageCollection(S2_SR)
        .filterBounds(aoi)
        .filterDate(start, end)
        .linkCollection(cs, ["cs_cdf"])
        .map(lambda img: img.updateMask(img.select("cs_cdf").gte(cs_threshold)))
    )


def add_indices(img):
    """Add NDVI (greenness) and NDMI (canopy moisture) bands."""
    ndvi = img.normalizedDifference(["B8", "B4"]).rename("NDVI")
    ndmi = img.normalizedDifference(["B8", "B11"]).rename("NDMI")
    return img.addBands([ndvi, ndmi])


def summer_composite(aoi, year, months=("06-01", "08-31")):
    """Median summer composite with NDVI and NDMI for one year."""
    col = s2_clean(aoi, f"{year}-{months[0]}", f"{year}-{months[1]}")
    return add_indices(col.median()).clip(aoi).set("year", year)
