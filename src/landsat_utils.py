import ee


def _mask_scale(img):
      qa = img.select("QA_PIXEL")
      # bits 1 dilated cloud, 3 cloud, 4 cloud shadow, 5 snow
      bad = qa.bitwiseAnd(0b111111).neq(0)
      lst = img.select("ST_B10").multiply(0.00341802).add(149.0).subtract(273.15).rename("LST")
      return lst.updateMask(bad.Not()).copyProperties(img, ["system:time_start"])


def lst_collection(aoi, start, end):
      cols = [
          ee.ImageCollection(f"LANDSAT/{s}/C02/T1_L2").filterBounds(aoi).filterDate(start, end)
          for s in ("LC08", "LC09")
      ]
      return cols[0].merge(cols[1]).map(_mask_scale)


def summer_lst(aoi, year, months=("06-01", "08-31")):
      """Median summer LST (°C) for one year."""
      col = lst_collection(aoi, f"{year}-{months[0]}", f"{year}-{months[1]}")
      return col.median().clip(aoi).set("year", year)