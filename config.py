"""Shared settings for the Climate Risk Atlas.

Change GEE_PROJECT to your own Google Cloud project id
(the one you use in the Earth Engine Code Editor).
"""

GEE_PROJECT = "striking-loop-456219-g7"  # <- your Earth Engine project id

# Study area: metropolitan core Nürnberg – Fürth – Erlangen incl. Nürnberger Reichswald
# [west, south, east, north] in degrees (WGS84). First version = simple box;
# later we replace it with official municipal boundaries (Bavarian open data).
AOI_BBOX = [10.85, 49.33, 11.40, 49.66]

# Projected CRS for area calculations and exports (UTM zone 32N)
CRS = "EPSG:32632"

# Summer window used for vegetation composites
SUMMER = ("06-01", "08-31")
