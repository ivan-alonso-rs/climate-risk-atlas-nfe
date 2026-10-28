"""Self-contained Leaflet maps (folium) for the website: PNG overlay + readable legend."""

import folium
from branca.element import MacroElement, Template

ESRI_IMAGERY = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
ESRI_GRAY = "https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}"

_LEGEND_CSS = ("position:absolute; bottom:24px; left:12px; z-index:9999; background:rgba(255,255,255,.92);"
               "padding:8px 10px; border-radius:6px; box-shadow:0 1px 4px rgba(0,0,0,.3);"
               "font:12px/1.4 system-ui,sans-serif; color:#222; max-width:240px;")


def _legend(inner_html):
    el = MacroElement()
    el._template = Template("{% macro html(this, kwargs) %}<div style='" + _LEGEND_CSS + "'>"
                            + inner_html + "</div>{% endmacro %}")
    return el


def gradient_legend(title, palette, vmin, vmax, low_label, high_label):
    bar = f"background:linear-gradient(to right,{','.join(palette)}); height:10px; border-radius:2px;"
    return (f"<b>{title}</b><div style='{bar} margin:4px 0 2px'></div>"
            f"<div style='display:flex; justify-content:space-between'><span>{vmin:+g}</span>"
            f"<span>0</span><span>{vmax:+g}</span></div>"
            f"<div style='display:flex; justify-content:space-between; color:#555'>"
            f"<span>{low_label}</span><span>{high_label}</span></div>")


def class_legend(title, colors, labels):
    rows = "".join(f"<div><span style='display:inline-block;width:12px;height:12px;background:{c};"
                   f"margin-right:6px;vertical-align:-2px'></span>{l}</div>" for c, l in zip(colors, labels))
    return f"<b>{title}</b>{rows}"


def overlay_map(png, bounds, layer_name, legend_html, out_html, opacity=0.8):
    """bounds = [[south, west], [north, east]] in WGS84; PNG must be in EPSG:3857."""
    (s, w), (n, e) = bounds
    m = folium.Map(location=[(s + n) / 2, (w + e) / 2], zoom_start=11,
                   tiles=None, control_scale=True)
    folium.TileLayer(ESRI_IMAGERY, attr="Esri, Maxar, Earthstar Geographics", name="Satellite").add_to(m)
    folium.TileLayer(ESRI_GRAY, attr="Esri, HERE, Garmin, © OpenStreetMap contributors",
                     name="Light map", max_zoom=16, show=False).add_to(m)
    folium.raster_layers.ImageOverlay(image=str(png), bounds=bounds, opacity=opacity,
                                      name=layer_name).add_to(m)
    folium.LayerControl(collapsed=True).add_to(m)
    m.get_root().add_child(_legend(legend_html))
    m.fit_bounds(bounds)
    m.save(out_html)
    return m
