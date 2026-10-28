# Climate Risk Atlas – Nürnberg · Fürth · Erlangen

Open-data case studies on climate risks in the Nürnberg–Fürth–Erlangen metropolitan core, built with satellite data (Sentinel-1/2, Landsat) and Google Earth Engine from Python.

**Author:** Iván Eduardo Alonso Ramírez · [LinkedIn](https://www.linkedin.com/in/ivan-alonso-rs/)

## Case studies

| # | Topic | Main data | Status |
|---|---|---|---|
| 1 | Forest health in the Nürnberger Reichswald (2017–2025) + wildfire risk | Sentinel-2 (NDVI, NDMI), DWD | planned |
| 2 | Urban heat islands and vulnerable population | Landsat 8/9 thermal, Zensus 2022 | planned |
| 3 | Watershed risk: heavy rain or drought (Regnitz/Pegnitz) | DEM (Bavarian open data), Sentinel-1 | planned |

Each case study answers: problem → data → method → result → recommendation for the municipality → limitations.

## Project structure

```
climate-risk-atlas-nfe/
├── config.py            # AOI, dates, GEE project id
├── environment.yml      # conda environment
├── notebooks/           # one notebook per step
├── src/                 # reusable functions (cloud masks, indices)
├── data/raw, processed  # local data (not pushed to GitHub)
├── figures/             # exported maps and charts
└── site/                # Quarto website (published on GitHub Pages)
```

## Setup

```bash
conda env create -f environment.yml
conda activate atlas
earthengine authenticate
jupyter lab
```

Then open `notebooks/00_setup_check.ipynb`.

## Disclaimer

Demonstration project based on open data. Not an official assessment by any public authority.
