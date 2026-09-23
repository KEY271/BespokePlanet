"""Build the Earth topography intermediate file of docs/dynamics/earth-topography.md.

Steps: fetch ETOPO 2022 (ice surface) at a 6 arc-minute stride through OPeNDAP
(cached under output/etopo2022/), aggregate to 0.5 degree cells (land fraction
and mean land height), and write core/data/earth_topography_0p5deg.{bin,json}.
The terrain the solver builds from it is shown by the visualizer (viz/).  Run with
`uv run --project ~/.local/share/llm-python --with netCDF4 python scripts/prepare_earth_topography.py`
(netCDF4 is only needed when the cache is missing).
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/"
              "60s_surface_elev_netcdf/ETOPO_2022_v1_60s_N90W180_surface.nc")
SOURCE_DOI = "10.25921/fd45-gt74"
STRIDE = 6  # 60 arc-second pixels -> 6 arc-minute samples
CACHE = ROOT / "output" / "etopo2022" / "ETOPO_2022_v1_surface_6min_stride.npz"
CELL_DEGREES = 0.5
DATA_BIN = ROOT / "core" / "data" / "earth_topography_0p5deg.bin"
DATA_JSON = ROOT / "core" / "data" / "earth_topography_0p5deg.json"

# (name, longitude east 0..360, latitude) used as landmarks in the docs and tests
LANDMARKS = [("Tibet", 90.0, 33.0), ("Andes (Altiplano)", 292.0, -18.0), ("Antarctica", 90.0, -80.0),
             ("Greenland", 320.0, 72.0), ("central Pacific", 200.0, 0.0)]


def fetch_samples() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (z[lat, lon], lat, lon) at the stride, from cache or OPeNDAP."""
    if CACHE.exists():
        with np.load(CACHE) as cached:
            return cached["z"], cached["lat"], cached["lon"]
    import netCDF4  # noqa: PLC0415  (optional dependency, only for the download)

    dataset = netCDF4.Dataset(SOURCE_URL)
    variable = dataset["z"]
    lat = np.asarray(dataset["lat"][::STRIDE], dtype=np.float64)
    lon = np.asarray(dataset["lon"][::STRIDE], dtype=np.float64)
    z = np.empty((lat.size, lon.size), dtype=np.float32)
    rows = 1800
    for start in range(0, variable.shape[0], rows):
        block = np.asarray(variable[start:start + rows:STRIDE, ::STRIDE])
        z[start // STRIDE:(start + rows) // STRIDE, :] = block
        print(f"fetched rows {start}-{start + rows}", file=sys.stderr, flush=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, z=z, lat=lat, lon=lon, stride=STRIDE, source=SOURCE_URL)
    return z, lat, lon


def to_model_longitudes(z: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lon360 = np.mod(lon, 360.0)
    order = np.argsort(lon360)
    return z[:, order], lon360[order]


def aggregate(z: np.ndarray, lat: np.ndarray, lon: np.ndarray):
    """Land fraction and mean land height on CELL_DEGREES cells (south to north, 0..360)."""
    nlat = int(round(180.0 / CELL_DEGREES))
    nlon = int(round(360.0 / CELL_DEGREES))
    lat_index = np.floor((lat + 90.0) / CELL_DEGREES).astype(int).clip(0, nlat - 1)
    lon_index = np.floor(lon / CELL_DEGREES).astype(int).clip(0, nlon - 1)
    land = (z > 0.0).astype(np.float64)
    height = np.maximum(z, 0.0).astype(np.float64)
    count = np.zeros((nlat, nlon))
    land_sum = np.zeros((nlat, nlon))
    height_sum = np.zeros((nlat, nlon))
    for j in range(lat.size):
        np.add.at(count[lat_index[j]], lon_index, 1.0)
        np.add.at(land_sum[lat_index[j]], lon_index, land[j])
        np.add.at(height_sum[lat_index[j]], lon_index, height[j])
    if count.min() <= 0:
        raise RuntimeError("some 0.5 degree cells received no samples")
    samples_per_cell = int(count.min())
    if count.max() != count.min():
        raise RuntimeError("uneven sample count per cell; check the stride")
    return land_sum / count, height_sum / count, samples_per_cell


def cell_centres():
    nlat = int(round(180.0 / CELL_DEGREES))
    nlon = int(round(360.0 / CELL_DEGREES))
    lat = -90.0 + (np.arange(nlat) + 0.5) * CELL_DEGREES
    lon = (np.arange(nlon) + 0.5) * CELL_DEGREES
    return lon, lat


def area_mean(field: np.ndarray, lat: np.ndarray) -> float:
    weight = np.cos(np.radians(lat))[:, None] * np.ones_like(field)
    return float((field * weight).sum() / weight.sum())


def nearest(field: np.ndarray, lon: np.ndarray, lat: np.ndarray, lon0: float, lat0: float) -> float:
    return float(field[np.abs(lat - lat0).argmin(), np.abs(lon - lon0).argmin()])


def write_intermediate(land_fraction: np.ndarray, land_height: np.ndarray, samples_per_cell: int) -> dict:
    DATA_BIN.parent.mkdir(parents=True, exist_ok=True)
    payload = np.concatenate([land_fraction.ravel(), land_height.ravel()]).astype("<f8").tobytes()
    DATA_BIN.write_bytes(payload)
    lon, lat = cell_centres()
    peak = np.unravel_index(land_height.argmax(), land_height.shape)
    info = {
        "source": "ETOPO 2022 v1, 60 arc-second, ice surface (NOAA NCEI)",
        "source_url": SOURCE_URL,
        "source_doi": SOURCE_DOI,
        "access": f"OPeNDAP, every {STRIDE}th pixel in latitude and longitude (6 arc-minute point samples)",
        "land_definition": "surface elevation z > 0 m",
        "cell_degrees": CELL_DEGREES,
        "samples_per_cell": samples_per_cell,
        "shape_lon_lat": [lon.size, lat.size],
        "layout": "real64 little-endian; land_fraction block then land_height block; "
                  "longitude varies fastest; latitude south to north; "
                  "cell centres lon=(i-1/2)*0.5 deg, lat=-90+(j-1/2)*0.5 deg",
        "land_height_units": "m, mean of max(z,0) over the cell (ocean counted as 0)",
        "generated": dt.date.today().isoformat(),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "global_land_fraction": area_mean(land_fraction, lat),
        "maximum_land_height_m": float(land_height.max()),
        "maximum_land_height_lon_lat": [float(lon[peak[1]]), float(lat[peak[0]])],
        "landmarks_m": {name: nearest(land_height, lon, lat, x, y) for name, x, y in LANDMARKS},
    }
    DATA_JSON.write_text(json.dumps(info, indent=2, ensure_ascii=False) + "\n")
    return info


def main() -> None:
    z, lat, lon = fetch_samples()
    if not np.isfinite(z).all() or (z < -20000).any():
        raise RuntimeError("fill values present in the ETOPO samples")
    z, lon = to_model_longitudes(z, lon)
    land_fraction, land_height, samples_per_cell = aggregate(z, lat, lon)
    info = write_intermediate(land_fraction, land_height, samples_per_cell)
    print(json.dumps({k: info[k] for k in ("global_land_fraction", "maximum_land_height_m",
                                           "maximum_land_height_lon_lat", "landmarks_m", "sha256")},
                     indent=2))
    print(f"wrote {DATA_BIN} ({DATA_BIN.stat().st_size} bytes), {DATA_JSON}")


if __name__ == "__main__":
    main()
