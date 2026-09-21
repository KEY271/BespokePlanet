"""Build the Earth topography intermediate file of docs/dynamics/earth-topography.md.

Steps: fetch ETOPO 2022 (ice surface) at a 6 arc-minute stride through OPeNDAP
(cached under output/etopo2022/), aggregate to 0.5 degree cells (land fraction
and mean land height), write core/data/earth_topography_0p5deg.{bin,json}, and
draw scripts/earth_topography_source.png for inspection.  Run with
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
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
SOURCE_URL = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/"
              "60s_surface_elev_netcdf/ETOPO_2022_v1_60s_N90W180_surface.nc")
SOURCE_DOI = "10.25921/fd45-gt74"
STRIDE = 6  # 60 arc-second pixels -> 6 arc-minute samples
CACHE = ROOT / "output" / "etopo2022" / "ETOPO_2022_v1_surface_6min_stride.npz"
CELL_DEGREES = 0.5
DATA_BIN = ROOT / "core" / "data" / "earth_topography_0p5deg.bin"
DATA_JSON = ROOT / "core" / "data" / "earth_topography_0p5deg.json"
FIGURE = ROOT / "scripts" / "earth_topography_source.png"

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


def draw(z: np.ndarray, sample_lon: np.ndarray, sample_lat: np.ndarray,
         land_fraction: np.ndarray, land_height: np.ndarray, info: dict) -> None:
    ocean = "#dfe7ee"
    ink, muted, grid = "#0b0b0b", "#898781", "#ffffff"
    land_cmap = LinearSegmentedColormap.from_list("land", ["#e9e2cf", "#c9a96e", "#8a5a2b", "#4a2c12"])
    fl_cmap = LinearSegmentedColormap.from_list("fl", [ocean, "#a9c4a4", "#4f7a45"])
    lon, lat = cell_centres()

    fig, axes = plt.subplots(3, 1, figsize=(10, 13.5), constrained_layout=True)
    fig.patch.set_facecolor("#fcfcfb")

    ax = axes[0]
    ax.set_facecolor(ocean)
    step = 2  # 12 arc-minute drawing resolution keeps the PNG small
    masked = np.ma.masked_less_equal(z[::step, ::step], 0.0)
    im = ax.pcolormesh(sample_lon[::step], sample_lat[::step], masked, cmap=land_cmap, vmin=0, vmax=5000,
                       shading="nearest", rasterized=True)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, extend="max")
    cb.set_label("surface elevation [m] (ocean: z <= 0, flat)")
    ax.set_title(f"ETOPO 2022 ice surface, 6' point samples   land (z > 0) fraction of area = "
                 f"{info['global_land_fraction']:.3f}", fontsize=10, loc="left")

    ax = axes[1]
    ax.set_facecolor(ocean)
    im = ax.pcolormesh(lon, lat, land_fraction, cmap=fl_cmap, vmin=0, vmax=1, shading="nearest", rasterized=True)
    ax.contour(lon, lat, land_fraction, levels=[0.5], colors=ink, linewidths=0.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("land fraction f(0)")
    ax.set_title(f"0.5 deg cells: land fraction ({info['samples_per_cell']} samples per cell); "
                 "line: f(0) = 0.5", fontsize=10, loc="left")

    ax = axes[2]
    ax.set_facecolor(ocean)
    masked = np.ma.masked_less_equal(land_height, 0.0)
    im = ax.pcolormesh(lon, lat, masked, cmap=land_cmap, vmin=0, vmax=5000, shading="nearest", rasterized=True)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, extend="max")
    cb.set_label("mean land height h(0) [m]")
    peak_lon, peak_lat = info["maximum_land_height_lon_lat"]
    ax.set_title(f"0.5 deg cells: mean of max(z, 0)   max {info['maximum_land_height_m']:.0f} m at "
                 f"({peak_lon:.1f}E, {peak_lat:.1f}N)", fontsize=10, loc="left")
    for name, x, y in LANDMARKS:
        ax.plot(x, y, marker="o", markersize=5, markerfacecolor="none", markeredgecolor="#2a78d6",
                markeredgewidth=1.2)
        ax.annotate(f"{name} {info['landmarks_m'][name]:.0f} m", (x, y), xytext=(6, 6),
                    textcoords="offset points", fontsize=7.5, color="#0b0b0b",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=1.0))

    for ax in axes:
        ax.set_xlim(0, 360)
        ax.set_ylim(-90, 90)
        ax.set_xticks(np.arange(0, 361, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
        ax.set_xlabel("longitude [deg]", color=muted)
        ax.set_ylabel("latitude [deg]", color=muted)
        ax.tick_params(colors=muted)
        ax.grid(True, color=grid, alpha=0.6, linewidth=0.5)
        ax.set_aspect("equal")
        for spine in ax.spines.values():
            spine.set_color("#c3c2b7")
    fig.savefig(FIGURE, dpi=130)


def main() -> None:
    z, lat, lon = fetch_samples()
    if not np.isfinite(z).all() or (z < -20000).any():
        raise RuntimeError("fill values present in the ETOPO samples")
    z, lon = to_model_longitudes(z, lon)
    land_fraction, land_height, samples_per_cell = aggregate(z, lat, lon)
    info = write_intermediate(land_fraction, land_height, samples_per_cell)
    draw(z, lon, lat, land_fraction, land_height, info)
    print(json.dumps({k: info[k] for k in ("global_land_fraction", "maximum_land_height_m",
                                           "maximum_land_height_lon_lat", "landmarks_m", "sha256")},
                     indent=2))
    print(f"wrote {DATA_BIN} ({DATA_BIN.stat().st_size} bytes), {DATA_JSON}, {FIGURE}")


if __name__ == "__main__":
    main()
