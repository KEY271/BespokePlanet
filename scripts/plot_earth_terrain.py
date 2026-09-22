"""Scan the smoothing width of docs/dynamics/earth-topography.md and draw the result.

Reproduces the Fortran generation in numpy on the octahedral Gaussian grid:
grid-cell average of the 0.5 degree land fraction (no smoothing), Gaussian-kernel
smoothing of the land height with width s = c * 180 / T, triangular truncation
of g z_s at T, and the diagnostics of the document (minimum, open-ocean ripple,
truncation and total RMS, landmarks) for c in {0.5, 0.75, 1.0, 1.25} and
T in {31, 63}.  Writes
scripts/earth_terrain_scan.png (diagnostics against c) and
scripts/earth_terrain.png (target and truncated z_s for the selected c).  Run with
`uv run --project ~/.local/share/llm-python python scripts/plot_earth_terrain.py`.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln, lpmv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[1]
DATA_BIN = ROOT / "core" / "data" / "earth_topography_0p5deg.bin"
GRAVITY = 9.80616
TRUNCATIONS = (31, 63)
SCALE_FACTORS = (0.5, 0.75, 1.0, 1.25)
WINDOW_FACTOR = 3.0
MIN_HEIGHT_CRITERION = -50.0  # m
OCEAN_RMS_CRITERION = 10.0  # m
OPEN_OCEAN_LAND_FRACTION = 0.01
LANDMARKS = [("Tibet", 90.0, 33.0), ("Andes", 292.0, -18.0), ("Antarctica", 90.0, -80.0),
             ("Greenland", 320.0, 72.0), ("Pacific", 200.0, 0.0)]


# --- intermediate file -----------------------------------------------------------------------

def load_intermediate():
    raw = np.frombuffer(DATA_BIN.read_bytes(), dtype="<f8")
    nlon, nlat = 720, 360
    if raw.size != 2 * nlon * nlat:
        raise RuntimeError(f"{DATA_BIN} has {raw.size} values; expected {2 * nlon * nlat}")
    land = raw[:nlon * nlat].reshape(nlat, nlon)
    height = raw[nlon * nlat:].reshape(nlat, nlon)
    lon = (np.arange(nlon) + 0.5) * 0.5
    lat = -90.0 + (np.arange(nlat) + 0.5) * 0.5
    return land, height, lon, lat


def unit_vectors(lon_deg, lat_deg):
    lon, lat = np.broadcast_arrays(np.radians(lon_deg), np.radians(lat_deg))
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


# --- octahedral Gaussian grid and spectral transform --------------------------------------

class Octahedral:
    def __init__(self, truncation: int):
        self.T = truncation
        L = 2 * (truncation + 1)
        self.mu, self.weights = leggauss(L)  # ascending mu: south to north
        self.lat = np.degrees(np.arcsin(self.mu))
        half = [20 + 4 * j for j in range(truncation + 1)]
        self.nlon = np.array(half + half[::-1])
        self.lons = [np.arange(n) * 360.0 / n for n in self.nlon]
        self.pnm = normalized_legendre_table(truncation, self.mu)  # [n, m, j]

    def analyze(self, field_rings):
        """Triangular truncation of a field given as a list of rings (values along longitude)."""
        T = self.T
        coefficients = np.zeros((T + 1, T + 1), dtype=complex)  # [n, m]
        for j, values in enumerate(field_rings):
            m_max = min(T, self.nlon[j] // 2 - 1)
            fourier = np.fft.rfft(values) / self.nlon[j]
            coefficients[:, :m_max + 1] += self.weights[j] * self.pnm[:, :m_max + 1, j] * fourier[None, :m_max + 1]
        return coefficients

    def synthesize_rings(self, coefficients):
        return [synthesize_ring(coefficients, self.pnm[:, :, j], self.lons[j], min(self.T, self.nlon[j] // 2 - 1))
                for j in range(len(self.nlon))]

    def global_mean(self, rings):
        return 0.5 * sum(w * np.mean(r) for w, r in zip(self.weights, rings))

    def nearest(self, rings, lon0, lat0):
        j = int(np.abs(self.lat - lat0).argmin())
        i = int(np.round(lon0 / (360.0 / self.nlon[j]))) % self.nlon[j]
        return float(rings[j][i])


def normalized_legendre_table(truncation, mu):
    """P_n^m(mu) normalized to integral over [-1,1] of P^2 = 1, shape [n, m, len(mu)]."""
    table = np.zeros((truncation + 1, truncation + 1, mu.size))
    for m in range(truncation + 1):
        for n in range(m, truncation + 1):
            log_norm = 0.5 * (np.log(2 * n + 1) - np.log(2.0) + gammaln(n - m + 1) - gammaln(n + m + 1))
            table[n, m, :] = lpmv(m, n, mu) * np.exp(log_norm)
    return table


def synthesize_ring(coefficients, pnm_j, lons_deg, m_max):
    fourier = np.einsum("nm,nm->m", coefficients[:, :m_max + 1], pnm_j[:, :m_max + 1])
    phase = np.exp(1j * np.outer(np.arange(m_max + 1), np.radians(lons_deg)))
    factor = np.where(np.arange(m_max + 1) == 0, 1.0, 2.0)
    return np.real(np.einsum("m,mk->k", factor * fourier, phase))


def synthesize_regular(coefficients, truncation, lon_deg, lat_deg):
    pnm = normalized_legendre_table(truncation, np.sin(np.radians(lat_deg)))
    return np.array([synthesize_ring(coefficients, pnm[:, :, j], lon_deg, truncation) for j in range(lat_deg.size)])


# --- smoothing and box averaging ---------------------------------------------------------------

class Smoother:
    def __init__(self, land, height, lon, lat):
        self.land, self.height, self.lon, self.lat = land, height, lon, lat
        self.vectors = unit_vectors(lon[None, :], lat[:, None])  # [lat, lon, 3]
        self.area = np.cos(np.radians(lat))

    def rows(self, lats_deg, lon_lists, s_deg):
        """Kernel-smoothed (land fraction, height) for rows of points at lats_deg with the given
        longitudes.  The smoothed land fraction only delimits the open ocean in the diagnostics;
        the model's land fraction is the box average."""
        s = np.radians(s_deg)
        cutoff = WINDOW_FACTOR * s
        land_rows, height_rows = [], []
        for lat0, lons in zip(lats_deg, lon_lists):
            band = np.abs(self.lat - lat0) <= WINDOW_FACTOR * s_deg
            cells = self.vectors[band].reshape(-1, 3)
            weights_area = np.repeat(self.area[band], self.lon.size)
            points = unit_vectors(lons, np.full(lons.shape, lat0))
            theta = np.arccos(np.clip(points @ cells.T, -1.0, 1.0))
            kernel = np.exp(-(theta / s) ** 2) * (theta <= cutoff) * weights_area[None, :]
            norm = kernel.sum(axis=1)
            land_rows.append(kernel @ self.land[band].ravel() / norm)
            height_rows.append(kernel @ self.height[band].ravel() / norm)
        return land_rows, height_rows

    def box_average(self, grid: Octahedral, field):
        """Intermediate cells of `field` area-averaged over each grid point's latitude band
        (between ring midpoints) and nearest-longitude sector; this is the model's f_L and
        the raw height of the total-error diagnostic."""
        edges = np.concatenate([[-90.0], 0.5 * (grid.lat[1:] + grid.lat[:-1]), [90.0]])
        ring_of_cell = np.clip(np.searchsorted(edges, self.lat) - 1, 0, grid.lat.size - 1)
        sums = [np.zeros(n) for n in grid.nlon]
        counts = [np.zeros(n) for n in grid.nlon]
        for jc, j in enumerate(ring_of_cell):
            n = grid.nlon[j]
            i = np.round(self.lon / (360.0 / n)).astype(int) % n
            np.add.at(sums[j], i, self.area[jc] * field[jc])
            np.add.at(counts[j], i, self.area[jc])
        return [s / c for s, c in zip(sums, counts)]


# --- diagnostics ----------------------------------------------------------------------------------

def rms(grid, rings):
    return float(np.sqrt(grid.global_mean([r ** 2 for r in rings])))


def evaluate(grid: Octahedral, smoother: Smoother, c: float, land_rings, raw_rings):
    s_deg = c * 180.0 / grid.T
    smoothed_land_rings, target_rings = smoother.rows(grid.lat, grid.lons, s_deg)
    coefficients = grid.analyze([GRAVITY * r for r in target_rings])
    truncated_rings = [r / GRAVITY for r in grid.synthesize_rings(coefficients)]
    flat = np.concatenate(truncated_rings)
    # open ocean: beyond the kernel's reach of any land, where z_s should be exactly zero
    ocean = np.concatenate(smoothed_land_rings) < OPEN_OCEAN_LAND_FRACTION
    # every ocean grid point of the model, including the coast the smoothing leaks height onto
    coastal = np.concatenate(land_rings) < OPEN_OCEAN_LAND_FRACTION
    ring_index = np.repeat(np.arange(grid.lat.size), grid.nlon)
    lon_flat = np.concatenate(grid.lons)
    k_min, k_max = int(flat.argmin()), int(flat.argmax())
    area = np.repeat(grid.weights / grid.nlon, grid.nlon)
    ocean_area = area[ocean]
    result = {
        "T": grid.T, "c": c, "s_deg": s_deg,
        "land_fraction": float(grid.global_mean(land_rings)),
        "min_m": float(flat[k_min]), "min_lon_lat": (float(lon_flat[k_min]), float(grid.lat[ring_index[k_min]])),
        "max_m": float(flat[k_max]), "max_lon_lat": (float(lon_flat[k_max]), float(grid.lat[ring_index[k_max]])),
        "truncation_rms_m": rms(grid, [a - b for a, b in zip(truncated_rings, target_rings)]),
        "total_rms_m": rms(grid, [a - b for a, b in zip(truncated_rings, raw_rings)]),
        "ocean_rms_m": float(np.sqrt(np.sum(ocean_area * flat[ocean] ** 2) / np.sum(ocean_area))),
        "ocean_min_m": float(flat[ocean].min()),
        "ocean_height_rms_m": float(np.sqrt(np.sum(area[coastal] * flat[coastal] ** 2) / np.sum(area[coastal]))),
        "ocean_height_max_m": float(flat[coastal].max()),
        "landmarks_m": {name: grid.nearest(truncated_rings, x, y) for name, x, y in LANDMARKS},
    }
    return result, coefficients


def passes(result):
    return result["min_m"] >= MIN_HEIGHT_CRITERION and result["ocean_rms_m"] <= OCEAN_RMS_CRITERION


def markdown_table(results):
    names = [name for name, *_ in LANDMARKS]
    lines = ["| T | c | s [deg] | <f_L> | min z_s [m] (lon, lat) | open-ocean RMS [m] | open-ocean min [m] | "
             "ocean z_s RMS [m] | ocean z_s max [m] | trunc. RMS [m] | total RMS [m] | max z_s [m] | "
             + " | ".join(names) + " | pass |",
             "|" + "---|" * (13 + len(names))]
    for r in results:
        lines.append(
            f"| {r['T']} | {r['c']:.2f} | {r['s_deg']:.2f} | {r['land_fraction']:.4f} | "
            f"{r['min_m']:.0f} ({r['min_lon_lat'][0]:.0f}, {r['min_lon_lat'][1]:.0f}) | {r['ocean_rms_m']:.1f} | "
            f"{r['ocean_min_m']:.0f} | {r['ocean_height_rms_m']:.0f} | {r['ocean_height_max_m']:.0f} | "
            f"{r['truncation_rms_m']:.1f} | {r['total_rms_m']:.0f} | {r['max_m']:.0f} | "
            + " | ".join(f"{r['landmarks_m'][n]:.0f}" for n in names) + f" | {'yes' if passes(r) else 'no'} |")
    return "\n".join(lines)


# --- figures ---------------------------------------------------------------------------------------

INK, MUTED, GRIDLINE, SURFACE = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"
SERIES = {31: "#2a78d6", 63: "#eb6834"}  # categorical slots 1 and 2
OCEAN = "#dfe7ee"
LAND_CMAP = LinearSegmentedColormap.from_list("land", ["#e9e2cf", "#c9a96e", "#8a5a2b", "#4a2c12"])


def style_axis(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRIDLINE, linewidth=0.6)
    ax.tick_params(colors=MUTED, labelsize=8)
    for name, spine in ax.spines.items():
        spine.set_visible(name in ("left", "bottom"))
        spine.set_color("#c3c2b7")


def draw_scan(results, selected_c, path):
    panels = [("minimum z_s [m]", "min_m", MIN_HEIGHT_CRITERION),
              ("open-ocean RMS of z_s [m]", "ocean_rms_m", OCEAN_RMS_CRITERION),
              ("truncation RMS [m]", "truncation_rms_m", None),
              ("total RMS vs raw 0.5 deg cells [m]", "total_rms_m", None),
              ("z_s RMS over ocean grid points (f_L < 0.01) [m]", "ocean_height_rms_m", None),
              ("max z_s over ocean grid points [m]", "ocean_height_max_m", None),
              ("Tibet z_s [m]", ("landmarks_m", "Tibet"), None),
              ("Andes z_s [m]", ("landmarks_m", "Andes"), None)]
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, (title, key, criterion) in zip(axes.ravel(), panels):
        style_axis(ax)
        for T in TRUNCATIONS:
            rows = [r for r in results if r["T"] == T]
            xs = [r["c"] for r in rows]
            ys = [r[key] if isinstance(key, str) else r[key[0]][key[1]] for r in rows]
            ax.plot(xs, ys, color=SERIES[T], linewidth=2, marker="o", markersize=7, label=f"T = {T}")
            ax.annotate(f"{ys[-1]:.0f}", (xs[-1], ys[-1]), xytext=(6, 0), textcoords="offset points",
                        fontsize=8, color=INK, va="center")
        if criterion is not None:
            ax.axhline(criterion, color=MUTED, linewidth=1, linestyle="--")
            ax.annotate(f"criterion {criterion:g}", (SCALE_FACTORS[0], criterion), xytext=(0, 4),
                        textcoords="offset points", fontsize=8, color=MUTED)
        ax.axvline(selected_c, color="#c3c2b7", linewidth=1)
        ax.set_title(title, fontsize=10, loc="left", color=INK)
        ax.set_xlabel("c   (s = c * 180 deg / T)", color=MUTED, fontsize=9)
        ax.set_xticks(SCALE_FACTORS)
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"Earth topography: smoothing-width scan (selected c = {selected_c:.2f}, vertical line)",
                 fontsize=11, color=INK, x=0.01, ha="left")
    fig.savefig(path, dpi=130)


def draw_maps(selected, smoother, path):
    lon = (np.arange(240) + 0.5) * 1.5
    lat = -90.0 + (np.arange(120) + 0.5) * 1.5
    fig, axes = plt.subplots(2, 2, figsize=(16, 8.6), constrained_layout=True)
    fig.patch.set_facecolor(SURFACE)
    for row, (T, result, coefficients) in enumerate(selected):
        target = np.array(smoother.rows(lat, [lon] * lat.size, result["s_deg"])[1])
        truncated = synthesize_regular(coefficients, T, lon, lat) / GRAVITY
        for col, (field, title) in enumerate([
                (target, f"T = {T}: target z_s after smoothing (s = {result['s_deg']:.1f} deg); line: 0.5 deg f(0) = 0.5"),
                (truncated, f"T = {T}: z_s after truncation   min {result['min_m']:.0f} m, max {result['max_m']:.0f} m, "
                            f"open-ocean RMS {result['ocean_rms_m']:.1f} m")]):
            ax = axes[row, col]
            ax.set_facecolor(OCEAN)
            masked = np.ma.masked_less(field, 50.0)
            im = ax.pcolormesh(lon, lat, masked, cmap=LAND_CMAP, vmin=0, vmax=4000, shading="nearest", rasterized=True)
            ax.contour(smoother.lon, smoother.lat, smoother.land, levels=[0.5], colors=INK, linewidths=0.5,
                       linestyles="--" if col else "-")
            if col == 1:
                ax.contour(lon, lat, field, levels=[-50, -20], colors=["#b03a2e", "#e08a7a"], linewidths=0.8)
            cb = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02, extend="max")
            cb.set_label("z_s [m]", color=MUTED)
            cb.ax.tick_params(colors=MUTED, labelsize=8)
            ax.set_title(title, fontsize=9.5, loc="left", color=INK)
            ax.set_xlim(0, 360)
            ax.set_ylim(-90, 90)
            ax.set_xticks(np.arange(0, 361, 60))
            ax.set_yticks(np.arange(-90, 91, 30))
            ax.tick_params(colors=MUTED, labelsize=8)
            ax.grid(True, color="#ffffff", alpha=0.6, linewidth=0.5)
            ax.set_aspect("equal")
    axes[1, 1].text(2, -88, "red: z_s = -20 / -50 m (Gibbs undershoot); dashed: 0.5 deg cell land fraction = 0.5",
                    fontsize=8, color="#333333", va="bottom",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))
    fig.savefig(path, dpi=130)


# --- main ---------------------------------------------------------------------------------------------

def main():
    started = time.time()
    land, height, lon, lat = load_intermediate()
    smoother = Smoother(land, height, lon, lat)
    results, coefficient_store = [], {}
    for T in TRUNCATIONS:
        grid = Octahedral(T)
        raw_rings = smoother.box_average(grid, smoother.height)
        land_rings = [np.clip(r, 0.0, 1.0) for r in smoother.box_average(grid, smoother.land)]
        for c in SCALE_FACTORS:
            result, coefficients = evaluate(grid, smoother, c, land_rings, raw_rings)
            results.append(result)
            coefficient_store[(T, c)] = coefficients
            print(f"T={T} c={c:.2f} done ({time.time() - started:.0f} s)", file=sys.stderr, flush=True)
    candidates = [c for c in SCALE_FACTORS if all(passes(r) for r in results if r["c"] == c)]
    selected_c = candidates[0] if candidates else SCALE_FACTORS[-1]
    print(markdown_table(results))
    print(f"\nselected c = {selected_c:.2f}" + ("" if candidates else " (no c meets both criteria; largest used)"))
    draw_scan(results, selected_c, ROOT / "scripts" / "earth_terrain_scan.png")
    selected = [(T, next(r for r in results if r["T"] == T and r["c"] == selected_c), coefficient_store[(T, selected_c)])
                for T in TRUNCATIONS]
    draw_maps(selected, smoother, ROOT / "scripts" / "earth_terrain.png")
    (ROOT / "scripts" / "earth_terrain_scan.json").write_text(json.dumps(results, indent=1) + "\n")
    print(f"figures written ({time.time() - started:.0f} s)", file=sys.stderr)


if __name__ == "__main__":
    main()
