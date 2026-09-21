"""Draw the default planet of docs/dynamics/topography.md.

Panels: the analytic land fraction f_L with the analytic surface height z_s,
and z_s after triangular spectral truncation at T31 and T63 (the field the
solver actually sees).  Run with
`uv run --project ~/.local/share/llm-python python scripts/plot_default_terrain.py`.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import gammaln, lpmv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

GRAVITY = 9.80616
COAST_WIDTH = 6.0  # w [deg]
BASE_HEIGHT = 300.0  # h_0 [m]

# (lambda_c, phi_c, a, b, theta)
CONTINENTS = {"A": (60.0, 45.0, 55.0, 28.0, 0.0),
              "B": (300.0, -10.0, 28.0, 40.0, 0.0),
              "C": (150.0, -30.0, 22.0, 16.0, 0.0)}
CAPS = {"S": (70.0, -1.0)}  # (phi_0, sign)
# (lambda_c, phi_c, theta, L, sigma, H)
RIDGES = {"1": (75.0, 35.0, 0.0, 50.0, 6.0, 2500.0),
          "2": (295.0, -10.0, 90.0, 60.0, 6.0, 2500.0)}


def local_frame(lon, lat, lon_c, lat_c, theta):
    dlon = ((lon - lon_c + 180.0) % 360.0) - 180.0
    x = np.cos(np.radians(lat)) * dlon
    y = lat - lat_c
    t = np.radians(theta)
    return x * np.cos(t) + y * np.sin(t), -x * np.sin(t) + y * np.cos(t)


def land_fraction(lon, lat):
    one_minus = np.ones_like(lon)
    for lon_c, lat_c, a, b, theta in CONTINENTS.values():
        xp, yp = local_frame(lon, lat, lon_c, lat_c, theta)
        d = np.sqrt((xp / a) ** 2 + (yp / b) ** 2)
        one_minus *= 1.0 - 0.5 * (1.0 - np.tanh((d - 1.0) * min(a, b) / COAST_WIDTH))
    for phi0, sign in CAPS.values():
        one_minus *= 1.0 - 0.5 * (1.0 + np.tanh((sign * lat - phi0) / COAST_WIDTH))
    return 1.0 - one_minus


def surface_height(lon, lat):
    ridges = np.zeros_like(lon)
    for lon_c, lat_c, theta, length, sigma, height in RIDGES.values():
        xp, yp = local_frame(lon, lat, lon_c, lat_c, theta)
        ridges += height * np.exp(-((yp / sigma) ** 2)) * 0.5 * (1.0 - np.tanh((np.abs(xp) - length / 2) / sigma))
    return land_fraction(lon, lat) * (BASE_HEIGHT + ridges)


def normalized_legendre(n, m, mu):
    log_norm = 0.5 * (np.log(2 * n + 1) - np.log(2.0) + gammaln(n - m + 1) - gammaln(n + m + 1))
    return lpmv(m, n, mu) * np.exp(log_norm)


def truncate(field_fn, truncation, lon_out, lat_out):
    """Triangular truncation via Gaussian quadrature analysis and synthesis on an output grid."""
    nlat = 2 * (truncation + 1)
    mu, weights = leggauss(nlat)
    nlon = 2 * nlat
    lon = np.arange(nlon) * 360.0 / nlon
    lon2, lat2 = np.meshgrid(lon, np.degrees(np.arcsin(mu)))
    fourier = np.fft.rfft(field_fn(lon2, lat2), axis=1) / nlon  # [lat, m]
    mu_out = np.sin(np.radians(lat_out))
    out = np.zeros((lat_out.size, lon_out.size))
    for m in range(truncation + 1):
        synth_m = np.zeros(lat_out.size, dtype=complex)
        for n in range(m, truncation + 1):
            coefficient = np.sum(weights * fourier[:, m] * normalized_legendre(n, m, mu))
            synth_m += coefficient * normalized_legendre(n, m, mu_out)
        phase = np.exp(1j * m * np.radians(lon_out))
        contribution = np.real(synth_m[:, None] * phase[None, :])
        out += contribution if m == 0 else 2.0 * contribution
    return out


def area_mean(field, lat):
    w = np.cos(np.radians(lat))[:, None] * np.ones_like(field)
    return float((field * w).sum() / w.sum())


def main():
    lon = np.linspace(0, 360, 721)
    lat = np.linspace(-90, 90, 361)
    lon2, lat2 = np.meshgrid(lon, lat)
    f_l = land_fraction(lon2, lat2)
    z_analytic = surface_height(lon2, lat2)
    z_t31 = truncate(surface_height, 31, lon, lat)
    z_t63 = truncate(surface_height, 63, lon, lat)

    ocean = "#dfe7ee"
    land_cmap = LinearSegmentedColormap.from_list("land", ["#e9e2cf", "#c9a96e", "#8a5a2b", "#4a2c12"])
    fig, axes = plt.subplots(3, 1, figsize=(10, 13.5), constrained_layout=True)
    fig.patch.set_facecolor("white")

    titles = [
        f"Analytic: land fraction f_L (fill) and surface height z_s (contours)   <f_L> = {area_mean(f_l, lat):.3f}",
        f"z_s after triangular truncation T31   min {z_t31.min():.0f} m, max {z_t31.max():.0f} m, "
        f"RMS diff {np.sqrt(area_mean((z_t31 - z_analytic) ** 2, lat)):.0f} m",
        f"z_s after triangular truncation T63   min {z_t63.min():.0f} m, max {z_t63.max():.0f} m, "
        f"RMS diff {np.sqrt(area_mean((z_t63 - z_analytic) ** 2, lat)):.0f} m",
    ]
    levels = np.arange(0, 3001, 250)

    ax = axes[0]
    ax.set_facecolor(ocean)
    fl_cmap = LinearSegmentedColormap.from_list("fl", [ocean, "#a9c4a4", "#4f7a45"])
    im = ax.pcolormesh(lon, lat, f_l, cmap=fl_cmap, vmin=0, vmax=1, shading="auto", rasterized=True)
    cs = ax.contour(lon, lat, z_analytic, levels=levels[1:], colors="#4a2c12", linewidths=0.7)
    ax.clabel(cs, levels=[500, 1500, 2500], fmt="%d m", fontsize=7)
    ax.contour(lon, lat, f_l, levels=[0.5], colors="#1f1f1f", linewidths=1.0)
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("land fraction f_L")
    for name, (lon_c, lat_c, *_rest) in CONTINENTS.items():
        ax.text(lon_c, lat_c, name, ha="center", va="center", fontsize=12, fontweight="bold", color="#1f1f1f")
    ax.text(180, -82, "S cap", ha="center", va="center", fontsize=10, fontweight="bold", color="#1f1f1f")
    for name, (lon_c, lat_c, *_rest) in RIDGES.items():
        ax.annotate(f"ridge {name}", (lon_c, lat_c), xytext=(lon_c + 18, lat_c - 14), fontsize=8,
                    arrowprops=dict(arrowstyle="-", color="#4a2c12", lw=0.7), color="#4a2c12")

    for ax, field in zip(axes[1:], [z_t31, z_t63]):
        ax.set_facecolor(ocean)
        masked = np.ma.masked_less(field, 50.0)
        im = ax.pcolormesh(lon, lat, masked, cmap=land_cmap, vmin=0, vmax=3000, shading="auto", rasterized=True)
        ax.contour(lon, lat, field, levels=[-50, -20], colors=["#b03a2e", "#e08a7a"], linewidths=0.7)
        ax.contour(lon, lat, f_l, levels=[0.5], colors="#1f1f1f", linewidths=0.6, linestyles="--")
        cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, extend="min")
        cb.set_label("z_s [m]")

    for ax, title in zip(axes, titles):
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlim(0, 360)
        ax.set_ylim(-90, 90)
        ax.set_xticks(np.arange(0, 361, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
        ax.set_xlabel("longitude [deg]")
        ax.set_ylabel("latitude [deg]")
        ax.grid(True, color="#ffffff", alpha=0.6, linewidth=0.5)
        ax.set_aspect("equal")
    axes[1].text(2, -88, "red: z_s = -20 / -50 m (Gibbs undershoot); dashed: analytic coast f_L = 0.5",
                 fontsize=7.5, color="#333333", va="bottom", bbox=dict(facecolor="white", edgecolor="none", alpha=0.8))

    fig.savefig("scripts/default_terrain.png", dpi=130)
    print(f"land fraction {area_mean(f_l, lat):.4f}")
    print(f"T31 min/max {z_t31.min():.1f} {z_t31.max():.1f}; T63 min/max {z_t63.min():.1f} {z_t63.max():.1f}")


if __name__ == "__main__":
    main()
