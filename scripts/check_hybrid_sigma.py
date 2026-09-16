#!/usr/bin/env python3
# /// script
# requires-python = ">=3.9"
# dependencies = ["matplotlib>=3.10,<4"]
# ///
"""Check the fixed N=12 hybrid-sigma A/B profile and plot the pressures.

The coordinate is ordered from the model top to the surface (eta=1).
Interface pressure is calculated directly from the prescribed coefficients:

    p_i(ps) = A_i + B_i ps

The figure shows both interface pressures and layer pressure thicknesses.  A
profile is strictly monotonic from the model top to the surface exactly when
all of its pressure thicknesses are positive.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_SURFACE_PRESSURES_HPA = (100.0, 200.0, 400.0, 600.0, 800.0, 1000.0)
REFERENCE_SURFACE_PRESSURE_PA = 100000.0
FIXED_A_PA = (
    100.0,
    300.0,
    1000.0,
    5000.0,
    10000.0,
    8000.0,
    8000.0,
    10000.0,
    12000.0,
    10000.0,
    7000.0,
    3000.0,
    0.0,
)
FIXED_B = (0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 1.0)


def make_coefficients(
    a_pa: Sequence[float],
    b: Sequence[float],
) -> tuple[list[float], list[float], list[float]]:
    """Return eta, A [hPa], and dimensionless B at all interfaces."""
    if len(a_pa) < 2:
        raise ValueError("at least two reference interfaces are required")
    if len(a_pa) != len(b):
        raise ValueError("A and B must have the same number of interfaces")
    if any(value < 0.0 for value in a_pa):
        raise ValueError("A must be nonnegative")
    if any(not 0.0 <= value <= 1.0 for value in b):
        raise ValueError("B must lie inside [0, 1]")
    if any(lower < upper for upper, lower in zip(b, b[1:])):
        raise ValueError("B must be nondecreasing")

    eta = [
        a / REFERENCE_SURFACE_PRESSURE_PA + weight
        for a, weight in zip(a_pa, b)
    ]
    a_hpa = [value / 100.0 for value in a_pa]
    return eta, a_hpa, list(b)


def interface_pressures(
    a_hpa: Sequence[float], b: Sequence[float], surface_pressure_hpa: float
) -> list[float]:
    return [a + weight * surface_pressure_hpa for a, weight in zip(a_hpa, b)]


def non_increasing_layers(pressures_hpa: Sequence[float]) -> list[int]:
    """Return k for every invalid interval p[k+1] <= p[k]."""
    return [
        k
        for k, (upper, lower) in enumerate(
            zip(pressures_hpa, pressures_hpa[1:])
        )
        if lower <= upper
    ]


def critical_surface_pressure_hpa(a_hpa: Sequence[float], b: Sequence[float]) -> float:
    """Find the lower ps bound implied by delta A + delta B * ps > 0.

    This profile has nondecreasing B, so each layer with delta B > 0 gives
    a lower bound on ps.  Layers with delta B == 0 must have delta A > 0.
    """
    lower_bounds: list[float] = []
    for k, (a0, a1, b0, b1) in enumerate(zip(a_hpa, a_hpa[1:], b, b[1:])):
        delta_a = a1 - a0
        delta_b = b1 - b0
        if delta_b < 0.0:
            raise ValueError(f"B decreases across layer {k}; no single lower bound exists")
        if delta_b == 0.0:
            if delta_a <= 0.0:
                raise ValueError(f"layer {k} is never strictly monotonic")
            continue
        lower_bounds.append(-delta_a / delta_b)
    return max(lower_bounds, default=float("-inf"))


def print_report(
    eta: Sequence[float],
    a_hpa: Sequence[float],
    b: Sequence[float],
    surface_pressures_hpa: Iterable[float],
) -> dict[float, list[float]]:
    n_layers = len(eta) - 1
    print(f"A/B coefficients (N={n_layers} gives {len(eta)} layer interfaces)")
    print("interface   eta       A [hPa]           B")
    for k, (eta_k, a_k, b_k) in enumerate(zip(eta, a_hpa, b)):
        print(f"{k:9d}  {eta_k:7.4f}  {a_k:12.6f}  {b_k:10.7f}")

    profiles: dict[float, list[float]] = {}
    print("\nInterface pressures [hPa], ordered from model top to surface")
    header = "ps [hPa]  " + "".join(f"{k:>10d}" for k in range(len(eta))) + "  result"
    print(header)
    for ps in surface_pressures_hpa:
        pressure = interface_pressures(a_hpa, b, ps)
        profiles[ps] = pressure
        invalid = non_increasing_layers(pressure)
        result = "PASS" if not invalid else f"FAIL at layers {invalid}"
        values = "".join(f"{value:10.3f}" for value in pressure)
        print(f"{ps:8.1f}  {values}  {result}")

    critical = critical_surface_pressure_hpa(a_hpa, b)
    print(
        "\nStrict monotonicity requires "
        f"ps > {critical:.6f} hPa for this discretized A/B profile."
    )
    return profiles


def plot_profiles(
    output: Path,
    eta: Sequence[float],
    profiles: dict[float, Sequence[float]],
    critical_ps_hpa: float,
) -> None:
    interfaces = list(range(len(eta)))
    layer_centers = [k + 0.5 for k in range(len(eta) - 1)]
    colors = ("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#000000")

    figure, (pressure_axis, thickness_axis) = plt.subplots(
        2,
        1,
        figsize=(10, 9),
        sharex=True,
        gridspec_kw={"height_ratios": (2.0, 1.0)},
    )
    all_delta_pressures: list[float] = []
    for index, (surface_pressure, pressure) in enumerate(profiles.items()):
        color = colors[index % len(colors)]
        invalid = non_increasing_layers(pressure)
        result = "PASS" if not invalid else "FAIL"
        label = rf"$p_s$={surface_pressure:g} hPa ({result})"
        pressure_axis.plot(
            interfaces,
            pressure,
            marker="o",
            linewidth=2,
            markersize=4,
            color=color,
            label=label,
        )

        delta_pressure = [lower - upper for upper, lower in zip(pressure, pressure[1:])]
        all_delta_pressures.extend(delta_pressure)
        thickness_axis.plot(
            layer_centers,
            delta_pressure,
            marker="o",
            linewidth=2,
            markersize=4,
            color=color,
        )

    pressure_axis.set_title(
        f"Hybrid-sigma coordinate (N={len(eta) - 1}, fixed A/B coefficients)\n"
        rf"strictly increasing only for $p_s > {critical_ps_hpa:.3f}$ hPa"
    )
    pressure_axis.set_ylabel("interface pressure [hPa]")
    pressure_axis.legend(loc="upper left", ncols=2, fontsize=9)
    pressure_axis.grid(True, alpha=0.3)

    thickness_axis.axhline(0.0, color="black", linewidth=1.2)
    delta_span = max(all_delta_pressures) - min(all_delta_pressures)
    y_min = min(0.0, min(all_delta_pressures)) - 0.08 * delta_span
    y_max = max(all_delta_pressures) + 0.08 * delta_span
    thickness_axis.set_ylim(y_min, y_max)
    thickness_axis.axhspan(
        y_min, 0.0, color="#d62728", alpha=0.08, label=r"invalid: $\Delta p \leq 0$"
    )
    thickness_axis.set_xlim(-0.2, len(eta) - 0.8)
    thickness_axis.set_xticks(interfaces)
    thickness_axis.set_xlabel(
        r"vertical index (interfaces at $k$, layers at $k+1/2$; top to surface)"
    )
    thickness_axis.set_ylabel(r"layer $\Delta p$ [hPa]")
    thickness_axis.legend(loc="lower left", fontsize=9)
    thickness_axis.grid(True, alpha=0.3)

    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("hybrid_sigma_pressure.png"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    eta, a_hpa, b = make_coefficients(FIXED_A_PA, FIXED_B)
    profiles = print_report(eta, a_hpa, b, DEFAULT_SURFACE_PRESSURES_HPA)
    critical = critical_surface_pressure_hpa(a_hpa, b)
    plot_profiles(args.output, eta, profiles, critical)
    print(f"\nPlot written to {args.output}")


if __name__ == "__main__":
    main()
