#!/usr/bin/env python3
"""Isolate the upper-atmosphere temperature/radiation system.

Only the temperatures of the four layers between 1 and 100 hPa are advanced.
The lower eight atmospheric layers and the surface temperature are held at the
reference state, but their upward longwave radiation is retained.  Dynamics,
convection, diffusion, surface exchange, and Rayleigh friction are excluded.
Ozone transmission always uses the vertical optical depth; solar zenith angle
changes the incident flux but does not lengthen the ozone path.

Two questions are tested:

1. Does the autonomous system with globally averaged shortwave absorption have
   a stable radiative-equilibrium fixed point?
2. Does a local equatorial column under the model's diurnal and seasonal solar
   forcing approach a stable annual cycle?

The time integrator reproduces the model's leapfrog/RAW update, including its
first-step midpoint startup.  The script has no third-party dependencies and
writes CSV summaries plus a compact SVG plot.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence


SECONDS_PER_DAY = 86_400.0
DAYS_PER_YEAR = 360
SECONDS_PER_YEAR = DAYS_PER_YEAR * SECONDS_PER_DAY

STEFAN_BOLTZMANN = 5.670374419e-8
LONGWAVE_SURFACE_OPTICAL_DEPTH = 1.0
DRY_AIR_SPECIFIC_HEAT = 1004.0
GRAVITY = 9.80616
SOLAR_CONSTANT = 1361.0
UV_SHORTWAVE_FRACTION = 0.02
AXIAL_TILT = math.radians(23.4)
ORBITAL_PERIOD = SECONDS_PER_YEAR
PLANETARY_ROTATION_RATE = 2.0 * math.pi * (
    1.0 / SECONDS_PER_DAY + 1.0 / ORBITAL_PERIOD
)

OZONE_SHORTWAVE_OPTICAL_DEPTH = 1.5
OZONE_LONGWAVE_OPTICAL_DEPTH = 0.005
OZONE_PRESSURE_LOWER_BOUND = 100.0
OZONE_PRESSURE_UPPER_BOUND = 10_000.0
OZONE_PEAK_PRESSURE = 1_000.0
OZONE_LOG_PRESSURE_WIDTH = math.log(3.0)

RAW_FILTER_EPSILON = 0.1
RAW_FILTER_ALPHA = 0.53

REFERENCE_SURFACE_PRESSURE = 100_000.0
A_HALF_PA = (
    100.0,
    300.0,
    1_000.0,
    5_000.0,
    10_000.0,
    8_000.0,
    8_000.0,
    10_000.0,
    12_000.0,
    10_000.0,
    7_000.0,
    3_000.0,
    0.0,
)
B_HALF = (0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.55, 0.7, 0.85, 1.0)
PRESSURE_HALF_PA = tuple(
    a + b * REFERENCE_SURFACE_PRESSURE for a, b in zip(A_HALF_PA, B_HALF)
)
UPPER_LEVEL_COUNT = 4


Vector = tuple[float, ...]
TendencyFunction = Callable[[Vector, float], Vector]


@dataclass(frozen=True)
class HistoryPoint:
    years: float
    temperatures: Vector
    residual_k_day: float


@dataclass(frozen=True)
class ExperimentResult:
    name: str
    forcing: str
    dt_seconds: float
    initial_scale: float
    converged: bool
    elapsed_years: float
    steps: int
    temperatures: Vector
    max_tendency_k_day: float
    time_level_difference_k: float
    cycle_mismatch_k: float | None
    history: tuple[HistoryPoint, ...]


def reference_temperature(eta: float) -> float:
    """Match dry_vertical_coordinate:jablonowski_mean_temperature."""
    base_temperature = 288.0
    lapse_rate = 0.005
    dry_air_gas_constant = 287.0
    tropopause_eta = 0.2
    stratospheric_adjustment = 4.8e5
    temperature = base_temperature * eta ** (
        dry_air_gas_constant * lapse_rate / GRAVITY
    )
    if eta < tropopause_eta:
        temperature += stratospheric_adjustment * (tropopause_eta - eta) ** 5
    return temperature


def reference_profile() -> Vector:
    return tuple(
        reference_temperature(
            0.5 * (top + bottom) / REFERENCE_SURFACE_PRESSURE
        )
        for top, bottom in zip(PRESSURE_HALF_PA, PRESSURE_HALF_PA[1:])
    )


def ozone_layer_optical_depth(
    pressure_top: float,
    pressure_bottom: float,
    total_optical_depth: float,
) -> float:
    """Distribute one total ozone optical depth over the prescribed profile."""
    pressure_lower = max(pressure_top, OZONE_PRESSURE_LOWER_BOUND)
    pressure_upper = min(pressure_bottom, OZONE_PRESSURE_UPPER_BOUND)
    if pressure_upper <= pressure_lower:
        return 0.0

    scale = math.sqrt(2.0) * OZONE_LOG_PRESSURE_WIDTH
    x_lower = math.log(pressure_lower / OZONE_PEAK_PRESSURE) / scale
    x_upper = math.log(pressure_upper / OZONE_PEAK_PRESSURE) / scale
    x_profile_lower = math.log(
        OZONE_PRESSURE_LOWER_BOUND / OZONE_PEAK_PRESSURE
    ) / scale
    x_profile_upper = math.log(
        OZONE_PRESSURE_UPPER_BOUND / OZONE_PEAK_PRESSURE
    ) / scale
    return total_optical_depth * (
        math.erf(x_upper) - math.erf(x_lower)
    ) / (math.erf(x_profile_upper) - math.erf(x_profile_lower))


LAYER_OZONE_SHORTWAVE_OPTICAL_DEPTH = tuple(
    ozone_layer_optical_depth(top, bottom, OZONE_SHORTWAVE_OPTICAL_DEPTH)
    for top, bottom in zip(PRESSURE_HALF_PA, PRESSURE_HALF_PA[1:])
)
LAYER_OZONE_LONGWAVE_OPTICAL_DEPTH = tuple(
    ozone_layer_optical_depth(top, bottom, OZONE_LONGWAVE_OPTICAL_DEPTH)
    for top, bottom in zip(PRESSURE_HALF_PA, PRESSURE_HALF_PA[1:])
)


def shortwave_absorption_at_cos_zenith(cos_zenith: float) -> Vector:
    """Return absorbed shortwave flux [W m-2] in each model layer."""
    absorbed = [0.0] * (len(PRESSURE_HALF_PA) - 1)
    if cos_zenith <= 0.0:
        return tuple(absorbed)
    # Only the UV fraction interacts with ozone.  The remaining shortwave
    # radiation is transparent to the atmosphere and therefore contributes no
    # atmospheric temperature tendency in this isolated experiment.
    downward = UV_SHORTWAVE_FRACTION * SOLAR_CONSTANT * cos_zenith
    for k, optical_depth in enumerate(LAYER_OZONE_SHORTWAVE_OPTICAL_DEPTH):
        transmission = math.exp(-optical_depth)
        absorbed[k] = downward * (1.0 - transmission)
        downward *= transmission
    return tuple(absorbed)


def global_mean_shortwave_absorption(quadrature_points: int = 20_000) -> Vector:
    """Area-average instantaneous absorption over the illuminated hemisphere.

    The cosine of the solar zenith angle is uniformly distributed on a sphere,
    with area weight d(mu)/2.  It changes the incident flux, while ozone always
    uses its vertical optical depth.  The result is independent of season and
    axial tilt and has TOA mean incoming shortwave S0/4.
    """
    sums = [0.0] * (len(PRESSURE_HALF_PA) - 1)
    for index in range(quadrature_points):
        mu = (index + 0.5) / quadrature_points
        layer_absorption = shortwave_absorption_at_cos_zenith(mu)
        for k, value in enumerate(layer_absorption):
            sums[k] += value
    weight = 0.5 / quadrature_points
    return tuple(value * weight for value in sums)


def solar_cos_zenith(
    sin_latitude: float, longitude: float, time_seconds: float
) -> float:
    """Match dry_radiation:shortwave_downward_flux geometry."""
    orbital_longitude = (2.0 * math.pi * time_seconds / ORBITAL_PERIOD) % (
        2.0 * math.pi
    )
    solar_right_ascension = math.atan2(
        math.cos(AXIAL_TILT) * math.sin(orbital_longitude),
        math.cos(orbital_longitude),
    )
    sin_declination = math.sin(AXIAL_TILT) * math.sin(orbital_longitude)
    cos_declination = math.sqrt(max(0.0, 1.0 - sin_declination**2))
    cos_latitude = math.sqrt(max(0.0, 1.0 - sin_latitude**2))
    hour_angle = (
        PLANETARY_ROTATION_RATE * time_seconds
        + longitude
        - solar_right_ascension
    )
    return max(
        0.0,
        sin_latitude * sin_declination
        + cos_latitude * cos_declination * math.cos(hour_angle),
    )


def longwave_temperature_tendency(
    temperature: Sequence[float], *, include_ozone: bool = True
) -> Vector:
    """Compute longwave tendency with optional additive ozone opacity."""
    number_of_levels = len(temperature)
    total_pressure_depth = PRESSURE_HALF_PA[-1] - PRESSURE_HALF_PA[0]
    transmission = []
    emission = []
    for k in range(number_of_levels):
        pressure_thickness = PRESSURE_HALF_PA[k + 1] - PRESSURE_HALF_PA[k]
        non_ozone_optical_depth = (
            LONGWAVE_SURFACE_OPTICAL_DEPTH
            * pressure_thickness
            / total_pressure_depth
        )
        ozone_optical_depth = (
            LAYER_OZONE_LONGWAVE_OPTICAL_DEPTH[k] if include_ozone else 0.0
        )
        layer_transmission = math.exp(
            -(non_ozone_optical_depth + ozone_optical_depth)
        )
        transmission.append(layer_transmission)
        emission.append(
            (1.0 - layer_transmission)
            * STEFAN_BOLTZMANN
            * temperature[k] ** 4
        )

    downward = [0.0] * (number_of_levels + 1)
    for k in range(number_of_levels):
        downward[k + 1] = transmission[k] * downward[k] + emission[k]

    surface_temperature = REFERENCE_TEMPERATURE[-1]
    upward = [0.0] * (number_of_levels + 1)
    upward[number_of_levels] = STEFAN_BOLTZMANN * surface_temperature**4
    for k in range(number_of_levels - 1, -1, -1):
        upward[k] = transmission[k] * upward[k + 1] + emission[k]

    net = [up - down for up, down in zip(upward, downward)]
    tendency = []
    for k in range(number_of_levels):
        pressure_thickness = PRESSURE_HALF_PA[k + 1] - PRESSURE_HALF_PA[k]
        tendency.append(
            GRAVITY
            / (DRY_AIR_SPECIFIC_HEAT * pressure_thickness)
            * (net[k + 1] - net[k])
        )
    return tuple(tendency)


REFERENCE_TEMPERATURE = reference_profile()


def upper_radiative_tendency(
    upper_temperature: Vector, shortwave_absorption: Sequence[float]
) -> Vector:
    complete_temperature = upper_temperature + REFERENCE_TEMPERATURE[UPPER_LEVEL_COUNT:]
    longwave = longwave_temperature_tendency(complete_temperature)
    result = []
    for k in range(UPPER_LEVEL_COUNT):
        pressure_thickness = PRESSURE_HALF_PA[k + 1] - PRESSURE_HALF_PA[k]
        shortwave_tendency = (
            GRAVITY
            / (DRY_AIR_SPECIFIC_HEAT * pressure_thickness)
            * shortwave_absorption[k]
        )
        result.append(longwave[k] + shortwave_tendency)
    return tuple(result)


def add_scaled(left: Vector, scale: float, right: Vector) -> Vector:
    return tuple(a + scale * b for a, b in zip(left, right))


def raw_filter(previous: Vector, current: Vector, candidate: Vector) -> tuple[Vector, Vector]:
    change = tuple(
        0.5 * RAW_FILTER_EPSILON * (old - 2.0 * now + new)
        for old, now, new in zip(previous, current, candidate)
    )
    filtered_current = tuple(
        now + RAW_FILTER_ALPHA * delta for now, delta in zip(current, change)
    )
    next_state = tuple(
        new - (1.0 - RAW_FILTER_ALPHA) * delta
        for new, delta in zip(candidate, change)
    )
    return filtered_current, next_state


def leapfrog_start(initial: Vector, dt: float, tendency: TendencyFunction) -> tuple[Vector, Vector]:
    """Reproduce dry_atmosphere's two-stage first advance."""
    half = add_scaled(initial, 0.5 * dt, tendency(initial, 0.0))
    candidate = add_scaled(initial, dt, tendency(half, 0.5 * dt))
    _, next_state = raw_filter(initial, half, candidate)
    return initial, next_state


def leapfrog_step(
    previous: Vector,
    current: Vector,
    time_seconds: float,
    dt: float,
    tendency: TendencyFunction,
) -> tuple[Vector, Vector]:
    candidate = add_scaled(previous, 2.0 * dt, tendency(current, time_seconds))
    return raw_filter(previous, current, candidate)


def run_fixed_point(
    initial_scale: float,
    dt: float,
    max_years: float,
    tolerance_k_day: float,
    shortwave_absorption: Vector,
) -> ExperimentResult:
    initial = tuple(
        initial_scale * value for value in REFERENCE_TEMPERATURE[:UPPER_LEVEL_COUNT]
    )

    def tendency(state: Vector, _time_seconds: float) -> Vector:
        return upper_radiative_tendency(state, shortwave_absorption)

    previous, current = leapfrog_start(initial, dt, tendency)
    steps = 1
    maximum_steps = round(max_years * SECONDS_PER_YEAR / dt)
    sample_interval = max(1, round(30.0 * SECONDS_PER_DAY / dt))
    required_consecutive_samples = 6
    consecutive_samples = 0
    history: list[HistoryPoint] = []

    while steps < maximum_steps:
        time_seconds = steps * dt
        previous, current = leapfrog_step(
            previous, current, time_seconds, dt, tendency
        )
        steps += 1
        if steps % sample_interval == 0 or steps == maximum_steps:
            residual = max(abs(value) for value in tendency(current, steps * dt))
            residual_k_day = residual * SECONDS_PER_DAY
            history.append(
                HistoryPoint(
                    steps * dt / SECONDS_PER_YEAR,
                    current,
                    residual_k_day,
                )
            )
            if residual_k_day < tolerance_k_day:
                consecutive_samples += 1
            else:
                consecutive_samples = 0
            if consecutive_samples >= required_consecutive_samples:
                break

    final_tendency = tendency(current, steps * dt)
    maximum_tendency = max(abs(value) for value in final_tendency) * SECONDS_PER_DAY
    time_level_difference = max(abs(a - b) for a, b in zip(previous, current))
    return ExperimentResult(
        name=f"fixed_scale_{initial_scale:g}_dt_{dt:g}",
        forcing="global-mean constant shortwave",
        dt_seconds=dt,
        initial_scale=initial_scale,
        converged=maximum_tendency < tolerance_k_day,
        elapsed_years=steps * dt / SECONDS_PER_YEAR,
        steps=steps,
        temperatures=current,
        max_tendency_k_day=maximum_tendency,
        time_level_difference_k=time_level_difference,
        cycle_mismatch_k=None,
        history=tuple(history),
    )


def run_periodic_cycle(
    dt: float,
    max_years: int,
    cycle_tolerance_k: float,
) -> ExperimentResult:
    initial = REFERENCE_TEMPERATURE[:UPPER_LEVEL_COUNT]

    def tendency(state: Vector, time_seconds: float) -> Vector:
        mu = solar_cos_zenith(0.0, 0.0, time_seconds)
        return upper_radiative_tendency(
            state, shortwave_absorption_at_cos_zenith(mu)
        )

    steps_per_year = round(SECONDS_PER_YEAR / dt)
    if not math.isclose(steps_per_year * dt, SECONDS_PER_YEAR, abs_tol=1.0e-9):
        raise ValueError("periodic experiment dt must divide the 360-day year")

    previous, current = leapfrog_start(initial, dt, tendency)
    steps = 1
    maximum_steps = max_years * steps_per_year
    prior_year_state: tuple[Vector, Vector] | None = None
    cycle_mismatch = math.inf
    consecutive_years = 0
    history: list[HistoryPoint] = []

    while steps < maximum_steps:
        time_seconds = steps * dt
        previous, current = leapfrog_step(
            previous, current, time_seconds, dt, tendency
        )
        steps += 1
        if steps % steps_per_year == 0:
            residual = max(abs(value) for value in tendency(current, steps * dt))
            history.append(
                HistoryPoint(
                    steps * dt / SECONDS_PER_YEAR,
                    current,
                    residual * SECONDS_PER_DAY,
                )
            )
            if prior_year_state is not None:
                old_previous, old_current = prior_year_state
                cycle_mismatch = max(
                    max(abs(a - b) for a, b in zip(previous, old_previous)),
                    max(abs(a - b) for a, b in zip(current, old_current)),
                )
                if cycle_mismatch < cycle_tolerance_k:
                    consecutive_years += 1
                else:
                    consecutive_years = 0
                if consecutive_years >= 2:
                    break
            prior_year_state = (previous, current)

    final_tendency = tendency(current, steps * dt)
    return ExperimentResult(
        name=f"periodic_equator_dt_{dt:g}",
        forcing="equatorial diurnal/seasonal shortwave",
        dt_seconds=dt,
        initial_scale=1.0,
        converged=cycle_mismatch < cycle_tolerance_k,
        elapsed_years=steps * dt / SECONDS_PER_YEAR,
        steps=steps,
        temperatures=current,
        max_tendency_k_day=max(abs(value) for value in final_tendency)
        * SECONDS_PER_DAY,
        time_level_difference_k=max(
            abs(a - b) for a, b in zip(previous, current)
        ),
        cycle_mismatch_k=cycle_mismatch,
        history=tuple(history),
    )


def run_self_checks(global_absorption: Vector) -> None:
    if PRESSURE_HALF_PA[: UPPER_LEVEL_COUNT + 1] != (
        100.0,
        300.0,
        1_000.0,
        5_000.0,
        10_000.0,
    ):
        raise AssertionError("upper-layer pressure selection changed")
    if not math.isclose(
        sum(LAYER_OZONE_SHORTWAVE_OPTICAL_DEPTH),
        OZONE_SHORTWAVE_OPTICAL_DEPTH,
        rel_tol=0.0,
        abs_tol=1.0e-14,
    ):
        raise AssertionError("ozone optical depth does not normalize to 1.5")
    if not math.isclose(
        sum(LAYER_OZONE_LONGWAVE_OPTICAL_DEPTH),
        OZONE_LONGWAVE_OPTICAL_DEPTH,
        rel_tol=0.0,
        abs_tol=1.0e-16,
    ):
        raise AssertionError("ozone longwave optical depth does not normalize to 0.005")
    if LAYER_OZONE_SHORTWAVE_OPTICAL_DEPTH[0] <= 0.0:
        raise AssertionError("the 1-3 hPa layer must absorb UV")

    global_mean_uv = UV_SHORTWAVE_FRACTION * SOLAR_CONSTANT / 4.0
    global_transmitted = global_mean_uv * math.exp(
        -OZONE_SHORTWAVE_OPTICAL_DEPTH
    )
    if not math.isclose(
        sum(global_absorption) + global_transmitted,
        global_mean_uv,
        rel_tol=0.0,
        abs_tol=2.0e-7,
    ):
        raise AssertionError("global-mean UV column does not conserve energy")
    if any(value <= 0.0 for value in REFERENCE_TEMPERATURE):
        raise AssertionError("reference temperature contains a nonpositive value")


def layer_midpoint_hpa(k: int) -> float:
    return math.sqrt(PRESSURE_HALF_PA[k] * PRESSURE_HALF_PA[k + 1]) / 100.0


def write_summary_csv(output: Path, results: Sequence[ExperimentResult]) -> None:
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "experiment",
                "forcing",
                "uv_shortwave_fraction",
                "ozone_uv_optical_depth",
                "ozone_longwave_optical_depth",
                "ozone_lower_bound_hpa",
                "ozone_optical_path",
                "dt_seconds",
                "initial_scale",
                "converged",
                "elapsed_years",
                "steps",
                "max_tendency_k_day",
                "time_level_difference_k",
                "cycle_mismatch_k",
                *[f"temperature_level_{k + 1}_k" for k in range(UPPER_LEVEL_COUNT)],
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.name,
                    result.forcing,
                    f"{UV_SHORTWAVE_FRACTION:.12g}",
                    f"{OZONE_SHORTWAVE_OPTICAL_DEPTH:.12g}",
                    f"{OZONE_LONGWAVE_OPTICAL_DEPTH:.12g}",
                    f"{OZONE_PRESSURE_LOWER_BOUND / 100.0:.12g}",
                    "vertical",
                    f"{result.dt_seconds:.12g}",
                    f"{result.initial_scale:.12g}",
                    str(result.converged).lower(),
                    f"{result.elapsed_years:.12g}",
                    result.steps,
                    f"{result.max_tendency_k_day:.12g}",
                    f"{result.time_level_difference_k:.12g}",
                    ""
                    if result.cycle_mismatch_k is None
                    else f"{result.cycle_mismatch_k:.12g}",
                    *[f"{value:.12g}" for value in result.temperatures],
                ]
            )


def write_profile_csv(
    output: Path,
    global_absorption: Vector,
    results: Sequence[ExperimentResult],
) -> None:
    fixed_results = [result for result in results if result.cycle_mismatch_k is None]
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "level",
                "pressure_top_hpa",
                "pressure_bottom_hpa",
                "pressure_midpoint_hpa",
                "reference_temperature_k",
                "ozone_uv_optical_depth",
                "ozone_longwave_optical_depth",
                "global_mean_shortwave_absorption_w_m-2",
                *[f"{result.name}_temperature_k" for result in fixed_results],
            ]
        )
        for k in range(UPPER_LEVEL_COUNT):
            writer.writerow(
                [
                    k + 1,
                    PRESSURE_HALF_PA[k] / 100.0,
                    PRESSURE_HALF_PA[k + 1] / 100.0,
                    layer_midpoint_hpa(k),
                    f"{REFERENCE_TEMPERATURE[k]:.12g}",
                    f"{LAYER_OZONE_SHORTWAVE_OPTICAL_DEPTH[k]:.12g}",
                    f"{LAYER_OZONE_LONGWAVE_OPTICAL_DEPTH[k]:.12g}",
                    f"{global_absorption[k]:.12g}",
                    *[f"{result.temperatures[k]:.12g}" for result in fixed_results],
                ]
            )


def write_history_csv(output: Path, results: Sequence[ExperimentResult]) -> None:
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "experiment",
                "elapsed_years",
                "max_tendency_k_day",
                *[f"temperature_level_{k + 1}_k" for k in range(UPPER_LEVEL_COUNT)],
            ]
        )
        for result in results:
            for point in result.history:
                writer.writerow(
                    [
                        result.name,
                        f"{point.years:.12g}",
                        f"{point.residual_k_day:.12g}",
                        *[f"{value:.12g}" for value in point.temperatures],
                    ]
                )


def write_budget_csv(
    output: Path,
    global_absorption: Vector,
    results: Sequence[ExperimentResult],
) -> None:
    """Write the final layer-by-layer ozone and non-ozone radiation budget.

    The ozone longwave term is the exact increment from adding ozone opacity to
    the otherwise unchanged longwave calculation.  This assigns the nonlinear
    interaction between the two opacities to the added ozone term.
    """
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "experiment",
                "level",
                "pressure_midpoint_hpa",
                "temperature_k",
                "uv_absorption_tendency_k_day",
                "non_ozone_longwave_tendency_k_day",
                "ozone_longwave_tendency_k_day",
                "net_radiative_tendency_k_day",
                "uv_absorption_w_m-2",
            ]
        )
        for result in results:
            complete_temperature = (
                result.temperatures + REFERENCE_TEMPERATURE[UPPER_LEVEL_COUNT:]
            )
            longwave_total = longwave_temperature_tendency(complete_temperature)
            longwave_non_ozone = longwave_temperature_tendency(
                complete_temperature, include_ozone=False
            )
            longwave_ozone = tuple(
                total - non_ozone
                for total, non_ozone in zip(
                    longwave_total, longwave_non_ozone
                )
            )
            if result.cycle_mismatch_k is None:
                shortwave_absorption = global_absorption
            else:
                mu = solar_cos_zenith(
                    0.0, 0.0, result.elapsed_years * SECONDS_PER_YEAR
                )
                shortwave_absorption = shortwave_absorption_at_cos_zenith(mu)
            for k in range(UPPER_LEVEL_COUNT):
                pressure_thickness = PRESSURE_HALF_PA[k + 1] - PRESSURE_HALF_PA[k]
                uv_absorption_tendency = (
                    GRAVITY
                    / (DRY_AIR_SPECIFIC_HEAT * pressure_thickness)
                    * shortwave_absorption[k]
                )
                writer.writerow(
                    [
                        result.name,
                        k + 1,
                        f"{layer_midpoint_hpa(k):.12g}",
                        f"{result.temperatures[k]:.12g}",
                        f"{uv_absorption_tendency * SECONDS_PER_DAY:.12g}",
                        f"{longwave_non_ozone[k] * SECONDS_PER_DAY:.12g}",
                        f"{longwave_ozone[k] * SECONDS_PER_DAY:.12g}",
                        f"{(longwave_total[k] + uv_absorption_tendency) * SECONDS_PER_DAY:.12g}",
                        f"{shortwave_absorption[k]:.12g}",
                    ]
                )


def svg_polyline(points: Iterable[tuple[float, float]], color: str) -> str:
    coordinates = " ".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return (
        f'<polyline points="{coordinates}" fill="none" stroke="{color}" '
        'stroke-width="2" stroke-linejoin="round"/>'
    )


def write_svg(output: Path, results: Sequence[ExperimentResult]) -> None:
    fixed_results = [
        result
        for result in results
        if result.cycle_mismatch_k is None and result.dt_seconds == 1200.0
    ]
    colors = ("#0072B2", "#009E73", "#D55E00", "#CC79A7")
    width, height = 900, 520
    left, right, top, bottom = 95, 845, 50, 445
    all_temperatures = [
        value
        for result in fixed_results
        for value in result.temperatures
    ] + list(REFERENCE_TEMPERATURE[:UPPER_LEVEL_COUNT])
    minimum_temperature = math.floor((min(all_temperatures) - 5.0) / 10.0) * 10.0
    maximum_temperature = math.ceil((max(all_temperatures) + 5.0) / 10.0) * 10.0
    log_pressure_min = math.log10(1.0)
    log_pressure_max = math.log10(100.0)

    def x_coordinate(temperature: float) -> float:
        return left + (temperature - minimum_temperature) / (
            maximum_temperature - minimum_temperature
        ) * (right - left)

    def y_coordinate(pressure_hpa: float) -> float:
        return top + (math.log10(pressure_hpa) - log_pressure_min) / (
            log_pressure_max - log_pressure_min
        ) * (bottom - top)

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="450" y="26" text-anchor="middle" font-family="sans-serif" font-size="18">Vertical path: UV tau={OZONE_SHORTWAVE_OPTICAL_DEPTH:g}, LW ozone tau={OZONE_LONGWAVE_OPTICAL_DEPTH:g}</text>',
        f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="black"/>',
    ]
    for temperature in range(
        int(minimum_temperature), int(maximum_temperature) + 1, 20
    ):
        x = x_coordinate(float(temperature))
        lines.extend(
            [
                f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#dddddd"/>',
                f'<text x="{x:.2f}" y="{bottom + 22}" text-anchor="middle" font-family="sans-serif" font-size="12">{temperature}</text>',
            ]
        )
    for pressure in (1.0, 3.0, 10.0, 30.0, 100.0):
        y = y_coordinate(pressure)
        lines.extend(
            [
                f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#dddddd"/>',
                f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="12">{pressure:g}</text>',
            ]
        )
    pressures = [layer_midpoint_hpa(k) for k in range(UPPER_LEVEL_COUNT)]
    reference_points = [
        (x_coordinate(value), y_coordinate(pressure))
        for value, pressure in zip(
            REFERENCE_TEMPERATURE[:UPPER_LEVEL_COUNT], pressures
        )
    ]
    lines.append(svg_polyline(reference_points, "#555555"))
    lines.append(
        '<text x="690" y="72" font-family="sans-serif" font-size="12" fill="#555555">reference</text>'
    )
    for index, result in enumerate(fixed_results):
        color = colors[index % len(colors)]
        points = [
            (x_coordinate(value), y_coordinate(pressure))
            for value, pressure in zip(result.temperatures, pressures)
        ]
        lines.append(svg_polyline(points, color))
        lines.append(
            f'<text x="690" y="{90 + 18 * index}" font-family="sans-serif" font-size="12" fill="{color}">initial x{result.initial_scale:g}</text>'
        )
    lines.extend(
        [
            f'<text x="{0.5 * (left + right):.1f}" y="490" text-anchor="middle" font-family="sans-serif" font-size="14">temperature [K]</text>',
            '<text x="22" y="250" text-anchor="middle" transform="rotate(-90 22 250)" font-family="sans-serif" font-size="14">pressure [hPa]</text>',
            '</svg>',
        ]
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def print_result(result: ExperimentResult) -> None:
    state = "PASS" if result.converged else "FAIL"
    temperatures = ", ".join(f"{value:.6f}" for value in result.temperatures)
    print(
        f"{state:4s}  {result.name:28s}  years={result.elapsed_years:8.3f}  "
        f"max|dT/dt|={result.max_tendency_k_day:.3e} K/day  "
        f"time-level diff={result.time_level_difference_k:.3e} K"
    )
    if result.cycle_mismatch_k is not None:
        print(f"      annual-cycle mismatch={result.cycle_mismatch_k:.3e} K")
    print(f"      T[1:4] = [{temperatures}] K")


def print_fixed_point_budget(
    result: ExperimentResult, global_absorption: Vector
) -> None:
    """Print the requested equilibrium budget in K/day."""
    complete_temperature = (
        result.temperatures + REFERENCE_TEMPERATURE[UPPER_LEVEL_COUNT:]
    )
    longwave_total = longwave_temperature_tendency(complete_temperature)
    longwave_non_ozone = longwave_temperature_tendency(
        complete_temperature, include_ozone=False
    )
    print("\nEquilibrium radiative budget [K/day]")
    print("level  pressure     T [K]   UV absorption   non-O3 LW       O3 LW          net")
    for k in range(UPPER_LEVEL_COUNT):
        pressure_thickness = PRESSURE_HALF_PA[k + 1] - PRESSURE_HALF_PA[k]
        uv_absorption = (
            GRAVITY
            / (DRY_AIR_SPECIFIC_HEAT * pressure_thickness)
            * global_absorption[k]
            * SECONDS_PER_DAY
        )
        non_ozone = longwave_non_ozone[k] * SECONDS_PER_DAY
        ozone = (longwave_total[k] - longwave_non_ozone[k]) * SECONDS_PER_DAY
        net = uv_absorption + non_ozone + ozone
        print(
            f"{k + 1:5d}  {PRESSURE_HALF_PA[k] / 100.0:3.0f}-{PRESSURE_HALF_PA[k + 1] / 100.0:<3.0f} hPa"
            f"  {result.temperatures[k]:8.3f}  {uv_absorption:+13.6f}"
            f"  {non_ozone:+13.6f}  {ozone:+13.6f}  {net:+12.3e}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "output"
        / "upper_radiative_equilibrium_uv02_tau15_vertical_pa1hpa_lw_tau0005",
    )
    parser.add_argument("--max-years", type=float, default=40.0)
    parser.add_argument("--periodic-max-years", type=int, default=40)
    parser.add_argument("--equilibrium-tolerance-k-day", type=float, default=1.0e-4)
    parser.add_argument("--cycle-tolerance-k", type=float, default=1.0e-3)
    parser.add_argument("--skip-periodic", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.max_years <= 0.0 or args.periodic_max_years <= 0:
        raise ValueError("experiment durations must be positive")
    if args.equilibrium_tolerance_k_day <= 0.0 or args.cycle_tolerance_k <= 0.0:
        raise ValueError("convergence tolerances must be positive")

    global_absorption = global_mean_shortwave_absorption()
    run_self_checks(global_absorption)

    print("Upper layers: 1-100 hPa (model levels 1-4)")
    print(f"Global-mean TOA shortwave: {SOLAR_CONSTANT / 4.0:.6f} W m-2")
    print(
        "Global-mean TOA UV shortwave: "
        f"{UV_SHORTWAVE_FRACTION * SOLAR_CONSTANT / 4.0:.6f} W m-2 "
        f"({100.0 * UV_SHORTWAVE_FRACTION:g}% of shortwave)"
    )
    print(f"Ozone UV optical depth: {OZONE_SHORTWAVE_OPTICAL_DEPTH:g}")
    print(f"Ozone longwave optical depth: {OZONE_LONGWAVE_OPTICAL_DEPTH:g}")
    print("Ozone optical path: vertical (no slant-path correction)")
    print(
        "Ozone absorption lower bound: "
        f"{OZONE_PRESSURE_LOWER_BOUND / 100.0:g} hPa"
    )
    print(
        "Global-mean upper-atmosphere shortwave absorption: "
        f"{sum(global_absorption[:UPPER_LEVEL_COUNT]):.6f} W m-2"
    )
    print("\nFixed-point experiments")
    results: list[ExperimentResult] = []
    for initial_scale in (0.8, 1.0, 1.2):
        result = run_fixed_point(
            initial_scale,
            1200.0,
            args.max_years,
            args.equilibrium_tolerance_k_day,
            global_absorption,
        )
        results.append(result)
        print_result(result)

    print("\nTime-step sensitivity")
    sensitivity = run_fixed_point(
        1.0,
        600.0,
        args.max_years,
        args.equilibrium_tolerance_k_day,
        global_absorption,
    )
    results.append(sensitivity)
    print_result(sensitivity)

    if not args.skip_periodic:
        print("\nPeriodic-forcing experiment")
        periodic = run_periodic_cycle(
            1200.0,
            args.periodic_max_years,
            args.cycle_tolerance_k,
        )
        results.append(periodic)
        print_result(periodic)

    reference_fixed = next(
        result
        for result in results
        if result.cycle_mismatch_k is None
        and result.dt_seconds == 1200.0
        and result.initial_scale == 1.0
    )
    print_fixed_point_budget(reference_fixed, global_absorption)

    output_directory = args.output_dir
    output_directory.mkdir(parents=True, exist_ok=True)
    write_summary_csv(output_directory / "summary.csv", results)
    write_profile_csv(
        output_directory / "equilibrium_profiles.csv", global_absorption, results
    )
    write_history_csv(output_directory / "history.csv", results)
    write_budget_csv(
        output_directory / "radiative_budget.csv", global_absorption, results
    )
    write_svg(output_directory / "equilibrium_profiles.svg", results)
    print(f"\nResults written to {output_directory}")

    fixed_results = [result for result in results if result.cycle_mismatch_k is None]
    initial_spread = max(
        max(abs(a - b) for a, b in zip(result.temperatures, reference_fixed.temperatures))
        for result in fixed_results
        if result.dt_seconds == 1200.0
    )
    dt_spread = max(
        abs(a - b)
        for a, b in zip(reference_fixed.temperatures, sensitivity.temperatures)
    )
    print(f"Maximum equilibrium spread across initial states: {initial_spread:.3e} K")
    print(f"Maximum equilibrium difference between dt=1200 and 600 s: {dt_spread:.3e} K")

    return 0 if all(result.converged for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
