#!/usr/bin/env python3
"""Dependency-free API server for BespokePlanet simulation output."""

from __future__ import annotations

import argparse
import json
import math
import mimetypes
import re
import sys
import threading
from array import array
from dataclasses import dataclass
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, unquote, urlparse


FRAME_RE = re.compile(r"^(zeta|delta|eta|u|v)_(\d{5})\.bin$")
DRY_LEVEL_FRAME_RE = re.compile(
    r"^(zeta|delta|temperature|u|v)_l(\d{2})_(\d{5})\.bin$"
)
DRY_SURFACE_FRAME_RE = re.compile(r"^surface_pressure_(\d{5})\.bin$")
RUN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
FIELD_DEFINITIONS: dict[str, dict[str, Any]] = {
    "zeta": {"label": "zeta", "symbol": "ζ", "unit": "s⁻¹", "signed": True},
    "delta": {"label": "delta", "symbol": "δ", "unit": "s⁻¹", "signed": True},
    "eta": {"label": "eta", "symbol": "η", "unit": "m", "signed": True},
    "temperature": {
        "label": "temperature",
        "symbol": "T",
        "unit": "K",
        "signed": False,
    },
    "surface_pressure": {
        "label": "surface pressure",
        "symbol": "pₛ",
        "unit": "Pa",
        "signed": False,
    },
    "u": {
        "label": "eastward wind u",
        "symbol": "u",
        "unit": "m s⁻¹",
        "signed": True,
    },
    "v": {
        "label": "northward wind v",
        "symbol": "v",
        "unit": "m s⁻¹",
        "signed": True,
    },
    "speed": {
        "label": "sqrt(u² + v²)",
        "symbol": "|u|",
        "unit": "m s⁻¹",
        "signed": False,
        "zero_based": True,
    },
}


class DataError(RuntimeError):
    """Raised when a run contains malformed or incomplete data."""


@dataclass(frozen=True)
class Run:
    name: str
    path: Path
    metadata: dict[str, Any]
    steps: tuple[int, ...]
    fields: tuple[str, ...]
    data_generation: int
    is_dry: bool = False
    level_count: int = 0


def _full_level_pressure(top: float, bottom: float) -> float:
    """Pressure represented by a layer-centred prognostic value."""
    if not (math.isfinite(top) and math.isfinite(bottom) and 0.0 < top < bottom):
        raise DataError("Dry-atmosphere half-level pressures must increase downward")
    thickness = bottom - top
    alpha = 1.0 - top * math.log(bottom / top) / thickness
    return bottom * math.exp(-alpha)


def _normalized_metadata(raw: dict[str, Any]) -> dict[str, Any]:
    """Normalize legacy dry metadata to the schema used by the visualizer."""
    metadata = dict(raw)
    grid = dict(metadata["grid"])
    nlon = [int(value) for value in grid["nlon"]]
    offsets = [0]
    for ring_size in nlon:
        offsets.append(offsets[-1] + ring_size)
    grid.setdefault("ring_offsets", offsets)
    grid.setdefault("point_count", offsets[-1])
    metadata["grid"] = grid

    equation = metadata.get("equation", "barotropic_vorticity")
    if equation != "dry_hydrostatic_atmosphere":
        metadata.setdefault("supports_conservation_diagnostics", True)
        return metadata

    metadata.setdefault(
        "simulation",
        {
            "duration_seconds": metadata.get("duration_seconds"),
            "time_step_seconds": metadata.get("time_step_seconds"),
            "number_of_steps": metadata.get("number_of_steps"),
            "snapshot_interval_steps": metadata.get("snapshot_interval_steps"),
            "number_of_snapshots": metadata.get("number_of_snapshots"),
        },
    )
    metadata.setdefault(
        "numerics",
        {
            "spectral_truncation": metadata.get("spectral_truncation"),
            "maximum_cfl": metadata.get("maximum_advective_cfl"),
        },
    )
    half_levels = [float(value) for value in metadata["reference_half_level_pressure_pa"]]
    level_count = int(metadata.get("number_of_levels", len(half_levels) - 1))
    if len(half_levels) != level_count + 1:
        raise DataError("Dry-atmosphere metadata has an invalid vertical level count")
    full_levels = [
        _full_level_pressure(half_levels[k], half_levels[k + 1])
        for k in range(level_count)
    ]
    default_level = min(
        range(1, level_count + 1),
        key=lambda level: abs(math.log(full_levels[level - 1] / 50000.0)),
    )
    vertical = dict(metadata.get("vertical_coordinate", {}))
    for obsolete_key in (
        "minimum_pressure_pa",
        "maximum_pressure_pa",
        "default_pressure_pa",
        "interpolation",
    ):
        vertical.pop(obsolete_key, None)
    vertical.update(
        {
            "type": "hybrid_sigma_pressure",
            "number_of_levels": level_count,
            "reference_half_level_pressure_pa": half_levels,
            "reference_full_level_pressure_pa": full_levels,
            "default_level": default_level,
            "level_order": "top_to_bottom",
        }
    )
    a_half = metadata.get("hybrid_a_half_pa")
    b_half = metadata.get("hybrid_b_half")
    if a_half is not None or b_half is not None:
        if a_half is None or b_half is None:
            raise DataError("Dry-atmosphere metadata must contain both hybrid A and B")
        vertical["a_half_pa"] = [float(value) for value in a_half]
        vertical["b_half"] = [float(value) for value in b_half]
        if len(vertical["a_half_pa"]) != level_count + 1 or len(vertical["b_half"]) != level_count + 1:
            raise DataError("Dry-atmosphere hybrid coefficients have an invalid length")
        vertical["pressure_is_column_dependent"] = True
    else:
        vertical["pressure_is_column_dependent"] = False
    metadata["vertical_coordinate"] = vertical
    metadata["supports_conservation_diagnostics"] = False
    return metadata


def _read_float64(path: Path, expected_count: int) -> array:
    expected_bytes = expected_count * 8
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise DataError(f"Could not read {path.name}: {exc}") from exc
    if len(raw) != expected_bytes:
        raise DataError(
            f"{path.name} is {len(raw)} bytes; expected {expected_bytes} bytes"
        )
    values = array("d")
    values.frombytes(raw)
    if sys.byteorder != "little":
        values.byteswap()
    return values


def _gauss_legendre_weights(nodes: list[float]) -> list[float]:
    """Recover Gaussian quadrature weights from the roots stored in metadata."""
    degree = len(nodes)
    weights: list[float] = []
    for x in nodes:
        p_nm2 = 1.0
        p_nm1 = x
        if degree == 1:
            p_n = p_nm1
            p_before = p_nm2
        else:
            for n in range(2, degree + 1):
                p_n = ((2 * n - 1) * x * p_nm1 - (n - 1) * p_nm2) / n
                p_nm2, p_nm1 = p_nm1, p_n
            p_before = p_nm2
        derivative = degree * (p_before - x * p_n) / (1.0 - x * x)
        weights.append(2.0 / ((1.0 - x * x) * derivative * derivative))
    return weights


class Repository:
    def __init__(self, output_root: Path) -> None:
        self.output_root = output_root.resolve()
        self._conservation_cache: dict[tuple[str, int], dict[str, Any]] = {}
        self._field_statistics_cache: dict[
            tuple[str, int, str, int | None], dict[str, float]
        ] = {}
        self._cache_lock = threading.Lock()

    def discover(self) -> list[Run]:
        if not self.output_root.is_dir():
            raise DataError(f"Output directory does not exist: {self.output_root}")
        runs: list[Run] = []
        for metadata_path in sorted(self.output_root.glob("*/metadata.json")):
            try:
                runs.append(self._load_run(metadata_path.parent.name))
            except (DataError, OSError, ValueError, KeyError, json.JSONDecodeError):
                continue
        return runs

    def get_run(self, name: str) -> Run:
        if not RUN_RE.fullmatch(name):
            raise KeyError(name)
        run = self._load_run(name)
        return run

    def _load_run(self, name: str) -> Run:
        run_path = (self.output_root / name).resolve()
        if run_path.parent != self.output_root or not run_path.is_dir():
            raise KeyError(name)
        metadata_path = run_path / "metadata.json"
        metadata = _normalized_metadata(
            json.loads(metadata_path.read_text(encoding="utf-8"))
        )
        grid = metadata["grid"]
        mu = grid["mu"]
        nlon = grid["nlon"]
        offsets = grid["ring_offsets"]
        if len(mu) != len(nlon) or len(offsets) != len(nlon) + 1:
            raise DataError(f"Grid arrays disagree in {metadata_path}")
        if offsets[0] != 0 or offsets[-1] != sum(nlon):
            raise DataError(f"Invalid ring offsets in {metadata_path}")
        if int(grid["point_count"]) != offsets[-1]:
            raise DataError(f"Invalid point count in {metadata_path}")

        is_dry = metadata.get("equation") == "dry_hydrostatic_atmosphere"
        level_count = int(metadata.get("vertical_coordinate", {}).get("number_of_levels", 0))
        if is_dry:
            return self._load_dry_run(
                name, run_path, metadata_path, metadata, level_count
            )

        step_sets = {field: set() for field in ("zeta", "delta", "eta", "u", "v")}
        for path in run_path.iterdir():
            match = FRAME_RE.match(path.name)
            if match:
                step_sets[match.group(1)].add(int(match.group(2)))

        fields: list[str] = []
        required_components: list[str] = []
        for field in ("zeta", "delta", "eta"):
            if step_sets[field]:
                fields.append(field)
                required_components.append(field)
        if step_sets["u"] and step_sets["v"]:
            fields.append("speed")
            required_components.extend(("u", "v"))
        if not fields:
            raise DataError(f"No supported fields found in {run_path}")

        steps = tuple(
            sorted(set.intersection(*(step_sets[field] for field in required_components)))
        )
        simulation = metadata.get("simulation", {})
        if "snapshot_interval_steps" in simulation:
            interval = int(simulation["snapshot_interval_steps"])
            final_step = int(simulation["number_of_steps"])
            if interval <= 0:
                raise DataError(f"Invalid snapshot interval in {metadata_path}")
            steps = tuple(
                step
                for step in steps
                if 0 <= step <= final_step
                and (step % interval == 0 or step == final_step)
            )
        if not steps:
            raise DataError(f"No complete frames found in {run_path}")
        return Run(
            name=name,
            path=run_path,
            metadata=metadata,
            steps=steps,
            fields=tuple(fields),
            data_generation=metadata_path.stat().st_mtime_ns,
        )

    @staticmethod
    def _load_dry_run(
        name: str,
        run_path: Path,
        metadata_path: Path,
        metadata: dict[str, Any],
        level_count: int,
    ) -> Run:
        if level_count <= 0 or level_count > 99:
            raise DataError(f"Invalid dry-atmosphere level count in {metadata_path}")
        level_steps = {
            field: {level: set() for level in range(1, level_count + 1)}
            for field in ("zeta", "delta", "temperature", "u", "v")
        }
        surface_steps: set[int] = set()
        for path in run_path.iterdir():
            level_match = DRY_LEVEL_FRAME_RE.match(path.name)
            if level_match:
                field, level_text, step_text = level_match.groups()
                level = int(level_text)
                if 1 <= level <= level_count:
                    level_steps[field][level].add(int(step_text))
                continue
            surface_match = DRY_SURFACE_FRAME_RE.match(path.name)
            if surface_match:
                surface_steps.add(int(surface_match.group(1)))

        def complete_level_steps(field: str) -> set[int]:
            return set.intersection(
                *(level_steps[field][level] for level in range(1, level_count + 1))
            )

        complete = {
            field: complete_level_steps(field)
            for field in ("zeta", "delta", "temperature", "u", "v")
        }
        fields: list[str] = []
        required_step_sets: list[set[int]] = []
        if surface_steps:
            fields.append("surface_pressure")
            required_step_sets.append(surface_steps)
        for field in ("temperature", "zeta", "delta", "u", "v"):
            if complete[field]:
                fields.append(field)
                required_step_sets.append(complete[field])
        if complete["u"] and complete["v"]:
            fields.append("speed")
        if not fields or not required_step_sets:
            raise DataError(f"No supported dry-atmosphere fields found in {run_path}")

        steps = tuple(sorted(set.intersection(*required_step_sets)))
        simulation = metadata["simulation"]
        interval = int(simulation["snapshot_interval_steps"])
        final_step = int(simulation["number_of_steps"])
        if interval <= 0:
            raise DataError(f"Invalid snapshot interval in {metadata_path}")
        steps = tuple(
            step
            for step in steps
            if 0 <= step <= final_step
            and (step % interval == 0 or step == final_step)
        )
        if not steps:
            raise DataError(f"No complete dry-atmosphere frames found in {run_path}")
        return Run(
            name=name,
            path=run_path,
            metadata=metadata,
            steps=steps,
            fields=tuple(fields),
            data_generation=metadata_path.stat().st_mtime_ns,
            is_dry=True,
            level_count=level_count,
        )

    @staticmethod
    def public_metadata(run: Run) -> dict[str, Any]:
        metadata = dict(run.metadata)
        metadata["available_steps"] = list(run.steps)
        metadata["available_frame_count"] = len(run.steps)
        metadata["data_generation"] = run.data_generation
        metadata["supports_streamlines"] = "speed" in run.fields
        metadata["available_fields"] = [
            {
                "id": field,
                **FIELD_DEFINITIONS[field],
                "uses_level": run.is_dry and field != "surface_pressure",
            }
            for field in run.fields
        ]
        return metadata

    def field_frame(
        self, run: Run, field: str, step: int, level: int | None = None
    ) -> tuple[bytes, dict[str, float]]:
        if step not in run.steps:
            raise KeyError(step)
        if field not in run.fields:
            raise KeyError(field)
        values = self._field_values(run, field, step, level)
        count = int(run.metadata["grid"]["point_count"])
        finite = [value for value in values if math.isfinite(value)]
        if len(finite) != count:
            raise DataError(f"Frame {step} contains non-finite values in {field}")
        ordered = sorted(finite)
        p005 = ordered[min(count - 1, max(0, math.floor(0.005 * count)))]
        p98 = ordered[min(count - 1, int(0.98 * (count - 1)))]
        p995 = ordered[min(count - 1, max(0, math.ceil(0.995 * count) - 1))]
        absolute_ordered = sorted(abs(value) for value in finite)
        p98_absolute = absolute_ordered[min(count - 1, int(0.98 * (count - 1)))]
        p995_absolute = absolute_ordered[
            min(count - 1, max(0, math.ceil(0.995 * count) - 1))
        ]
        stats = {
            "minimum": ordered[0],
            "maximum": ordered[-1],
            "p005": p005,
            "p98": p98,
            "p995": p995,
            "maximum_absolute": max(abs(ordered[0]), abs(ordered[-1])),
            "p98_absolute": p98_absolute,
            "p995_absolute": p995_absolute,
            "mean": math.fsum(ordered) / count,
        }
        if sys.byteorder != "little":
            values.byteswap()
        return values.tobytes(), stats

    @staticmethod
    def _field_values(
        run: Run, field: str, step: int, level: int | None = None
    ) -> array:
        count = int(run.metadata["grid"]["point_count"])
        if run.is_dry:
            if field == "surface_pressure":
                return array(
                    "f",
                    _read_float64(
                        run.path / f"surface_pressure_{step:05d}.bin", count
                    ),
                )
            selected_level = int(
                run.metadata["vertical_coordinate"]["default_level"]
                if level is None
                else level
            )
            if not 1 <= selected_level <= run.level_count:
                raise DataError(
                    f"Level must be between 1 and {run.level_count} for this run"
                )
            if field == "speed":
                east = _read_float64(
                    run.path / f"u_l{selected_level:02d}_{step:05d}.bin", count
                )
                north = _read_float64(
                    run.path / f"v_l{selected_level:02d}_{step:05d}.bin", count
                )
                return array(
                    "f", (math.hypot(u_value, v_value) for u_value, v_value in zip(east, north))
                )
            source = _read_float64(
                run.path / f"{field}_l{selected_level:02d}_{step:05d}.bin", count
            )
            return array("f", source)
        if field == "speed":
            u = _read_float64(run.path / f"u_{step:05d}.bin", count)
            v = _read_float64(run.path / f"v_{step:05d}.bin", count)
            return array("f", (math.hypot(east, north) for east, north in zip(u, v)))
        source = _read_float64(run.path / f"{field}_{step:05d}.bin", count)
        return array("f", source)

    def speed_frame(self, run: Run, step: int) -> tuple[bytes, dict[str, float]]:
        """Backward-compatible alias for clients using the original speed API."""
        return self.field_frame(run, "speed", step)

    def wind_frame(
        self, run: Run, step: int, level: int | None = None
    ) -> tuple[bytes, float]:
        """Return u followed by v as two contiguous float32 grid arrays."""
        if step not in run.steps:
            raise KeyError(step)
        if "speed" not in run.fields:
            raise KeyError("wind")
        eastward = self._field_values(run, "u", step, level)
        northward = self._field_values(run, "v", step, level)
        maximum_speed = 0.0
        for east, north in zip(eastward, northward):
            speed = math.hypot(east, north)
            if not math.isfinite(speed):
                raise DataError(f"Frame {step} contains non-finite wind values")
            maximum_speed = max(maximum_speed, speed)
        if sys.byteorder != "little":
            eastward.byteswap()
            northward.byteswap()
        return eastward.tobytes() + northward.tobytes(), maximum_speed

    def field_statistics(
        self, run: Run, field: str, level: int | None = None
    ) -> dict[str, float]:
        if field not in run.fields:
            raise KeyError(field)
        level_key = None if level is None else int(level)
        key = (run.name, run.data_generation, field, level_key)
        with self._cache_lock:
            cached = self._field_statistics_cache.get(key)
        if cached is not None:
            return cached

        minimum = math.inf
        maximum = -math.inf
        for step in run.steps:
            for value in self._field_values(run, field, step, level):
                if not math.isfinite(value):
                    raise DataError(f"Frame {step} contains non-finite values in {field}")
                minimum = min(minimum, value)
                maximum = max(maximum, value)
        result = {
            "minimum": minimum,
            "maximum": maximum,
            "maximum_absolute": max(abs(minimum), abs(maximum)),
        }
        with self._cache_lock:
            for stale_key in list(self._field_statistics_cache):
                if stale_key[0] == run.name and stale_key[1] != run.data_generation:
                    del self._field_statistics_cache[stale_key]
            self._field_statistics_cache[key] = result
        return result

    def conservation(self, run: Run) -> dict[str, Any]:
        if run.is_dry:
            raise DataError(
                "Conservation diagnostics are not yet defined for dry-atmosphere output"
            )
        cache_key = (run.name, run.data_generation)
        with self._cache_lock:
            cached = self._conservation_cache.get(cache_key)
        if cached is not None:
            return cached

        metadata = run.metadata
        grid = metadata["grid"]
        nlon = [int(value) for value in grid["nlon"]]
        mu = [float(value) for value in grid["mu"]]
        offsets = [int(value) for value in grid["ring_offsets"]]
        weights = _gauss_legendre_weights(mu)
        point_count = int(grid["point_count"])
        radius = float(metadata["physical_constants"]["earth_radius_m"])
        omega = float(metadata["physical_constants"]["rotation_rate_rad_s"])
        dt = float(metadata["simulation"]["time_step_seconds"])

        series = {
            "mean_relative_vorticity": [],
            "mean_kinetic_energy": [],
            "mean_absolute_enstrophy": [],
            "mean_relative_axial_angular_momentum": [],
        }
        times: list[float] = []
        for step in run.steps:
            zeta = _read_float64(run.path / f"zeta_{step:05d}.bin", point_count)
            u = _read_float64(run.path / f"u_{step:05d}.bin", point_count)
            v = _read_float64(run.path / f"v_{step:05d}.bin", point_count)
            circulation = energy = enstrophy = angular_momentum = 0.0
            for j, ring_size in enumerate(nlon):
                start, end = offsets[j], offsets[j + 1]
                point_weight = weights[j] / (2.0 * ring_size)
                coriolis = 2.0 * omega * mu[j]
                cos_latitude = math.sqrt(max(0.0, 1.0 - mu[j] * mu[j]))
                zeta_sum = energy_sum = enstrophy_sum = momentum_sum = 0.0
                for index in range(start, end):
                    z = zeta[index]
                    east = u[index]
                    north = v[index]
                    zeta_sum += z
                    energy_sum += 0.5 * (east * east + north * north)
                    enstrophy_sum += 0.5 * (z + coriolis) ** 2
                    momentum_sum += east * radius * cos_latitude
                circulation += point_weight * zeta_sum
                energy += point_weight * energy_sum
                enstrophy += point_weight * enstrophy_sum
                angular_momentum += point_weight * momentum_sum
            times.append(step * dt)
            series["mean_relative_vorticity"].append(circulation)
            series["mean_kinetic_energy"].append(energy)
            series["mean_absolute_enstrophy"].append(enstrophy)
            series["mean_relative_axial_angular_momentum"].append(angular_momentum)

        result = {
            "times_seconds": times,
            "steps": list(run.steps),
            "metrics": [
                {
                    "id": "mean_kinetic_energy",
                    "label": "平均運動エネルギー",
                    "unit": "m² s⁻²",
                    "values": series["mean_kinetic_energy"],
                },
                {
                    "id": "mean_absolute_enstrophy",
                    "label": "平均絶対エンストロフィー",
                    "unit": "s⁻²",
                    "values": series["mean_absolute_enstrophy"],
                },
                {
                    "id": "mean_relative_vorticity",
                    "label": "平均相対渦度",
                    "unit": "s⁻¹",
                    "values": series["mean_relative_vorticity"],
                },
                {
                    "id": "mean_relative_axial_angular_momentum",
                    "label": "平均相対軸角運動量",
                    "unit": "m² s⁻¹",
                    "values": series["mean_relative_axial_angular_momentum"],
                },
            ],
            "normalization": "global spherical mean",
        }
        with self._cache_lock:
            for stale_key in list(self._conservation_cache):
                if stale_key[0] == run.name and stale_key != cache_key:
                    del self._conservation_cache[stale_key]
            self._conservation_cache[cache_key] = result
        return result


class AppHandler(BaseHTTPRequestHandler):
    repository: Repository
    static_root: Path
    three_build_root: Path

    def log_message(self, fmt: str, *args: object) -> None:
        print(f"[{self.log_date_time_string()}] {fmt % args}")

    def do_GET(self) -> None:  # noqa: N802
        try:
            self._route_get()
        except KeyError:
            self._json_error(HTTPStatus.NOT_FOUND, "Requested run or frame was not found")
        except DataError as exc:
            self._json_error(HTTPStatus.UNPROCESSABLE_ENTITY, str(exc))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            self._json_error(HTTPStatus.INTERNAL_SERVER_ERROR, str(exc))

    def _route_get(self) -> None:
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        query = parse_qs(parsed.query)
        level: int | None = None
        if "level" in query:
            try:
                if len(query["level"]) != 1:
                    raise ValueError
                level = int(query["level"][0])
            except ValueError as exc:
                raise DataError("level must be one integer") from exc
        if path == "/api/runs":
            runs = self.repository.discover()
            self._send_json(
                {
                    "runs": [
                        {
                            "name": run.name,
                            "case_name": run.metadata.get("case_name", run.name),
                            "frame_count": len(run.steps),
                            "point_count": run.metadata["grid"]["point_count"],
                            "equation": run.metadata.get("equation", "barotropic_vorticity"),
                            "available_fields": list(run.fields),
                        }
                        for run in runs
                    ]
                }
            )
            return
        parts = [part for part in path.split("/") if part]
        if len(parts) >= 3 and parts[:2] == ["api", "runs"]:
            run = self.repository.get_run(parts[2])
            if len(parts) == 4 and parts[3] == "metadata":
                self._send_json(self.repository.public_metadata(run))
                return
            if len(parts) == 4 and parts[3] == "conservation":
                self._send_json(self.repository.conservation(run))
                return
            if len(parts) == 6 and parts[3] == "fields" and parts[5] == "statistics":
                self._send_json(
                    self.repository.field_statistics(run, parts[4], level)
                )
                return
            if len(parts) == 6 and parts[3] == "fields":
                field = parts[4]
                try:
                    step = int(parts[5])
                except ValueError as exc:
                    raise KeyError(parts[5]) from exc
                payload, stats = self.repository.field_frame(
                    run, field, step, level
                )
                self._send_field(payload, stats, level=level)
                return
            if len(parts) == 5 and parts[3] == "frame":
                try:
                    step = int(parts[4])
                except ValueError as exc:
                    raise KeyError(parts[4]) from exc
                payload, stats = self.repository.speed_frame(run, step)
                self._send_field(payload, stats, legacy_speed_headers=True)
                return
            if len(parts) == 5 and parts[3] == "wind":
                try:
                    step = int(parts[4])
                except ValueError as exc:
                    raise KeyError(parts[4]) from exc
                payload, maximum_speed = self.repository.wind_frame(run, step, level)
                self._send_wind(payload, maximum_speed, level)
                return
        if path.startswith("/vendor/"):
            vendor_name = path.removeprefix("/vendor/")
            vendor_path = (self.three_build_root / vendor_name).resolve()
            if (
                vendor_path.parent != self.three_build_root.resolve()
                or not vendor_name.startswith("three.")
                or not vendor_name.endswith(".js")
            ):
                raise KeyError(vendor_name)
            if not vendor_path.is_file():
                raise DataError("Three.js is not installed. Run: cd viz && npm install")
            self._send_file(vendor_path, "text/javascript")
            return
        relative = "index.html" if path == "/" else path.lstrip("/")
        candidate = (self.static_root / relative).resolve()
        if candidate != self.static_root and self.static_root not in candidate.parents:
            raise KeyError(relative)
        if not candidate.is_file():
            raise KeyError(relative)
        self._send_file(candidate)

    def _send_file(self, path: Path, content_type: str | None = None) -> None:
        payload = path.read_bytes()
        mime = content_type or mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", f"{mime}; charset=utf-8" if mime.startswith("text/") else mime)
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def _send_field(
        self,
        payload: bytes,
        stats: dict[str, float],
        legacy_speed_headers: bool = False,
        level: int | None = None,
    ) -> None:
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        for name, value in stats.items():
            header = "-".join(part.capitalize() for part in name.split("_"))
            self.send_header(f"X-Field-{header}", f"{value:.9g}")
        if level is not None:
            self.send_header("X-Field-Level", str(level))
        if legacy_speed_headers:
            self.send_header("X-Speed-Min", f"{stats['minimum']:.9g}")
            self.send_header("X-Speed-Max", f"{stats['maximum']:.9g}")
            self.send_header("X-Speed-P98", f"{stats['p98']:.9g}")
            self.send_header("X-Speed-Mean", f"{stats['mean']:.9g}")
        self.end_headers()
        self.wfile.write(payload)

    def _send_wind(
        self, payload: bytes, maximum_speed: float, level: int | None
    ) -> None:
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Wind-Layout", "u-then-v-float32-le")
        self.send_header("X-Wind-Maximum-Speed", f"{maximum_speed:.9g}")
        if level is not None:
            self.send_header("X-Field-Level", str(level))
        self.end_headers()
        self.wfile.write(payload)

    def _send_json(self, value: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        payload = json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(payload)

    def _json_error(self, status: HTTPStatus, message: str) -> None:
        self._send_json({"error": message}, status)


def parse_args() -> argparse.Namespace:
    script_root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--output", type=Path, default=script_root.parent / "output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script_root = Path(__file__).resolve().parent
    handler = type(
        "ConfiguredAppHandler",
        (AppHandler,),
        {
            "repository": Repository(args.output),
            "static_root": script_root / "static",
            "three_build_root": script_root / "node_modules" / "three" / "build",
        },
    )
    server = ThreadingHTTPServer((args.host, args.port), handler)
    print(f"BespokePlanet visualizer: http://{args.host}:{args.port}")
    print(f"Reading output from: {args.output.resolve()}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
