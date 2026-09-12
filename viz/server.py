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
from urllib.parse import unquote, urlparse


FRAME_RE = re.compile(r"^(zeta|u|v)_(\d{5})\.bin$")
RUN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


class DataError(RuntimeError):
    """Raised when a run contains malformed or incomplete data."""


@dataclass(frozen=True)
class Run:
    name: str
    path: Path
    metadata: dict[str, Any]
    steps: tuple[int, ...]


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
        self._conservation_cache: dict[str, dict[str, Any]] = {}
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
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
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

        step_sets = {field: set() for field in ("zeta", "u", "v")}
        for path in run_path.iterdir():
            match = FRAME_RE.match(path.name)
            if match:
                step_sets[match.group(1)].add(int(match.group(2)))
        steps = tuple(sorted(set.intersection(*step_sets.values())))
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
        return Run(name=name, path=run_path, metadata=metadata, steps=steps)

    @staticmethod
    def public_metadata(run: Run) -> dict[str, Any]:
        metadata = dict(run.metadata)
        metadata["available_steps"] = list(run.steps)
        metadata["available_frame_count"] = len(run.steps)
        return metadata

    def speed_frame(self, run: Run, step: int) -> tuple[bytes, dict[str, float]]:
        if step not in run.steps:
            raise KeyError(step)
        count = int(run.metadata["grid"]["point_count"])
        u = _read_float64(run.path / f"u_{step:05d}.bin", count)
        v = _read_float64(run.path / f"v_{step:05d}.bin", count)
        speeds = array("f", (math.hypot(east, north) for east, north in zip(u, v)))
        finite = [value for value in speeds if math.isfinite(value)]
        if len(finite) != count:
            raise DataError(f"Frame {step} contains non-finite velocity")
        ordered = sorted(finite)
        p98 = ordered[min(count - 1, int(0.98 * (count - 1)))]
        stats = {
            "minimum": ordered[0],
            "maximum": ordered[-1],
            "p98": p98,
            "mean": math.fsum(ordered) / count,
        }
        if sys.byteorder != "little":
            speeds.byteswap()
        return speeds.tobytes(), stats

    def conservation(self, run: Run) -> dict[str, Any]:
        with self._cache_lock:
            cached = self._conservation_cache.get(run.name)
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
            self._conservation_cache[run.name] = result
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
        path = unquote(urlparse(self.path).path)
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
            if len(parts) == 5 and parts[3] == "frame":
                try:
                    step = int(parts[4])
                except ValueError as exc:
                    raise KeyError(parts[4]) from exc
                payload, stats = self.repository.speed_frame(run, step)
                self.send_response(HTTPStatus.OK)
                self.send_header("Content-Type", "application/octet-stream")
                self.send_header("Content-Length", str(len(payload)))
                self.send_header("Cache-Control", "private, max-age=3600")
                self.send_header("X-Speed-Min", f"{stats['minimum']:.9g}")
                self.send_header("X-Speed-Max", f"{stats['maximum']:.9g}")
                self.send_header("X-Speed-P98", f"{stats['p98']:.9g}")
                self.send_header("X-Speed-Mean", f"{stats['mean']:.9g}")
                self.end_headers()
                self.wfile.write(payload)
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
