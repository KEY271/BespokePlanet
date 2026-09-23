#!/usr/bin/env python3
"""Dependency-free data server for the BespokePlanet land-sea visualizer.

Only the current land-sea output (moist atmosphere with land/sea tiles, sea
ice and snow) is served. The server validates and forwards the raw float64
grids; every climate diagnostic is computed in the browser (static/climate.js)
so that the land threshold can be changed interactively.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import mimetypes
import re
import sys
import threading
from array import array
from dataclasses import dataclass, field
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, unquote, urlparse


RUN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
FIELD_RE = re.compile(r"^[a-z][a-z0-9_]*$")
MONTHLY_RE = re.compile(r"^monthly_([a-z0-9_]+)_m(\d{4})\.bin$")
YEARLY_SURFACE_RE = re.compile(r"^yearly_([a-z0-9_]+)_y(\d{4})\.bin$")
YEARLY_LEVEL_RE = re.compile(r"^yearly_([a-z0-9_]+)_y(\d{4})_l(\d{2})\.bin$")
YEARLY_TIME_RE = re.compile(r"^yearly_time_y(\d{4})\.json$")
ZONAL_PREFIXES = ("zonal_", "eddy_")
STATIC_FIELDS = ("land_fraction", "surface_height", "ocean_q_flux")
REQUIRED_STATIC = ("land_fraction", "surface_height")
# Monthly fields read by the final-year analysis; a run without all of them is
# not the current land-sea output and is not offered.
REQUIRED_MONTHLY = (
    "surface_temperature", "land_temperature", "surface_air_temperature", "ocean_temperature",
    "precipitation", "evaporation", "cloud_cover", "surface_pressure",
    "surface_water", "sea_ice_fraction", "sea_ice_volume",
    "sea_ice_thickness", "sea_ice_temperature", "snow_fraction",
    "snow_water", "zonal_v",
)
MAXIMUM_INDICES = 1000


class DataError(RuntimeError):
    """Raised when a run contains malformed or incomplete data."""


@dataclass
class Run:
    name: str
    path: Path
    metadata: dict[str, Any]
    point_count: int
    nlat: int
    level_count: int
    months: tuple[int, ...]
    monthly: dict[str, set[int]]
    yearly_surface: dict[str, set[int]]
    yearly_level: dict[str, set[int]]
    years: tuple[int, ...]
    static: tuple[str, ...]
    generation: int
    year_times: dict[int, float] = field(default_factory=dict)

    @property
    def terrain(self) -> str:
        source = str(self.metadata.get("topography", {}).get("source", ""))
        return "earth" if "ETOPO" in source else "analytic"


def _grid_shape(metadata: dict[str, Any]) -> tuple[int, int]:
    grid = metadata["grid"]
    mu = [float(value) for value in grid["mu"]]
    nlon = [int(value) for value in grid["nlon"]]
    if not mu or len(mu) != len(nlon) or any(n <= 0 for n in nlon):
        raise DataError("格子の配列の長さが合いません")
    if any(not -1.0 < value < 1.0 for value in mu) or any(b <= a for a, b in zip(mu, mu[1:])):
        raise DataError("格子の mu が南から北へ増加していません")
    offsets = grid.get("ring_offsets")
    if offsets is not None:
        expected = [0]
        for n in nlon:
            expected.append(expected[-1] + n)
        if [int(value) for value in offsets] != expected:
            raise DataError("緯度リングのオフセットが不正です")
    return sum(nlon), len(nlon)


def _is_land_sea(metadata: dict[str, Any]) -> bool:
    return (
        int(metadata.get("schema_version", 0)) >= 3
        and metadata.get("equation") == "moist_hydrostatic_atmosphere"
        and bool(metadata.get("sea_ice", {}).get("enabled"))
        and bool(metadata.get("snow", {}).get("enabled"))
        and "static_land_fraction" in metadata.get("output", {})
    )


def _read_float64(path: Path, expected_count: int) -> array:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise DataError(f"{path.name} を読めません: {exc}") from exc
    if len(raw) != expected_count * 8:
        raise DataError(f"{path.name} は {len(raw)} バイトですが、{expected_count * 8} バイトのはずです")
    values = array("d")
    values.frombytes(raw)
    if sys.byteorder != "little":
        values.byteswap()
    if not all(map(math.isfinite, values)):
        raise DataError(f"{path.name} に有限でない値があります")
    return values


def parse_indices(text: str) -> list[int]:
    """Parse "1,2,5-8" into a list of positive integers, keeping the order."""
    indices: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            raise DataError("index が空です")
        if "-" in part:
            first, _, last = part.partition("-")
            if not (first.isdigit() and last.isdigit()) or int(last) < int(first):
                raise DataError(f"index の範囲が不正です: {part}")
            indices.extend(range(int(first), int(last) + 1))
        elif part.isdigit():
            indices.append(int(part))
        else:
            raise DataError(f"index が不正です: {part}")
        if len(indices) > MAXIMUM_INDICES:
            raise DataError("index が多すぎます")
    if any(index <= 0 for index in indices):
        raise DataError("index は 1 から始まります")
    return indices


class Repository:
    def __init__(self, output_root: Path) -> None:
        self.output_root = output_root.resolve()
        self._cache: dict[str, tuple[tuple[int, int], Run]] = {}
        self._lock = threading.Lock()

    def discover(self) -> list[Run]:
        if not self.output_root.is_dir():
            raise DataError(f"出力ディレクトリがありません: {self.output_root}")
        runs = []
        for metadata_path in sorted(self.output_root.glob("*/metadata.json")):
            try:
                runs.append(self.get_run(metadata_path.parent.name))
            except (DataError, KeyError, OSError, ValueError, TypeError, json.JSONDecodeError):
                continue
        return runs

    def get_run(self, name: str) -> Run:
        if not RUN_RE.fullmatch(name):
            raise KeyError(name)
        path = (self.output_root / name).resolve()
        if path.parent != self.output_root or not (path / "metadata.json").is_file():
            raise KeyError(name)
        stamp = ((path / "metadata.json").stat().st_mtime_ns, path.stat().st_mtime_ns)
        with self._lock:
            cached = self._cache.get(name)
        if cached and cached[0] == stamp:
            return cached[1]
        run = self._load_run(name, path)
        with self._lock:
            self._cache[name] = (stamp, run)
        return run

    @staticmethod
    def _load_run(name: str, path: Path) -> Run:
        metadata_path = path / "metadata.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not _is_land_sea(metadata):
            raise DataError(f"{name} は最新形式の陸海ケースではありません")
        point_count, nlat = _grid_shape(metadata)
        a_half = metadata.get("hybrid_a_half_pa")
        b_half = metadata.get("hybrid_b_half")
        half_pressure = metadata.get("reference_half_level_pressure_pa")
        if not a_half or not b_half or not half_pressure or not len(a_half) == len(b_half) == len(half_pressure):
            raise DataError("hybrid 座標の係数がありません")
        level_count = len(a_half) - 1
        for key in ("days_per_month", "months_per_year"):
            if int(metadata["calendar"][key]) <= 0:
                raise DataError("暦の設定が不正です")

        monthly: dict[str, set[int]] = {}
        yearly_surface: dict[str, set[int]] = {}
        yearly_level_steps: dict[str, dict[int, set[int]]] = {}
        year_times: dict[int, float] = {}
        generation = metadata_path.stat().st_mtime_ns
        for entry in path.iterdir():
            filename = entry.name
            if match := MONTHLY_RE.fullmatch(filename):
                monthly.setdefault(match[1], set()).add(int(match[2]))
            elif match := YEARLY_LEVEL_RE.fullmatch(filename):
                if not match[1].endswith("_spectral") and 1 <= int(match[3]) <= level_count:
                    yearly_level_steps.setdefault(match[1], {}).setdefault(int(match[2]), set()).add(int(match[3]))
            elif match := YEARLY_SURFACE_RE.fullmatch(filename):
                if not match[1].endswith("_spectral"):
                    yearly_surface.setdefault(match[1], set()).add(int(match[2]))
            elif match := YEARLY_TIME_RE.fullmatch(filename):
                try:
                    timestamp = float(json.loads(entry.read_text(encoding="utf-8"))["time_seconds"])
                except (OSError, ValueError, KeyError, TypeError):
                    continue
                if math.isfinite(timestamp) and timestamp >= 0:
                    year_times[int(match[1])] = timestamp
        static = tuple(field for field in STATIC_FIELDS if (path / f"{field}.bin").is_file())
        if any(field not in static for field in REQUIRED_STATIC):
            raise DataError("陸面率または地表高度のファイルがありません")
        if not (path / "daily_global.csv").is_file():
            raise DataError("daily_global.csv がありません")
        if any(field not in monthly for field in REQUIRED_MONTHLY):
            raise DataError("解析に使う月平均の場がそろっていません")
        months = tuple(sorted(set.intersection(*(monthly[field] for field in REQUIRED_MONTHLY))))
        if len(months) < int(metadata["calendar"]["months_per_year"]):
            raise DataError("月平均の出力が1年分に足りません")
        yearly_level = {
            field: {year for year, levels in steps.items() if len(levels) == level_count}
            for field, steps in yearly_level_steps.items()
        }
        yearly_level = {field: years for field, years in yearly_level.items() if years}
        year_sets = [years for years in list(yearly_surface.values()) + list(yearly_level.values())]
        years = tuple(sorted(set.union(*year_sets))) if year_sets else ()
        generation = max(generation, path.stat().st_mtime_ns)
        return Run(
            name=name, path=path, metadata=metadata, point_count=point_count, nlat=nlat,
            level_count=level_count, months=months, monthly=monthly,
            yearly_surface=yearly_surface, yearly_level=yearly_level, years=years,
            static=static, generation=generation, year_times=year_times,
        )

    @staticmethod
    def summary(run: Run) -> dict[str, Any]:
        return {
            "name": run.name,
            "case_name": run.metadata.get("case_name", run.name),
            "terrain": run.terrain,
            "truncation": run.metadata.get("numerics", {}).get("spectral_truncation"),
            "point_count": run.point_count,
            "month_count": len(run.months),
            "year_count": len(run.years),
        }

    @classmethod
    def public_metadata(cls, run: Run) -> dict[str, Any]:
        grid = dict(run.metadata["grid"])
        offsets = [0]
        for n in grid["nlon"]:
            offsets.append(offsets[-1] + int(n))
        grid["ring_offsets"] = offsets
        grid["point_count"] = run.point_count
        grid_fields = sorted(name for name in run.monthly if not name.startswith(ZONAL_PREFIXES))
        return {
            **cls.summary(run),
            "metadata": run.metadata,
            "grid": grid,
            "level_count": run.level_count,
            "months": list(run.months),
            "monthly_fields": {name: sorted(run.monthly[name]) for name in grid_fields},
            "zonal_fields": sorted(name for name in run.monthly if name.startswith(ZONAL_PREFIXES)),
            "years": [
                {"index": year, "time_seconds": run.year_times.get(year)} for year in run.years
            ],
            "yearly_surface_fields": {name: sorted(years) for name, years in sorted(run.yearly_surface.items())},
            "yearly_level_fields": {name: sorted(years) for name, years in sorted(run.yearly_level.items())},
            "static_fields": list(run.static),
            "generation": run.generation,
        }

    def grid_values(
        self, run: Run, sampling: str, field: str, indices: list[int] | None, level: int | None
    ) -> bytes:
        if not FIELD_RE.fullmatch(field):
            raise KeyError(field)
        paths: list[Path] = []
        if sampling == "static":
            if field not in run.static or indices or level is not None:
                raise KeyError(field)
            paths.append(run.path / f"{field}.bin")
        elif sampling == "monthly":
            if field.startswith(ZONAL_PREFIXES) or field not in run.monthly or level is not None:
                raise KeyError(field)
            if not indices:
                raise DataError("index を指定してください")
            for index in indices:
                if index not in run.monthly[field]:
                    raise KeyError(index)
                paths.append(run.path / f"monthly_{field}_m{index:04d}.bin")
        elif sampling == "yearly":
            if not indices:
                raise DataError("index を指定してください")
            if field in run.yearly_level:
                if level is None or not 1 <= level <= run.level_count:
                    raise DataError(f"level は 1 から {run.level_count} の範囲で指定してください")
                available = run.yearly_level[field]
                suffix = f"_l{level:02d}"
            elif field in run.yearly_surface:
                if level is not None:
                    raise DataError(f"{field} には層がありません")
                available = run.yearly_surface[field]
                suffix = ""
            else:
                raise KeyError(field)
            for index in indices:
                if index not in available:
                    raise KeyError(index)
                paths.append(run.path / f"yearly_{field}_y{index:04d}{suffix}.bin")
        else:
            raise KeyError(sampling)
        return self._concatenate(paths, run.point_count)

    def zonal_values(self, run: Run, field: str, indices: list[int] | None) -> bytes:
        if not FIELD_RE.fullmatch(field) or not field.startswith(ZONAL_PREFIXES) or field not in run.monthly:
            raise KeyError(field)
        if not indices:
            raise DataError("index を指定してください")
        paths = []
        for index in indices:
            if index not in run.monthly[field]:
                raise KeyError(index)
            paths.append(run.path / f"monthly_{field}_m{index:04d}.bin")
        return self._concatenate(paths, run.nlat * run.level_count)

    @staticmethod
    def _concatenate(paths: list[Path], count: int) -> bytes:
        payload = array("d")
        for path in paths:
            payload.extend(_read_float64(path, count))
        if sys.byteorder != "little":
            payload.byteswap()
        return payload.tobytes()

    @staticmethod
    def daily(run: Run) -> dict[str, Any]:
        path = run.path / "daily_global.csv"
        with path.open(encoding="utf-8", newline="") as source:
            reader = csv.reader(source)
            try:
                header = [name.strip() for name in next(reader)]
            except StopIteration as exc:
                raise DataError("daily_global.csv が空です") from exc
            columns: list[list[float | None]] = [[] for _ in header]
            for row in reader:
                if not row:
                    continue
                if len(row) != len(header):
                    raise DataError("daily_global.csv に列数の合わない行があります")
                for column, text in zip(columns, row):
                    try:
                        value = float(text)
                    except ValueError:
                        value = math.nan
                    column.append(value if math.isfinite(value) else None)
        return {"columns": header, "data": dict(zip(header, columns)), "row_count": len(columns[0]) if columns else 0}


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
            self._send_json({"error": "指定したケース・場・index が見つかりません"}, HTTPStatus.NOT_FOUND)
        except DataError as exc:
            self._send_json({"error": str(exc)}, HTTPStatus.UNPROCESSABLE_ENTITY)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            self._send_json({"error": str(exc)}, HTTPStatus.INTERNAL_SERVER_ERROR)

    def _route_get(self) -> None:
        parsed = urlparse(self.path)
        path = unquote(parsed.path)
        query = parse_qs(parsed.query)
        indices = parse_indices(query["index"][0]) if "index" in query else None
        level = None
        if "level" in query:
            text = query["level"][0]
            if len(query["level"]) != 1 or not text.isdigit():
                raise DataError("level は整数を1つ指定してください")
            level = int(text)

        if path == "/api/runs":
            self._send_json({"runs": [Repository.summary(run) for run in self.repository.discover()]})
            return
        parts = [part for part in path.split("/") if part]
        if len(parts) >= 4 and parts[:2] == ["api", "runs"]:
            run = self.repository.get_run(parts[2])
            if parts[3:] == ["metadata"]:
                self._send_json(Repository.public_metadata(run))
                return
            if parts[3:] == ["daily"]:
                self._send_json(Repository.daily(run))
                return
            if len(parts) == 6 and parts[3] == "grid":
                payload = self.repository.grid_values(run, parts[4], parts[5], indices, level)
                self._send_binary(payload, run.point_count)
                return
            if len(parts) == 5 and parts[3] == "zonal":
                payload = self.repository.zonal_values(run, parts[4], indices)
                self._send_binary(payload, run.nlat * run.level_count)
                return
            raise KeyError(path)
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
                raise DataError("Three.js がインストールされていません。cd viz && npm install を実行してください")
            self._send_file(vendor_path, "text/javascript")
            return
        relative = "index.html" if path == "/" else path.lstrip("/")
        candidate = (self.static_root / relative).resolve()
        if self.static_root.resolve() not in candidate.parents or not candidate.is_file():
            raise KeyError(relative)
        self._send_file(candidate)

    def _send_file(self, path: Path, content_type: str | None = None) -> None:
        payload = path.read_bytes()
        mime = content_type or mimetypes.guess_type(path.name)[0] or "application/octet-stream"
        if path.suffix in {".js", ".mjs"}:
            mime = "text/javascript"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", f"{mime}; charset=utf-8" if mime.startswith("text/") else mime)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(payload)

    def _send_binary(self, payload: bytes, values_per_record: int) -> None:
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Layout", "float64-le")
        self.send_header("X-Values-Per-Record", str(values_per_record))
        self.end_headers()
        self.wfile.write(payload)

    def _send_json(self, value: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        payload = json.dumps(value, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)


def make_handler(output: Path) -> type[AppHandler]:
    script_root = Path(__file__).resolve().parent
    return type(
        "ConfiguredAppHandler",
        (AppHandler,),
        {
            "repository": Repository(output),
            "static_root": script_root / "static",
            "three_build_root": script_root / "node_modules" / "three" / "build",
        },
    )


def parse_args() -> argparse.Namespace:
    script_root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--output", type=Path, default=script_root.parent / "output")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    server = ThreadingHTTPServer((args.host, args.port), make_handler(args.output))
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
