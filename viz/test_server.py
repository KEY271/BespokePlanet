import json
import math
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from array import array
from http.server import ThreadingHTTPServer
from pathlib import Path

from server import REQUIRED_MONTHLY, DataError, Repository, make_handler, parse_indices

MU = [-0.5, 0.5]
NLON = [3, 3]
LEVELS = 2


def write_values(path: Path, values) -> None:
    path.write_bytes(array("d", values).tobytes())


def land_sea_metadata(**overrides):
    metadata = {
        "schema_version": 3,
        "case_name": "tiny_land_sea",
        "equation": "moist_hydrostatic_atmosphere",
        "simulation": {"start_calendar_time": "0001-04-01 00:00:00", "time_step_seconds": 1200},
        "calendar": {"solar_day_seconds": 86400, "days_per_month": 30, "months_per_year": 12, "days_per_year": 360},
        "grid": {"type": "octahedral_gaussian", "mu": MU, "nlon": NLON},
        "hybrid_a_half_pa": [0, 100, 0],
        "hybrid_b_half": [0, 0.5, 1],
        "reference_half_level_pressure_pa": [0, 50000, 100000],
        "sea_ice": {"enabled": True},
        "snow": {"enabled": True, "masking_water_equivalent_kg_m-2": 50},
        "topography": {"source": "analytic"},
        "output": {"static_land_fraction": "land_fraction.bin"},
    }
    metadata.update(overrides)
    return metadata


class LandSeaFixture(unittest.TestCase):
    months = 13

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.run_path = self.make_run("tiny", land_sea_metadata())
        self.repository = Repository(self.root)

    def tearDown(self):
        self.temporary.cleanup()

    def make_run(self, name, metadata):
        path = self.root / name
        path.mkdir()
        (path / "metadata.json").write_text(json.dumps(metadata))
        points = sum(NLON)
        write_values(path / "land_fraction.bin", [1, 0.5, 0, 0, 0.2, 0.9])
        write_values(path / "surface_height.bin", [100, 50, 0, 0, 10, 2000])
        write_values(path / "ocean_q_flux.bin", [0] * points)
        for month in range(1, self.months + 1):
            for field in REQUIRED_MONTHLY:
                count = len(MU) * LEVELS if field.startswith("zonal_") else points
                write_values(path / f"monthly_{field}_m{month:04d}.bin", [month + 0.25 * i for i in range(count)])
        for year in (1, 2):
            write_values(path / f"yearly_surface_temperature_y{year:04d}.bin", [280 + year] * points)
            write_values(path / f"yearly_log_surface_pressure_spectral_y{year:04d}.bin", [0] * 4)
            (path / f"yearly_time_y{year:04d}.json").write_text(json.dumps({"time_seconds": (year - 1) * 31104000.0}))
            for level in range(1, LEVELS + 1):
                write_values(path / f"yearly_u_y{year:04d}_l{level:02d}.bin", [10 * level + year] * points)
                write_values(path / f"yearly_temperature_spectral_y{year:04d}_l{level:02d}.bin", [0] * 4)
        (path / "daily_global.csv").write_text("time_seconds,simulation_day,mean_surface_temperature_k\n0,0,288.5\n86400,1,NaN\n")
        return path


class RepositoryTest(LandSeaFixture):
    def test_discovers_only_current_land_sea_runs(self):
        self.make_run("old_moist", land_sea_metadata(snow={"enabled": False}))
        self.make_run("dry", land_sea_metadata(equation="dry_hydrostatic_atmosphere"))
        broken = self.make_run("broken", land_sea_metadata())
        (broken / "monthly_zonal_v_m0001.bin").unlink()
        for month in range(2, self.months + 1):
            (broken / f"monthly_zonal_v_m{month:04d}.bin").unlink()
        self.assertEqual([run.name for run in self.repository.discover()], ["tiny"])

    def test_requires_one_complete_year(self):
        short = self.make_run("short", land_sea_metadata())
        for month in range(12, self.months + 1):
            (short / f"monthly_precipitation_m{month:04d}.bin").unlink()
        with self.assertRaises(DataError):
            self.repository.get_run("short")

    def test_public_metadata_lists_fields(self):
        run = self.repository.get_run("tiny")
        public = Repository.public_metadata(run)
        self.assertEqual(public["months"], list(range(1, self.months + 1)))
        self.assertEqual(public["grid"]["ring_offsets"], [0, 3, 6])
        self.assertEqual(public["grid"]["point_count"], 6)
        self.assertIn("precipitation", public["monthly_fields"])
        self.assertNotIn("zonal_v", public["monthly_fields"])
        self.assertEqual(public["zonal_fields"], ["zonal_v"])
        self.assertEqual(public["yearly_level_fields"], {"u": [1, 2]})
        self.assertEqual(sorted(public["yearly_surface_fields"]), ["surface_temperature"])
        self.assertEqual(public["years"], [{"index": 1, "time_seconds": 0.0}, {"index": 2, "time_seconds": 31104000.0}])
        self.assertEqual(public["terrain"], "analytic")

    def test_grid_values_concatenate_months(self):
        run = self.repository.get_run("tiny")
        values = array("d")
        values.frombytes(self.repository.grid_values(run, "monthly", "precipitation", [2, 3], None))
        self.assertEqual(list(values), [2 + 0.25 * i for i in range(6)] + [3 + 0.25 * i for i in range(6)])

    def test_level_fields_require_a_level(self):
        run = self.repository.get_run("tiny")
        with self.assertRaises(DataError):
            self.repository.grid_values(run, "yearly", "u", [1], None)
        values = array("d")
        values.frombytes(self.repository.grid_values(run, "yearly", "u", [2], 2))
        self.assertEqual(list(values), [22.0] * 6)
        with self.assertRaises(DataError):
            self.repository.grid_values(run, "yearly", "surface_temperature", [1], 1)

    def test_rejects_unknown_fields_and_indices(self):
        run = self.repository.get_run("tiny")
        with self.assertRaises(KeyError):
            self.repository.grid_values(run, "monthly", "precipitation", [99], None)
        with self.assertRaises(KeyError):
            self.repository.grid_values(run, "monthly", "../metadata", [1], None)
        with self.assertRaises(KeyError):
            self.repository.grid_values(run, "yearly", "temperature_spectral", [1], 1)
        with self.assertRaises(KeyError):
            self.repository.grid_values(run, "monthly", "zonal_v", [1], None)
        with self.assertRaises(KeyError):
            self.repository.get_run("..")

    def test_zonal_values_have_level_latitude_layout(self):
        run = self.repository.get_run("tiny")
        values = array("d")
        values.frombytes(self.repository.zonal_values(run, "zonal_v", [1]))
        self.assertEqual(len(values), len(MU) * LEVELS)

    def test_non_finite_and_truncated_files_are_rejected(self):
        write_values(self.run_path / "monthly_precipitation_m0002.bin", [math.nan] * 6)
        write_values(self.run_path / "monthly_precipitation_m0003.bin", [1.0] * 5)
        run = self.repository.get_run("tiny")
        with self.assertRaises(DataError):
            self.repository.grid_values(run, "monthly", "precipitation", [2], None)
        with self.assertRaises(DataError):
            self.repository.grid_values(run, "monthly", "precipitation", [3], None)

    def test_daily_columns_turn_non_finite_values_into_null(self):
        daily = Repository.daily(self.repository.get_run("tiny"))
        self.assertEqual(daily["row_count"], 2)
        self.assertEqual(daily["data"]["mean_surface_temperature_k"], [288.5, None])

    def test_parse_indices(self):
        self.assertEqual(parse_indices("1,3-5,2"), [1, 3, 4, 5, 2])
        for text in ("", "0", "3-1", "a", "1,,2"):
            with self.assertRaises(DataError):
                parse_indices(text)


class HttpTest(LandSeaFixture):
    def setUp(self):
        super().setUp()
        self.server = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(self.root))
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.base = f"http://127.0.0.1:{self.server.server_address[1]}"

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        super().tearDown()

    def get(self, path):
        with urllib.request.urlopen(self.base + path) as response:
            return response.status, response.headers, response.read()

    def test_routes(self):
        _, _, body = self.get("/api/runs")
        self.assertEqual([run["name"] for run in json.loads(body)["runs"]], ["tiny"])
        _, _, body = self.get("/api/runs/tiny/metadata")
        self.assertEqual(json.loads(body)["level_count"], LEVELS)
        _, headers, body = self.get("/api/runs/tiny/grid/monthly/precipitation?index=1-12")
        self.assertEqual(len(body), 12 * 6 * 8)
        self.assertEqual(headers["X-Values-Per-Record"], "6")
        _, _, body = self.get("/api/runs/tiny/grid/static/land_fraction")
        self.assertEqual(len(body), 6 * 8)
        _, _, body = self.get("/api/runs/tiny/zonal/zonal_v?index=1,2")
        self.assertEqual(len(body), 2 * len(MU) * LEVELS * 8)
        _, _, body = self.get("/api/runs/tiny/daily")
        self.assertEqual(json.loads(body)["columns"][0], "time_seconds")
        status, _, body = self.get("/")
        self.assertEqual(status, 200)
        self.assertIn(b"app.js", body)

    def test_errors(self):
        for path, status in (
            ("/api/runs/missing/metadata", 404),
            ("/api/runs/tiny/grid/monthly/precipitation?index=0", 422),
            ("/api/runs/tiny/grid/yearly/u?index=1", 422),
            ("/api/runs/tiny/grid/monthly/nope?index=1", 404),
            ("/../server.py", 404),
        ):
            with self.subTest(path=path):
                with self.assertRaises(urllib.error.HTTPError) as caught:
                    self.get(path)
                self.assertEqual(caught.exception.code, status)


if __name__ == "__main__":
    unittest.main()
