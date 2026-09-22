import json
import csv
import math
import os
import tempfile
import unittest
from array import array
from pathlib import Path

from server import DataError, Repository, _gauss_legendre_weights


@unittest.skipUnless(os.environ.get("BESPOKE_SEA_ICE_OUTPUT"), "requires the check_sea_ice Fortran output fixture")
class FortranSeaIceOutputTest(unittest.TestCase):
    def test_generated_output(self):
        path = Path(os.environ["BESPOKE_SEA_ICE_OUTPUT"]).resolve()
        repository = Repository(path.parent)
        for name in (path.name, path.name + ".yearly"):
            run = repository.get_run(name)
            self.assertTrue(run.metadata["sea_ice"]["enabled"])
            for field in run.fields:
                payload, stats = repository.field_frame(run, field, 1)
                self.assertEqual(len(payload), 4 * run.metadata["grid"]["point_count"])
                self.assertTrue(math.isfinite(stats["mean"]))
            payload, _ = repository.field_frame(run, "ocean_temperature", 1)
            self.assertTrue(math.isnan(array("f", payload)[0]))
            self.assertGreaterEqual(min(v for v in array("f", payload) if math.isfinite(v)), 271.35 - 1e-4)
        with (path / "daily_global.csv").open() as source:
            row = next(csv.DictReader(source))
        self.assertNotIn(None, row)
        self.assertTrue(all(v is not None and math.isfinite(float(v)) for v in row.values()))
        self.assertLess(float(row["maximum_ice_energy_residual_j_m-2"]), 2e-5)
        self.assertLess(float(row["maximum_ice_projection_energy_residual_j_m-2"]), 2e-5)
        self.assertGreater(int(row["ice_checked_cells"]), 0)


class RadiationRepositoryTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.run_path = self.root / "ice"
        self.run_path.mkdir()
        metadata = {
            "equation": "moist_hydrostatic_atmosphere",
            "simulation": {"time_step_seconds": 1200},
            "calendar": {"solar_day_seconds": 86400, "days_per_month": 30, "months_per_year": 12},
            "grid": {"mu": [-0.5, 0.5], "nlon": [2, 2]},
            "output": {"monthly_surface_temperature": "monthly_surface_temperature_m{month:04d}.bin"},
            "surface_tiles": {},
        }
        (self.run_path / "metadata.json").write_text(json.dumps(metadata))
        self.write("land_fraction.bin", [1, 0.5, 0, 0])
        for prefix, suffix in (("monthly", "m0001"), ("yearly", "y0001")):
            for field, values in {
                "surface_temperature": [280, 270, 270, 280],
                "ocean_temperature": [0, 271.35, 271.35, 280],
                "sea_ice_fraction": [0, 0.5, 1, 0],
                "sea_ice_volume": [0, 0.25, 2, 0],
                "sea_ice_thickness": [0, 0.5, 2, 0],
                "sea_ice_temperature": [0, 260, 270, 0],
            }.items():
                self.write(f"{prefix}_{field}_{suffix}.bin", values)
        self.repository = Repository(self.root)

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, filename, values):
        (self.run_path / filename).write_bytes(array("d", values).tobytes())

    def test_monthly_and_yearly_discovery(self):
        self.assertEqual([run.name for run in self.repository.discover()], ["ice", "ice.yearly"])
        monthly = self.repository.get_run("ice")
        self.assertEqual(monthly.steps, (1,))
        metadata = self.repository.public_metadata(monthly)
        self.assertFalse(metadata["supports_streamlines"])
        self.assertFalse(any(field["uses_level"] for field in metadata["available_fields"]))
        self.assertEqual(metadata["frame_times_seconds"]["1"], 15 * 86400)
        yearly = self.repository.get_run("ice.yearly")
        self.assertEqual(yearly.steps, (1,))
        self.assertEqual(yearly.metadata["frame_times_seconds"]["1"], 0)

    def test_ice_mask_and_color_statistics(self):
        run = self.repository.get_run("ice")
        payload, stats = self.repository.field_frame(run, "sea_ice_temperature", 1)
        values = array("f", payload)
        self.assertTrue(math.isnan(values[0]))
        self.assertTrue(math.isnan(values[3]))
        self.assertEqual(list(values[1:3]), [260, 270])
        self.assertEqual(stats["minimum"], 260)
        self.assertEqual(stats["mean"], 265)
        self.assertEqual(self.repository.field_statistics(run, "sea_ice_temperature")["minimum"], 260)

    def test_all_ice_missing_is_renderable(self):
        self.write("monthly_sea_ice_fraction_m0001.bin", [0, 0, 0, 0])
        run = self.repository.get_run("ice")
        payload, stats = self.repository.field_frame(run, "sea_ice_thickness", 1)
        self.assertTrue(all(math.isnan(value) for value in array("f", payload)))
        self.assertEqual(stats["maximum"], 0)
        self.assertEqual(self.repository.field_statistics(run, "sea_ice_thickness")["maximum"], 0)

    def test_corrupt_source_is_not_hidden_by_mask(self):
        self.write("monthly_sea_ice_temperature_m0001.bin", [math.nan, 260, 270, 0])
        with self.assertRaises(DataError):
            self.repository.field_frame(self.repository.get_run("ice"), "sea_ice_temperature", 1)

    def test_ocean_temperature_masks_only_land(self):
        payload, stats = self.repository.field_frame(self.repository.get_run("ice.yearly"), "ocean_temperature", 1)
        values = array("f", payload)
        self.assertTrue(math.isnan(values[0]))
        self.assertTrue(all(math.isfinite(value) for value in values[1:]))
        self.assertAlmostEqual(stats["minimum"], 271.35, places=4)

    def test_legacy_ground_temperature_without_land_mask(self):
        metadata_path = self.run_path / "metadata.json"
        metadata = json.loads(metadata_path.read_text())
        metadata.pop("surface_tiles")
        metadata["ground"] = {}
        metadata_path.write_text(json.dumps(metadata))
        (self.run_path / "land_fraction.bin").unlink()
        self.write("monthly_deep_temperature_m0001.bin", [280, 275, 270, 265])
        payload, stats = self.repository.field_frame(self.repository.get_run("ice"), "deep_temperature", 1)
        self.assertEqual(list(array("f", payload)), [280, 275, 270, 265])
        self.assertEqual(stats["minimum"], 265)


class RepositoryTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
        self.root = root
        run = root / "tiny"
        run.mkdir()
        mu = [-1 / math.sqrt(3), 1 / math.sqrt(3)]
        metadata = {
            "case_name": "tiny",
            "simulation": {
                "time_step_seconds": 10.0,
                "number_of_steps": 4,
                "snapshot_interval_steps": 2,
            },
            "physical_constants": {"earth_radius_m": 2.0, "rotation_rate_rad_s": 0.25},
            "grid": {"mu": mu, "nlon": [2, 2], "ring_offsets": [0, 2, 4], "point_count": 4},
        }
        (run / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
        for step in (0, 1, 2, 4):
            for field, values in {
                "zeta": [0.0, 0.0, 0.0, 0.0],
                "delta": [-2.0, -1.0, 1.0, 2.0],
                "eta": [10.0, 20.0, 30.0, 40.0],
                "u": [3.0, 3.0, 3.0, 3.0],
                "v": [4.0, 4.0, 4.0, 4.0],
            }.items():
                (run / f"{field}_{step:05d}.bin").write_bytes(array("d", values).tobytes())
        self.repository = Repository(root)

    def tearDown(self):
        self.temporary.cleanup()

    def test_quadrature_weights_integrate_constant(self):
        nodes = [-1 / math.sqrt(3), 1 / math.sqrt(3)]
        self.assertAlmostEqual(sum(_gauss_legendre_weights(nodes)), 2.0)

    def test_discovers_complete_frames_and_encodes_speed(self):
        run = self.repository.get_run("tiny")
        self.assertEqual(run.steps, (0, 2, 4))
        payload, stats = self.repository.speed_frame(run, 0)
        self.assertEqual(len(payload), 4 * 4)
        self.assertAlmostEqual(stats["maximum"], 5.0)

    def test_wind_frame_packs_u_then_v(self):
        run = self.repository.get_run("tiny")
        payload, maximum_speed = self.repository.wind_frame(run, 0)
        values = list(array("f", payload))
        self.assertEqual(values[:4], [3.0] * 4)
        self.assertEqual(values[4:], [4.0] * 4)
        self.assertAlmostEqual(maximum_speed, 5.0)

    def test_exposes_and_encodes_all_supported_fields(self):
        run = self.repository.get_run("tiny")
        self.assertEqual(run.fields, ("zeta", "delta", "eta", "speed"))
        metadata = self.repository.public_metadata(run)
        self.assertIsInstance(metadata["data_generation"], int)
        self.assertEqual(
            [field["id"] for field in metadata["available_fields"]],
            ["zeta", "delta", "eta", "speed"],
        )
        payload, stats = self.repository.field_frame(run, "delta", 0)
        self.assertEqual(list(array("f", payload)), [-2.0, -1.0, 1.0, 2.0])
        self.assertEqual(stats["minimum"], -2.0)
        self.assertEqual(stats["maximum_absolute"], 2.0)
        self.assertEqual(stats["p995_absolute"], 2.0)

    def test_field_statistics_use_the_whole_run(self):
        run = self.repository.get_run("tiny")
        stats = self.repository.field_statistics(run, "eta")
        self.assertEqual(stats, {"minimum": 10.0, "maximum": 40.0, "maximum_absolute": 40.0})

    def test_rejects_unavailable_field(self):
        run = self.repository.get_run("tiny")
        with self.assertRaises(KeyError):
            self.repository.field_frame(run, "temperature", 0)

    def test_barotropic_run_only_exposes_zeta_and_speed(self):
        run_path = self.root / "barotropic"
        run_path.mkdir()
        metadata = json.loads((self.root / "tiny" / "metadata.json").read_text())
        (run_path / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
        for field in ("zeta", "u", "v"):
            source = self.root / "tiny" / f"{field}_00000.bin"
            (run_path / source.name).write_bytes(source.read_bytes())
        run = self.repository.get_run("barotropic")
        self.assertEqual(run.fields, ("zeta", "speed"))

    def test_conservation_for_uniform_flow(self):
        run = self.repository.get_run("tiny")
        data = self.repository.conservation(run)
        metrics = {item["id"]: item["values"] for item in data["metrics"]}
        self.assertEqual(data["times_seconds"], [0.0, 20.0, 40.0])
        self.assertAlmostEqual(metrics["mean_kinetic_energy"][0], 12.5)
        self.assertAlmostEqual(metrics["mean_relative_vorticity"][0], 0.0)
        self.assertEqual(metrics["mean_kinetic_energy"][0], metrics["mean_kinetic_energy"][1])

    def test_dry_run_exposes_level_fields_and_normalizes_legacy_grid(self):
        run_path = self._write_dry_run()
        run = self.repository.get_run(run_path.name)
        self.assertTrue(run.is_dry)
        self.assertEqual(run.steps, (0, 2))
        self.assertEqual(
            run.fields,
            ("surface_pressure", "temperature", "zeta", "delta", "u", "v", "speed"),
        )
        metadata = self.repository.public_metadata(run)
        self.assertEqual(metadata["grid"]["ring_offsets"], [0, 2, 4])
        self.assertEqual(metadata["grid"]["point_count"], 4)
        self.assertEqual(metadata["simulation"]["time_step_seconds"], 10.0)
        self.assertFalse(metadata["supports_conservation_diagnostics"])
        fields = {field["id"]: field for field in metadata["available_fields"]}
        self.assertTrue(fields["temperature"]["uses_level"])
        self.assertFalse(fields["surface_pressure"]["uses_level"])
        vertical = metadata["vertical_coordinate"]
        self.assertEqual(vertical["default_level"], 2)
        self.assertEqual(len(vertical["reference_full_level_pressure_pa"]), 3)
        self.assertNotIn("interpolation", vertical)

    def test_dry_field_is_read_from_requested_layer(self):
        run = self.repository.get_run(self._write_dry_run().name)
        payload, _ = self.repository.field_frame(run, "zeta", 0, 2)
        self.assertEqual(list(array("f", payload)), [1.0] * 4)

        speed_payload, _ = self.repository.field_frame(run, "speed", 0, 2)
        for value in array("f", speed_payload):
            self.assertAlmostEqual(value, 5.0, places=5)

        wind_payload, maximum_speed = self.repository.wind_frame(run, 0, 2)
        wind = list(array("f", wind_payload))
        self.assertEqual(wind[:4], [3.0] * 4)
        self.assertEqual(wind[4:], [4.0] * 4)
        self.assertAlmostEqual(maximum_speed, 5.0)

    def test_dry_field_rejects_layer_outside_run(self):
        run = self.repository.get_run(self._write_dry_run().name)
        with self.assertRaisesRegex(DataError, "Level must be between 1 and 3"):
            self.repository.field_frame(run, "temperature", 0, 0)

    def test_dry_field_statistics_are_computed_per_layer(self):
        run = self.repository.get_run(self._write_dry_run().name)
        self.assertEqual(
            self.repository.field_statistics(run, "zeta", 1)["maximum"], 0.0
        )
        self.assertEqual(
            self.repository.field_statistics(run, "zeta", 3)["maximum"], 2.0
        )

    def test_surface_pressure_does_not_require_a_layer(self):
        run = self.repository.get_run(self._write_dry_run().name)
        payload, _ = self.repository.field_frame(run, "surface_pressure", 0)
        self.assertEqual(list(array("f", payload)), [100000.0, 80000.0, 100000.0, 80000.0])

    def _write_dry_run(self) -> Path:
        run_path = self.root / "dry"
        if run_path.exists():
            return run_path
        run_path.mkdir()
        surface_pressure = [100000.0, 80000.0, 100000.0, 80000.0]
        a_half = [10000.0, 0.0, 0.0, 0.0]
        b_half = [0.0, 0.4, 0.7, 1.0]
        reference_half = [
            a + b * 100000.0 for a, b in zip(a_half, b_half)
        ]
        metadata = {
            "equation": "dry_hydrostatic_atmosphere",
            "time_step_seconds": 10.0,
            "number_of_steps": 2,
            "snapshot_interval_steps": 2,
            "number_of_levels": 3,
            "reference_half_level_pressure_pa": reference_half,
            "hybrid_a_half_pa": a_half,
            "hybrid_b_half": b_half,
            "grid": {
                "mu": [-1 / math.sqrt(3), 1 / math.sqrt(3)],
                "nlon": [2, 2],
            },
        }
        (run_path / "metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
        for step in (0, 2):
            (run_path / f"surface_pressure_{step:05d}.bin").write_bytes(
                array("d", surface_pressure).tobytes()
            )
            full_pressure_by_point = []
            for ps in surface_pressure:
                half = [a + b * ps for a, b in zip(a_half, b_half)]
                full_pressure_by_point.append(
                    [
                        self._full_pressure(half[level], half[level + 1])
                        for level in range(3)
                    ]
                )
            for level in range(3):
                temperature = [
                    math.log(full_pressure_by_point[point][level])
                    for point in range(4)
                ]
                values_by_field = {
                    "temperature": temperature,
                    "zeta": [float(level)] * 4,
                    "delta": [-float(level)] * 4,
                    "u": [3.0] * 4,
                    "v": [4.0] * 4,
                }
                for field, values in values_by_field.items():
                    (run_path / f"{field}_l{level + 1:02d}_{step:05d}.bin").write_bytes(
                        array("d", values).tobytes()
                    )
        return run_path

    @staticmethod
    def _full_pressure(top: float, bottom: float) -> float:
        alpha = 1.0 - top * math.log(bottom / top) / (bottom - top)
        return bottom * math.exp(-alpha)


if __name__ == "__main__":
    unittest.main()
