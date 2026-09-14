import json
import math
import tempfile
import unittest
from array import array
from pathlib import Path

from server import DataError, Repository, _gauss_legendre_weights


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
