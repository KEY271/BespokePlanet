import json
import math
import tempfile
import unittest
from array import array
from pathlib import Path

from server import Repository, _gauss_legendre_weights


class RepositoryTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        root = Path(self.temporary.name)
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

    def test_conservation_for_uniform_flow(self):
        run = self.repository.get_run("tiny")
        data = self.repository.conservation(run)
        metrics = {item["id"]: item["values"] for item in data["metrics"]}
        self.assertEqual(data["times_seconds"], [0.0, 20.0, 40.0])
        self.assertAlmostEqual(metrics["mean_kinetic_energy"][0], 12.5)
        self.assertAlmostEqual(metrics["mean_relative_vorticity"][0], 0.0)
        self.assertEqual(metrics["mean_kinetic_energy"][0], metrics["mean_kinetic_energy"][1])


if __name__ == "__main__":
    unittest.main()
