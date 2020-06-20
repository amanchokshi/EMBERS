import sys
import json
import unittest
from pathlib import Path

sys.path.append("../../code/sat_ephemeris")
import plot_ephemeris as pe


class Test_plot_ephemeris(unittest.TestCase):

    json_dir = "../../data/tests/ephem_json"
    out_dir = "./"
    sat_id = 25986

    def test_plot_ephem(self):

        pe.plot_ephem(self.json_dir, self.out_dir, self.sat_id)
        out_file = Path(f"./{self.sat_id}.png")

        self.assertTrue(out_file.is_file())
        out_file.unlink()


if __name__ == "__main__":
    unittest.main()
