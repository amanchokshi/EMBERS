import sys
import unittest
import numpy as np
from pathlib import Path

sys.path.append("../../code/beam_pointings")
import download_metafits as dm


class Test_sort_pointings(unittest.TestCase):

    gps_time = 1253874336
    dm.cerberus_metafit(gps_time, "./")

    def test_cerberus_metafits(self):
        metafits = Path(f"{self.gps_time}.metafits")
        self.assertTrue(metafits.is_file())

        if metafits.is_file():
            metafits.unlink()


if __name__ == "__main__":
    unittest.main()
