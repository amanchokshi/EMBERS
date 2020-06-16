import sys
import unittest
from pathlib import Path
sys.path.append('../../code/beam_pointings')
import dead_dipoles as dd


class Test_dead_dipoles(unittest.TestCase):

    metafits     = Path('../../data/tests/metafits')
    
    tiles = [
            'HexS6', 'HexS7',
            'HexS8', 'HexS9',
            'HexS10', 'HexS12',
            'HexS29', 'HexS30',
            'HexS31', 'HexS32',
            'HexS33', 'HexS34', 
            'HexS35', 'HexS36']

        
    def test_org_pointing_json(self):
        flags = dd.find_flag(self.metafits, self.tiles)
        self.assertEqual(flags["HexS6X"][0], 0)
        self.assertEqual(flags["HexS33Y"][0], 9)
        

if __name__ == "__main__":
    unittest.main()


