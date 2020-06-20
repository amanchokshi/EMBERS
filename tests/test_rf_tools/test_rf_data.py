import pytest
from embers.rf_tools import rf_data as rfd

#power, times = rfd.read_data(rf_file='../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
#print(type(times).__name__)
print(len(rfd.tile_pairs(rfd.tile_names)))


def test_read_data_power_shape():
    power, _ = rfd.read_data(rf_file='../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
    assert power.shape == (16655, 112)

def test_read_data_times_shape():
    _, times = rfd.read_data(rf_file='../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
    assert list(times.shape)[0] == 16655

def test_read_data_power_times_shape():
    power, times = rfd.read_data(rf_file='../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
    assert list(times.shape)[0] == list(power.shape)[0]

def test_read_data_return_type():
    _, times = rfd.read_data(rf_file='../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
    assert type(times).__name__ == 'ndarray'

def test_read_data_time_arrow():
    _, times = rfd.read_data(rf_file='../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
    assert times[0] <= times[-1]

def test_tile_names_first():
    tiles = rfd.tile_names()
    assert tiles[0] == 'rf0XX'

def test_tile_names_last():
    tiles = rfd.tile_names()
    assert tiles[-1] == 'S36YY'

def test_tile_names_length():
    tiles = rfd.tile_names()
    assert len(tiles) == 32

def test_tile_pairs_first():
    pairs = rfd.tile_pairs(rfd.tile_names)
    assert pairs[0] == ('rf0XX', 'S06XX')

def test_tile_pairs_last():
    pairs = rfd.tile_pairs(rfd.tile_names)
    assert pairs[-1] == ('rf1YY', 'S36YY')

def test_tile_pairs_length():
    pairs = rfd.tile_pairs(rfd.tile_names)
    assert len(pairs) == 56




