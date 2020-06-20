import pytest
from embers.rf_tools import rf_data as rfd

def test_rf_data_power_shape():
    power, _ = rfd.rf_data(rf_file='rf_file='../data/')

