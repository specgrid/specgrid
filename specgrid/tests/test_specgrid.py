import pytest


@pytest.skipif(True, 'needs to be implemented')
def test_specgrid_parameters()
    specgrid = BaseSpectralGrid('test_specgrid.h5')
    assert specgrid.teff is None
    assert specgrid.logg is None
    assert specgrid.feh is None
    
