import pytest
from specgrid.model_star import assemble_observation
import numpy.testing as nptesting

@pytest.mark.parametrize("plugin_names",
                         [['doppler', 'rotation', 'resolution'],
                         ['rotation', 'resolution', 'doppler']]
                         )
def test_simple_assemble_modelstar(test_specgrid, h5_test_data, plugin_names):
    model_star = assemble_observation(test_specgrid, plugin_names=plugin_names)
    assert model_star.all_plugins[0].__class__.__name__ == 'RotationalBroadening'
    assert model_star.all_plugins[-1].__class__.__name__ == 'InstrumentConvolve'
    assert model_star.param_names == ['teff', 'logg', 'feh', 'vrot', 'vrad', 'R']
    nptesting.assert_allclose(model_star().flux.value,
                              h5_test_data['test_data/simple_model_star_flux'])
