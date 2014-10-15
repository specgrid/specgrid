import pytest
from specgrid.model_star import assemble_observation, ModelStar
from specgrid.plugins import DopplerShift
import numpy.testing as nptesting

@pytest.mark.parametrize("plugin_names",
                         [['doppler', 'rotation', 'resolution'],
                         ['rotation', 'resolution', 'doppler']]
                         )
def test_simple_assemble_observation(test_specgrid, h5_test_data, plugin_names):
    observation = assemble_observation(test_specgrid, plugin_names=plugin_names)
    assert observation.all_plugins[0].__class__.__name__ == 'RotationalBroadening'
    assert observation.all_plugins[-1].__class__.__name__ == 'InstrumentConvolve'
    assert observation.param_names == ['teff', 'logg', 'feh', 'vrot', 'vrad', 'R']
    nptesting.assert_allclose(observation().flux.value,
                              h5_test_data['test_data/simple_model_star_flux'])


def test_simple_modelstar(test_specgrid):
    doppler = DopplerShift()
    model_star = ModelStar(test_specgrid, [doppler])
    assert model_star.param_names == ['teff', 'logg', 'feh', 'vrad']