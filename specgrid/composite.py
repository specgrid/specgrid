from collections import OrderedDict
from astropy.units import Quantity


class ModelStar(object):

    param2model = OrderedDict()

    def __init__(self, models_list):
        self.models_list = models_list
        self.param2model = OrderedDict()
        for model in models_list:
            self.param2model.update(OrderedDict([(param, model)
                                          for param in model.parameters]))
        self.parameters = self.param2model.keys()

    def __getattr__(self, item):
        if item in self.param2model:
            return getattr(self.param2model[item], item)
        else:
            return super(ModelStar, self).__getattribute__(item)

    def __setattr__(self, item, value):
        if item in self.param2model:
            return setattr(self.param2model[item], item, value)
        else:
            super(ModelStar, self).__setattr__(item, value)

    def __call__(self):
        spectrum = self.models_list[0]()
        for model in self.models_list[1:]:
            spectrum = model(spectrum)
        return spectrum
    
    def evaluate(self, **kwargs):
        for kwarg in kwargs:
            current_value = getattr(self, kwarg)
            new_value = kwargs[kwarg]
            if hasattr(current_value, 'unit'):
                new_value = Quantity(new_value, current_value.unit)
            setattr(self, kwarg, new_value)
        return self()

def assemble_model_star(spectral_grid, spectrum=None, normalize_pol=None, **kwargs):

    for key in kwargs:
        pass


    convolve = InstrumentConvolve(R=self.spectral_parameters['R'])
    #rot = RotationalBroadening()
    doppler = DopplerShift()
    interp = Interpolate(self.to_spectrum_1d())
    one_chip_size = len(self.table['x'])
    parts = [slice(None, one_chip_size),
             slice(one_chip_size, 2*one_chip_size),
             slice(2*one_chip_size, None)]
    norm = NormalizeParts(self.to_spectrum_1d(), parts=parts,
                          npol=self.spectral_parameters['npol'])

    model_star = ModelStar([spectral_grid, doppler, convolve, interp, norm])

