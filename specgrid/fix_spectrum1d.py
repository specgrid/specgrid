import specutils
from astropy import units as u

class Spectrum1D(specutils.Spectrum1D):

    @property
    def flux(self):
        return self.data * u.Unit(self.unit)

#    @property
#    def uncertainty(self):
#        return self._uncertainty

#    @uncertainty.setter
#    def uncertainty_setter(self, value):
#        self._uncertainty = value

