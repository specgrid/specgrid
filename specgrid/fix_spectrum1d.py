import specutils
from astropy import units as u

class Spectrum1D(specutils.Spectrum1D):

    @property
    def flux(self):
        return u.Quantity(self.data, self.unit)



    def uncertainty_getter(self):
        return self._uncertainty

    def uncertainty_setter(self, value):
        self._uncertainty = value

    uncertainty = property(uncertainty_getter, uncertainty_setter)