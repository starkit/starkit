from astropy import modeling

class GridConnector(modeling.Model):
    inputs = tuple('teff', 'logg', 'feh', 'radius')
    outputs = ('wavelength', 'flux')

    def __init__(self, grid, **kwargs):
        super(GridConnector, self).__init__(**kwargs)
        self.grid = grid


    def evaluate(self, teff, logg, feh, radius):
        return self.grid.evaluate(teff, logg, feh, radius)

