import astropy.units as u
import h5py

from .loader import FileLoader


class ExoMolFileLoader(FileLoader):
    '''
    File loader for ExoMol file. It works for opacities file in the `TauREx.h5` format.
    These file can be downloaded from `ExoMol website <http://exomol.com/data/data-types/opacity>`_.

    Notes
    -----
    This class is implicitly called by the model when initialised if the matched input is used.

    Examples
    ---------
    Let's build the Rosseland method using an Exomol file as input

    >>> from rapoc import Rosseland
    >>> exomolFile = 'example.TauREx.h5'
    >>> model = Rosseland(input_data=exomolFile)

    Now the model is ready to be used

    '''

    def _open(self):
        return h5py.File(self.filename, 'r')

    def _close(self, opened_file):
        opened_file.close()

    def _read_molecule_name(self, opened_file):
        return opened_file['mol_name'][()][0]

    def _get_mol_mass(self, mol):
        from molmass import Formula
        f = Formula(mol)
        self.mol_mass = (f.mass * u.u / u.ct).to(u.g / u.ct)
        return self.mol_mass

    def _read_pressure_grid(self, opened_file):
        pressure_units = opened_file['p'].attrs['units']
        pressure_grid = opened_file['p'][:] * u.Unit(pressure_units)
        return pressure_grid.si

    def _read_temperature_grid(self, opened_file):
        temperature_units = opened_file['t'].attrs['units']
        if temperature_units == 'kelvin': temperature_units = 'K'
        temperature_grid = opened_file['t'][:] * u.Unit(temperature_units)
        return temperature_grid.si

    def _read_wavenumber_grid(self, opened_file):
        return opened_file['bin_edges'][:] / u.cm

    def _read_ktable(self, opened_file):
        xsec_grid_units = opened_file['xsecarr'].attrs['units']
        xsec_grid_raw = opened_file['xsecarr'][:]

        if xsec_grid_units == 'cm^2/molecule':  # it's molecule, not mol!!!
            xsec_grid_raw *= u.cm ** 2 / u.ct  # I use counts as molecule, to help astropy.units
        else:
            raise AttributeError('units to be handle by hand')
        xsec_grid = xsec_grid_raw / self.mol_mass
        return xsec_grid.si
