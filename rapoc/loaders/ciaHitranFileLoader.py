import astropy.units as u
import numpy as np
from taurex.cia.hitrancia import HitranCIA

from .loader import FileLoader


class CiaHitranFileLoader(FileLoader):
    '''
    File loader for ExoMol file. It works for opacities file in the `TauREx.h5` format.
    These file can be downloaded from `ExoMol website <http://exomol.com/data/data-types/opacity>`_.

    Parameters
    ----------
    input_data: str
        data file name
    '''

    def read_content(self):
        file_content = HitranCIA(self.filename)
        mol = None
        mol_mass = None
        temperature_grid = file_content.temperatureGrid * u.K
        wavenumber_grid = file_content.wavenumberGrid / u.cm
        ktable = file_content._xsec_grid[np.newaxis, ...]
        pressure_grid = None
        return mol, mol_mass, pressure_grid, temperature_grid, wavenumber_grid, ktable

    def _get_mol_mass(self, mol):
        from molmass import Formula
        f = Formula(mol)
        self.mol_mass = (f.mass * u.u / u.ct).to(u.g / u.ct)
        return self.mol_mass
