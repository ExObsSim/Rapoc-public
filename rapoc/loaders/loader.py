from abc import abstractmethod


class FileLoader:
    '''
    Base file loaders class.

    '''

    def __init__(self, filename):
        '''
        Parameters
        ----------
        filename: str
            data file name
        '''
        self.filename = filename

    def read_content(self):
        '''
        Reads the file content and returns the needed valued for the opacity models.

        Returns
        -------
        str
            molecule name
        astropy.units.Quantity:
            molecular mass
        astropy.units.Quantity:
            data pressure grid in si units
        astropy.units.Quantity:
            data temperature grid in si units
        astropy.units.Quantity:
            data wavenumber grid
        astropy.units.Quantity:
            data opacities grid in si units
        '''
        opened_file = self._open()
        mol = self._read_molecule_name(opened_file)
        mol_mass = self._get_mol_mass(mol)
        pressure_grid = self._read_pressure_grid(opened_file)
        temperature_grid = self._read_temperature_grid(opened_file)
        wavenumber_grid = self._read_wavenumber_grid(opened_file)
        ktable = self._read_ktable(opened_file)
        self._close(opened_file)
        return mol, mol_mass, pressure_grid, temperature_grid, wavenumber_grid, ktable

    @abstractmethod
    def _open(self):
        raise NotImplementedError

    @abstractmethod
    def _close(self, opened_file):
        raise NotImplementedError

    @abstractmethod
    def _read_molecule_name(self, opened_file):
        raise NotImplementedError

    @abstractmethod
    def _get_mol_mass(self, opened_file):
        raise NotImplementedError

    @abstractmethod
    def _read_pressure_grid(self, opened_file):
        raise NotImplementedError

    @abstractmethod
    def _read_temperature_grid(self, opened_file):
        raise NotImplementedError

    @abstractmethod
    def _read_wavenumber_grid(self, opened_file):
        raise NotImplementedError

    @abstractmethod
    def _read_ktable(self, opened_file):
        raise NotImplementedError
