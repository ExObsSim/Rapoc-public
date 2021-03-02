from astropy import units as u


class DictLoader:
    '''
    Dict loader class. The dictionary should contain the following information under certain keys:

    .. csv-table::
        :header: information, default units, supported keys

        molecular name, ,"mol"
        pressure grid, :math:`Pa`,"p, P, pressure, pressure_grid"
        temperature grid, :math:`K`,"t, T, temperature, temperature_grid"
        wavenumber grid, :math:`1/cm`,"wn, wavenumber, wavenumbers, wavenumbers_grid, wavenumber_grid"
        ktable, :math:`m^2/kg`, "ktable"

    If the input data has no units attached, the defaults units are assume.
    If they have units, a simple conversion is performed by the code to match the default values.


    Notes
    -----
    This class is implicitly called by the model when initialised if the matched input is used.

    Examples
    ---------
    First we prepare a data dictionary. Please be aware that this are not representative of any molecule.

    >>> import numpy as np
    >>> data_dict = {'mol': 'None',
    >>>              'mol_mass':42,
    >>>              'pressure': np.array([1.00000000e+00, 1.00000000e+01, 1.00000000e+02, 1.00000000e+03,
    >>>                                   1.00000000e+04, 1.00000000e+05, 1.00000000e+06, 1.00000000e+07]),
    >>>              'temperature': np.array([500, 1000, 1500, 2000, 3000]),
    >>>              'wavenumber': np.array([100000, 1000])}
    >>> ktab = np.ones((data_dict['pressure'].size,data_dict['temperature'].size,data_dict['wavenumber'].size,))
    >>> data_dict['ktable'] = ktab

    Let's build the Planck method using the dictionary as input

    >>> from rapoc import Planck
    >>> model = Planck(input_data=data_dict)

    Now the model is ready to be used
    '''

    def __init__(self, input_data):
        '''

        Parameters
        ----------
        input_data: dict
            input_data dictionary
        '''
        self.input_data = input_data

    def read_content(self):
        '''
        Reads the dict content and returns the needed valued for the opacity models.

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
        global pressure_grid, temperature_grid, wavenumber_grid

        try:
            mol = self.input_data['mol']
        except KeyError:
            raise KeyError('molecule name not found')

        try:
            mol_mass = self.input_data['mol_mass']
        except KeyError:
            from molmass import Formula
            f = Formula(mol)
            mol_mass = (f.mass * u.u / u.ct).to(u.g / u.ct)

        if any(key in self.input_data.keys() for key in ['p', 'P', 'pressure', 'pressure_grid']):
            for press_key in ['p', 'P', 'pressure', 'pressure_grid']:
                try:
                    pressure_grid = self.input_data[press_key]
                except KeyError:
                    pass
            if isinstance(pressure_grid, u.Quantity):
                pressure_grid = pressure_grid.to(u.Pa)
            else:
                pressure_grid *= u.Pa
        else:
            raise KeyError('Pressure grid not found')

        if any(key in self.input_data.keys() for key in ['t', 'T', 'temperature', 'temperature_grid']):
            for temp_key in ['t', 'T', 'temperature', 'temperature_grid']:
                try:
                    temperature_grid = self.input_data[temp_key]
                except KeyError:
                    pass
            if isinstance(temperature_grid, u.Quantity):
                temperature_grid = temperature_grid.to(u.K)
            else:
                temperature_grid *= u.K
        else:
            raise KeyError('Temperature grid not found')

        if any(key in self.input_data.keys()
               for key in ['wn', 'wavenumber', 'wavenumbers', 'wavenumbers_grid', 'wavenumber_grid']):
            for wn_key in ['wn', 'wavenumber', 'wavenumbers', 'wavenumbers_grid', 'wavenumber_grid']:
                try:
                    wavenumber_grid = self.input_data[wn_key]
                except KeyError:
                    pass
            if isinstance(wavenumber_grid, u.Quantity):
                wavenumber_grid = wavenumber_grid.to(1 / u.cm)
            else:
                wavenumber_grid /= u.cm
        else:
            raise KeyError('Wavenumber grid not found')

        try:
            ktable = self.input_data['ktable']
            if isinstance(ktable, u.Quantity):
                ktable = ktable.to(u.m ** 2 / u.kg)
            else:
                ktable *= u.m ** 2 / u.kg
        except KeyError:
            raise KeyError('ktable not found')

        return mol, mol_mass, pressure_grid, temperature_grid, wavenumber_grid, ktable
