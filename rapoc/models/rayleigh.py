import astropy.units as u
import numpy as np
import astropy.constants as const


class Rayleigh:
    '''
    The Rayleigh class estimates the Rayleigh scattering for the indicated atom using the
    equations from Chapter 10 of David R. Lide, ed., CRC Handbook of Chemistry and Physics, Internet Version 2005,
    `hbcponline <http://hbcponline.com/faces/contents/ContentsSearch.xhtml>`_, CRC Press, Boca Raton, FL, 2005.
    These values won't depend on pressure or temperature and therefore only a wavenumbers grid is needed to produce them.
    The wavenumbers grid can be given by the user or loaded from an already initialised
    :class:`~rapoc.models.model.Model` class (as :class:`~rapoc.Rosseland` or :class:`~rapoc.PLanck`).
    The latter solution might be of help in using the Rayleigh scattering along with other opacities.
    '''

    def __init__(self, atom, wavenumber_grid=None, model=None):
        '''
        Parameters
        ----------
        atom: str
            name of the considered atom
        wavenumber_grid: list or numpy.array or astropy.units.Quantity
            data wavenumber grid
        model: :class:`~rapoc.models.model.Model`
            built model to use to load the wavenumbers grid.


        Examples
        ---------
        First we prepare a data wavenumbers grid, then we can produce the Rayleigh data:

        >>> rayleigh = Rayleigh('H', wavenumber_grid=[100000, 10000, 1000])


        if we already have a molecular model loaded, as example built from a datafile,
        we can use it to initialise the Rayleigh class:

        >>> input_data = Rosseland(input_data=exomol_file)
        >>> rayleigh = Rayleigh(atom='Na', model=input_data)


        Now the Rayleigh data are ready to be used as input data for any RAPOC :class:`~rapoc.models.model.Model` class:

        >>> pl = Planck(input_data=rayleigh)
        >>> rs = Rosseland(input_data=rayleigh)

        '''
        self.atom = atom
        self.atom_mass = self._get_atommass()
        if model:
            self.wavenumber_grid = model.wavenumber_grid
            self.wavelength_grid = model.wavelength_grid
        else:
            self.wavenumber_grid, self.wavelength_grid = self._get_wavegrid(wavenumber_grid)
        self.pressure_grid, self.temperature_grid = np.array([1e-10, 1e100]) * u.Pa, np.array([1, 1e6]) * u.K
        self.opacities = self.compute_opacities()

    def _get_atommass(self):
        from molmass import Formula
        f = Formula(self.atom)
        return (f.mass * u.u / u.ct).to(u.g / u.ct)

    def _get_wavegrid(self, wavenumber_grid):
        if isinstance(wavenumber_grid, u.Quantity):
            wavenumber_grid = wavenumber_grid.to(1 / u.cm)
        else:
            wavenumber_grid /= u.cm
        wavelength_grid = 1 / wavenumber_grid.to(1 / u.um)

        return wavenumber_grid, wavelength_grid

    def compute_opacities(self):
        """
        It computes the opacities for Rayleigh scattering as described in David R. Lide, ed.,
        CRC Handbook of Chemistry and Physics, Internet Version 2005,
        `hbcponline <http://hbcponline.com/faces/contents/ContentsSearch.xhtml>`_, CRC Press, Boca Raton, FL, 2005.
        chapter 10 table 1:

            .. math::
                \\alpha(\\lambda) = \\frac{128 \\pi^5}{3 \\lambda^4} \\frac{a^2}{m}

        where :math:`a` is depending on the atom polarizabiility listed in table 2 chapter 10 of CRC handbook,
        and :math:`m` is the atom mass.

        Returns
        -------
        astropy.units.Quantity
            returns the estimated opacity sampled at the indicated wavenumbers.
        """

        a = a_table[self.atom] * (1e-24 * u.cm ** 3)
        alpha = 128 * np.pi ** 5 / (3 * self.wavelength_grid ** 4) * a ** 2
        alpha /= self.atom_mass * u.ct
        return alpha.si


    def read_content(self):
        '''
        Reads the class content and returns the needed valued for the opacity models.

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
        opac = np.ones((self.pressure_grid.size, self.temperature_grid.size, self.wavenumber_grid.size,))
        opac *= u.m ** 2 / u.kg
        opac[0, 0] = self.opacities
        return self.atom, self.atom_mass, self.pressure_grid, self.temperature_grid, self.wavenumber_grid, opac, 'rayleigh'


# from table 2 sec 10 of CRC
a_table = {
    'H': 0.666793,
    'He': 0.20496,
    "Li": 24.3,
    "Be": 5.6,
    "B": 3.03,
    "C": 1.76,
    "N": 1.1,
    "O": 0.802,
    "F": 0.557,
    "Ne": 0.3956,
    "Na": 24.11,
    "Mg": 10.6,
    "Al": 6.8,
    "Si": 5.38,
    "P": 3.63,
    "S": 2.9,
    "Cl": 2.18,
    "Ar": 1.6411,
    "K": 43.4,
    "Ca": 22.8,
    "Sc": 17.8,
    "Ti": 14.6,
    "V": 12.4,
    "Cr": 11.6,
    "Mn": 9.4,
    "Fe": 8.4,
    "Co": 7.5,
    "Ni": 6.8,
    "Cu": 6.2,
    "Zn": 5.75,
    "Ga": 8.12,
    "Ge": 6.07,
    "As": 4.31,
    "Se": 3.77,
    "Br": 3.05,
    "Kr": 2.4844,
    "Rb": 47.3,
    "Sr": 27.6,
    "Y": 22.7,
    "Zr": 17.9,
    "Nb": 15.7,
    "Mo": 12.8,
    "Tc": 11.4,
    "Ru": 9.6,
    "Rh": 8.6,
    "Pd": 4.8,
    "Ag": 7.2,
    "Cd": 7.36,
    "In": 10.2,
    "Sn": 7.7,
    "Sb": 6.6,
    "Te": 5.5,
    "I": 5.35,
    "Xe": 4.044,
    "Cs": 59.42,
    "Ba": 39.7,
    "La": 31.1,
    "Ce": 29.6,
    "Pr": 28.2,
    "Nd": 31.4,
    "Pm": 30.1,
    "Sm": 28.8,
    "Eu": 27.7,
    "Gd": 23.5,
    "Tb": 25.5,
    "Dy": 24.5,
    "Ho": 23.6,
    "Er": 22.7,
    "Tm": 21.8,
    "Yb": 21,
    "Lu": 21.9,
    "Hf": 16.2,
    "Ta": 13.1,
    "W": 11.1,
    "Re": 9.7,
    "Os": 8.5,
    "Ir": 7.6,
    "Pt": 6.5,
    "Au": 5.8,
    "Hg": 5.02,
    "Tl": 7.6,
    "Pb": 6.8,
    "Bi": 7.4,
    "Po": 6.8,
    "At": 6,
    "Rn": 5.3,
    "Fr": 47.1,
    "Ra": 38.3,
    "Ac": 32.1,
    "Th": 32.1,
    "Pa": 25.4,
    "U": 24.9,
    "Np": 24.8,
    "Pu": 24.5,
    "Am": 23.3,
    "Cm": 23,
    "Bk": 22.7,
    "Cf": 20.5,
    "Es": 19.7,
    "Fm": 23.8,
    "Md": 18.2,
    "No": 17.5,
}
