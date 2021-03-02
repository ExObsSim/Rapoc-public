import astropy.constants as const
import numpy as np

from .model import Model


class Planck(Model):
    '''
    Planck model class. This class implements the Planck Opacity model:

    .. math::
         k = \\frac{\\int_{\\nu_1}^{\\nu_2} k_{\\nu} B_\\nu(T) d\\nu}
         {\\int_{\\nu_1}^{\\nu_2} B_\\nu(T) d\\nu}

    where :math:`\\nu_1` and :math:`\\nu_2` are the investigated frequency band edges,
    :math:`k_{\\nu}` is the input data opacity at a given frequency,
    and :math:`B_\\nu(T)` is the black body radiation energy distribution computed at a certain temperature:

    .. math::
        B_\\nu(T) = \\frac{2h \\nu^3}{c^2} \\frac{1}{e^{\\frac{h\\nu}{k_B T}-1}}


    Attributes
    ----------
    input_data: str or dict
        data file name or data dictionary
    model_name: str
        model name
    band: tuple
        investigated band. Is initialised only after the use of :func:`~self.map`
    selected_grid: astropy.units.Quantity
        investigated grid. Is initialised only after the use of :func:`~self.map`
    mol: str
        molecule name
    mol_mass: astropy.units.Quantity
        molecular mass
    pressure_grid: astropy.units.Quantity
        data pressure grid in `si` units
    temperature_grid: astropy.units.Quantity
        data temperature grid in `si` units
    wavenumber_grid: astropy.units.Quantity
        data wavenumber grid
    ktable: astropy.units.Quantity
        data opacities grid in `si` units
    frequency_grid: astropy.units.Quantity
        data frequency grid
    wavelength_grid: astropy.units.Quantity
        data wavelength grid

    '''

    def __init__(self, input_data):
        '''
        Parameters
        ----------
        input_data: str or dict
            data file name or input_data dictionary.
            If the input is a `str`, the correspondent loader is used (:class:`~rapoc.loaders.ExoMolFileLoader`)
            if the file format is supported.
            If is `dict`, then :class:`~rapoc.loaders.DictLoader` is used.

        Raises
        -------
        IOError: if data format is not supported

        Examples
        ---------
        First we prepare a data dictionary:

        >>> import numpy as np
        >>> data_dict = {'mol': 'H2O',
        >>>              'pressure': np.array([1.00000000e+00, 1.00000000e+01, 1.00000000e+02, 1.00000000e+03,
        >>>                                   1.00000000e+04, 1.00000000e+05, 1.00000000e+06, 1.00000000e+07]),
        >>>              'temperature': np.array([500, 1000, 1500, 2000, 3000]),
        >>>              'wavenumber': np.array([100000, 1000])}
        >>> ktab = np.ones((data_dict['pressure'].size,data_dict['temperature'].size,data_dict['wavenumber'].size,))
        >>> data_dict['ktable'] = ktab

        Let's build the Planck method using an Exomol file as input

        >>> from rapoc import Planck
        >>> dict_data = 'example.TauREx.h5'
        >>> model = Planck(input_data=data_dict)

        Now the model is ready to be used
        '''
        super().__init__(input_data)
        self.model_name = 'Planck'

    def opacity_model(self, ktable, nu, T_input):
        """
        This function computes the Planck Opacity model:

        .. math::
             k = \\frac{\\int_{\\nu_1}^{\\nu_2} k_{\\nu} B_\\nu(T) d\\nu}
             {\\int_{\\nu_1}^{\\nu_2} B_\\nu(T) d\\nu}

        where :math:`\\nu_1` and :math:`\\nu_2` are the investigated frequency band edges,
        :math:`k_{\\nu}` is the input data opacity at a given frequency,
        and :math:`B_\\nu(T)` is the black body radiation energy distribution computed at a certain temperature:

        .. math::
            B_\\nu(T) = \\frac{2h \\nu^3}{c^2} \\frac{1}{e^{\\frac{h\\nu}{k_B T}-1}}

        Parameters
        ----------
        ktable: np.array
            opacity array. Has the same dimension of the frequency grid `nu`
        nu: np.array
            frequency grid
        T_input: float
            temperature

        Returns
        -------
        float
            mean opacity computed from the model
        """
        return self._planck(ktable, nu, T_input).si

    def _bb_nu(self, nu, T):
        first_term = const.h * nu ** 3 / const.c ** 2
        exp_term = np.exp(const.h * nu / const.k_B / T)
        second_term = 1 / (exp_term - 1)
        return 2 * first_term * second_term

    def _planck(self, xsec, nu, T):
        top = np.trapz(xsec * self._bb_nu(nu, T), x=nu)
        bottom = np.trapz(self._bb_nu(nu, T), x=nu)
        return top / bottom
