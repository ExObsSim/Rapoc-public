import warnings

import astropy.constants as const
import astropy.units as u
import numpy as np

from .model import Model


class Rosseland(Model):
    '''
    Rosseland model class. This class implements the Rosseland Opacity model:

    .. math::
         \\frac{1}{k} = \\frac{\\int_{\\nu_1}^{\\nu_2} k_{\\nu}^{-1} u(\\nu, T) d\\nu}
         {\\int_{\\nu_1}^{\\nu_2} u(\\nu, T) d\\nu}

    where :math:`\\nu_1` and :math:`\\nu_2` are the investigated frequency band edges,
    :math:`k_{\\nu}` is the input data opacity at a given frequency,
    and :math:`u(\\nu,T)` is the black body temperature derivative:

    .. math::
        u(\\nu, T) = \\frac{\\partial B(\\nu, T)}{\\partial T} =
        2 \\frac{h^2 \\nu^4}{c^2 k_b} \\frac{1}{T^2} \\frac{e^{\\frac{h \\nu}{k_b T}}}{(e^{\\frac{h \\nu}{k_b T}} -1)^2}


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
        Let's build the Rosseland method using an Exomol file as input

        >>> from rapoc import Rosseland
        >>> exomolFile = 'example.TauREx.h5'
        >>> model = Rosseland(input_data=exomolFile)

        Now the model is ready to be used
        '''
        super().__init__(input_data)
        self.model_name = 'Rosseland'

    def opacity_model(self, ktable, nu, T_input):
        """
        This function computes the Rosseland Opacity model:

        .. math::
             \\frac{1}{k} = \\frac{\\int_{\\nu_1}^{\\nu_2} k_{\\nu}^{-1} u(\\nu, T) d\\nu}
             {\\int_{\\nu_1}^{\\nu_2} u(\\nu, T) d\\nu}

        where :math:`\\nu_1` and :math:`\\nu_2` are the investigated frequency band edges,
        :math:`k_{\\nu}` is the input data opacity at a given frequency,
        and :math:`u(\\nu,T)` is the black body temperature derivative:

        .. math::
            u(\\nu, T) = \\frac{\\partial B(\\nu, T)}{\\partial T} =
            2 \\frac{h^2 \\nu^4}{c^2 k_b} \\frac{1}{T^2} \\frac{e^{\\frac{h \\nu}{k_b T}}}{(e^{\\frac{h \\nu}{k_b T}} -1)^2}


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
        return self._ross(ktable, nu, T_input).si

    def _u_nu(self, nu, T):
        first_term = const.h ** 2 * nu ** 4 / const.c ** 2 / const.k_B
        second_term = 1 / T ** 2
        exp_term = np.exp(const.h * nu / const.k_B / T)
        third_term = exp_term / (exp_term - 1) ** 2
        return 2 * first_term * second_term * third_term

    def _ross(self, xsec, nu, T):
        i = np.where(xsec == 0)[0]
        if i.size == xsec.size:
            warnings.warn('Ktable is zero. A zero result has been forced', RuntimeWarning)
            return 0.0 * u.m ** 2 / u.kg  # This is physically motivated
        else:
            xsec[i] = 1e-300 * xsec.unit  # for numerical reasons we replace the zeros with a very small value.
            warnings.warn('Ktable contains zeros. They have been replaced with 1e-300', RuntimeWarning)

        top = np.trapz(1 / xsec * self._u_nu(nu, T), x=nu)
        bottom = np.trapz(self._u_nu(nu, T), x=nu)
        #     top = np.sum( 1/xsec * u_nu(nu, T))
        #     bottom = np.sum( u_nu(nu, T))
        return bottom / top
