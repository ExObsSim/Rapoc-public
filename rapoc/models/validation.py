import warnings

import numpy as np
from astropy import units as u


def _validate_inputs_units(p_input, t_input):
    '''
    It checks the units in input pressure an temperature, and return them into `SI`.

    Parameters
    ----------
    p_input: float or np.array
            pressure to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `Pa` are assumed as units.
            Can either be a single value or an array.
    t_input: float or np.array
            temperature to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `K` are assumed as units.
            Can either be a single value or an array.

    Returns
    -------
    float or np.array
        p_input
    float or np.array
        t_input

    Raises
    ------
    astropy.units.UnitConversionError
        if invalid units are used
    '''

    if not hasattr(p_input, 'unit'):
        p_input *= u.Pa
    else:
        try:
            p_input = p_input.to(u.Pa).si
        except u.UnitConversionError:
            raise u.UnitConversionError('invalid pressure unit')
    if not hasattr(t_input, 'unit'):
        t_input *= u.K
    else:
        try:
            t_input = t_input.to(u.K)
        except u.UnitConversionError:
            raise u.UnitConversionError('invalid temperature unit')
    return p_input, t_input


def _validate_inputs_size(p_input, t_input):
    '''
    It checks the sizes of input pressure an temperature. If are just float, it returns arrays instead.

    Parameters
    ----------
    p_input: float or np.array
            pressure to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `Pa` are assumed as units.
            Can either be a single value or an array.
    t_input: float or np.array
            temperature to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `K` are assumed as units.
            Can either be a single value or an array.

    Returns
    -------
    float or np.array
        p_input
    float or np.array
        t_input
    '''

    if p_input.size == 1:
        p_input = [p_input.value] * p_input.unit
    if t_input.size == 1:
        t_input = [t_input.value] * t_input.unit
    return p_input, t_input


def _validate_inputs_boundaries(p_input, t_input, pressure_grid, temperature_grid):
    '''
    It checks the input pressure an temperature are included the data grid boundaries.

    Parameters
    ----------
    p_input: float or np.array
            pressure to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `Pa` are assumed as units.
            Can either be a single value or an array.
    t_input: float or np.array
            temperature to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `K` are assumed as units.
            Can either be a single value or an array.
    pressure_grid: astropy.units.Quantity
        data pressure grid in `si` units
    temperature_grid: astropy.units.Quantity
        data temperature grid in `si` units

    Raises
    ------
    ValueError
        if input pressure or temperature are outside the data grid.
    '''

    for p in p_input:
        if p > pressure_grid.max() or p < pressure_grid.min():
            print(p, pressure_grid.max(), pressure_grid.min())
            raise ValueError('Input pressure out of pressure data grid')
    for t in t_input:
        if t > temperature_grid.max() or t < temperature_grid.min():
            raise ValueError('Input temperature out of temperature data grid')


def validate_inputs(p_input, t_input, pressure_grid, temperature_grid):
    '''
    Function to validate the input data before the model estimation.

    Parameters
    ----------
    p_input: float or np.array
            pressure to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `Pa` are assumed as units.
            Can either be a single value or an array.
    t_input: float or np.array
            temperature to use for the estimation.
            This should be a :class:`astropy.units.Quantity` with specified units, if not `K` are assumed as units.
            Can either be a single value or an array.
    pressure_grid: astropy.units.Quantity
        data pressure grid in `si` units
    temperature_grid: astropy.units.Quantity
        data temperature grid in `si` units

    Returns
    -------
    np.array
        pressure to use for the estimation.
    np.array
        temperature to use for the estimation.

    Raises
    -------
        ValueError
            if input temperature of pressure are out of the data grid
        u.UnitConversionError:
            if the input units cannot be converted to `si`
    '''
    p_input, t_input = _validate_inputs_units(p_input, t_input)
    p_input, t_input = _validate_inputs_size(p_input, t_input)
    _validate_inputs_boundaries(p_input, t_input, pressure_grid, temperature_grid)
    return p_input, t_input


def _validate_output_size(estimate):
    '''
    If the input estimate is an array of a single element, it return single number instead.

    Parameters
    ----------
    estimate:  astropy.units.Quantity
        estimated opacity

    Returns
    -------
    astropy.units.Quantity
    '''

    try:
        if estimate.shape[0] == 1 and estimate.shape[1] == 1:
            estimate = estimate[0, 0]
    except IndexError:
        if estimate.shape == (1, 1):
            estimate = estimate[0, 0]
        if estimate.shape == (1,):
            estimate = estimate[0]
    return estimate


def _validate_output_units(estimate, units):
    '''
    It checks the estimate units and convert them into the desired units, if possible.

    Parameters
    ----------
    estimate:  astropy.units.Quantity
        estimated opacity
    units: str or astropy.units, optional
        indicates the output units system.
        If `si` the opacity is returned as :math:`m^2/kg`, if `cgs` is returned as :math:`cm^2/g`.
        Instead of a string, you can also indicate the units using astropy.units. Units is `si` by default.

    Returns
    -------
    astropy.units.Quantity

    Raises
    ------
    astropy.units.UnitConversionError
        if invalid units are indicated
    '''
    try:
        if units == 'si':
            return estimate.si
        elif units == 'cgs':
            return estimate.cgs
        else:
            return estimate.to(units)
    except u.UnitConversionError:
        raise u.UnitConversionError('invalid output units')


def _validate_output_values(estimate):
    '''
    It checks for wrong numbers inside for estimate.
    The code replaces these values with 0.

    Parameters
    ----------
    estimate:  astropy.units.Quantity
        estimated opacity

    Returns
    -------
    astropy.units.Quantity

    Raises
    ------
    ValueError
        if NaNs are inside estimate
    '''

    if estimate.size > 1:
        idx = np.where(estimate.value < 0)
        if idx[0].size != 0:
            warnings.warn('Found negative value in output. Converted to 0.', RuntimeWarning)
            estimate[idx] = 0

        idx = np.where(estimate.value < 1e-200)  # such small values are due to replaced zeros in Rosselan model
        if idx[0].size != 0:
            warnings.warn('Found close to zero value in output. Converted to 0.', RuntimeWarning)
            estimate[idx] = 0

    else:
        if estimate.value < 0:
            warnings.warn('Found negative value in output. Converted to 0.', RuntimeWarning)
            estimate = 0 * estimate.unit

        if estimate.value < 1e-200:  # such small values are due to replaced zeros in Rosselan model
            warnings.warn('Found close to zero value in output. Converted to 0.', RuntimeWarning)
            estimate = 0 * estimate.unit

    if np.isnan(estimate).any():
        raise ValueError('output is Nan')
    return estimate


def validate_output(estimate, units='si'):
    '''
    Function to validate the estimation output.

    Parameters
    ----------
    estimate:  astropy.units.Quantity
            estimated opacity
    units: str or astropy.units, optional
            indicates the output units system.
            If `si` the opacity is returned as :math:`m^2/kg`, if `cgs` is returned as :math:`cm^2/g`.
            Instead of a string, you can also indicate the units using astropy.units. Units is `si` by default.

    Returns
    -------
    astropy.units.Quantity
            returns the estimated opacity for the pressures and temperatures indicated.
            If only one pressure and temperature are indicated, it returns a single Quantity.
            If `n` pressures and `m` temperatures are indicated, it return a Quantity array of dimensions (n,m)

    Raises
    -------
        ValueError
            if output is nan
        u.UnitConversionError:
            if the output units are invalid
    '''

    estimate = _validate_output_units(estimate, units)
    estimate = _validate_output_size(estimate)
    estimate = _validate_output_values(estimate)
    return estimate
