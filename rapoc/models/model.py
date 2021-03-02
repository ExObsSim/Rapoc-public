import warnings

import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from rapoc.models.validation import validate_inputs, validate_output


class Model:
    """
    Base model class.

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

    """

    def __init__(self, input_data):
        """
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
        For the sake of semplicity let's assume we use an ExoMol file as input

        >>> from rapoc.models import Model
        >>> exomolFile = 'example.TauREx.h5'
        >>> model = Model(input_data=exomolFile)
        """

        self.input_data = input_data
        self.model_name = 'model'
        self.band = None
        self._map = None
        self.selected_grid = None
        loaded = self._loader_selector()
        self.mol, self.mol_mass, self.pressure_grid, self.temperature_grid, self.wavenumber_grid, self.ktable = loaded.read_content()
        self.frequency_grid = (self.wavenumber_grid * const.c).si
        self.wavelength_grid = 1 / self.wavenumber_grid.to(1 / u.um)

    def opacity_model(self, ktable, nu, T_input):
        """
        Computes the mean opacity in the investigated range.

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
        raise NotImplementedError

    def extract_ktable(self, P_input, T_input, band, units='si'):
        """
        Returns the input_data ktable for given pressure, temperature and band.

        Parameters
        ----------
        P_input: float or np.array
            pressure to use for the estimation.
            This should be a `astropy.units.Quantity` with specified units, if not `Pa` are assumed as units.
            Can either be a single value or an array.
        T_input: float or np.array
            temperature to use for the estimation.
            This should be a `astropy.units.Quantity` with specified units, if not `K` are assumed as units.
            Can either be a single value or an array.
        band: tuple
            this should be a tuple of `astropy.units.Quantity` indicating the edges of the band to use
            in the estimation.
            The units used reveal the intention to work in wavelengths, wavenumbers or frequencies.
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
        list
            list of wavenumber (or wavelength or frequency) index relative to the investigated band

        """
        P_input, T_input = validate_inputs(P_input, T_input, self.pressure_grid, self.temperature_grid)
        idw = self._wl_range_parser(band)
        _extract = np.zeros([P_input.size, T_input.size, idw.size]) * u.m ** 2 / u.kg

        for i, p in enumerate(P_input):
            for j, t in enumerate(T_input):
                idp = np.argmin(abs(self.pressure_grid - p.to(u.Pa)))
                idt = np.argmin(abs(self.temperature_grid - t.to(u.K)))
                _extract[i, j] = self.ktable[idp][idt][idw]
        return validate_output(_extract, units), idw

    def estimate(self, P_input, T_input, band, mode='closest', units='si'):
        """
        It estimates the model opacity in the indicated pressure and temperature grid for the indicated band,
        using the :func:`~rapoc.models.model.Model.opacity_model`.

        Parameters
        ----------
        P_input: float or np.array
            pressure to use for the estimation.
            This should be a `astropy.units.Quantity` with specified units, if not `Pa` are assumed as units.
            Can either be a single value or an array.
        T_input: float or np.array
            temperature to use for the estimation.
            This should be a `astropy.units.Quantity` with specified units, if not `K` are assumed as units.
            Can either be a single value or an array.
        band: tuple
            this should be a tuple of `astropy.units.Quantity` indicating the edges of the band to use
            in the estimation.
            The units used reveal the intention to work in wavelengths, wavenumbers or frequencies.
        mode: str, optional
            indicates the estimation mode desired. If `closest` the output will be the opacity computed from the data at
            the nearest pressure on the data grid to the indicated one, and at the nearest temperature on the data grid
            to the indicated one. If {'linear', 'cubic'} it's used the :func:`scipy.interpolate.griddata`
            with the indicated `kind`. 
            To interpolate a map of opacities over the data grid is computed first using :func:`~rapoc.models.model.Model.map`. 
            Mode is `closest` by default.
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
        ------
        NotImplementedError:
            if the indicated mode is not supported
        astropy.units.UnitConversionError:
            if the input units cannot be converted to `si`
        AttributeError:
            if the band cannot be interpreted as wavelength, wavenumber or frequency
        """

        P_input, T_input = validate_inputs(P_input, T_input, self.pressure_grid, self.temperature_grid)

        if mode == 'closest':
            _estimate = np.zeros([P_input.size, T_input.size]) * u.m ** 2 / u.kg
            idw = self._wl_range_parser(band)
            frequency = self.frequency_grid[idw]
            for i, p in enumerate(P_input):
                for j, t in enumerate(T_input):
                    idp = np.argmin(abs(self.pressure_grid - p.to(u.Pa)))
                    idt = np.argmin(abs(self.temperature_grid - t.to(u.K)))
                    _estimate[i, j] = self.opacity_model(self.ktable[idp][idt][idw], frequency,
                                                         self.temperature_grid[idt])

        elif mode in ['linear', 'cubic']:
            from scipy import interpolate
            _map = self.map(band).to(u.m ** 2 / u.kg).value

            points, map_values = [], []
            for i in range(self.pressure_grid.size):
                for j in range(self.temperature_grid.size):
                    points.append((self.pressure_grid[i].value, self.temperature_grid[j].value))
                    map_values.append(_map[i, j])
            _estimate = np.zeros([P_input.size, T_input.size]) * u.m ** 2 / u.kg
            for i, p in enumerate(P_input.value):
                for j, t in enumerate(T_input.value):
                    _estimate[i, j] = interpolate.griddata(points, np.array(map_values), (p, t),
                                                           method=mode, rescale=True) * u.m ** 2 / u.kg
        else:
            raise NotImplementedError('The indicated mode is not supported.')

        return validate_output(_estimate, units)

    def estimate_plot(self, P_input, T_input, band, mode='closest', force_map=False,
                      fig=None, ax=None, yscale='log', xscale='linear', grid='wavelength'):
        """
        Produces a plot of the estimates (produced with :func:`~rapoc.models.model.Model.estimate`),
        comparing the raw data with the result.
        If a single estimate is to be plotted, this method produces a plot of raw opacities vs wavelength
        from the ktable (in grey) and the mean opacity estimated.
        If multiple estimates are be plotted, it produces a 3D plot,  with the surface of mean opacity vs
        temperature and pressure from the data grid (using :func:`~rapoc.models.model.Model.map_plot`) 
        and the interpolated data superimposed.

        Parameters
        ----------
        P_input: float or np.array
            pressure to use for the estimation.
            This should be a `astropy.units.Quantity` with specified units, if not `Pa` are assumed as units.
            Can either be a single value or an array.
        T_input: float or np.array
            temperature to use for the estimation.
            This should be a `astropy.units.Quantity` with specified units, if not `K` are assumed as units.
            Can either be a single value or an array.
        band: tuple
            this should be a tuple of `astropy.units.Quantity` indicating the edges of the band to use
            in the estimation.
            The units used reveal the intention to work in wavelengths, wavenumbers or frequencies.
        mode: str, optional
            indicates the estimation mode desired. If `closest` the output will be the opacity computed from the data at
            the nearest pressure on the data grid to the indicated one, and at the nearest temperature on the data grid
            to the indicated one. If {'linear', 'cubic'} it's used the :func:`scipy.interpolate.griddata`
            with the indicated `kind`.
            Mode is `closest` by default.
        force_map: bool
            If True a 3D map plot is generate even for a single estimate. Default is `False`.
        fig: matplotlib.pyplot.figure, optional
            indicates the figure to use. If `None` a new figure is produced. Default is `None`
        ax: matplotlib.pyplot.axis, optional
            indicates the axis of the figure to use. If `None` a new axis is produced. Default is `None`
        yscale: str
            y-axis scale to use. Default is `log`
        xscale: str
            x-axis scale to use. Default is `linear`
        grid: str
            x-axis grid format. {`wavelength`, `frequency`, `wavenumber`} are supported. Default is `wavelength`.

        Returns
        -------
        matplotlib.pyplot.figure, optional
            figure containing the plot
        matplotlib.pyplot.axis, optional
            axis containing the plot

        Raises
        ------
        KeyError:
            if the indicated grid is not supported
        """
        r = self.estimate(P_input, T_input, band, mode=mode)
        P_input, T_input = validate_inputs(P_input, T_input, self.pressure_grid, self.temperature_grid)
        if r.size == 1 and not force_map:
            idw = self._wl_range_parser(band)
            if grid == 'wavelength':
                grid_plot = self.wavelength_grid[idw]
            elif grid == 'frequency':
                grid_plot = self.frequency_grid[idw]
            elif grid == 'wavenumber':
                grid_plot = self.wavenumber_grid[idw]
            else:
                raise KeyError('not supported grid')

            if (fig is None) and (ax is None):
                fig, ax = plt.subplots(1, 1)
                ax.set_xlabel(r'{} [${}$]'.format(grid, grid_plot.unit))
                ax.set_ylabel(r'${}$'.format(self.ktable.unit))
                ax.set_title(self.model_name)

            idp = np.argmin(abs(self.pressure_grid - P_input.to(u.Pa)))
            idt = np.argmin(abs(self.temperature_grid - T_input.to(u.K)))
            ax.plot(grid_plot, self.ktable[idp][idt][idw], label='ktable', c='grey', alpha=0.2)
            ax.axhline(r.value, lw=1, label=self.model_name)

            ax.set_yscale(yscale)
            ax.set_xscale(xscale)
            return fig, ax

        else:
            fig, ax = self.map_plot(band, fig=fig, ax=ax)
            T_input, P_input = np.meshgrid(T_input, P_input)
            ax.scatter(T_input, np.log10(P_input.to(u.Pa).value), r, )
            return fig, ax

    def map(self, band):
        """
        Returns the mean opacity map for the indicated model.

        Parameters
        ----------
        band: tuple
            this should be a tuple of `astropy.units.Quantity` indicating the edges of the band to use
            in the estimation.
            The units used reveal the intention to work in wavelengths, wavenumbers or frequencies.

        Returns
        -------
        np.array
            mean opacity map in `si` units

        Notes
        -----
        If the map has been already built for the indicated band, then the method returns the previous map.
        This speeds up the use of the code in external functions.
        """
        if self.band is not None:
            band_check = any(x != y for x, y in zip(band, self.band))
        else:
            band_check = True
        if band_check:
            warnings.warn('Computing map.', UserWarning)
            self.band = band
            idw = self._wl_range_parser(band)
            frequency = self.frequency_grid[idw]
            _map = np.zeros([self.pressure_grid.size, self.temperature_grid.size]) * u.m ** 2 / u.kg
            for i in range(self.pressure_grid.size):
                for j in range(self.temperature_grid.size):
                    _map[i, j] = self.opacity_model(self.ktable[i][j][idw], frequency, self.temperature_grid[j])
            _map = validate_output(_map)
            self._map = _map
        return self._map

    def map_plot(self, band, fig=None, ax=None):
        """
        Produces a 3D-plot of the model mean opacity map built with :func:`~rapoc.models.model.Model.map`.

        Parameters
        ----------
        band: tuple
            this should be a tuple of `astropy.units.Quantity` indicating the edges of the band to use
            in the estimation.
            The units used reveal the intention to work in wavelengths, wavenumbers or frequencies.
        fig: matplotlib.pyplot.figure, optional
            indicates the figure to use. If `None` a new figure is produced. Default is `None`
        ax: matplotlib.pyplot.axis, optional
            indicates the axis of the figure to use. If `None` a new axis is produced. Default is `None`

        Returns
        -------
        matplotlib.pyplot.figure, optional
            figure containing the plot
        matplotlib.pyplot.axis, optinal
            axis containing the plot

        """
        from matplotlib import cm
        t, p = np.meshgrid(self.temperature_grid, self.pressure_grid)

        if (fig is None) and (ax is None):
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_ylabel('log(Pressure [Pa])')
            ax.set_xlabel('Temperature [K]')
            ax.set_zlabel(r'${}$'.format(self.ktable.unit))
        ax.set_title(self.model_name)
        if not hasattr(ax, 'plot_surface'):
            raise AttributeError('input axes should have 3D projection. '
                                 'Please refer to https://matplotlib.org/3.1.1/gallery/mplot3d/subplot3d.html')
        ax.plot_surface(t, np.log10(p.value), self.map(band).to(u.m ** 2 / u.kg), cmap=cm.coolwarm, alpha=0.7)

        return fig, ax

    def _loader_selector(self):

        if isinstance(self.input_data, str):
            if 'TauREx.h5' in self.input_data:
                from rapoc.loaders import ExoMolFileLoader
                return ExoMolFileLoader(filename=self.input_data)
            else:
                raise IOError('file extension not supported')
        elif isinstance(self.input_data, dict):
            from rapoc.loaders import DictLoader
            return DictLoader(input_data=self.input_data)
        else:
            raise IOError('Data format not supported')

    def _wl_range_parser(self, investigated_range):
        try:
            investigated_range[0].to(u.um)
            lookup = self.wavelength_grid
        except u.UnitConversionError:
            try:
                investigated_range[0].to(1 / u.um)
                lookup = self.wavenumber_grid
            except u.UnitConversionError:
                try:
                    investigated_range[0].to(u.Hz)
                    lookup = self.frequency_grid
                except u.UnitConversionError:
                    raise AttributeError('wrong range format')

        if investigated_range[0] > investigated_range[1]:
            investigated_range = reversed(investigated_range)
        idx = np.where(np.logical_and(lookup < investigated_range[1], lookup > investigated_range[0]))[0]
        self.selected_grid = lookup[idx]
        return idx
