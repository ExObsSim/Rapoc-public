import unittest

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from rapoc.models import Rayleigh
from test_model import data_dict
from rapoc.models import Model


class RayleighTest(unittest.TestCase):

    def test_raw_input(self):

        model = Rayleigh('H', wavenumber_grid=data_dict['wavenumber'])

        np.testing.assert_almost_equal(1.007941, model.atom_mass.to(u.u / u.ct).value, decimal=5)


        np.testing.assert_equal(data_dict['wavenumber'], model.wavenumber_grid.value)

    def test_model_input(self):
        dict_model = Model(input_data=data_dict)
        model = Rayleigh(atom='H', model=dict_model)

        np.testing.assert_almost_equal(1.007941, model.atom_mass.to(u.u / u.ct).value, decimal=5)

        np.testing.assert_equal(data_dict['wavenumber'], model.wavenumber_grid.value)


    def test_inject_model(self):
        dict_model = Model(input_data=data_dict)
        model = Rayleigh(atom='H', model=dict_model)
        ray_model = Model(input_data=model)

        np.testing.assert_equal(data_dict['wavenumber'], ray_model.wavenumber_grid.value)

    def test_compute_opacity(self):
        model = Rayleigh(atom='Na', wavenumber_grid=[100000/u.cm])
        np.testing.assert_almost_equal(19.76, model.opacities.cgs.value, decimal=0)

        model = Rayleigh(atom='K', wavenumber_grid=[2000.00/u.cm])
        np.testing.assert_almost_equal(6, model.opacities.cgs.value * 1e6, decimal=0)

    def test_rosseland_planck(self):
        from rapoc.models import Rosseland, Planck
        dict_model = Model(input_data=data_dict)
        rayleigh = Rayleigh(atom='Na', model=dict_model)

        pl = Planck(input_data=rayleigh)

        est = pl.estimate(T_input=1000*u.K, band=(1000/u.cm, 100000/u.cm) )
        print(est)

        rs = Rosseland(input_data=rayleigh)

        est = rs.estimate(T_input=1000*u.K, band=(1000/u.cm, 100000/u.cm))
        print(est)


    def test_extract_opacity(self):
        from input import exomol_file
        from rapoc.models import Rosseland, Planck
        input_data = Model(input_data=exomol_file)
        rayleigh = Rayleigh(atom='H', model=input_data)
        pl = Planck(input_data=rayleigh)
        rs = Rosseland(input_data=rayleigh)

        T = 1000.0 *u.K
        wl = (1 * u.um, 10 * u.um)

        k_ext, _ = pl.extract_opacities(T_input=T, band=wl)
        print(k_ext)
        #self.assertListEqual(list(rayleigh.opacities[0,0]), list(k_ext.value))

        k_ext, _ = rs.extract_opacities(T_input=T, band=wl)
        #self.assertListEqual(list(rayleigh.opacities[0,0]), list(k_ext.value))
        print(k_ext)


    def test_estimate(self):
        from input import exomol_file
        from rapoc.models import Rosseland, Planck
        input_data = Model(input_data=exomol_file)
        rayleigh = Rayleigh(atom='H', model=input_data)
        pl = Planck(input_data=rayleigh)
        rs = Rosseland(input_data=rayleigh)

        T = 1000.0 *u.K
        wl = (1 * u.um, 10 * u.um)

        est= pl.estimate(T_input=T, band=wl)
        print(est)
        #self.assertListEqual(list(rayleigh.opacities[0,0]), list(k_ext.value))

        est = rs.estimate(T_input=T, band=wl)
        #self.assertListEqual(list(rayleigh.opacities[0,0]), list(k_ext.value))
        print(est)