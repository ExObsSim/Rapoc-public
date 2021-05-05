import unittest

import astropy.units as u
import numpy as np

from input import exomol_file
from rapoc.loaders import ExoMolFileLoader


class ExoMolLoaderTest(unittest.TestCase):
    loaded = ExoMolFileLoader(filename=exomol_file)
    mol, mol_mass, pressure_grid, temperature_grid, wavenumber_grid, ktable, force_t = loaded.read_content()

    def test_mol(self):
        self.assertEqual(self.mol, 'H2O')

    def test_mol_mass(self):
        from molmass import Formula
        f = Formula("H2O")
        mol_mass = (f.mass * u.u / u.ct).to(u.g / u.ct)
        self.assertEqual(self.mol_mass, mol_mass)

    def test_pressure_grid(self):
        self.assertEqual(self.pressure_grid.unit, u.Pa)
        p_content = [1.00000000e+00, 2.15443469e+00, 4.64158883e+00, 1.00000000e+01,
                     2.15443469e+01, 4.64158883e+01, 1.00000000e+02, 2.15443469e+02,
                     4.64158883e+02, 1.00000000e+03, 2.15443469e+03, 4.64158883e+03,
                     1.00000000e+04, 2.15443469e+04, 4.64158883e+04, 1.00000000e+05,
                     2.15443469e+05, 4.64158883e+05, 1.00000000e+06, 2.15443469e+06,
                     4.64158883e+06, 1.00000000e+07]
        np.testing.assert_almost_equal(self.pressure_grid.value, np.array(p_content), verbose=True)

    def test_temperature_grid(self):
        self.assertEqual(self.temperature_grid.unit, u.K)
        t_content = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
                     1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800,
                     3000, 3200, 3400]
        self.assertListEqual(list(self.temperature_grid.value), t_content)

    def test_wavenumber_grid(self):
        self.assertEqual(self.wavenumber_grid.unit, 1 / u.cm)

    def test_ktable(self):
        self.assertEqual(self.ktable.unit, u.m ** 2 / u.kg)

# class HitranLoaderTest(unittest.TestCase):
#     loaded = CiaHitranFileLoader(input_data=hitran_file)
#     mol, mol_mass, pressure_grid, temperature_grid, wavenumber_grid, opacities = loaded.read_content()
#
#     def test_mol(self):
#         print(self.mol)
#         self.assertEqual(self.mol, 'H2O')
#
#     def test_mol_mass(self):
#         from molmass import Formula
#         print(self.mol_mass)
#         f = Formula("H2O")
#         mol_mass = (f.mass * u.u / u.ct).to(u.g / u.ct)
#         self.assertEqual(self.mol_mass, mol_mass)
#
#     def test_pressure_grid(self):
#         self.assertEqual(self.pressure_grid, None)
#
#     def test_temperature_grid(self):
#         self.assertEqual(self.temperature_grid.unit, u.K)
#         print(self.temperature_grid)
#         print(self.temperature_grid.shape)
#
#     def test_wavenumber_grid(self):
#         print(self.wavenumber_grid)
#         print(self.wavenumber_grid.shape)
#         self.assertEqual(self.wavenumber_grid.unit, 1 / u.cm)
#
#     def test_ktable(self):
#         print(self.opacities)
#         print(self.opacities.shape)
#         self.assertEqual(self.opacities.unit, u.m ** 2 / u.kg)
