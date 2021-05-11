import unittest

import astropy.units as u
import numpy as np

from input import exomol_file
from rapoc.models import Model
from rapoc.models.utils.validation import validate_output

data_dict = {'mol': 'H2O',
             'pressure': np.array([1.00000000e+00, 2.15443469e+00, 4.64158883e+00, 1.00000000e+01,
                                   2.15443469e+01, 4.64158883e+01, 1.00000000e+02, 2.15443469e+02,
                                   4.64158883e+02, 1.00000000e+03, 2.15443469e+03, 4.64158883e+03,
                                   1.00000000e+04, 2.15443469e+04, 4.64158883e+04, 1.00000000e+05,
                                   2.15443469e+05, 4.64158883e+05, 1.00000000e+06, 2.15443469e+06,
                                   4.64158883e+06, 1.00000000e+07]),
             'temperature': np.array([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
                                      1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800,
                                      3000, 3200, 3400]),
             'wavenumber': np.array([100000, 1000])}
ktab = np.ones((data_dict['pressure'].size, data_dict['temperature'].size, data_dict['wavenumber'].size,))
data_dict['opacities'] = ktab


class ModelTest(unittest.TestCase):
    file_model = Model(input_data=exomol_file)
    dict_model = Model(input_data=data_dict)

    def test_mol(self):
        self.assertEqual(self.file_model.mol, 'H2O')
        self.assertEqual(self.dict_model.mol, 'H2O')

    def test_mol_mass(self):
        from molmass import Formula
        f = Formula("H2O")
        mol_mass = (f.mass * u.u / u.ct).to(u.g / u.ct)
        self.assertEqual(self.file_model.mol_mass, mol_mass)
        self.assertEqual(self.dict_model.mol_mass, mol_mass)

    def test_pressure_grid(self):
        self.assertEqual(self.file_model.pressure_grid.unit, u.Pa)
        self.assertEqual(self.dict_model.pressure_grid.unit, u.Pa)

        p_content = [1.00000000e+00, 2.15443469e+00, 4.64158883e+00, 1.00000000e+01,
                     2.15443469e+01, 4.64158883e+01, 1.00000000e+02, 2.15443469e+02,
                     4.64158883e+02, 1.00000000e+03, 2.15443469e+03, 4.64158883e+03,
                     1.00000000e+04, 2.15443469e+04, 4.64158883e+04, 1.00000000e+05,
                     2.15443469e+05, 4.64158883e+05, 1.00000000e+06, 2.15443469e+06,
                     4.64158883e+06, 1.00000000e+07]
        np.testing.assert_almost_equal(self.file_model.pressure_grid.value, np.array(p_content), verbose=True)
        np.testing.assert_almost_equal(self.dict_model.pressure_grid.value, np.array(p_content), verbose=True)

    def test_temperature_grid(self):
        self.assertEqual(self.file_model.temperature_grid.unit, u.K)
        self.assertEqual(self.dict_model.temperature_grid.unit, u.K)

        t_content = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
                     1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800,
                     3000, 3200, 3400]
        self.assertListEqual(list(self.file_model.temperature_grid.value), t_content)
        self.assertListEqual(list(self.dict_model.temperature_grid.value), t_content)

    def test_wavenumber_grid(self):
        self.assertEqual(self.file_model.wavenumber_grid.unit, 1 / u.cm)
        self.assertEqual(self.file_model.frequency_grid.unit, 1 / u.s)
        self.assertEqual(self.file_model.wavelength_grid.unit, u.um)

        self.assertEqual(self.dict_model.wavenumber_grid.unit, 1 / u.cm)
        self.assertEqual(self.dict_model.frequency_grid.unit, 1 / u.s)
        self.assertEqual(self.dict_model.wavelength_grid.unit, u.um)

    def test_opacity(self):
        self.assertEqual(self.file_model.opacities.unit, u.m ** 2 / u.kg)
        self.assertEqual(self.dict_model.opacities.unit, u.m ** 2 / u.kg)


class OutputValidationTest(unittest.TestCase):

    def test_validate_outputs_units(self):
        r = 0.033092609716294547 * u.m ** 2 / u.kg

        r = validate_output(r, units='si')
        self.assertAlmostEqual(r.value, 0.033092609716294547)

        r = validate_output(r, units=u.m ** 2 / u.kg)
        self.assertAlmostEqual(r.value, 0.033092609716294547)

        r = validate_output(r, units='cgs')
        self.assertAlmostEqual(r.value, 0.33092609716294547)

        with self.assertRaises(ValueError):
            r = validate_output(r, units=u.m / u.kg)

        with self.assertRaises(ValueError):
            r = validate_output(r, units='cgcc')

    def test_output_values(self):
        r = -1 * u.m ** 2 / u.kg
        r = validate_output(r)
        self.assertEqual(r.value, 0.0)

        r = np.array([1, 0, -1]) * u.m ** 2 / u.kg
        r = validate_output(r)
        self.assertEqual(list(r.value), [1, 0, 0])

        r = 1e-250 * u.m ** 2 / u.kg
        r = validate_output(r)
        self.assertEqual(r.value, 0.0)

        r = np.array([1e-250, 0, -1]) * u.m ** 2 / u.kg
        r = validate_output(r)
        self.assertEqual(list(r.value), [0, 0, 0])

        with self.assertRaises(ValueError):
            r = [1, np.nan] * u.m ** 2 / u.kg
            r = validate_output(r)

        with self.assertRaises(ValueError):
            r = np.nan * u.m ** 2 / u.kg
            r = validate_output(r)
