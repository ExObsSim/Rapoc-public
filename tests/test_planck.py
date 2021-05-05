import unittest

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np

from input import exomol_file
from rapoc.models import Planck


class PlanckTest(unittest.TestCase):
    model = Planck(input_data=exomol_file)

    def test_extract_ktable(self):
        k_ext, _ = self.model.extract_opacities(P_input=1 * u.bar, T_input=1000 * u.K, band=(0.3 * u.um, 50 * u.um))

        idw = self.model._wl_range_parser((0.3 * u.um, 50 * u.um))
        idp = np.argmin(abs(self.model.pressure_grid - (1 * u.bar).to(u.Pa)))
        idt = np.argmin(abs(self.model.temperature_grid - (1000 * u.K).to(u.K)))
        k_true = self.model.opacities[idp][idt][idw]
        self.assertListEqual(list(k_true.value), list(k_ext.value))

        P_input = [1, 10] * u.bar
        T_input = [1000, 2000] * u.K
        k_ext, _ = self.model.extract_opacities(P_input=P_input, T_input=T_input, band=(0.3 * u.um, 50 * u.um))

        for i, p in enumerate(P_input):
            for j, t in enumerate(T_input):
                k_true, _ = self.model.extract_opacities(P_input=p, T_input=t, band=(0.3 * u.um, 50 * u.um))
                self.assertListEqual(list(k_true.value), list(k_ext[i, j].value))

    def test_estimate(self):
        r = self.model.estimate(P_input=1 * u.bar, T_input=1000 * u.K, band=(0.3 * u.um, 50 * u.um), mode='closest')
        self.assertAlmostEqual(r.value, 23.0784696099494)

        fig, ax = self.model.estimate_plot(P_input=1 * u.bar, T_input=1000 * u.K, band=(0.3 * u.um, 50 * u.um),
                                           mode='closest')
        plt.show()

        r = self.model.estimate(P_input=[1, 10] * u.bar, T_input=[1000, 2000] * u.K, band=(0.3 * u.um, 50 * u.um),
                                mode='closest')

        fig, ax = self.model.estimate_plot(P_input=[1, 10] * u.bar, T_input=[1000, 2000] * u.K,
                                           band=(0.3 * u.um, 50 * u.um), mode='closest')

    def test_interpolation(self):
        r = self.model.estimate(P_input=1 * u.bar, T_input=1000 * u.K, band=(0.3 * u.um, 50 * u.um), mode='linear',
                                )
        self.assertAlmostEqual(np.float_(r.value), 23.0784696099494, places=5)

        # r = self.model.estimate(P_input=1 * u.bar, T_input=1000 * u.K, band=(0.3 * u.um, 50 * u.um), mode='cubic',
        #                         )
        # self.assertAlmostEqual(np.float_(r.value), 23.0784696099494, places=5)
        r = self.model.estimate(P_input=1 * u.bar, T_input=1000 * u.K, band=(0.3 * u.um, 50 * u.um), mode='loglinear',
                                )
        self.assertAlmostEqual(np.float_(r.value), 23.0784696099494, places=5)

        fig, ax = self.model.estimate_plot(P_input=1 * u.bar, T_input=1000 * u.K, band=(0.3 * u.um, 50 * u.um),
                                           mode='linear', )
        plt.show()

        r = self.model.estimate(P_input=[1, 10] * u.bar, T_input=[1000, 2000] * u.K, band=(0.3 * u.um, 50 * u.um),
                                mode='linear', )
        fig, ax = self.model.estimate_plot(P_input=[1, 10] * u.bar, T_input=[1000, 2000] * u.K,
                                           band=(0.3 * u.um, 50 * u.um), mode='linear', )
        plt.show()
