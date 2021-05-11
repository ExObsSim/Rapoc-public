import unittest

import astropy.units as u
from input import exomol_file, dace_file
from rapoc.models import Rosseland, Planck

class DataTest(unittest.TestCase):

    def test_rosseland(self):
        ross_exo = Rosseland(input_data=exomol_file)
        ross_dace = Rosseland(input_data=dace_file)

        T = 1500*u.K
        P = 1*u.bar
        band = (0.3*u.um, 0.7*u.um)

        r_e = ross_exo.estimate(P_input=P, T_input=T, band=band)
        r_d = ross_dace.estimate(P_input=P, T_input=T, band=band)
        print (1-r_e/r_d)

        T = 1700*u.K
        r_e = ross_exo.estimate(P_input=P, T_input=T, band=band)
        r_d = ross_dace.estimate(P_input=P, T_input=T, band=band)
        print (1-r_e/r_d)

    def test_planck(self):
        planck_exo = Planck(input_data=exomol_file)
        plack_dace = Planck(input_data=dace_file)

        T = 1500 * u.K
        P = 1 * u.bar
        band = (0.3 * u.um, 0.7 * u.um)

        r_e = planck_exo.estimate(P_input=P, T_input=T, band=band)
        r_d = plack_dace.estimate(P_input=P, T_input=T, band=band)
        print(1-r_e/r_d)

        T = 1700 * u.K
        r_e = planck_exo.estimate(P_input=P, T_input=T, band=band)
        r_d = plack_dace.estimate(P_input=P, T_input=T, band=band)
        print(1-r_e/r_d)