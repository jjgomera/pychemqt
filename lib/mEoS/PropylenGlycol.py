#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from unittest import TestCase

from lib import unidades
from lib.meos import MEoS


# TODO: Add plot with extrapolation lines and range of validity, to show in
# docummentation. Maybe I've delete the script...


class PropylenGlycol(MEoS):
    """Multiparameter equation of state for propylene glycol"""
    name = "Propylene glycol"
    CASNumber = "57-55-6"
    formula = "CH3-CH(OH)-CH2OH"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(339.3811132)
    Tc = unidades.Temperature(674)
    Pc = unidades.Pressure(7291.8, "kPa")
    M = 76.09442  # g/mol
    Tt = unidades.Temperature(242.8)
    Tb = unidades.Temperature(461.224)
    f_acent = 0.72
    momentoDipolar = unidades.DipoleMoment(None)
    id = 266

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [1.45359225002898, 2.58396332560320],
           "ao_exp": [5, 28],
           "titao": [1000/Tc, 1330/Tc]}

    eisenbach = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propylen glycol of "
                    "Eisenbach et al. (2021).",
        "__doi__": {"autor": "Eisenbach, T., Scholz, C., Span, R., "
                             "Cristancho, D., Lemmon, E.W., Thol, M.",
                    "title": "Speed-of-Sound Measurements and a Fundamental "
                             "Equation of State for Propylene Glycol",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023105",
                    "doi":  "10.1063/5.0050021"},

        "R": 8.314462618,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 680, "Pmax": 350000, "rhomax": 14.4,

        "nr1":  [0.046611538, 2.0273992, -2.6048664, -0.58592792, 0.2967405,
                 0.053863656],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1.0, 0.14, 0.92, 1.254, 0.425, 0.688],

        "nr2": [-0.078280924, -0.76968025, 0.13016359, -0.015287585,
                -0.010000015, -0.1500221],
        "d2": [1, 3, 2, 7, 1, 2],
        "t2": [1.6, 2.23, 1.55, 0.9, 5.0, 3.0],
        "c2": [2, 2, 1, 1, 3, 2],
        "gamma2": [1]*6,

        "nr3": [-0.24426526, -0.00356737, -0.27150835, 1.2948298, -1.7031454,
                1.7600461, -1.0654478],
        "d3": [2, 2, 1, 1, 1, 1, 1],
        "t3": [1.1, 1.0, 1.5, 2.44, 2.37, 1.77, 2.28],
        "alfa3": [18.7, 18.7, 1.86, 0.63, 0.83, 1.278, 0.45],
        "beta3": [685, 1230, 2.28, 0.13, 0.07, 1.09, 0.13],
        "gamma3": [1.09, 1.04, 1.05, 1.5, 1.43, 1.13, 2.11],
        "epsilon3": [0.789, 0.99, 0.981, 1.004, 0.698, 0.808, 0.81],
        }

    eq = eisenbach,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-10.12, 3.15, -5.60, -0.337, -2.39],
        "t": [1.0, 1.5, 2.6, 4.0, 5.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.46, 2.06, 0.743, -1.905, 1.536],
        "t": [0.21, 0.43, 2.7, 3.7, 4.7]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.0507, -6.8362, -19.835, -10.097, -55.772, -144.55],
        "t": [0.32, 0.9, 2.5, 4.2, 5.7, 12.0]}


class Test(TestCase):

    def test_eisenbach(self):
        # Table 8, Pag 12
        st = PropylenGlycol(T=400, rhom=0.001)
        self.assertEqual(round(st.P.MPa, 10), 0.0033220476)
        self.assertEqual(round(st.w, 5), 214.57321)
        self.assertEqual(round(st.cpM.JmolK, 5), 158.29676)
        self.assertEqual(round(st.hM.Jmol, 3), 43040.574)
        self.assertEqual(round(st.sM.JmolK, 5), 119.95947)
        self.assertEqual(round(st.aM.Jmol, 4), -8265.2622)

        st = PropylenGlycol(T=400, rhom=13)
        self.assertEqual(round(st.P.MPa, 6), 61.287909)
        self.assertEqual(round(st.w, 4), 1467.8267)
        self.assertEqual(round(st.cpM.JmolK, 5), 227.48403)
        self.assertEqual(round(st.hM.Jmol, 3), -11727.204)
        self.assertEqual(round(st.sM.JmolK, 6), -38.665809)
        self.assertEqual(round(st.aM.Jmol, 5), -975.33462)

        st = PropylenGlycol(T=500, rhom=0.01)
        self.assertEqual(round(st.P.MPa, 9), 0.041324887)
        self.assertEqual(round(st.w, 5), 237.58157)
        self.assertEqual(round(st.cpM.JmolK, 5), 197.10614)
        self.assertEqual(round(st.hM.Jmol, 3), 60808.024)
        self.assertEqual(round(st.sM.JmolK, 5), 138.55317)
        self.assertEqual(round(st.aM.Jmol, 3), -12601.051)

        st = PropylenGlycol(T=500, rhom=13)
        self.assertEqual(round(st.P.MPa, 5), 196.80128)
        self.assertEqual(round(st.w, 4), 1633.6809)
        self.assertEqual(round(st.cpM.JmolK, 5), 252.37935)
        self.assertEqual(round(st.hM.Jmol, 3), 19730.659)
        self.assertEqual(round(st.sM.JmolK, 7), 8.1675599)
        self.assertEqual(round(st.aM.Jmol, 5), 508.31841)

        st = PropylenGlycol(T=680, rhom=13)
        self.assertEqual(round(st.P.MPa, 5), 430.54903)
        self.assertEqual(round(st.w, 4), 1814.6100)
        self.assertEqual(round(st.cpM.JmolK, 5), 280.86419)
        self.assertEqual(round(st.hM.Jmol, 3), 80765.269)
        self.assertEqual(round(st.sM.JmolK, 6), 81.467593)
        self.assertEqual(round(st.aM.Jmol, 4), -7751.8504)


if __name__ == "__main__":
    st = PropylenGlycol(T=400, rhom=0.001)
    print(st.w-214.57321)

