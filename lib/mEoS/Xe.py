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

from lib.meos import MEoS
from lib import unidades


class Xe(MEoS):
    """Multiparameter equation of state for xenon"""
    name = "xenon"
    CASNumber = "7440-63-3"
    formula = "Xe"
    synonym = ""
    _refPropName = "XENON"
    _coolPropName = "Xenon"
    rhoc = unidades.Density(1102.8612)
    Tc = unidades.Temperature(289.733)
    Pc = unidades.Pressure(5842.0, "kPa")
    M = 131.293  # g/mol
    Tt = unidades.Temperature(161.405)
    Tb = unidades.Temperature(165.05)
    f_acent = 0.00363
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    # id = 994
    id = 1

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-3.8227178129, 3.8416395351],
           "ao_exp": [], "titao": []}

    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [], "ao_hyp": [], "hyp": []}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for xenon of Lemmon "
                    "and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 700000.0, "rhomax": 28.78,
        "Pmin": 81.77, "rhomin": 22.59,

        "nr1": [0.83115, -2.3553, 0.53904, 0.014382, 0.066309, 0.00019649],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.14996, -0.035319, -0.15929, -0.027521, -0.023305, 0.0086941],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    eq = lemmon,

    _surface = {"sigma": [-0.11538, 0.16598], "exp": [1.0512, 1.098]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [10.122], "expt1": [0], "expd1": [1],
                   "a2": [31.97, 46.97, -82.51, -948.4],
                   "expt2": [0, 1, 0], "expd2": [2, 2, 2.7]}
    _melting = {"eq": 1, "Tref": 1, "Pref": 101.325,
                "Tmin": Tt, "Tmax": 800.0,
                "a1": [-2573.936225, 0.7983277028], "exp1": [0, 1.589165],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 81.750,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-13.9, 14.0], "exp2": [1.06, 3.1],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.60231e1, 0.14989e1, -0.74906, -0.12194e1, -0.44905],
        "exp": [1., 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.13570e2, -0.47545e2, 0.63876e2, -0.39983e2, 0.12701e2],
        "exp": [0.56, 0.8, 1.0, 1.3, 1.6]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.0026, -6.056, -0.60339e2, 0.48838e3, -0.81974e3, 0.47287e3],
        "exp": [0.435, 1.4, 4.4, 6.2, 7.0, 8.6]}

    visco0 = {"eq": 2, "omega": 2,
              "collision": [128.829355170398, -824.923907889772, 2218.37801659791,
                            -3223.99202732053, 2718.40030222947, -1324.88234523685,
                            356.464839471621, -40.0927287567597, 0],
              "__name__": "NIST12",
              "__doc__": """Coefficients from NIST Thermophysical Properties of Pure Fluids Database, NIST12, Version 3.0, National Institute of Standards and Technology, Boulder, CO, 1992.""",
              "ek": 300, "sigma": 0.3297,
              "n_chapman": 0.305864975918623,
              "t_chapman": 0.0,
              "F": [0.768059558541217, -.585958377425158, 2.984837805288,
                    26.32847824613],
              "E": [-10.78336030151, 50.05660460723, 11.1406641168716,
                    -779.716643301403, 6.15104211699e-2, 10.7552268985402,
                    70.1937254720167],
              "rhoc": 5.3593311454524}

    _viscosity = visco0,

    thermo0 = {"eq": 3, "critical": 0,
               "__name__": "NIST12",
               "__doc__": """Coefficients from NIST Thermophysical Properties of Pure Fluids Database, NIST12""",

               "ek": 300, "sigma": 0.3297,
               "Nchapman": 0.305864975918623,
               "tchapman": 0,
               "b": [1.52442313680e-2, -9.05313615496e-2, 0.220032138191832,
                     -.278004805199205, 0.189554114709829, -6.36328719931e-2,
                     9.29951868906e-3, 0, 0],
               "F": [2.64173335524e-4, 4.73502202366e-5, -.6198732951154,
                     1525.9253243],
               "E": [-18.50657092152, 222.4871694717, 11.0124644286886,
                     -3621.41559218313, 5.14892242754e-3, 16.2049998648212,
                     -11.4853001847611],
               "rhoc": 5.15587382303351,
               "ff": 1.7124,
               "rm": 0.00000003669}

    _thermal = thermo0,


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = Xe(T=291, rhom=8)
        self.assertEqual(round(st.P.kPa, 3), 5986.014)
        self.assertEqual(round(st.hM.kJkmol, 3), 9193.668)
        self.assertEqual(round(st.sM.kJkmolK, 3), 36.895)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 28.692)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 3063.309)
        self.assertEqual(round(st.w, 3), 125.648)
