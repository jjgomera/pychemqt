#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


class nC7(MEoS):
    """Multiparameter equation of state for n-heptane"""
    name = "heptane"
    CASNumber = "142-82-5"
    formula = "CH3-(CH2)5-CH3"
    synonym = ""
    _refPropName = "HEPTANE"
    _coolPropName = "n-Heptane"
    rhoc = unidades.Density(231.9977)
    Tc = unidades.Temperature(540.13)
    Pc = unidades.Pressure(2736.0, "kPa")
    M = 100.202  # g/mol
    Tt = unidades.Temperature(182.55)
    Tb = unidades.Temperature(371.53)
    f_acent = 0.349
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 11
    _Tr = unidades.Temperature(525.389862)
    _rhor = unidades.Density(235.977855)
    _w = 0.350780196

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [15.063786601, -97.345252349],
           "ao_sinh": [13.7266, 43.55610], "sinh": [169.789/Tc, 1760.46/Tc],
           "ao_cosh": [30.4707], "cosh": [836.195/Tc]}

    CP1 = {"ao": 4,
           "ao_sinh": [13.7266, 43.55610], "sinh": [169.789, 1760.46],
           "ao_cosh": [30.4707], "cosh": [836.195]}

    CP3 = {"ao": 1.157528,
           "an": [0.070489617, -2.3419686e-5, -1.4768221e-9, -2.0117611e-12],
           "pow": [1, 2, 3, 4]}

    CP4 = {"ao": 30.4029/8.3159524*4.184,
           "ao_sinh": [2.5273083e8/8.3159524*4.184], "sinh": [1.6693200e3],
           "ao_cosh": [3.9046536e7/8.3159524*4.184], "cosh": [7.86001e2]}

    CP5 = {"ao": 4,
           "ao_exp": [15.29994054, 31.86604737, 14.10640675],
           "exp": [401.5547607, 1813.365387, 5041.869289]}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for heptane of Span "
                    "and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "OTO",
        "M": 100.204, "Tc": 540.13, "rhoc": 232/100.204,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 7.75,

        "nr1": [0.10543747645262e1, -0.26500681506144e1, 0.81730047827543,
                -0.30451391253428, 0.12253868710800, 0.27266472743928e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.49865825681670, -0.71432815084176e-3, -0.54236895525450,
                -0.13801821610756, -0.61595287380011e-2, 0.48602510393022e-3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-heptane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 7.75,

        "nr1": [0.10543747645262e1, -0.26500681506144e1, 0.81730047827543,
                -0.30451391253428, 0.122538687108, 0.27266472743928e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.49865825681670, -0.71432815084176e-3, -0.54236895525450,
                -0.13801821610756, -0.61595287380011e-2, 0.48602510393022e-3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heptane of Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP3,
        "ref": "NBP",

        "Tmin": 273.0, "Tmax": 500.0, "Pmax": 510000.0, "rhomax": 7.3348901,

        "nr1": [-0.52030538102, 0.338196304523, -0.491117643215e-2,
                0.200594802481, -0.260824422526e-1, -0.191516844204e1,
                0.364407895089, -0.142523250539, -.160069782510, .578283584822,
                0.476898816887, 0.937511885529e-1, -0.442185898133,
                0.553661375084e-1, -0.303420126133e-1, 0.138649129298e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.520305381020, -0.338196304523, 0.491117643215e-2,
                0.256518106995e1, -0.528051955217e1, 0.266827442122e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1]*6}

    starling = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heptane of Starling "
                    "(1973)",
        "__doi__": {"autor": "Starling, K.E.",
                    "title": "Fluid Thermodynamic Properties for Light "
                             "Petroleum Systems",
                    "ref": "Gulf Publishing Company, 1973.",
                    "doi": ""},

        "R": 8.3159524,
        "cp": CP4,
        "ref": "NBP",

        "Tmin": 255.37, "Tmax": 644.0, "Pmax": 55000.0, "rhomax": 7.2015722,

        "nr1": [0.153471579811e1, 0.521386289098, -0.107860953728e1,
                -0.902616154206, 0.117182735038, -0.986768914864e-4,
                0.287014205217, -0.359887681359, -0.860848441514e-2,
                0.952855119365e-2, 0.227922178775e-3],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.153471579811e1, -0.397448776976],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.51794447]*2}

    # This eq don't work
    ratanapisit = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for pentane of Ratanapisit (1999)",
        "__doi__": {"autor": "Ratanapisit, J., Ely, J.F.",
                    "title": "Application of New, Modified BWR Equations of "
                             "State to the Corresponding-States Prediction of "
                             "Natural Gas Properties",
                    "ref": "Int. J. Thermophys., 20(6) (1999) 1721-1735",
                    "doi": "10.1023/A:1022610013596"},

        "R": 8.31451,
        "Tc": 540.15, "Pc": 2738.8, "rhoc": 2.341,

        "cp": CP5,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 70000.0, "rhomax": 11.2,

        "b": [None, -9.53769631187e-3, 9.72551866385e-1, -2.60081304889e1,
              5.20865382062e3, -1.07729056282e6, -6.20474297014e-4,
              2.08733258744, - 1.37572781583e3, 6.95627225584e4,
              1.90615930406e-4, -5.61551412281e-1, 2.73983005070e2,
              6.28902715950e-2, -1.11012478028e1, 6.22600247144e2,
              1.57273923084, -6.63204129629e-2, -1.79732347053e1,
              1.24881866033, 3.81777590060e5, - 3.56280298214e7,
              1.7565835641e4, 4.54695406896e9, 2.05985406654e3,
              8.72406003683e5, 5.62265877351e1, -3.20150071052e7,
              3.57524917645, 3.27649699126e3, -1.15729200586e-1,
              3.93007045330e1, 3.88225605345e3]}

    eq = shortSpan, GERG, polt, starling   # ratanapisit
    _PR = [0.0074, -19.6376]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [-0.022144, 0.074781], "exp": [1.0324, 1.1636]}
    _dielectric = {
        "eq": 1,
        "a": [34.96, 0.035], "b": [162.24, 308.90], "c": [-37446, -39684],
        "Au": 29.84, "D": 2}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.76995e1, 0.14238e1, -0.20583e2, -0.50767e1, 0.19986e2],
        "t": [1.0, 1.5, 3.4, 4.7, 3.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.26395e1, 0.21806e2, -0.28896e2, 0.12609e2, 0.40749],
        "t": [0.322, 0.504, 0.651, 0.816, 6.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.24571, -6.3004, -19.144, -96.97, 216.43, -279.53],
        "t": [0.097, 0.646, 2.56, 6.6, 9.3, 10.7]}

    visco0 = {"__name__": "Michailidou (2014)",
              "__doi__": {
                  "autor": "Michailidou, E.K., Assael, M.J., Huber, M.L., "
                           "Abdulagatov, I.M., Perkins, R.A.",
                  "title": "Reference Correlation of the Viscosity of "
                           "n-Heptane from the Triple Point to 600 K and up "
                           "to 248 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023103",
                  "doi": "10.1063/1.4875930"},

              "eq": 1, "omega": 1,

              "ek": 426.118, "sigma": 0.61362,
              "n_chapman": 0.021357, "M": 100.202,
              "collision": [0.33974, -0.49396, 0, 0.0805],

              "Tref_virial": 426.118,
              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.0125,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],

              "Tref_res": 540.13, "rhoref_res": 232,
              "nr": [22.15, -15.0087, 3.71791],
              "tr": [-0.5, -0.5, -0.5],
              "dr": [5/3, 8/3, 11/3],

              "nr_num": [77.72818],
              "tr_num": [-0.5],
              "dr_num": [5/3],
              "nr_den": [9.73449, 9.519, -6.34076, 1, -2.51909],
              "tr_den": [0, -1, 0, 0, -1],
              "dr_den": [0, 0, 1, 2, 1]}

    visco1 = {"__name__": "Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 540.13,
              "no": [19.6036, -59.7839, 50.7528],
              "to": [0, 0.25, 0.5],

              "a": [3.76297e-5, 0.0, -4.40242e-5],
              "b": [1.38068e-4, 0.0, -9.11096e-5],
              "c": [9.93871e-5, -6.36533e-6, 0.0],
              "A": [-3.76786e-9, 1.92500e-9, 0.0],
              "B": [0.0, 9.75463e-9, 2.71874e-9],
              "C": [-1.24466e-6, 8.83261e-7, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"__name__": "Assael (2013)",
               "__doi__": {
                   "autor": "Assael, M.J., Bogdanou, I., Mylona, S.K., Huber, "
                            "M.L. Huber, Perkins, R.A., Vesovic, V.",
                   "title": "Reference Correlation of the Thermal "
                            "Conductivity of n-Heptane from the Triple Point "
                            "to 600 K and up to 250 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 42(2) (2013) 023101",
                   "doi": "10.1063/1.4794091"},

               "eq": 1,

               "Toref": Tc, "koref": 1e-3,
               "no_num": [-1.83367, 16.2572, -39.0996, 47.8594, 15.1925,
                          -3.39115],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [0.250611, -0.320871, 1.0],
               "to_den": [0, 1, 2],

               "Tref": Tc, "rhoref_res": 232, "kref_res": 1,
               "nr": [5.17785e-2, -9.24052e-2, 5.11484e-2, -7.76896e-3,
                      1.21637e-4, -7.72433e-3, 2.18899e-2, 1.7172500e-3,
                      -7.91642e-3, 1.83379e-3],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.245e-9, "gam0": 0.0586, "qd": 0.8e-9, "Tcref": 810.2}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_shortSpan(self):
        """Table III, Pag 46"""
        st = nC7(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.1651)
        self.assertEqual(round(st.P.MPa, 3), 7.957)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.7079)

        st2 = nC7(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 211.90)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.31964)

    def test_Michailidou(self):
        """Table 7, Pag 10"""
        self.assertEqual(round(nC7(T=250, rho=0).mu.muPas, 4), 4.9717)
        self.assertEqual(round(nC7(T=400, rho=0).mu.muPas, 4), 7.8361)
        self.assertEqual(round(nC7(T=550, rho=0).mu.muPas, 4), 10.7394)
        # self.assertEqual(round(nC7(T=250, rho=720).mu.muPas, 2), 725.61)
        self.assertEqual(round(nC7(T=400, rho=600).mu.muPas, 2), 175.94)
        self.assertEqual(round(nC7(T=550, rho=500).mu.muPas, 3), 95.102)

    def test_Assael(self):
        """Table 4, Pag 8"""
        self.assertEqual(round(nC7(T=250, rho=720).k.mWmK, 2), 137.08)
        self.assertEqual(round(nC7(T=400, rho=2).k.mWmK, 3), 21.794)
        self.assertEqual(round(nC7(T=400, rho=650).k.mWmK, 2), 120.74)
        # Using viscosity value given in paper work
        # self.assertEqual(round(nC7(T=535, rho=100).k.mWmK, 3), 51.655)
        self.assertEqual(round(nC7(T=535, rho=100).k.mWmK, 3), 51.503)
