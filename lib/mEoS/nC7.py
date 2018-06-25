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


class nC7(MEoS):
    """Multiparameter equation of state for n-heptane"""
    name = "heptane"
    CASNumber = "142-82-5"
    formula = "CH3-(CH2)5-CH3"
    synonym = ""
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
           "ao_exp": [], "titao": [],
           "ao_hyp": [13.7266, 30.4707, 43.55610, 0],
           "hyp": [0.314348398, 1.54813656, 3.259326458, 0]}

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [13.7266, 30.4707, 43.55610, 0],
           "hyp": [169.789, 836.195, 1760.46, 0]}

    CP3 = {"ao": 1.157528,
           "an": [0.070489617, -2.3419686e-5, -1.4768221e-9, -2.0117611e-12],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 30.4029/8.3159524*4.184,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [2.5273083e8/8.3159524*4.184, 3.9046536e7/8.3159524*4.184,
                      0, 0],
           "hyp": [1.6693200e3, 7.86001e2, 0, 0]}

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
        "Pmin": 0.17549e-3, "rhomin": 7.7457,

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
        "Pmin": 0.61166, "rhomin": 55.497,

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
        "Pmin": 0.17549e-3, "rhomin": 7.7457,

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
        "Pmin": 0.17549e-3, "rhomin": 7.7457,

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

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for pentane of Ratanapisit (1999)",
        "__doi__": {"autor": "Ratanapisit, J., Ely, J.F.",
                    "title": "Application of New, Modified BWR Equations of "
                             "State to the Corresponding-States Prediction of "
                             "Natural Gas Properties",
                    "ref": "Int. J. Thermophys., 20(6) (1999) 1721-1735",
                    "doi": "10.1023/A:1022610013596"},

        "R": 8.31434,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 70000.0, "rhomax": 11.2,
        "Pmin": 0.0000815, "rhomin": 10.558,

        "b": [None, -9.53769631187e-3, 9.72551866385e-1, -2.60081304889e1,
              5.20865382062e3, -1.07729056282e6, -6.20474297014e-4,
              2.08733258744, - 1.37572781583e3, 6.95627225584e4,
              1.90615930406e-4, -5.61551412281e-1, 2.73983005070e2,
              6.28902715950e-2, -1.11012478028e1, 6.22600247144e2,
              1.57273923084, -6.63204129629e-2, -17.9732347053, 1.24881866033,
              3.81777590060e5, - 3.56280298214e7, 1.7565835641e4,
              4.54695406896e9, 2.05985406654e3, 8.72406003683e5,
              56.2265877351, -3.20150071052e7, 3.57524917645, 3.27649699126e3,
              -1.15729200586e-1, 3.93007045330e1, 3.88225605345e3]}

    # eq = shortSpan, GERG, polt, starling, MBWR
    eq = shortSpan, GERG, polt, starling

    _surface = {"sigma": [0.07765, -0.02599], "exp": [1.319, 1.6]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.10924],  "expt0": [-1.], "expd0": [1.],
                   "a1": [34.96, 0.035], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [162.24, 308.9, -37446, -39684],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.76995e1, 0.14238e1, -0.20583e2, -0.50767e1, 0.19986e2],
        "exp": [1.0, 1.5, 3.4, 4.7, 3.6]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.26395e1, 0.21806e2, -0.28896e2, 0.12609e2, 0.40749],
        "exp": [0.322, 0.504, 0.651, 0.816, 6.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.24571, -6.3004, -19.144, -96.97, 216.43, -279.53],
        "exp": [0.097, 0.646, 2.56, 6.6, 9.3, 10.7]}

    visco0 = {"__name__": u"Quiñones-Cisneros (2006)",
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

    _viscosity = visco0,

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

    _thermal = thermo0,


class Test(TestCase):

    def test_shortSpan(self):
        # Table III, Pag 46
        st = nC7(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.1651)
        self.assertEqual(round(st.P.MPa, 3), 7.957)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.7079)

        st2 = nC7(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 211.90)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.31964)

    def test_Assael(self):
        # Table 4, Pag 8
        self.assertEqual(round(nC7(T=250, rho=720).k.mWmK, 2), 137.08)
        self.assertEqual(round(nC7(T=400, rho=2).k.mWmK, 3), 21.794)
        self.assertEqual(round(nC7(T=400, rho=650).k.mWmK, 2), 120.74)
        # Using viscosity value given in paper work
        # self.assertEqual(round(nC7(T=535, rho=100).k.mWmK, 3), 51.655)
        self.assertEqual(round(nC7(T=535, rho=100).k.mWmK, 3), 51.266)
