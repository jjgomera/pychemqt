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
from lib.mEoS import C3


class nC5(MEoS):
    """Multiparameter equation of state for n-pentane"""
    name = "pentane"
    CASNumber = "109-66-0"
    formula = "CH3-(CH2)3-CH3"
    synonym = "R-601"
    _refPropName = "PENTANE"
    _coolPropName = "n-Pentane"
    rhoc = unidades.Density(232.)
    Tc = unidades.Temperature(469.7)
    Pc = unidades.Pressure(3370.0, "kPa")
    M = 72.14878  # g/mol
    Tt = unidades.Temperature(143.47)
    Tb = unidades.Temperature(309.21)
    f_acent = 0.251
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 8
    _Tr = unidades.Temperature(449.271155)
    _rhor = unidades.Density(233.873368)
    _w = 0.247058753

    Fi1 = {"ao_log": [1, 3.0],
           "ao_sinh": [8.95043, 33.4032], "sinh": [0.380391739, 3.777411113],
           "ao_cosh": [21.836], "cosh": [1.789520971]}

    CP0 = {"ao": 4,
           "ao_sinh": [8.95043, 33.4032], "sinh": [178.67, 1774.25],
           "ao_cosh": [21.836], "cosh": [840.538]}

    CP1 = {"ao": 10.288132,
           "an": [-0.2695377e-1, 0.20951065e-3, -0.27910773e-6, 0.12266269e-9],
           "pow": [1, 2, 3, 4]}

    f = 4.184/8.3159524
    CP2 = {"ao": 22.5012*f,
           "ao_sinh": [2.057417e8*f], "sinh": [1.71958e3],
           "ao_cosh": [2.972927e7*f], "cosh": [8.02069e2]}

    CP3 = {"ao": 4,
           "ao_exp": [9.751560716, 22.71445741, 11.65392685],
           "exp": [404.8796661, 1785.491483, 4504.430788]}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for pentane of Span "
                    "and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP0,
        "ref": "OTO",
        "M": 72.15, "Tc": 469.7, "rhoc": 232/72.15,

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 11.2,

        "nr1": [0.10968643e1, -0.29988888e1, 0.99516887, -0.16170709,
                0.11334460, 0.26760595e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40979882, -0.40876423e-1, -0.38169482, -0.10931957,
                -0.32073223e-1, 0.16877016e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Kunz and "
                    "Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi": "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 143.47, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 10.57,

        "nr1": [0.10968643098001e1, -0.29988888298061e1, 0.99516886799212,
                -0.16170708558539, 0.11334460072775, 0.26760595150748e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40979881986931, -0.40876423083075e-1, -0.38169482469447,
                -0.10931956843993, -0.32073223327990e-1, 0.16877016216975e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 238.0, "Tmax": 573.0, "Pmax": 30000.0, "rhomax": 9.410819,

        "nr1": [-0.117648900900e1, 0.163499095773e1, -0.366669005817,
                0.724947274043, -0.221919300269e1, 0.188671490348e1,
                -0.195774652096e1, 0.308440851184, 0.437424419722,
                -0.625853472351, .382868807091, -0.119467393955, .218631441082,
                0.485668874195e-1, -0.132198161379, 0.213549844850e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.117648900900e1, -0.163499095773e1, 0.366669005817,
                -0.363660829618e-2, 0.633672105685, -0.705792643982],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.968832]*6}

    starling = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Starling "
                    "(1973)",
        "__doi__": {"autor": "Starling, K.E.",
                    "title": "Fluid Thermodynamic Properties for Light "
                             "Petroleum Systems",
                    "ref": "Gulf Publishing Company, 1973.",
                    "doi": ""},

        "R": 8.3159524,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 177.0, "Tmax": 589.0, "Pmax": 55000.0, "rhomax": 10.2534,

        "nr1": [0.175873733594e1, 0.485604047435, -0.111896446456e1,
                -0.685918143315, 0.368714111378e-1, -0.167498784887e-2,
                0.327765295239, -0.352742092747, -0.999487301826e-1,
                0.781999120830e-2, 0.221577806386e-2],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.175873733594e1, -0.411653507564],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.46812392]*2}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L., Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,

        "nr1": [2.20261753, 1.07797592, -3.82130221, 1.06627357e-1,
                3.07513215e-4, -2.84309667e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-7.28441220e-2, -4.60943732e-1, 8.39360011e-2, -1.50650444e-2,
                -0.203771872, -7.90244277e-3, -5.68993564e-2, -2.99387974e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

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
        "Tc": 469.65, "Pc": 3364.56, "rhoc": 3.2155,

        "cp": CP3,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 70000.0, "rhomax": 11.2,

        "b": [None, -7.41533782499e-2, 7.54044021950, -1.93328401588e2,
              3.39428034054e4, -5.12571561595e6, 1.51195406963e-3,
              -7.12225059892, 4.12664185793e3, 8.40258305443e5,
              -4.68416651753e-4, 3.03565637672, -1.42146321204e3,
              -1.10170659283e-1, -9.80664356304, 1.10979804446e3,
              2.98029604130, -1.41484307201e-1, -3.39208006239e1,
              2.08782048763, 5.38055429992e5, -6.40401885304e8,
              -1.19676622034e5, 1.71973349582e10, -3.06383363882e3,
              1.43168348944e6, 1.41452433419e1, -2.52955687564e7,
              -3.85316416299, 2.65416349789e3, 4.76643876980e-3,
              -8.37595968663, -1.35160880503e3]}

    eq = shortSpan, GERG, polt, starling, sun, ratanapisit
    _PR = [-0.0752, -17.5]

    _surface = {
        "__doi__": {"autor": "Mulero, A., Cachadiña, I., Bautista, D.",
                    "title": "Recommended Correlations for the Surface "
                             "Tension of n-Alkanes",
                    "ref": "J. Phys. Chem. Ref. Data 50(2) (2021) 023104",
                    "doi": "10.1063/5.0048675"},
        "sigma": [-0.20756, 0.40732, -0.14796],
        "exp": [1.0995, 1.1576, 1.2289]}
    _dielectric = {
        "eq": 1,
        "a": [25.39, 0.025], "b": [78.39, 54.15], "c": [-12480, -4800.0],
        "Au": 29.84, "D": 2}

    _melting = {
        "eq": 2,
        "__doi__": {
            "autor": "Reeves, L.E., Scott, G.J., Babb, S.E. Jr.",
            "title": "Melting Curves of Pressure-Transmitting fluids",
            "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
            "doi": "10.1063/1.1725068"},

        "Tmin": Tt, "Tmax": 2000.0,
        "Tref": Tt, "Pref": 0.076321,
        "a2": [6600e5], "exp2": [1.649]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.73918e1, 0.31102e1, -0.22415e1, -0.31585e1, -0.90451],
        "t": [1., 1.5, 1.74, 3.75, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.10178e1, 0.42703, 0.11334e1, 0.41518, -0.47950e-1],
        "t": [0.27, 0.44, 0.6, 4.0, 5.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.9389, -6.2784, -19.941, -16.709, -36.543, -127.99],
        "t": [0.4, 1.18, 3.2, 6.6, 7.0, 15.0]}

    trnECS = {"__name__": "Huber (2018)",

              "__doi__": {
                  "autor": "Huber, M.L.",
                  "title": "Models for Viscosity, Thermal Conductivity, and "
                           "Surface Tension of Selected Pure Fluids as "
                           "Implemented in REFPROP v10.0",
                  "ref": "NISTIR 8209",
                  "doi": "10.6028/NIST.IR.8209"},

              "eq": "ecs",
              "ref": C3,
              "visco": "visco1",

              "ek": 349.44, "sigma": 0.5790, "omega": 6,
              "n_chapman": 26.692e-3, "Fc": 1,

              "psi": [0.455019, 0.677221, -0.277823, 0.0372505],
              "psi_d": [0, 1, 2, 3]}

    visco0 = {"__name__": "Quiñones-Cisneros (2006)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Deiters, U.K.",
                  "title": "Generalization of the Friction Theory for "
                           "Viscosity Modeling",
                  "ref": "J. Phys. Chem. B, 110(25) (2006) 12820-12834",
                  "doi": "10.1021/jp0618577"},

              "eq": 4, "omega": 0,

              "Toref": 469.7,
              "no": [17.6805, -55.6942, 48.7177],
              "to": [0, 0.25, 0.5],

              "a": [1.08193e-5, -4.71699e-5, 0.0],
              "b": [1.21502e-4, -9.84766e-5, 0.0],
              "c": [5.08307e-5, -1.07e-5, 0.0],
              "A": [-2.10025e-10, -1.56583e-9, 0.0],
              "B": [1.98521e-8, 2.05972e-9, 0.0],
              "C": [-1.18487e-7, 1.69571e-7, 0.0]}

    _viscosity = (trnECS, visco0)

    thermo0 = {"__name__": "Vassiliou (2015)",
               "__doi__": {
                   "autor": "Vassiliou, C.-M., Assael, M.J., Huber, M.L., "
                            "Perkins, R.A.",
                   "title": "Reference Correlation of the Thermal Conductivity"
                            " of Cyclopentane, iso-pentane, and n-Pentane",
                   "ref": "J. Phys. Chem. Ref. Data 44(3) (2015) 033102",
                   "doi": "10.1063/1.4927095"},

               "eq": 1,

               "Toref": 469.7, "koref": 1e-3,
               "no_num": [-3.96685, 35.3805, 5.11554, -108.585, 179.573,
                          39.2128],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [2.71636, -5.76265, 6.77885, -0.59135, 1],
               "to_den": [0, 1, 2, 3, 4],

               "Tref_res": 469.7, "rhoref_res": 232, "kref_res": 1e-3,
               "nr": [7.76054e-1, 1.17655e2, -1.33101e2, 5.34026e1, -6.8793,
                      7.97696, -7.85888e1, 9.16089e1, -3.70431e1, 5.0962],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.227e-9, "gam0": 0.058, "qd": 0.668e-9, "Tcref": 704.55}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_shortSpan(self):
        """Table III, Pag 46"""
        st = nC5(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 3.2053)
        self.assertEqual(round(st.P.MPa, 3), 13.454)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.6052)

        st2 = nC5(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 213.42)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.34915)

    def test_Huber(self):
        """Table 7, pag 266"""
        self.assertEqual(round(
            # nC5(T=422.7, rhom=6.548).mu.muPas, 5), 80.03148)
            nC5(T=422.7, rhom=6.548).mu.muPas, 5), 80.16498)

    def test_Vassiliou(self):
        """Section 3.3.2, Pag 14"""
        # Viscosity value different to used in paper
        # self.assertEqual(round(nC5(T=460, P=3.3e6).k.mWmK, 3), 71.300)
        self.assertEqual(round(nC5(T=460, P=3.3e6).k.mWmK, 3), 71.394)
