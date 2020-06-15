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


class iC5(MEoS):
    """Multiparameter equation of state for isopentane"""
    name = "isopentane"
    CASNumber = "78-78-4"
    formula = "(CH3)2-CH-CH2-CH3"
    synonym = "R-601a"
    _refPropName = "IPENTANE"
    _coolPropName = "Isopentane"
    rhoc = unidades.Density(235.99865938)
    Tc = unidades.Temperature(460.35)
    Pc = unidades.Pressure(3378.0, "kPa")
    M = 72.14878  # g/mol
    Tt = unidades.Temperature(112.65)
    Tb = unidades.Temperature(300.98)
    f_acent = 0.2274
    momentoDipolar = unidades.DipoleMoment(0.11, "Debye")
    id = 7

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [2.5822330405, 1.1609103419],
           "ao_exp": [7.4056, 9.5772, 15.765, 12.119],
           "titao": [442/Tc, 1109/Tc, 2069/Tc, 4193/Tc]}

    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [15.449907693, -101.298172792],
           "ao_sinh": [11.7618, 33.1688],
           "sinh": [0.635392636, 4.169371131],
           "ao_cosh": [20.1101],
           "cosh": [1.977271641]}

    f = 72.151/8.3143
    CP3 = {"ao": 0.396504*f,
           "an": [0.260678e-2*f, 0.93677e-5*f, -0.158286e-7*f,  0.76525e-11*f],
           "pow": [1, 2, 3, 4]}

    f = 4.184/8.3159524
    CP4 = {"ao": 21.3861*f,
           "ao_sinh": [2.1524504e8*f], "sinh": [1.70158e3],
           "ao_cosh": [2.8330244e7*f], "cosh": [7.75899e2]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isopentane of "
                    "Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 1000000.0, "rhomax": 13.3,

        "nr1": [1.0963, -3.0402, 1.0317, -0.15410, 0.11535, 0.00029809],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.39571, -0.045881, -0.35804, -0.10107, -0.035484, 0.018156],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isopentane of Kunz and "
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

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 1000000.0, "rhomax": 13.3,

        "nr1": [0.11017531966644e1, -0.30082368531980e1, 0.99411904271336,
                -0.14008636562629, 0.11193995351286, 0.29548042541230e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.36370108598133, -0.48236083488293e-1, -0.35100280270615,
                -0.10185043812047, -0.35242601785454e-1, 0.19756797599888e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isopentane of Polt "
                    "(1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP3,
        "ref": "NBP",

        "Tmin": 200.0, "Tmax": 553.0, "Pmax": 7500.0, "rhomax": 5.2252,

        "nr1": [-0.143819012123e1, 0.138298276836e1, -0.203328695121,
                0.619304204378, -0.311353942178e1, 0.316914412369e1,
                -0.218812895934e1, 0.211230723299, 0.765790344231,
                -0.851773312153, 0.706192861166, -0.165802139239,
                0.781356542750e-1, 0.106516957202, -0.205642736936,
                0.360787537633e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.143819012123e1, -0.138298276836e1, 0.203328695121,
                -0.213463476736e1, 0.547491842897e1, -0.335666356499e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.002528]*6}

    starling = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isopentane of Starling "
                    "(1973)",
        "__doi__": {"autor": "Starling, K.E.",
                    "title": "Fluid Thermodynamic Properties for Light "
                             "Petroleum Systems",
                    "ref": "Gulf Publishing Company, 1973.",
                    "doi": ""},

        "R": 8.3159524,
        "cp": CP4,
        "ref": "NBP",

        "Tmin": 199.82, "Tmax": 589.0, "Pmax": 55000.0, "rhomax": 9.9258626,

        "nr1": [0.179378842786e1, 0.258488286720, -0.812072482201,
                -0.753941018871, 0.565338153509e-1, -0.115706201242e-2,
                0.406090628523, -0.469700474204, -0.967480812300e-1,
                0.958936263943e-2, 0.197520012548e-2],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.179378842786e1, -0.431019031876],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2, 2],
        "gamma2": [0.48056842]*2}

    eq = lemmon, GERG, polt, starling
    _PR = [-0.0840, -17.8184]

    _surface = {"sigma": [0.051], "exp": [1.209]}
    _dielectric = {
        "eq": 1,
        "a": [25.31, 0.025], "b": [108.9, 63.68], "c": [-15447, -5449.3],
        "Au": 73.69, "D": 2}

    _melting = {
            "eq": 2,
            "__doi__": {
                "autor": "Reeves, L.E., Scott, G.J., Babb, S.E. Jr.",
                "title": "Melting Curves of Pressure-Transmitting fluids",
                "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                "doi": "10.1063/1.1725068"},

            "Tmin": Tt, "Tmax": 2000.0,
            "Tref": Tt, "Pref": 8.952745e-5,
            "a2": [5916e5], "exp2": [1.563]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.72392e1, 0.22635e1, -0.18237e1, -0.29997e1, -0.27752e1],
        "t": [1., 1.5, 2.02, 4.24, 16.1]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.18367e2, -0.30283e2, 0.13557e2, -0.90533, 0.20927e1],
        "t": [1.21, 1.41, 1.65, 0.09, 0.164]}
    _vapor_Density = {
        "eq": 2,
        "n": [-38.825, 79.040, -48.791, -21.603, -57.218, -151.64],
        "t": [0.565, 0.66, 0.77, 3.25, 7.3, 16.6]}

    thermo0 = {"__name__": "Vassiliou (2015)",
               "__doi__": {
                   "autor": "Vassiliou, C.-M., Assael, M.J., Huber, M.L., "
                            "Perkins, R.A.",
                   "title": "Reference Correlation of the Thermal Conductivity"
                            " of Cyclopentane, iso-pentane, and n-Pentane",
                   "ref": "J. Phys. Chem. Ref. Data 44(3) (2015) 033102",
                   "doi": "10.1063/1.4927095"},

               "eq": 1,

               "Toref": 460.35, "koref": 1e-3,
               "no_num": [0.773049, -15.9754, 218.987, -329.556, 281.075,
                          53.326],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [5.10467, -8.12044, 8.11607, -0.294969, 1],
               "to_den": [0, 1, 2, 3, 4],

               "Tref_res": 460.35, "rhoref_res": 236, "kref_res": 1e-3,
               "nr": [-1.17507e1, -1.61346e1, 5.27254e1, -2.7494e1, 4.54817,
                      5.14003, 5.58445e1, -9.51474e1, 4.75268e1, -7.29296],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02,
               "Xio": 0.227e-9, "gam0": 0.058, "qd": 0.664e-9, "Tcref": 690.53}

    _thermal = thermo0,


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = iC5(T=462, rhom=3)
        self.assertEqual(round(st.P.kPa, 3), 3458.617)
        self.assertEqual(round(st.hM.kJkmol, 3), 37318.534)
        self.assertEqual(round(st.sM.kJkmolK, 3), 94.634)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 190.631)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 4660.943)
        self.assertEqual(round(st.w, 3), 96.324)

    def test_Vassiliou(self):
        # Section 3.2.2, Pag 11
        # Viscosity value different to used in paper
        # self.assertEqual(round(iC5(T=460, P=3.5e6).k.mWmK, 3), 59.649)
        self.assertEqual(round(iC5(T=460, P=3.5e6).k.mWmK, 3), 59.388)
