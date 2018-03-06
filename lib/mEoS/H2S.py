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


class H2S(MEoS):
    """Multiparameter equation of state for hydrogen sulfide"""
    name = "hydrogen sulfide"
    CASNumber = "7783-06-4"
    formula = "H2S"
    synonym = ""
    rhoc = unidades.Density(347.2841672)
    Tc = unidades.Temperature(373.1)
    Pc = unidades.Pressure(9000.0, "kPa")
    M = 34.08088  # g/mol
    Tt = unidades.Temperature(187.7)
    Tb = unidades.Temperature(212.85)
    f_acent = 0.1005
    momentoDipolar = unidades.DipoleMoment(0.97, "Debye")
    id = 50

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1.5],
           "ao_pow": [-4.0740770957, 3.7632137341, -0.002753352822675789],
           "ao_exp": [1.1364, 1.9721],
           "titao": [1823/Tc, 3965/Tc]}

    Fi2 = {"ao_log": [1, 3],
           "pow": [0, 1], "ao_pow": [9.336197742, -16.266508995],
           "ao_exp": [], "titao": [],
           "ao_hyp": [3.11942, 1.00243, 0, 0],
           "hyp": [4.914580541, 2.27065398, 0, 0]}

    Fi3 = {"ao_log": [1, 3.],
           "pow": [0, 1], "ao_pow": [7.881037, -3.20986],
           "ao_exp": [0.9767422, 2.151898], "titao": [4.506266, 10.15526],
           "ao_hyp": [], "hyp": []}

    CP3 = {"ao": 4.1012105,
           "an": [-0.16720073e-2, 0.75303152e-5, -0.62421053e-8, 0.18098453e-11],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 7.9468/8.3159524*4.184,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [-1.5769761e4/8.3159524*4.184, 2.0329947e6/8.3159524*4.184,
                      1.3861204e7/8.3159524*4.184, -3.5044957e6/8.3159524*4.184],
           "hyp": [4.33801e2, 8.43792e2, 1.48143e3, 1.10223e3]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for hydrogen sulfide of"
                    " Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 "
                             "Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 760.0, "Pmax": 170000.0, "rhomax": 29.12,
        "Pmin": 23.3, "rhomin": 29.12,

        "nr1": [0.87641, -2.0367, 0.21634, -0.050199, 0.066994, 0.00019076],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.20227, -0.0045348, -0.22230, -0.034714, -0.014885, 0.0074154],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    sakoda = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of Sakoda and Uematsu (2004)",
        "__doi__": {"autor": "Sakoda, N., Uematsu, M.",
                    "title": "A Thermodynamic Property Model for Fluid Phase Hydrogen Sulfide",
                    "ref": "Int. J. Thermophys., 25(3):709-737, 2004",
                    "doi": "10.1023/B:IJOT.0000034234.06341.8a"},
        "R": 8.314472,
        "cp": Fi3,
        "ref": {"Tref": 298.15, "Pref": 100., "ho": 0, "so": 0},

        "Tmin": 187.67, "Tmax": 760.0, "Pmax": 170000.0, "rhomax": 29.13,
        "Pmin": 23.3, "rhomin": 29.12,

        "nr1": [0.1545780, -0.1717693e1, -0.1595211e1, 0.2046589e1,
                -0.1690358e1, 0.9483623, -0.6800772e-1, 0.4372273e-2,
                0.3788552e-4, -0.3680980e-4, 0.8710726e-5],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4, 8, 9, 10],
        "t1": [0.241, 0.705, 1, 0.626, 1.15, 1.63, 0.21, 3.08, 0.827, 3.05, 3.05],

        "nr2": [0.6886876, 0.2751922e1, -0.1492558e1, 0.9202832, -0.2103469,
                0.1084359e-2, 0.3754723e-1, -0.5885793e-1, -0.2329265e-1,
                -0.1272600e-3, -0.1336824e-1, 0.1053057e-1],
        "d2": [1, 1, 1, 2, 5, 1, 4, 4, 3, 8, 2, 3],
        "t2": [0.11, 1.07, 1.95, 0.142, 2.13, 4.92, 1.75, 3.97, 11.8, 10,
               9.83, 14.2],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4],
        "gamma2": [1]*12}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of Polt et al. (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP3,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 760.0, "Pmax": 142000.0, "rhomax": 29.1,
        "Pmin": 23.85, "rhomin": 29.07,

        "nr1": [0.135782366339e1, -0.153224981014e1, 0.329107661253,
                0.195802782279e1, -0.301125182071e1, -0.126614059078e1,
                0.129960331548e1, -0.185645977138, -0.160919744092e1,
                0.234395817019e1, -0.378573094883, 0.758423219040,
                -0.973372615169, -0.120786235447, 0.209004959689,
                -0.919656385346e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.135782366339e1, 0.153224981014e1, -0.329107661253,
                0.891427552242, -0.204776100441e1, 0.101366381241e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.9873538]*6}

    starling = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of Starling (1973)",
        "__doi__": {"autor": "Starling, K.E.",
                    "title": "Fluid Thermodynamic Properties for Light Petroleum Systems",
                    "ref": "Gulf Publishing Company, 1973.",
                    "doi": ""},
        "R": 8.3159524,
        "cp": CP4,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 589.0, "Pmax": 55000.0, "rhomax": 29.578,
        "Pmin": 23.85, "rhomin": 29.07,

        "nr1": [0.110928333109e1, 0.188834546108, -0.930906931583,
                -0.411249591635, 0.140676923412e-1, -0.169077883177e-4,
                0.510265859853, -0.572402742986, -0.828859606622e-3,
                0.971664064871e-2, 0.140700425434e-4],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.110928333109e1, -0.26913741657],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.48524558]*2}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Kunz and "
                    "Wagner (2008).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for "
                             "Natural Gases and Other Mixtures: An Expansion "
                             "of GERG-2004",
                    "ref": "J. Chem.Eng. Data 57(11) (2012) 3032-3091",
                    "doi":  "10.1021/je300655b"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 760.0, "Pmax": 170000.0, "rhomax": 29.12,
        "Pmin": 23.3, "rhomin": 29.12,

        "nr1":  [0.87641, -2.0367, 0.21634, -0.050199, 0.066994, 0.19076e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.20227, -0.45348e-2, -0.2223, -0.034714, -.014885, .74154e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = lemmon, sakoda, polt, starling, GERG

    _surface = {"sigma": [0.078557], "exp": [1.2074]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.5884, 2.1582, -1.6054, -2.3870],
        "exp": [1, 1.5, 2, 4.8]}
    _liquid_Density = {
        "eq": 2,
        "ao": [11.833, -17.019, 7.8047],
        "exp": [1.63, 1.95, 2.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.9156, -7.7093, -19.543, -49.418],
        "exp": [0.49, 1.69, 4, 8]}

    visco0 = {"eq": 2, "omega": 1,
              "__name__": "NIST",
              "__doi__": {"autor": "",
                          "title": "Coefficients are taken from NIST14, Version 9.08",
                          "ref": "",
                          "doi": ""},

              "ek": 301.1, "sigma": 0.36237,
              "n_chapman": 0.1558117/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-12.328630418994, 782.29421491, 11.840322553,
                    -10401.582791, -0.0482407464, 69.709031672, 256.31792390],
              "rhoc": 10.2}

    visco1 = {"eq": 4, "omega": 1,
              "__name__": "Schmidt (2007)",
              "__doi__": {"autor": "Schmidt, K.A.G., Carroll, J.J., Quinones-Cisneros, S.E., and Kvamme, B.",
                          "title": "Hydrogen Sulfide Viscosity Modeling",
                          "ref": "Energy Fuels, 2008, 22 (5), pp 3424–3434",
                          "doi": "10.1021/ef700701h"},

              "Tref": 373.1, "muref": 1.0,
              "ek": 301.1, "sigma": 0.36237, "n_chapman": 0,
              "n_ideal": [4.36694e1, -12.1530e1, 9.35279e1],
              "t_ideal": [0, 0.25, 0.5],

              "a": [5.46919e-5, -7.32295e-6, -7.35622e-6],
              "b": [4.56159e-5, -1.82572e-5, -6.59654e-6],
              "c": [-4.33882e-6, 6.13716e-6, 0.0],
              "A": [6.67324e-9, -2.16365e-9, 0.0],
              "B": [-1.53973e-9, 2.17652e-9, 0.0],
              "C": [3.54228e-7, -4.76258e-8, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 1,
               "__name__": "NIST14",
               "__doi__": {"autor": "",
                           "title": "Coefficients are taken from NIST14 V9.08",
                           "ref": "",
                           "doi": ""},

               "Tref": 301.1, "kref": 1e-3,
               "no": [1.35558587, -0.1402163, 1],
               "co": [0, -1, -96],

               "Trefb": 373.4, "rhorefb": 10.2, "krefb": 1e-3,
               "nb": [21.7827447865, 10.8880930411, -7.45794247629,
                      3.65609005216, 1.89362258187, -1.10975687736],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
               "gam0": 0.0496, "qd": 0.3211e-9, "Tcref": 559.65}

    _thermal = thermo0,


class Test(TestCase):
    def test_shortLemmon(self):
        # Table 10, Pag 842
        st = H2S(T=375, rhom=10)
        self.assertEqual(round(st.P.kPa, 3), 9289.914)
        self.assertEqual(round(st.hM.kJkmol, 3), 15571.227)
        self.assertEqual(round(st.sM.kJkmolK, 3), 49.880)
        self.assertEqual(round(st.cvM.kJkmolK, 3), 42.574)
        self.assertEqual(round(st.cpM.kJkmolK, 3), 2645.392)
        self.assertEqual(round(st.w, 3), 248.512)
