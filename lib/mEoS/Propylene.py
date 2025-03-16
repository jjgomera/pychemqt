#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


class Propylene(MEoS):
    """Multiparameter equation of state for propylene"""
    name = "propylene"
    CASNumber = "115-07-1"
    formula = "CH2=CH-CH3"
    synonym = "R-1270"
    _refPropName = "PROPYLEN"
    _coolPropName = "Propylene"
    rhoc = unidades.Density(230.08)
    Tc = unidades.Temperature(364.211)
    Pc = unidades.Pressure(4555.0, "kPa")
    M = 42.07974  # g/mol
    Tt = unidades.Temperature(87.953)
    Tb = unidades.Temperature(225.531)
    f_acent = 0.146
    momentoDipolar = unidades.DipoleMoment(0.366, "Debye")
    id = 23

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-5.1823279651, 4.3639902765],
           "ao_exp": [1.544, 4.013, 8.923, 6.020],
           "titao": [324/Tc, 973/Tc, 1932/Tc, 4317/Tc]}

    Fi2 = {"ao_log": [1, 3.07317535],
           "pow": [0, 1],
           "ao_pow": [9.48120502357782, -4.47976952867319],
           "ao_exp": [1.7018443, 3.61342025, 8.83689058, 6.27183616],
           "titao": [1.01164134251849, 2.75278088800174, 5.16557061703243,
                     11.68984352477]}

    CP1 = {"ao": 0.65591381,
           "an": [0.44359641e-1, -.36650786e-4, 0.16822223e-7,
                  -.32651013e-11, 0.33747826e4],
           "pow": [1, 2, 3, 4, -2],
           "ao_exp": [-4.7032420], "exp": [615.8]}

    lemmon = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propylene of Lemmon et "
                    "al. (2013)",
        "__doi__": {
            "autor": "Lemmon, E.W., Overhoff, U., McLinden, M.O., Wagner, W.",
            "title": "A Reference Equation of State for the Thermodynamic "
                     "Properties of Propene for Temperatures from the Melting "
                     "Line to 575 K and Pressures up to 1000 MPa",
            "ref": "to be submitted to J. Phys. Chem. Ref. Data",
            "doi": ""},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 1000000.0, "rhomax": 23.1,

        "nr1": [0.4341002e-1, 0.1136592e1, -0.8528611, 0.5216669, -0.1382953e1,
                0.1214347],
        "d1": [4, 1, 1, 2, 2, 3],
        "t1": [1.0, 0.205, 0.56, 0.676, 1.0, 0.5],

        "nr2": [-0.5984662, -0.1391883e1, -0.1008434e1, 0.1961249, -0.3606930,
                -0.2407175e-2],
        "d2": [1, 1, 3, 2, 2, 8],
        "t2": [1.0, 1.94, 2.0, 1.0, 2.66, 0.83],
        "c2": [1, 2, 2, 1, 2, 1],
        "gamma2": [1]*6,

        "nr3": [0.7432121, 0.1475162, -0.2503391e-1, -0.2734409, 0.6378889e-2,
                0.1502940e-1, -0.3162971e-1, -0.4107194e-1, -0.1190241e1],
        "d3": [1, 1, 2, 3, 3, 2, 1, 2, 3],
        "t3": [1.6, 2.5, 3.0, 2.5, 2.72, 4.0, 4.0, 1.0, 4.0],
        "alfa3": [1.07, 0.66, 1.2, 1.12, 1.47, 1.93, 3.3, 15.4, 6],
        "beta3": [0.77, 0.83, 0.607, 0.4, 0.66, 0.07, 3.1, 387, 41],
        "gamma3": [1.21, 1.08, 0.83, 0.56, 1.22, 1.81, 1.54, 1.12, 1.4],
        "epsilon3": [0.78, 0.82, 1.94, 0.69, 1.96, 1.3, 0.38, 0.91, 0.7]}

    overhoff = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propylene of Overhoff "
                    "(2006).",
        "__doi__": {
            "autor": "Overhoff, U.",
            "title": "Development of a New Equation of State for the Fluid "
                     "Region of Propene for Temperatures from the Melting "
                     "Line to 575 K with Pressures to 1000 MPa as well as "
                     "Software for the computation of thermodynamic "
                     "properties of fluids",
            "ref": "PhD. Dissertation, Ruhr University, Bochum Germany, 2006",
            "doi": ""},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 1000000.0, "rhomax": 23.4,

        "nr1": [0.11167427541961e1, -0.76114879497376, -0.18654354344883e1,
                0.41500701892893e-1, 0.10706545719025e-1, 0.17481482892991e-1],
        "d1": [1, 1, 1, 3, 4, 4],
        "t1": [0.125, 0.625, 1.25, 0, 0.25, 1.25],

        "nr2": [0.56509607629258, 0.99156795771235, -0.16341922173416,
                -0.37037920319844e-1, -0.80058345775777e-1, 0.17004662808796,
                0.81351262137108e-1, -0.23817885171378, 0.12962562859214e-1,
                0.22577442976798e2, -0.43611886043491e2, 0.21944325628071e2,
                -0.66234078215924, -0.22258580712469e1, 0.29538388307646e1,
                -0.10257185828694e1, 0.20521625234481e-1, -0.36462809205891e-1,
                0.17625833164005e-1],
        "d2": [2, 3, 3, 3, 4, 4, 5, 5, 6, 1, 1, 1, 1, 2, 2, 2, 5, 6, 1],
        "t2": [2.25, 1.25, 2.125, 2.75, 0.125, 2, 1.125, 1.5, 1.375, 3.5, 3.75,
               4, 5, 3, 3.5, 4.5, 4.75, 3.25, 3, ],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3],
        "gamma2": [1]*19,

        "nr3": [0.31819374579431, -0.32648950998998, -0.37684374593786e2,
                0.72265437094447e2, -0.34814669335983e2, -0.39854778355193e1,
                0.37313453915501],
        "d3": [2, 2, 1, 1, 1, 2, 2],
        "t3": [3, 4, 2, 3, 4, 1, 1],
        "alfa3": [10, 10, 11, 11, 11, 25, 30],
        "beta3": [150, 150, 225, 225, 225, 300, 350],
        "gamma3": [1.13, 1.13, 1.19, 1.19, 1.19, 1.19, 1.19],
        "epsilon3": [0.85, 0.85, 1, 1, 1, 1, 1]}

    angus = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propene of Angus (1980)",
        "__doi__": {"autor": "Angus, S., Armstrong, B., de Reuck, K.M.",
                    "title": "International Thermodynamic Tables of the Fluid "
                             "State-7: Propylene (Propene)",
                    "ref": "IUPAC Chemical Data Series nº25, Pergamon Press, "
                           "1980",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP1,
        "ref": "IIR",
        "M": 42.0804, "Tc": 365.57, "Pc": 4664.6, "rhoc": 5.3086,

        "Tmin": 100.0, "Tmax": 600.0, "Pmax": 200000.0, "rhomax": 9.73,

        "nr1": [0.631922681460, 0.102655250604, -0.70798923e-2, 0.18624829,
                -0.1292611017e1, -0.5410160974e-1, 0.5069017035,
                -0.10606146125e1, 0.763136083, -0.850733053e-1, 0.438262575,
                0.2316495716e-1, 0.25503741325e-1, -0.57327581,
                -0.1141334722e-1, 0.2502895522, -0.468392547833e-1,
                0.325228355714e-2],
        "d1": [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 6, 7],
        "t1": [3, 4, 5, 1, 2, 3, 0, 1, 2, 2, 3, 0, 1, 3, -1, 3, 3, 3],

        "nr2": [-0.63192268146, -0.102655250604, 0.70798923e-2, -0.63192268146,
                -0.102655250604, -0.11049992895, -0.31596134073,
                -0.51327625302e-1, -0.4918627871e-1, -0.17109208434e-1,
                -0.1492467645e-1, -0.42773021085e-2, -0.8554604217e-3,
                -0.14257673695e-3],
        "d2": [0, 0, 0, 2, 2, 2, 4, 4, 6, 6, 8, 8, 10, 12],
        "t2": [3, 4, 5, 3, 4, 5, 3, 4, 3, 4, 3, 4, 4, 4],
        "c2": [2]*14,
        "gamma2": [1]*14}

    eq = lemmon, overhoff, angus

    # Unimplemented mBWR equation
    # Bender, E.
    # Equations of state for Ethylene and propylene
    # Cryogenics (1975) 667-673
    # 10.1016/0011-2275(75)90100-9

    _PR = [-0.1815, -16.3103]

    _surface = {"sigma": [0.05268], "exp": [1.186]}

    _melting = {
        "eq": 2,
        "__doi__": {
            "autor": "Reeves, L.E., Scott, G.J., Babb, S.E. Jr.",
            "title": "Melting Curves of Pressure-Transmitting fluids",
            "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
            "doi": "10.1063/1.1725068"},

        "Tmin": Tt, "Tmax": 2000.0,
        "Tref": Tt, "Pref": 0.00074864,
        "a2": [3196e5], "exp2": [2.821]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-6.75625, 2.02700, -1.35883, -2.74671, -0.936445],
        "t": [1.0, 1.5, 1.9, 4.3, 15.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.405430, 2.02481, 0.304022, 0.179159],
        "t": [0.195, 0.47, 2.25, 8.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-1.59841, -4.73840, -10.8886, -31.0312, -56.9431, -143.544],
        "t": [0.309, 0.853, 2.37, 5.2, 10., 20.]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",

              "ref": C3,
              "visco": "visco1",

              "ek": 298.9, "sigma": 0.4678, "omega": 5,

              "psi": [1.33962, -0.256307, 4.68211e-2], "psi_d": [0, 1, 2],
              "fint": [1.09939e-3, 3.72539e-7], "fint_t": [0, 1],
              "chi": [1.3521, -0.123177], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5e-10, "Tcref": 1.5*Tc}

    _viscosity = (trnECS, )

    thermo0 = {"__name__": "Assael (2016)",
               "__doi__": {
                   "autor": "Assael, M.J., Koutian, A., Huber, M.L., Perkins, "
                            "R.A.",
                   "title": "Reference Correlations of the Thermal "
                            "Conductivity of Ethene and Propene",
                   "ref": "J. Phys. Chem. Ref. Data 45(3) (2016) 033104",
                   "doi": "10.1063/1.4958984"},

               "eq": 1,

               "Toref": 364.211, "koref": 1e-3,
               "no_num": [-1.37218, 17.3386, -3.27682, 9.34452, 12.88,
                          -1.5705],
               "to_num": [0, 1, 2, 3, 4, 5],
               "no_den": [1.39367, -1.04648, 1],
               "to_den": [0, 1, 2],

               "Tref_res": 364.211, "rhoref_res": 229.63, "kref_res": 1e-3,
               "nr": [0.271511e1, -0.363839e2, 0.106159e3, -0.616755e2,
                      0.105424e2, 0.994697e1, 0.242705e2, -0.659429e2,
                      0.379916e2, -0.569120e1],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.02, "Xio": 0.198e-9,
               "gam0": 0.057, "qd": 0.43e-9, "Tcref": 546.32}

    _thermal = thermo0, trnECS


class Test(TestCase):
    """Testing"""

    def test_Assael(self):
        """Table 9, pag 12"""
        # Tiny error about inconsistency in ecs viscosity method between the
        # eq used for conformal state solver (Overhoff) and the eq used in
        # critical enhancement (Lemmon)
        kw = {"eq": "overhoff"}

        self.assertEqual(round(Propylene(T=200, rho=0, **kw).k.mWmK, 2), 8.75)
        self.assertEqual(round(Propylene(T=300, rho=0, **kw).k.mWmK, 2), 17.55)
        self.assertEqual(round(Propylene(T=400, rho=0, **kw).k.mWmK, 2), 29.18)
        self.assertEqual(round(Propylene(T=500, rho=0, **kw).k.mWmK, 2), 42.64)
        self.assertEqual(round(Propylene(T=200, P=1e5, **kw).k.mWmK, 1), 152.3)
        self.assertEqual(round(Propylene(T=300, P=1e5, **kw).k.mWmK, 2), 17.64)
        self.assertEqual(round(Propylene(T=400, P=1e5, **kw).k.mWmK, 2), 29.26)
        self.assertEqual(round(Propylene(T=500, P=1e5, **kw).k.mWmK, 2), 42.71)
        self.assertEqual(round(
            Propylene(T=200, P=2.5e7, **kw).k.mWmK, 1), 171.9)
        self.assertEqual(round(
            Propylene(T=300, P=2.5e7, **kw).k.mWmK, 1), 126.7)
        self.assertEqual(round(
            Propylene(T=400, P=2.5e7, **kw).k.mWmK, 2), 98.92)
        self.assertEqual(round(
            Propylene(T=500, P=2.5e7, **kw).k.mWmK, 2), 80.31)
        self.assertEqual(round(Propylene(T=200, P=5e7, **kw).k.mWmK, 1), 190.9)
        self.assertEqual(round(Propylene(T=300, P=5e7, **kw).k.mWmK, 1), 145.1)
        self.assertEqual(round(Propylene(T=400, P=5e7, **kw).k.mWmK, 1), 122.2)
        self.assertEqual(round(Propylene(T=500, P=5e7, **kw).k.mWmK, 1), 106.5)

        # Critical enhancement point, section 3.2.4, pag 12
        self.assertEqual(round(
            Propylene(T=350, rho=385, **kw).k.mWmK, 2), 81.48)
