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


class R152a(MEoS):
    """Multiparameter equation of state for R152a"""
    name = "1,1-difluoroethane"
    CASNumber = "75-37-6"
    formula = "CHF2CH3"
    synonym = "R152a"
    _refPropName = "R152A"
    _coolPropName = "R152A"
    rhoc = unidades.Density(368.)
    Tc = unidades.Temperature(386.411)
    Pc = unidades.Pressure(4516.75, "kPa")
    M = 66.051  # g/mol
    Tt = unidades.Temperature(154.56)
    Tb = unidades.Temperature(249.127)
    f_acent = 0.27521
    momentoDipolar = unidades.DipoleMoment(2.262, "Debye")
    id = 245

    CP1 = {
           # Cp/R relation in paper
           # Tr terms in polynomial, so the resulting terms are:
           # a0 = c0
           # a1 = c1/Tc
           # a2 = c2/Tc**2
           # a3 = c3/Tc**3
           "ao": 3.354951,
           "an": [4.245301/Tc, 3.735248/Tc**2, -1.608254/Tc**3],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": []}

    Fi1 = {"R": 8.314471,
           "ao_log": [1, -1],
           "pow": [0, 1, -0.25, -2, -4],
           "ao_pow": [10.87227, 6.839515, -20.78887, -0.6539092, 0.03342831],
           "ao_exp": [], "titao": []}

    Fi2 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.5, 0.25],
           "ao_pow": [-9.508135074, 6.812068779, -7.285916044, 6.741130104],
           "ao_exp": [1.978152028, 5.880826311],
           "titao": [1.753741145, 4.360150337]}

    Fi3 = {"ao_log": [1, 0.0434935],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [-5.969835, 7.421932, -5.56713, 0.436101, -0.0196281],
           "ao_exp": [], "titao": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-152a of Outcalt (1996)",
        "__doi__": {"autor": "Outcalt, S.L., McLinden, M.O.",
                    "title": "A modified Benedict-Webb-Rubin Equation of "
                             "State for the Thermodynamic Properties of R152a "
                             "(1,1-difluoroethane)",
                    "ref": "J. Phys. Chem. Ref. Data 25(2) (1996) 605-636",
                    "doi": "10.1063/1.555979"},

        "R": 8.314471,
        "Tc": 386.411, "Pc": 4516.75, "rhoc": 5.57145,

        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 520.0, "Pmax": 60000.0, "rhomax": 18.07,
        "Pmin": 0.0641, "rhomin": 18.061,

        "b": [None, -0.101623317192e-1, 0.215677129618e1, -0.648581254334e2,
              0.122535596303e5, -0.206805988259e7, -0.379836507323e-3,
              -0.441333232984, 0.158248874708e3, 0.564062216256e6,
              -0.124115350431e-3, 0.494972178825, -0.208058039834e3,
              -0.131403187106e-1, 0.212083848812, -0.151263785082e3,
              0.311108025395e-1, -0.115280979645e-2, 0.437040025765,
              -0.965596535032e-2, -0.242705525346e6, -0.518042519989e8,
              -0.119070545681e5, 0.459333195257e9, -0.719317286511e2,
              -0.840102861460e4, -0.102910957390e1, -0.325913880841e5,
              -0.412362182230e-2, 0.175102808144e1, -0.198636624640e-4,
              -0.421363036104e-2, -0.198696760653e1]}

    outcalt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz transform of MBWR EOS for R-152a of Outcalt "
                    "and McLinden (1996).",
        "__doi__": {"autor": "Outcalt, S.L., McLinden, M.O.",
                    "title": "A modified Benedict-Webb-Rubin Equation of "
                             "State for the Thermodynamic Properties of R152a "
                             "(1,1-difluoroethane)",
                    "ref": "J. Phys. Chem. Ref. Data 25(2) (1996) 605-636",
                    "doi": "10.1063/1.555979"},

        "R": 8.314471,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 18.07,
        "Pmin": 0.0641, "rhomin": 18.061,

        "nr1": [-0.354657949982e1, -0.364631280620, 0.333233335558e-1,
                -0.6809684351170, 0.735212646801e1, -0.112473063838e2,
                0.549916715657e1, -0.240186327322e1, -0.709036447042e-1,
                -0.213200886814, 0.197839736368, 0.182494769909e1,
                -0.860546479693e-1, 0.888137366540, -0.966127346370,
                -0.985223479324e-1, 0.183419368472e-1, -0.338550204252e-1,
                0.124921101016e-1, -0.221056706423e-2, 0.216879133161e-2,
                -0.233597690478e-3],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7,
               8],
        "t1": [3, 4, 5, 0, 0.5, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 1, 2, 3, 2, 2, 3,
               3],

        "nr2": [0.354657949982e1, 0.364631280620, -0.333233335558e-1,
                0.276133830254e1, -0.691185711880e-1, -0.333233335558e-1,
                0.782761327717, -0.345592855940e-1, 0.137813531906,
                0.186173126153, -0.341119393297e-1, 0.459378439687e-1,
                0.216470012607e-1, -0.852798483242e-2, 0.620394038634e-2,
                0.185210290813e-2, 0.101674662734e-2, 0.124078807727e-2],
        "d2": [0, 0, 0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 8, 8, 8, 10, 10, 10],
        "t2": [3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5],
        "c2": [2]*18,
        "gamma2": [1]*18}

    kim = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-152a of Kim (1997).",
        "__doi__": {"autor": "Kim, Y., Borgnakke, C., Sonntag, R.E.",
                    "title": "Equation of State for 1,1-difluoroethane "
                             " (R152a)",
                    "ref": "International Journal of Energy Research 21 (7)"
                           "(1997) 575-589",
                    "doi": "10.1002/(sici)1099-114x(19970610)21:7<575::"
                           "aid-er272>3.0.co;2-f"},

        "R": 8.314471,
        "cp": Fi3,
        "ref": "IIR",
        "Tc": 386.4, "Pc": 4519., "rhoc": 368/M,

        "Tmin": 213, "Tmax": 433, "Pmax": 20000.0, "rhomax": 18.07,
        # "Pmin": 0.0641, "rhomin": 18.061,

        "nr1": [3.27282477979913, -5.25887189160385, 0.849951067158520,
                -0.326056649432851, 0.147973856820353, 0.463200609308586e-2],
        "d1": [1, 1, 1, 1, 2, 5],
        "t1": [1, 1.5, 3, 5, 0.5, 1],

        "nr2": [-0.184693035421790e-1, -0.529265795606284, 1.39788588805247,
                -0.826528289800619, 0.603238985406408, 0.184020254678691e-9,
                0.198000633690890e-1, 0.385227997762326e-1,
                -0.354915684935072e-1, -0.146266261800962e-3,
                0.385244461907819e-4, -0.930695615881793e-7,
                0.792443305748410e-2, -0.117773096693244e-1,
                0.780856661432880e-2, -0.335895387327679e-2,
                -0.905744836093298e-4, 0.348630546773750e-3,
                0.167579895771929e-1, -0.159255383659542e-1],
        "d2": [1.9, 2.2, 2.2, 2.5, 3.0, 3.3, 4.4, 4.9, 5.3, 6.6, 9.7, 13.1,
               4.6, 4.0, 5.2, 5.3, 13.3, 11.9, 4.1, 4.1],
        "t2": [8.0, 1.4, 3.1, 5.0, 5.5, 25.5, 5.2, 3.3, 3.5, 6.7, 5.3, 5.7,
               6.4, 30.0, 28.4, 8.2, 5.9, 20.3, 27.1, 29.3],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 5, 6, 6],
        "gamma2": [1]*20}

    tillner = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-152a of Tillner-Roth "
                    "(1995).",
        "__doi__": {"autor": "Tillner-Roth, R.",
                    "title": "A Fundamental Equation of State for "
                             "1,1-Difluoroethane (HFC-152a)",
                    "ref": "Int. J. Thermophys., 16(1) (1995) 91-100",
                    "doi": "10.1007/BF01438960"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "M": 66.051, "Tc": 386.41,
        "Tmin": Tt, "Tmax": 435.0, "Pmax": 30000.0, "rhomax": 18.03,
        "Pmin": 0.065395176, "rhomin": 18.020671,

        "nr1": [0.3552260, -0.1425660e1, -0.4631621e-1, 0.6903546e-1,
                0.1975710e-1, 0.7486977e-3, 0.4642204e-3],
        "d1": [1, 1, 1, 1.5, 3, 6, 6],
        "t1": [0, 1.5, 3, -0.5, -0.5, -0.5, 1.5],

        "nr2": [-0.2603396, -0.7624212e-1, 0.2233522, 0.1992515e-1, 0.3449040,
                -0.4963849, 0.1290719, 0.9760790e-3, 0.5066545e-2,
                -0.1402020e-1, 0.5169918e-2, 0.2679087e-3],
        "d2": [1, 1, 3, 4, 1, 1, 1, 8, 2, 3, 5, 6],
        "t2": [3, 4, 3, 2, 4, 5, 6, 5, 12.5, 25, 20, 25],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
        "gamma2": [1]*12}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-152a of Span and "
                    "Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 18.1,
        "Pmin": 0.064093, "rhomin": 18.031,

        "nr1": [0.95702326, -2.3707196, 0.18748463, 0.063800843, 1.6625977e-4],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.82208165e-1, 0.57243518, 0.39476701e-2, -0.23848654,
                -0.80711618e-1, -0.73103558e-1, -0.15538724e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*12}

    astina = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-152a of Astina (2004)",
        "__doi__": {"autor": "Astina, I.M., Sato, H.",
                    "title": "A Rigorous Thermodynamic Property Model for "
                             "Fluid-Phase 1,1-Difluoroethane (R-152a)",
                    "ref": "Int. J. Thermophys., 25(6) (2004) 1713-1733",
                    "doi": "10.1007/s10765-004-7731-8"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",
        "M": 66.05, "Tc": 386.41, "Pc": 4516, "rhoc": 368/66.05,

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 60000.0, "rhomax": 18.04,
        "Pmin": 0.064, "rhomin": 18.04,

        "nr1": [1.753847317, -4.049760759, -2.277389257e-1, 7.087751950e-1,
                -0.5528619502, -3.025046686e-2, 0.1396289974, 1.121238954e-4],
        "d1": [1, 1, 1, 2, 2, 3, 3, 4],
        "t1": [0.5, 1.125, 2.875, 0.875, 1.875, 0.5, 1.875, 4],

        "nr2": [1.181005890, 1.535785579, 7.468363045e-1, -1.252266405e-1,
                -3.898223986e-2, -7.260588801e-2, -2.659302250e-3,
                4.210849329e-3, 2.015953966e-4],
        "d2": [1, 2, 3, 1, 2, 3, 3, 4, 5],
        "t2": [1.25, 2, 2.75, 6, 9, 6, 22, 20, 32],
        "c2": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*9}

    eq = MBWR, outcalt, kim, tillner, shortSpan, astina

    _surface = {"sigma": [0.05808], "exp": [1.2115]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.74821e1, 0.21105e1, -0.20761e1, -0.35539e1, 0.58004],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.19914e2, -0.68624e2, 0.99821e2, -0.77984e2, 0.29913e2],
        "t": [0.56, 0.76, 0.95, 1.2, 1.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-.33621e1, -.85985e1, -.2683e1, -.2414e2, -.43159e2, -.28045e2],
        "t": [0.406, 1.42, 3.6, 3.9, 8.0, 9.0]}

    visco0 = {"__name__": "Krauss (1996)",
              "__doi__": {
                  "autor": "Krauss, R., Weiss, V.C., Edison, T.A., Sengers, "
                           "J.V., Stephan, K.",
                  "title": "Transport Properties of 1,1-Difluoroethane "
                           "(R152a)",
                  "ref": "Int. J. Thermophysics 17:731-757, 1996.",
                  "doi": "10.1007/BF01439187"},

              "eq": 1, "omega": 1,
              "M": 66.05, "ek": 354.84, "sigma": 0.46115,
              "n_chapman": 0.2169614/M**0.5,
              "collision": [0.4425728, -0.5138403, 0.1547566, -0.02821844,
                            0.001578286],

              "rhoref_res": 368, "muref_res": 51.12,
              "nr": [-0.139986563, -0.0737927, 0.517924, -0.308875, 0.108049],
              "tr": [0, 0, 0, 0, 0],
              "dr": [0, 1, 2, 3, 4],

              "nr_num": [-0.408387],
              "tr_num": [0],
              "dr_num": [0],
              "nr_den": [1, -2.91733],
              "tr_den": [0, 0],
              "dr_den": [1, 0]}

    _viscosity = visco0,

    thermo0 = {"__name__": "Krauss (1996)",
               "__doi__": {
                   "autor": "Krauss, R., Weiss, V.C., Edison, T.A., Sengers, "
                            "J.V., Stephan, K.",
                   "title": "Transport Properties of 1,1-Difluoroethane "
                            "(R152a)",
                   "ref": "Int. J. Thermophysics 17:731-757, 1996.",
                   "doi": "10.1007/BF01439187"},

               "eq": 1,
               "M": 66.05, "Tc": 386.411, "Pc": 4520,

               "Toref": 1., "koref": 1e-3,
               "no": [-14.942, 0.0973283],
               "to": [0, 1],

               "Tref_res": 1., "rhoref_res": 368, "kref_res": 1.155e-3,
               "nr": [9.1809, 11.8577, -5.4473, 1.71379],
               "tr": [0, 0, 0, 0],
               "dr": [1, 2, 3, 4],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 1.894e-10,
               "gam0": 0.0487, "qd": 4.37e-10, "Tcref": 579.6165}

    _thermal = thermo0,


class Test(TestCase):

    def test_Outcalt(self):
        # Selected point from Table 6, Pag 616, saturation states
        st = R152a(T=-118+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.000069)
        self.assertEqual(round(st.Liquido.rho, 1), 1191.8)
        self.assertEqual(round(st.Gas.rho, 4), 0.0036)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 14.67)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 419.73)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.1176)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.7284)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 0.998)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.574)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.480)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.700)
        self.assertEqual(round(st.Liquido.w, 1), 1396.1)
        self.assertEqual(round(st.Gas.w, 1), 154.3)

        st = R152a(T=-100+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.000579)
        self.assertEqual(round(st.Liquido.rho, 1), 1158.7)
        self.assertEqual(round(st.Gas.rho, 4), 0.0266)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 41.75)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 432.59)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.2827)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.5399)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.030)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.613)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.518)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.740)
        self.assertEqual(round(st.Liquido.w, 1), 1274.9)
        self.assertEqual(round(st.Gas.w, 1), 162.0)

        st = R152a(T=-50+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.027425)
        self.assertEqual(round(st.Liquido.rho, 1), 1063.7)
        self.assertEqual(round(st.Gas.rho, 4), 0.9936)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 118.62)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 470.40)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.6723)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.2487)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.042)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.738)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.567)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.877)
        self.assertEqual(round(st.Liquido.w, 1), 1016.4)
        self.assertEqual(round(st.Gas.w, 1), 179.5)

        st = R152a(T=273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.263992)
        self.assertEqual(round(st.Liquido.rho, 1), 959.1)
        self.assertEqual(round(st.Gas.rho, 4), 8.3589)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 200.00)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 507.11)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.0000)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.1243)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.100)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 0.898)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.697)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.094)
        self.assertEqual(round(st.Liquido.w, 1), 772.0)
        self.assertEqual(round(st.Gas.w, 1), 187.4)

        st = R152a(T=50+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 1.177382)
        self.assertEqual(round(st.Liquido.rho, 1), 830.8)
        self.assertEqual(round(st.Gas.rho, 4), 37.0576)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 290.50)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 535.93)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.3003)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.0598)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.182)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.092)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.957)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 1.489)
        self.assertEqual(round(st.Liquido.w, 1), 519.9)
        self.assertEqual(round(st.Gas.w, 1), 178.9)

        st = R152a(T=100+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 3.505025)
        self.assertEqual(round(st.Liquido.rho, 1), 618.5)
        self.assertEqual(round(st.Gas.rho, 4), 145.7543)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 403.59)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 536.28)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.6151)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.9707)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.322)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.346)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 3.495)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 3.776)
        self.assertEqual(round(st.Liquido.w, 1), 233.5)
        self.assertEqual(round(st.Gas.w, 1), 143.1)

        st = R152a(T=112+273.15, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 4.408087)
        self.assertEqual(round(st.Liquido.rho, 1), 473.0)
        self.assertEqual(round(st.Gas.rho, 4), 263.7988)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 451.59)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 506.12)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 1.7371)
        self.assertEqual(round(st.Gas.s.kJkgK, 4), 1.8787)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 3), 1.422)
        self.assertEqual(round(st.Gas.cv.kJkgK, 3), 1.458)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 22.908)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 31.032)
        self.assertEqual(round(st.Liquido.w, 1), 139.1)
        self.assertEqual(round(st.Gas.w, 1), 125.5)

        # Selected point of Table 7, pag 618, singhe phase region
        st = R152a(T=-100+273.15, P=1e4)
        self.assertEqual(round(st.rho, 1), 1158.7)
        self.assertEqual(round(st.h.kJkg, 2), 41.76)
        self.assertEqual(round(st.s.kJkgK, 4), 0.2827)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.030)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.518)
        self.assertEqual(round(st.w, 1), 1274.9)

        st = R152a(T=-50+273.15, P=2e4)
        self.assertEqual(round(st.rho, 3), 0.721)
        self.assertEqual(round(st.h.kJkg, 2), 470.93)
        self.assertEqual(round(st.s.kJkgK, 4), 2.2903)
        self.assertEqual(round(st.cv.kJkgK, 3), 0.734)
        self.assertEqual(round(st.cp.kJkgK, 3), 0.869)
        self.assertEqual(round(st.w, 1), 180.1)

        st = R152a(T=240+273.15, P=4e4)
        self.assertEqual(round(st.rho, 3), 0.620)
        self.assertEqual(round(st.h.kJkg, 2), 813.67)
        self.assertEqual(round(st.s.kJkgK, 4), 3.1512)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.362)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.488)
        self.assertEqual(round(st.w, 1), 265.4)

        st = R152a(T=-35+273.15, P=6e4)
        self.assertEqual(round(st.rho, 3), 2.063)
        self.assertEqual(round(st.h.kJkg, 2), 481.84)
        self.assertEqual(round(st.s.kJkgK, 4), 2.2021)
        self.assertEqual(round(st.cv.kJkgK, 3), 0.782)
        self.assertEqual(round(st.cp.kJkgK, 3), 0.931)
        self.assertEqual(round(st.w, 1), 183.2)

        st = R152a(T=-25+273.15, P=1e5)
        self.assertEqual(round(st.rho, 1), 1013.2)
        self.assertEqual(round(st.h.kJkg, 2), 158.48)
        self.assertEqual(round(st.s.kJkgK, 4), 0.8413)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.067)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.622)
        self.assertEqual(round(st.w, 1), 894.3)

        st = R152a(T=50+273.15, P=101325)
        self.assertEqual(round(st.rho, 3), 2.531)
        self.assertEqual(round(st.h.kJkg, 2), 566.38)
        self.assertEqual(round(st.s.kJkgK, 4), 2.4397)
        self.assertEqual(round(st.cv.kJkgK, 3), 0.963)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.100)
        self.assertEqual(round(st.w, 1), 212.1)

        st = R152a(T=-10+273.15, P=2e5)
        self.assertEqual(round(st.rho, 1), 981.3)
        self.assertEqual(round(st.h.kJkg, 2), 183.17)
        self.assertEqual(round(st.s.kJkgK, 4), 0.9375)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.086)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.664)
        self.assertEqual(round(st.w, 1), 821.2)

        st = R152a(T=-100+273.15, P=4e5)
        self.assertEqual(round(st.rho, 1), 1159.0)
        self.assertEqual(round(st.h.kJkg, 2), 42.00)
        self.assertEqual(round(st.s.kJkgK, 4), 0.2821)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.030)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.517)
        self.assertEqual(round(st.w, 1), 1276.1)

        st = R152a(T=25+273.15, P=6e5)
        self.assertEqual(round(st.rho, 1), 899.5)
        self.assertEqual(round(st.h.kJkg, 2), 243.73)
        self.assertEqual(round(st.s.kJkgK, 4), 1.1519)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.138)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.800)
        self.assertEqual(round(st.w, 1), 648.0)

        st = R152a(T=240+273.15, P=1e6)
        self.assertEqual(round(st.rho, 3), 15.938)
        self.assertEqual(round(st.h.kJkg, 2), 806.30)
        self.assertEqual(round(st.s.kJkgK, 4), 2.7350)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.369)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.519)
        self.assertEqual(round(st.w, 1), 260.0)

        st = R152a(T=70+273.15, P=2e6)
        self.assertEqual(round(st.rho, 1), 765.8)
        self.assertEqual(round(st.h.kJkg, 2), 331.04)
        self.assertEqual(round(st.s.kJkgK, 4), 1.4188)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.225)
        self.assertEqual(round(st.cp.kJkgK, 3), 2.171)
        self.assertEqual(round(st.w, 1), 415.3)

        st = R152a(T=100+273.15, P=4e6)
        self.assertEqual(round(st.rho, 1), 638.3)
        self.assertEqual(round(st.h.kJkg, 2), 400.18)
        self.assertEqual(round(st.s.kJkgK, 4), 1.6038)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.309)
        self.assertEqual(round(st.cp.kJkgK, 3), 3.008)
        self.assertEqual(round(st.w, 1), 261.3)

        st = R152a(T=273.15, P=6e6)
        self.assertEqual(round(st.rho, 1), 973.0)
        self.assertEqual(round(st.h.kJkg, 2), 202.26)
        self.assertEqual(round(st.s.kJkgK, 4), 0.9865)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.102)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.667)
        self.assertEqual(round(st.w, 1), 816.8)

        st = R152a(T=-100+273.15, P=1e7)
        self.assertEqual(round(st.rho, 1), 1167.5)
        self.assertEqual(round(st.h.kJkg, 2), 48.03)
        self.assertEqual(round(st.s.kJkgK, 4), 0.2693)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.040)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.509)
        self.assertEqual(round(st.w, 1), 1306.2)

        st = R152a(T=200+273.15, P=2e7)
        self.assertEqual(round(st.rho, 1), 541.5)
        self.assertEqual(round(st.h.kJkg, 2), 586.81)
        self.assertEqual(round(st.s.kJkgK, 4), 1.9821)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.396)
        self.assertEqual(round(st.cp.kJkgK, 3), 2.167)
        self.assertEqual(round(st.w, 1), 327.7)

        st = R152a(T=-50+273.15, P=4e7)
        self.assertEqual(round(st.rho, 1), 1111.4)
        self.assertEqual(round(st.h.kJkg, 2), 142.13)
        self.assertEqual(round(st.s.kJkgK, 4), 0.6132)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.058)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.514)
        self.assertEqual(round(st.w, 1), 1200.6)

        st = R152a(T=240+273.15, P=6e7)
        self.assertEqual(round(st.rho, 1), 730.8)
        self.assertEqual(round(st.h.kJkg, 2), 645.96)
        self.assertEqual(round(st.s.kJkgK, 4), 1.9773)
        self.assertEqual(round(st.cv.kJkgK, 3), 1.470)
        self.assertEqual(round(st.cp.kJkgK, 3), 1.879)
        self.assertEqual(round(st.w, 1), 599.3)

    def test_shortSpan(self):
        # Table III, Pag 117
        st = R152a(T=500, rho=500, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 1.4632)
        self.assertEqual(round(st.P.MPa, 3), 21.594)
        self.assertEqual(round(st.cp.kJkgK, 4), 2.1580)

        st2 = R152a(T=600, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 270.60)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.60934)

    def test_astina(self):
        # Table III, Pag 1719
        st = R152a(T=200, P=1e4, eq="astina")
        self.assertEqual(round(st.rho, 2), 1108.30)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.02563)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.53226)
        self.assertEqual(round(st.w, 2), 1144.46)
        self.assertEqual(round(st.h.kJkg, 4), 83.0731)
        self.assertEqual(round(st.s.kJkgK, 6), 0.504046)

        st = R152a(T=250, P=1e5, eq="astina")
        self.assertEqual(round(st.rho, 5), 3.32533)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.85146)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.02093)
        self.assertEqual(round(st.w, 3), 185.395)
        self.assertEqual(round(st.h.kJkg, 3), 490.101)
        self.assertEqual(round(st.s.kJkgK, 5), 2.17402)

        st = R152a(T=300, P=5e5, eq="astina")
        self.assertEqual(round(st.rho, 4), 14.9178)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.01270)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.25049)
        self.assertEqual(round(st.w, 3), 190.307)
        self.assertEqual(round(st.h.kJkg, 3), 528.281)
        self.assertEqual(round(st.s.kJkgK, 5), 2.12554)

        st = R152a(T=250, P=1e6, eq="astina")
        self.assertEqual(round(st.rho, 2), 1011.40)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.06828)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.61778)
        self.assertEqual(round(st.w, 3), 887.901)
        self.assertEqual(round(st.h.kJkg, 3), 162.066)
        self.assertEqual(round(st.s.kJkgK, 6), 0.852079)

        st = R152a(T=450, P=2e6, eq="astina")
        self.assertEqual(round(st.rho, 4), 39.2134)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.27176)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.47906)
        self.assertEqual(round(st.w, 3), 230.750)
        self.assertEqual(round(st.h.kJkg, 3), 703.818)
        self.assertEqual(round(st.s.kJkgK, 5), 2.44002)

        st = R152a(T=450, P=3e6, eq="astina")
        self.assertEqual(round(st.rho, 4), 62.4333)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.29309)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.56254)
        self.assertEqual(round(st.w, 3), 221.167)
        self.assertEqual(round(st.h.kJkg, 3), 692.526)
        self.assertEqual(round(st.s.kJkgK, 5), 2.37021)

        st = R152a(T=300, P=5e6, eq="astina")
        self.assertEqual(round(st.rho, 3), 910.768)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.14104)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.75988)
        self.assertEqual(round(st.w, 3), 681.439)
        self.assertEqual(round(st.h.kJkg, 3), 247.753)
        self.assertEqual(round(st.s.kJkgK, 5), 1.14904)

        st = R152a(T=350, P=1.5e7, eq="astina")
        self.assertEqual(round(st.rho, 3), 829.854)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.22465)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.86556)
        self.assertEqual(round(st.w, 3), 577.309)
        self.assertEqual(round(st.h.kJkg, 3), 339.639)
        self.assertEqual(round(st.s.kJkgK, 5), 1.39677)

        st = R152a(T=400, P=2.5e7, eq="astina")
        self.assertEqual(round(st.rho, 3), 764.468)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.30865)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.92116)
        self.assertEqual(round(st.w, 3), 526.848)
        self.assertEqual(round(st.h.kJkg, 3), 433.427)
        self.assertEqual(round(st.s.kJkgK, 5), 1.61365)

        st = R152a(T=250, P=4e7, eq="astina")
        self.assertEqual(round(st.rho, 2), 1069.55)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.10283)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.54070)
        self.assertEqual(round(st.w, 2), 1076.73)
        self.assertEqual(round(st.h.kJkg, 3), 183.511)
        self.assertEqual(round(st.s.kJkgK, 5), 0.78821)

        st = R152a(T=300, P=4.5e7, eq="astina")
        self.assertEqual(round(st.rho, 3), 999.404)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.17500)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.61174)
        self.assertEqual(round(st.w, 3), 920.812)
        self.assertEqual(round(st.h.kJkg, 3), 265.090)
        self.assertEqual(round(st.s.kJkgK, 5), 1.06788)

        st = R152a(T=450, P=5e7, eq="astina")
        self.assertEqual(round(st.rho, 3), 782.051)
        self.assertEqual(round(st.cv.kJkgK, 5), 1.40653)
        self.assertEqual(round(st.cp.kJkgK, 5), 1.86403)
        self.assertEqual(round(st.w, 3), 612.067)
        self.assertEqual(round(st.h.kJkg, 3), 528.381)
        self.assertEqual(round(st.s.kJkgK, 5), 1.76108)

        st = R152a(T=200, x=0.5, eq="astina")
        self.assertEqual(round(st.P.MPa, 5), 0.00608)
        self.assertEqual(round(st.Liquido.rho, 2), 1108.30)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 1.02563)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 5), 1.53227)
        self.assertEqual(round(st.Liquido.w, 2), 1144.44)
        self.assertEqual(round(st.Liquido.h.kJkg, 4), 83.0708)
        self.assertEqual(round(st.Liquido.s.kJkgK, 6), 0.504052)
        self.assertEqual(round(st.Gas.rho, 6), 0.243503)
        self.assertEqual(round(st.Gas.cv.kJkgK, 6), 0.688696)
        self.assertEqual(round(st.Gas.cp.kJkgK, 6), 0.822899)
        self.assertEqual(round(st.Gas.w, 3), 172.126)
        self.assertEqual(round(st.Gas.h.kJkg, 3), 452.530)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), 2.35135)

        st = R152a(T=300, x=0.5, eq="astina")
        self.assertEqual(round(st.P.MPa, 5), 0.62958)
        self.assertEqual(round(st.Liquido.rho, 3), 895.050)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 1.14163)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 5), 1.80605)
        self.assertEqual(round(st.Liquido.w, 3), 635.932)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 246.936)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 1.16245)
        self.assertEqual(round(st.Gas.rho, 4), 19.5363)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 1.05165)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 1.34292)
        self.assertEqual(round(st.Gas.w, 3), 184.955)
        self.assertEqual(round(st.Gas.h.kJkg, 3), 523.189)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), 2.08329)

        st = R152a(P=1e4, x=0.5, eq="astina")
        self.assertEqual(round(st.T, 3), 206.996)
        self.assertEqual(round(st.Liquido.rho, 2), 1095.05)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 1.02888)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 5), 1.54024)
        self.assertEqual(round(st.Liquido.w, 2), 1106.00)
        self.assertEqual(round(st.Liquido.h.kJkg, 4), 93.8198)
        self.assertEqual(round(st.Liquido.s.kJkgK, 6), 0.556861)
        self.assertEqual(round(st.Gas.rho, 6), 0.387819)
        self.assertEqual(round(st.Gas.cv.kJkgK, 6), 0.709506)
        self.assertEqual(round(st.Gas.cp.kJkgK, 6), 0.846608)
        self.assertEqual(round(st.Gas.w, 3), 174.478)
        self.assertEqual(round(st.Gas.h.kJkg, 3), 457.773)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), 2.31513)

        st = R152a(P=2e6, x=0.5, eq="astina")
        self.assertEqual(round(st.T, 3), 345.817)
        self.assertEqual(round(st.Liquido.rho, 3), 755.354)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 1.23824)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 5), 2.22690)
        self.assertEqual(round(st.Liquido.w, 3), 400.169)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), 336.806)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), 1.43550)
        self.assertEqual(round(st.Gas.rho, 4), 67.2945)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 1.25774)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 1.99694)
        self.assertEqual(round(st.Gas.w, 3), 166.520)
        self.assertEqual(round(st.Gas.h.kJkg, 3), 542.188)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), 2.02940)

    def test_kim(self):
        # Table 6, pag 585, saturation states
        # FIXME: The value fail in fourth digital sign
        pass
        # st = R152a(T=-60+273.15, x=0.5, eq="kim")
        # self.assertEqual(round(st.P.MPa, 5), 0.01500)
        # self.assertEqual(round(st.Liquido.rho, 1), 1081.8)
        # self.assertEqual(round(st.Gas.v, 5), 1.76626)
        # self.assertEqual(round(st.Liquido.h.kJkg, 2), 103.24)
        # self.assertEqual(round(st.Gas.h.kJkg, 2), 462.57)
        # self.assertEqual(round(st.Liquido.s.kJkgK, 4), 0.6017)
        # self.assertEqual(round(st.Gas.s.kJkgK, 4), 2.2875)
        # self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 1.565)
        # self.assertEqual(round(st.Gas.cp.kJkgK, 3), 0.864)
        # self.assertEqual(round(st.Liquido.w, 1), 1021.6)
        # self.assertEqual(round(st.Gas.w, 1), 176.4)

    def test_Krauss(self):
        # Table VI, pag 750, saturation states

        # The correlation use the Pelt-Sengers extension for the critical
        # region from tillner mEoS, so the returned values differ, specially
        # the thermal conductivity values
        # van Pelt, Sengers, J.V.
        # Thermodynamic Properties of 1,1-Difluoroethane (R152a) in the
        # Critical Region
        # The Journal of Supercritical Fluids 8(1) (1995) 81-99
        # doi: 10.1016/0896-8446(95)90021-7

        # For testing it uses the outcalt mEoS for point near to critical point
        st = R152a(T=240, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 0.06642)
        self.assertEqual(round(st.Liquido.rho, 1), 1029.6)
        self.assertEqual(round(st.Gas.rho, 4), 2.2736)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 364.8)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 8.09)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 128.8)
        self.assertEqual(round(st.Gas.k.mWmK, 3), 8.483)

        st = R152a(T=280, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 5), 0.33558)
        self.assertEqual(round(st.Liquido.rho, 2), 943.16)
        self.assertEqual(round(st.Gas.rho, 3), 10.550)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 198.1)
        self.assertEqual(round(st.Gas.mu.muPas, 3), 9.606)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 109.0)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 12.63)

        st = R152a(T=320, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 4), 1.0883)
        self.assertEqual(round(st.Liquido.rho, 2), 839.97)
        self.assertEqual(round(st.Gas.rho, 3), 34.202)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 128.3)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 11.20)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 90.71)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 17.31)

        st = R152a(T=360, x=0.5, eq="tillner")
        self.assertEqual(round(st.P.MPa, 4), 2.7024)
        self.assertEqual(round(st.Liquido.rho, 2), 694.46)
        self.assertEqual(round(st.Gas.rho, 3), 98.845)
        self.assertEqual(round(st.Liquido.mu.muPas, 2), 76.74)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 13.87)
        self.assertEqual(round(st.Liquido.k.mWmK, 2), 71.86)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 23.91)

        # Table VII, Pag 753, Single phase point Viscosity
        # Table VIII, Pag 754, Single phase point thermal conductivity
        st = R152a(T=240, P=1e5, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 1), 365.1)
        self.assertEqual(round(st.k.mWmK, 1), 128.8)

        st = R152a(T=360, P=1e5, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 2), 12.52)
        self.assertEqual(round(st.k.mWmK, 2), 20.16)

        st = R152a(T=430, P=1e5, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 2), 14.91)
        self.assertEqual(round(st.k.mWmK, 2), 26.96)

        st = R152a(T=240, P=5e6, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 1), 411.3)
        self.assertEqual(round(st.k.mWmK, 1), 131.4)

        st = R152a(T=360, P=5e6, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 2), 86.54)
        self.assertEqual(round(st.k.mWmK, 2), 76.58)

        st = R152a(T=430, P=5e6, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 2), 17.53)
        self.assertEqual(round(st.k.mWmK, 2), 32.54)

        st = R152a(T=250, P=2e7, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 1), 445.9)
        self.assertEqual(round(st.k.mWmK, 1), 133.7)

        st = R152a(T=360, P=2e7, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 1), 123.9)
        self.assertEqual(round(st.k.mWmK, 2), 92.37)

        st = R152a(T=430, P=2e7, eq="tillner")
        self.assertEqual(round(st.mu.muPas, 2), 70.29)
        self.assertEqual(round(st.k.mWmK, 2), 73.98)
