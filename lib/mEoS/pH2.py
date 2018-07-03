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

from scipy import exp, log

from lib.meos import MEoS
from lib import unidades


class pH2(MEoS):
    """Multiparameter equation of state for hydrogen (para)"""
    name = "parahydrogen"
    CASNumber = "1333-74-0p"
    formula = "H2"
    synonym = "R-702p"
    _refPropName = "PARAHYD"
    _coolPropName = "ParaHydrogen"
    rhoc = unidades.Density(31.32274344)
    Tc = unidades.Temperature(32.938)
    Pc = unidades.Pressure(1285.8, "kPa")
    M = 2.01588  # g/mol
    Tt = unidades.Temperature(13.8033)
    Tb = unidades.Temperature(20.271)
    f_acent = -0.219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-1.4485891134, 1.884521239],
           "ao_exp": [4.30256, 13.0289, -47.7365, 50.0013, -18.6261, 0.993973,
                      0.536078],
           "titao": [15.1496751472, 25.0925982148, 29.4735563787,
                     35.4059141417, 40.724998482, 163.7925799988,
                     309.2173173842]}

    CP1 = {"ao": 2.4995169,
           "an": [-0.11125185e-2, 0.27491461e-3, -0.10005269e-4, 0.22695404e-8,
                  -0.21031029e-12],
           "pow": [1, 1.5, 2, 3, 4],
           "ao_exp": [0.12353388e2, -0.17777676e2, 0.64309174e1, 0.73347521e1],
           "exp": [598, 778, 1101, 6207],
           "ao_hyp": [], "hyp": []}

    leachman = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for parahydrogen of Leachman "
                    "et al. (2007)",
        "__doi__": {"autor": "Leachman, J.W., Jacobsen, R.T, Penoncello, S.G.,"
                             " Lemmon, E.W.",
                    "title": "Fundamental equations of state for Parahydrogen,"
                             " Normal Hydrogen, and Orthohydrogen",
                    "ref": "J. Phys. Chem. Ref. Data, 38(3) (2009) 721-748",
                    "doi": "10.1063/1.3160306"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 2000000.0, "rhomax": 104.0,
        "Pmin": 7.041, "rhomin": 38.185,

        "nr1": [-7.33375, .01, 2.60375, 4.66279, 0.682390, -1.47078, 0.135801],
        "d1": [1, 4, 1, 1, 2, 2, 3],
        "t1": [0.6855, 1, 1, 0.489, 0.774, 1.133, 1.386],

        "nr2": [-1.05327, 0.328239],
        "d2": [1, 3],
        "t2": [1.619, 1.162],
        "c2": [1, 1],
        "gamma2": [1]*2,

        "nr3": [-0.0577833, 0.0449743, 0.0703464, -0.0401766, 0.119510],
        "d3": [2, 1, 3, 1, 1],
        "t3": [3.96, 5.276, 0.99, 6.791, 3.19],
        "alfa3": [1.7437, 0.5516, 0.0634, 2.1341, 1.777],
        "beta3": [0.194, 0.2019, 0.0301, 0.2383, 0.3253],
        "gamma3": [0.8048, 1.5248, 0.6648, 0.6832, 1.493],
        "epsilon3": [1.5487, 0.1785, 1.28, 0.6319, 1.7104],
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for parahydrogen of Younglove "
                    "(1982)",
        "__doi__": {"autor": "Younglove, B.A.",
                    "title": "Thermophysical Properties of Fluids. I. Argon, "
                             "Ethylene, Parahydrogen, Nitrogen, Nitrogen "
                             "Trifluoride, and Oxygen",
                    "ref": "J. Phys. Chem. Ref. Data, 11(Suppl. 1) (1982)",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP1,
        "ref": "IIR",
        "Tc": 32.938, "Pc": 1.28377, "rhoc": 15.556,

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 44.0,
        "Pmin": 7.042, "rhomin": 38.21,

        "b": [None, 0.4675528393416e-3, 0.4289274251454e-1, -0.5164085596504,
              0.2961790279801e1, -0.3027194968412e2, 0.1908100320379e-4,
              -0.1339776859288e-2, 0.3056473115421, 0.5161197159532e2,
              0.1999981550224e-6, 0.2896367059356e-3, -0.2257803939041e-1,
              -0.2287392761826e-5, 0.2446261478645e-4, -0.1718181601119e-2,
              -0.5465142603459e-6, 0.4051941401315e-8, 0.1157595123961e-5,
              -0.1269162728389e-7, -0.4983023605519e2, -0.1606676092098e3,
              -0.1926799185310, 0.9319894638928e1, -0.3222596554434e-3,
              0.1206839307669e-2, -0.3841588197470e-6, -0.4036157453608e-4,
              -0.1250868123513e-9, 0.1976107321888e-8, -0.2411883474011e-12,
              -0.4127551498251e-12, 0.8917972883610e-11]}

    # eq = leachman, MBWR
    eq = leachman,

    _surface = {"sigma": [0.005314], "exp": [1.06]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [2.0297, 0.0069], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [0.181, 0.021, -7.4],
                   "expt2": [0, 1, 0], "expd2": [2, 2, 3]}
    _melting = {"Tmin": Tt, "Tmax": 400.0}
    _sublimation = {"eq": 2, "Tref": 1, "Pref": 0.133332237,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [4.009857354, -90.77568949], "exp1": [0, -1],
                    "a2": [], "exp2": [], "a3": [2.48983094], "exp3": [1]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.487767e1, 0.103359e1, 0.826680, -0.129412],
        "exp": [1.0, 1.5, 2.65, 7.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.13509, 0.40739e1, -0.53985e1, 0.55230e1, -0.23643e1],
        "exp": [0.15, 0.44, 0.7, 0.99, 1.31]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.57545e1, 0.38153e1, -0.12293e2, 0.15095e2, -0.17295e2, -0.34190e2],
        "exp": [0.53, 0.7, 1.7, 2.4, 3.3, 10]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "McCarty (1972)",
              "__doi__": {"autor": "McCarty, R.D. and Weber, L.A.",
                          "title": "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 psia",
                          "ref": "Natl. Bur. Stand., Tech. Note 617, 1972.",
                          "doi": ""}}

    _viscosity = visco0,

    def _visco0(self, rho, T, fase):
        DELV = lambda rho1, T1, rho2, T2: DILV(T1)+EXCESV(rho1, T2)-DILV(T2)-EXCESV(rho2, T2)

        def EXVDIL(rho, T):
            A = exp(5.7694+log(rho.gcc)+0.65e2*rho.gcc**1.5-6e-6*exp(127.2*rho.gcc))
            B = 10+7.2*((rho.gcc/0.07)**6-(rho.gcc/0.07)**1.5)-17.63*exp(-58.75*(rho.gcc/0.07)**3)
            return A*exp(B/T)*0.1

        def DILV(T):
            b = [-0.1841091042788e2, 0.3185762039455e2, -0.2308233586574e2,
                 0.9129812714730e1, -0.2163626387630e1, 0.3175128582601,
                 -0.2773173035271e-1, 0.1347359367871e-2, -0.2775671778154e-4]
            suma = 0
            for i, b in enumerate(b):
                suma += b*T**((-3.+i)/3)
            return suma*100

        def EXCESV(rho, T):
            c = [-0.1099981128e2, 0.1895876508e2, -0.3813005056e3,
                 0.5950473265e2, 0.1099399458e1, 0.8987269839e1,
                 0.1231422148e4, 0.311]
            R2 = rho.gcc**0.5*(rho.gcc-c[7])/c[7]
            A = c[0]+c[1]*R2+c[2]*rho.gcc**0.1+c[3]*R2/T**2+c[4]*rho.gcc**0.1/T**1.5+c[5]/T+c[6]*R2/T
            B = c[0]+c[5]/T
            return 0.1*(exp(A)-exp(B))

        if T > 100:
            n = DILV(100)+EXVDIL(rho, 100)+DELV(rho, T, rho, 100)
        else:
            n = DILV(T)+EXVDIL(rho, T)
        return unidades.Viscosity(n, "muPas")

    @classmethod
    def _Melting_Pressure(cls, T=None):
        Tref = 1
        Pref = 1000
        Tita = T/Tref
#        a = [-0.265289115e2, 0.248578596, -0.21272389e2, 0.125746643]
#        expo = [0, 0.1764739e1, 0, 0.1955e1]
        if T > 22:
            suma = -0.265289115e2+0.248578596*Tita**0.1764739e1
        else:
            suma = -0.21272389e2+0.125746643*Tita**0.1955e1

        return unidades.Pressure(suma*Pref, "kPa")

    thermo0 = {"eq": 1,
               "__name__": "Assael (2011)",
               "__doi__": {"autor": " Assael, M.J., Assael. J.-A.M., Huber, M.L., Perkins, R.A. and Takata, Y.",
                           "title": "Correlation of the Thermal Conductivity of Normal and Parahydrogen from the Triple Point to 1000 K and up to 100 MPa",
                           "ref": "J. Phys. Chem. Ref. Data 40, 033101 (2011)",
                           "doi": "10.1063/1.3606499"},
               "__test__": """
                   >>> st=pH2(T=298.15, rho=0)
                   >>> print "%0.5g" % st.k.mWmK
                   192.38
                   >>> st=pH2(T=298.15, rho=0.80844)
                   >>> print "%0.5g" % st.k.mWmK
                   192.81
                   >>> st=pH2(T=298.15, rho=14.4813)
                   >>> print "%0.5g" % st.k.mWmK
                   207.85
                   >>> st=pH2(T=35, rho=0)
                   >>> print "%0.5g" % st.k.mWmK
                   27.222
                   >>> st=pH2(T=35, rho=30)
                   >>> print "%0.5g" % st.k.mWmK
                   70.335
                   >>> st=pH2(T=35, rho=30)
                   >>> print "%0.5g" % st.k.mWmK
                   68.611
                   >>> st=pH2(T=18, rho=0)
                   >>> print "%0.5g" % st.k.mWmK
                   13.643
                   >>> st=pH2(T=18, rho=75)
                   >>> print "%0.5g" % st.k.mWmK
                   100.52
                   """, # Table 4, Pag 8

               "Tref": 1.0, "kref": 1e-3,
               "no": [-1.24500e3, 9.41806e3, -3.05098e2, 6.88449,
                      -5.58871e-2, 2.79243e-4, -4.06944e-7, 3.42309e-10],
               "co": [0, 1, 2, 3, 4, 5, 6, 7],
               "noden": [1.42304e4, -5.88749e2, 1.45983e1, -1.34830e-1,
                         6.19047e-4, -9.21777e-7, 7.83099e-10],
               "coden": [0, 1, 2, 3, 4, 5, 6],

               "Trefb": 32.938, "rhorefb": 15.538, "krefb": 1.,
               "nb": [0.265975e-1, -0.133826e-2, 0.130219e-1, -0.567678e-2,
                      -0.92338e-4, -0.121727e-2, 0.366663e-2, 0.388715e-2,
                      -0.921055e-2, 0.400723e-2],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.15e-9, "gam0": 0.052, "qd": 0.5e-9, "Tcref": 49.407}

    thermo1 = {"eq": 0,
               "method": "_thermo1",
               "__name__": "McCarty (1972)",
               "__doi__": {"autor": "McCarty, R.D. and Weber, L.A.",
                           "title": "Thermophysical properties of parahydrogen from the freezing liquid line to 5000 R for pressures to 10,000 psia",
                           "ref": "Natl. Bur. Stand., Tech. Note 617, 1972.",
                           "doi": ""}}

    def _thermo1(self, rho, T, fase):
        # TODO:
        return unidades.ThermalConductivity(0)

    # _thermal = thermo0, thermo1


class Test(TestCase):

    def test_leachman(self):
        # Selected point from Table 13, Pag 745, saturation states
        st = pH2(T=13.8033, x=0.5)
        self.assertEqual(round(st.P.kPa, 3), 7.0410)
        self.assertEqual(round(st.Liquido.rho, 3), 76.977)
        self.assertEqual(round(st.Gas.rho, 5), 0.12555)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), -53.741)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 396.31)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), -3.0840)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 29.521)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 5.1313)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 6.2265)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 6.9241)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 10.534)
        self.assertEqual(round(st.Liquido.w, 1), 1263.1)
        self.assertEqual(round(st.Gas.w, 2), 305.65)

        st = pH2(T=20, x=0.5)
        self.assertEqual(round(st.P.kPa, 3), 93.414)
        self.assertEqual(round(st.Liquido.rho, 3), 71.135)
        self.assertEqual(round(st.Gas.rho, 4), 1.2440)
        self.assertEqual(round(st.Liquido.h.kJkg, 4), -2.6915)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 444.54)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.12814)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 22.234)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 5.6371)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 6.4499)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 9.5688)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 11.920)
        self.assertEqual(round(st.Liquido.w, 1), 1118.6)
        self.assertEqual(round(st.Gas.w, 2), 353.50)

        st = pH2(T=30, x=0.5)
        self.assertEqual(round(st.P.kPa, 2), 823.19)
        self.assertEqual(round(st.Liquido.rho, 3), 53.976)
        self.assertEqual(round(st.Gas.rho, 3), 10.871)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 144.24)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 435.71)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 5.2108)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 14.926)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 6.4715)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 7.6246)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 26.649)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 32.583)
        self.assertEqual(round(st.Liquido.w, 2), 693.04)
        self.assertEqual(round(st.Gas.w, 2), 377.20)

        st = pH2(T=32, x=0.5)
        self.assertEqual(round(st.P.kPa, 1), 1120.3)
        self.assertEqual(round(st.Liquido.rho, 3), 45.901)
        self.assertEqual(round(st.Gas.rho, 3), 17.492)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), 204.08)
        self.assertEqual(round(st.Gas.h.kJkg, 2), 392.26)
        self.assertEqual(round(st.Liquido.s.kJkgK, 4), 6.9452)
        self.assertEqual(round(st.Gas.s.kJkgK, 3), 12.826)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 7.0097)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 8.3364)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 68.189)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 92.392)
        self.assertEqual(round(st.Liquido.w, 2), 522.59)
        self.assertEqual(round(st.Gas.w, 2), 372.03)
