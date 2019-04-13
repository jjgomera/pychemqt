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

from lib import unidades
from lib.meos import MEoS


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

    # Used the correlation in Leachaman et al. paper
    # Younglove use a fractional correlation depend of temperature not easy
    # of implement here
    CP1 = {"ao": 2.5,
           "an": [], "pow": [],
           "ao_exp": [4.30256, 13.0289, -47.7365, 50.0013, -18.6261, 0.993973,
                      0.536078],
           "exp": [499, 826.5, 970.8, 1166.2, 1341.4, 5395, 10185]}

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

    younglove = {
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
        "ref": {"Tref": 300, "Pref": 101.325, "ho": 8466.7, "so": 130.5},
        "Tc": 32.938, "Pc": 1.28377, "rhoc": 15.556,

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 104.0,

        "gamma": -0.0041,
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

    eq = leachman, younglove

    _surface = {"sigma": [0.005314], "exp": [1.06]}
    _dielectric = {
        "eq": 1,
        "a": [2.0297, 0.0069], "b": [0.181, 0.021], "c": [-7.4, 0],
        "Au": 0, "D": 2}

    _melting = {
        "__doi__": younglove["__doi__"],
        "Tmin": Tt, "Tmax": 400.0}

    @classmethod
    def _Melting_Pressure(cls, T):
        """The correlation has two different equation so hardcoded here"""
        if T > 22:
            P = -26.5289115 + 0.248578596*T**1.764739
        else:
            P = -21.272389 + 0.125746643*T**1.955

        return unidades.Pressure(P, "MPa")

    _sublimation = {
        "eq": 2,
        "__doi__": {
            "autor": "McCarty, R.D., Hord, J., Roder, H.M.",
            "title": "Selected Properties of Hydrogen (Engineering Design "
                     "Data)",
            "ref": "NBS Monograph 168, NBS 1981.",
            "doi": ""},

        "Tmin": Tt, "Tmax": Tt,
        # The returned pressure is in mmHg
        "Tref": 1, "Pref": 101325/760,
        "a0": 4.009857354,
        "a1": [-90.77568949], "exp1": [-1],
        "a3": [2.48983094], "exp3": [1]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.487767e1, 0.103359e1, 0.826680, -0.129412],
        "t": [1.0, 1.5, 2.65, 7.4]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.13509, 0.40739e1, -0.53985e1, 0.55230e1, -0.23643e1],
        "t": [0.15, 0.44, 0.7, 0.99, 1.31]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.57545e1, 0.38153e1, -0.12293e2, 0.15095e2, -0.17295e2,
              -0.34190e2],
        "t": [0.53, 0.7, 1.7, 2.4, 3.3, 10]}

    visco0 = {"__name__": "McCarty (1972)",
              "__doi__": {
                  "autor": "McCarty, R.D. and Weber, L.A.",
                  "title": "Thermophysical Properties of Parahydrogen from "
                           "the Freezing Liquid Line to 5000 R for Pressures "
                           "to 10000 psia",
                  "ref": "NBS Technical Note 617",
                  "doi": ""},

              "eq": 0,
              "method": "_visco0"}

    def _visco0(self, rho, T, fase):
        def muo(T):
            b = [-0.1841091042788e2, 0.3185762039455e2, -0.2308233586574e2,
                 0.9129812714730e1, -0.2163626387630e1, 0.3175128582601,
                 -0.2773173035271e-1, 0.1347359367871e-2, -0.2775671778154e-4]
            suma = 0
            for i, b in enumerate(b):
                suma += b*T**((-3.+i)/3)
            return suma*100

        def mu1(rho, T):
            A = exp(5.7694 + log(rho.gcc) + 0.65e2*rho.gcc**1.5 -
                    6e-6*exp(127.2*rho.gcc))
            B = 10 + 7.2*((rho.gcc/0.07)**6-(rho.gcc/0.07)**1.5) - \
                17.63*exp(-58.75*(rho.gcc/0.07)**3)
            return A*exp(B/T)*0.1

        def mu2(rho, T):
            c = [-0.1324266117873e2, 0.1895048470537e2, 0.2184151514282e2,
                 0.9771827164811e5, -0.1157010275059e4, 0.1911147702539e3,
                 -0.3186427506942e4, 0.0705565000000]
            R2 = rho.gcc**0.5*(rho.gcc-c[7])/c[7]
            A = c[0] + c[1]*R2 + c[2]*rho.gcc**0.1 + c[3]*R2/T**2 + \
                c[4]*rho.gcc**0.1/T**1.5 + c[5]/T + c[6]*R2/T
            B = c[0]+c[5]/T
            return 0.1*(exp(A)-exp(B))

        def mur(rho1, T1, rho2, T2):
            return muo(T1) + mu2(rho1, T2) - muo(T2) - mu2(rho2, T2)

        if T > 100:
            mu = muo(100) + mu1(rho, 100) + mur(rho, T, rho, 100)
        else:
            mu = muo(T)+mu1(rho, T)
        return unidades.Viscosity(mu, "muPas")

    _viscosity = visco0,

    thermo0 = {"__name__": "Assael (2011)",
               "__doi__": {
                   "autor": "Assael, M.J., Assael. J.-A.M., Huber, M.L., "
                            "Perkins, R.A., Takata, Y.",
                   "title": "Correlation of the Thermal Conductivity of "
                            "Normal and Parahydrogen from the Triple Point to "
                            "1000 K and up to 100 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 40(3) (2011) 033101",
                   "doi": "10.1063/1.3606499"},

               "eq": 1,

               "Toref": 32.938, "koref": 1,
               "no_num": [-1.245, 3.10212e2, -3.31004e2, 2.46016e2, -6.57810e1,
                          1.08260e1, -5.19659e-1, 1.43979e-2],
               "to_num": [0, 1, 2, 3, 4, 5, 6, 7],
               "no_den": [1.42304e4, -1.93922e4, 1.58379e4, -4.81812e3,
                          7.28639e2, -3.57365e1, 1],
               "to_den": [0, 1, 2, 3, 4, 5, 6],

               "Tref_res": 32.938, "rhoref_res": 31.323, "kref_res": 1.,
               "nr": [2.65975e-2, -1.33826e-3, 1.30219e-2, -5.67678e-3,
                      -9.2338e-5, -1.21727e-3, 3.66663e-3, 3.88715e-3,
                      -9.21055e-3, 4.00723e-3],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.15e-9, "gam0": 0.052, "qd": 0.5e-9, "Tcref": 49.407}

    _thermal = thermo0,


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

    def test_younglove(self):
        # Selected point from Appendix H. Pag 97
        st = pH2(T=40, P=5e3, eq="younglove")
        self.assertEqual(round(st.rho, 5), 0.03033)
        self.assertEqual(round(st.rhoM, 5), 0.01505)
        self.assertEqual(round(st.uM.Jmol, 1), 499.0)
        self.assertEqual(round(st.hM.Jmol, 1), 831.3)
        self.assertEqual(round(st.sM.JmolK, 1), 100.7)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.49)
        self.assertEqual(round(st.cpM.JmolK, 2), 20.83)
        self.assertEqual(round(st.w, 1), 524.0)

        st = pH2(T=110, P=1e4, eq="younglove")
        self.assertEqual(round(st.rho, 5), 0.02204)
        self.assertEqual(round(st.rhoM, 5), 0.01093)
        self.assertEqual(round(st.uM.Jmol, 0), 1569)
        self.assertEqual(round(st.hM.Jmol, 0), 2483)
        self.assertEqual(round(st.sM.JmolK, 1), 118.1)
        self.assertEqual(round(st.cvM.JmolK, 2), 20.37)
        self.assertEqual(round(st.cpM.JmolK, 2), 28.69)
        self.assertEqual(round(st.w, 1), 799.4)

        st = pH2(T=400, P=2e4, eq="younglove")
        self.assertEqual(round(st.rho, 5), 0.01212)
        self.assertEqual(round(st.rhoM, 6), 0.006013)
        self.assertEqual(round(st.uM.Jmol, 0), 8093)
        self.assertEqual(round(st.hM.Jmol, 0), 11419)
        self.assertEqual(round(st.sM.JmolK, 1), 152.5)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.02)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.33)
        self.assertEqual(round(st.w, 0), 1518)

        st = pH2(T=14, P=4e4, eq="younglove")
        self.assertEqual(round(st.rho, 2), 76.89)
        self.assertEqual(round(st.rhoM, 2), 38.14)
        self.assertEqual(round(st.uM.Jmol, 1), -619.7)
        self.assertEqual(round(st.hM.Jmol, 1), -618.6)
        self.assertEqual(round(st.sM.JmolK, 2), 10.12)
        self.assertEqual(round(st.cvM.JmolK, 2), 10.60)
        self.assertEqual(round(st.cpM.JmolK, 2), 15.16)
        self.assertEqual(round(st.w, 0), 1366)

        st = pH2(T=20, P=6e4, eq="younglove", rho0=0.7)
        self.assertEqual(round(st.rho, 4), 0.7705)
        self.assertEqual(round(st.rhoM, 4), 0.3822)
        self.assertEqual(round(st.uM.Jmol, 1), 236.2)
        self.assertEqual(round(st.hM.Jmol, 1), 393.2)
        self.assertEqual(round(st.sM.JmolK, 2), 64.89)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.93)
        self.assertEqual(round(st.cpM.JmolK, 2), 22.83)
        self.assertEqual(round(st.w, 1), 359.6)

        st = pH2(T=80, P=8e4, eq="younglove")
        self.assertEqual(round(st.rho, 4), 0.2427)
        self.assertEqual(round(st.rhoM, 4), 0.1204)
        self.assertEqual(round(st.uM.Jmol, 0), 1031)
        self.assertEqual(round(st.hM.Jmol, 0), 1695)
        self.assertEqual(round(st.sM.JmolK, 2), 92.50)
        self.assertEqual(round(st.cvM.JmolK, 2), 15.31)
        self.assertEqual(round(st.cpM.JmolK, 2), 23.71)
        self.assertEqual(round(st.w, 1), 713.9)

        st = pH2(T=20, P=1e5, eq="younglove")
        self.assertEqual(round(st.rho, 2), 71.12)
        self.assertEqual(round(st.rhoM, 2), 35.28)
        self.assertEqual(round(st.uM.Jmol, 1), -524.3)
        self.assertEqual(round(st.hM.Jmol, 1), -521.5)
        self.assertEqual(round(st.sM.JmolK, 2), 15.75)
        self.assertEqual(round(st.cvM.JmolK, 2), 11.32)
        self.assertEqual(round(st.cpM.JmolK, 2), 19.11)
        self.assertEqual(round(st.w, 0), 1111)

        st = pH2(T=300, P=101325, eq="younglove")
        self.assertEqual(round(st.rho, 5), 0.08184)
        self.assertEqual(round(st.rhoM, 5), 0.04060)
        self.assertEqual(round(st.uM.Jmol, 0), 5971)
        self.assertEqual(round(st.hM.Jmol, 0), 8467)
        self.assertEqual(round(st.sM.JmolK, 1), 130.5)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.61)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.93)
        self.assertEqual(round(st.w, 0), 1310)

        st = pH2(T=25, P=2e5, eq="younglove", rho0=2)
        self.assertEqual(round(st.rho, 3), 2.187)
        self.assertEqual(round(st.rhoM, 3), 1.085)
        self.assertEqual(round(st.uM.Jmol, 1), 277.2)
        self.assertEqual(round(st.hM.Jmol, 1), 461.6)
        self.assertEqual(round(st.sM.JmolK, 2), 58.77)
        self.assertEqual(round(st.cvM.JmolK, 2), 13.13)
        self.assertEqual(round(st.cpM.JmolK, 2), 25.15)
        self.assertEqual(round(st.w, 1), 391.4)

        st = pH2(T=400, P=3e5, eq="younglove")
        self.assertEqual(round(st.rho, 4), 0.1816)
        self.assertEqual(round(st.rhoM, 5), 0.09007)
        self.assertEqual(round(st.uM.Jmol, 0), 8093)
        self.assertEqual(round(st.hM.Jmol, 0), 11423)
        self.assertEqual(round(st.sM.JmolK, 1), 130.0)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.02)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.34)
        self.assertEqual(round(st.w, 0), 1520)

        st = pH2(T=25, P=4e5, eq="younglove")
        self.assertEqual(round(st.rho, 2), 64.65)
        self.assertEqual(round(st.rhoM, 2), 32.07)
        self.assertEqual(round(st.uM.Jmol, 1), -415.3)
        self.assertEqual(round(st.hM.Jmol, 1), -402.9)
        self.assertEqual(round(st.sM.JmolK, 2), 20.61)
        self.assertEqual(round(st.cvM.JmolK, 2), 12.64)
        self.assertEqual(round(st.cpM.JmolK, 2), 27.19)
        self.assertEqual(round(st.w, 0), 927)

        st = pH2(T=14, P=5e5, eq="younglove")
        self.assertEqual(round(st.rho, 2), 77.24)
        self.assertEqual(round(st.rhoM, 2), 38.32)
        self.assertEqual(round(st.uM.Jmol, 1), -621.5)
        self.assertEqual(round(st.hM.Jmol, 1), -608.4)
        self.assertEqual(round(st.sM.JmolK, 2), 9.99)
        self.assertEqual(round(st.cvM.JmolK, 2), 10.90)
        self.assertEqual(round(st.cpM.JmolK, 2), 15.10)
        self.assertEqual(round(st.w, 0), 1356)

        st = pH2(T=30, P=6e5, eq="younglove", rho0=6)
        self.assertEqual(round(st.rho, 3), 6.368)
        self.assertEqual(round(st.rhoM, 3), 3.159)
        self.assertEqual(round(st.uM.Jmol, 1), 278.9)
        self.assertEqual(round(st.hM.Jmol, 1), 468.8)
        self.assertEqual(round(st.sM.JmolK, 2), 51.47)
        self.assertEqual(round(st.cvM.JmolK, 2), 13.79)
        self.assertEqual(round(st.cpM.JmolK, 2), 34.16)
        self.assertEqual(round(st.w, 1), 406.1)

        st = pH2(T=400, P=8e5, eq="younglove")
        self.assertEqual(round(st.rho, 4), 0.4830)
        self.assertEqual(round(st.rhoM, 4), 0.2396)
        self.assertEqual(round(st.uM.Jmol, 0), 8092)
        self.assertEqual(round(st.hM.Jmol, 0), 11430)
        self.assertEqual(round(st.sM.JmolK, 1), 121.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.03)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.35)
        self.assertEqual(round(st.w, 0), 1523)

        st = pH2(T=28.5, P=1e6, eq="younglove")
        self.assertEqual(round(st.rho, 2), 59.48)
        self.assertEqual(round(st.rhoM, 2), 29.51)
        self.assertEqual(round(st.uM.Jmol, 1), -326.1)
        self.assertEqual(round(st.hM.Jmol, 1), -292.2)
        self.assertEqual(round(st.sM.JmolK, 2), 24.02)
        self.assertEqual(round(st.cvM.JmolK, 2), 13.15)
        self.assertEqual(round(st.cpM.JmolK, 2), 35.43)
        self.assertEqual(round(st.w, 1), 811.5)

        st = pH2(T=85, P=1.5e6, eq="younglove")
        self.assertEqual(round(st.rho, 3), 4.334)
        self.assertEqual(round(st.rhoM, 3), 2.150)
        self.assertEqual(round(st.uM.Jmol, 0), 1058)
        self.assertEqual(round(st.hM.Jmol, 0), 1755)
        self.assertEqual(round(st.sM.JmolK, 2), 68.99)
        self.assertEqual(round(st.cvM.JmolK, 2), 16.22)
        self.assertEqual(round(st.cpM.JmolK, 2), 25.84)
        self.assertEqual(round(st.w, 1), 738.8)

        st = pH2(T=120, P=2e6, eq="younglove")
        self.assertEqual(round(st.rho, 3), 4.009)
        self.assertEqual(round(st.rhoM, 3), 1.989)
        self.assertEqual(round(st.uM.Jmol, 0), 1733)
        self.assertEqual(round(st.hM.Jmol, 0), 2739)
        self.assertEqual(round(st.sM.JmolK, 2), 76.23)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.96)
        self.assertEqual(round(st.cpM.JmolK, 2), 31.04)
        self.assertEqual(round(st.w, 1), 843.7)

        st = pH2(T=340, P=3e6, eq="younglove")
        self.assertEqual(round(st.rho, 3), 2.104)
        self.assertEqual(round(st.rhoM, 3), 1.044)
        self.assertEqual(round(st.uM.Jmol, 0), 6816)
        self.assertEqual(round(st.hM.Jmol, 0), 9690)
        self.assertEqual(round(st.sM.JmolK, 1), 106.0)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.29)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.67)
        self.assertEqual(round(st.w, 0), 1421)

        st = pH2(T=80, P=4e6, eq="younglove")
        self.assertEqual(round(st.rho, 2), 12.56)
        self.assertEqual(round(st.rhoM, 3), 6.230)
        self.assertEqual(round(st.uM.Jmol, 1), 878.8)
        self.assertEqual(round(st.hM.Jmol, 0), 1521)
        self.assertEqual(round(st.sM.JmolK, 2), 58.17)
        self.assertEqual(round(st.cvM.JmolK, 2), 15.61)
        self.assertEqual(round(st.cpM.JmolK, 2), 27.86)
        self.assertEqual(round(st.w, 1), 750.5)

        st = pH2(T=16, P=5e6, eq="younglove")
        self.assertEqual(round(st.rho, 2), 79.16)
        self.assertEqual(round(st.rhoM, 2), 39.27)
        self.assertEqual(round(st.uM.Jmol, 1), -607.8)
        self.assertEqual(round(st.hM.Jmol, 1), -480.4)
        self.assertEqual(round(st.sM.JmolK, 2), 10.79)
        self.assertEqual(round(st.cvM.JmolK, 2), 10.22)
        self.assertEqual(round(st.cpM.JmolK, 2), 13.25)
        self.assertEqual(round(st.w, 0), 1348)

        st = pH2(T=400, P=1e7, eq="younglove")
        self.assertEqual(round(st.rho, 3), 5.783)
        self.assertEqual(round(st.rhoM, 3), 2.869)
        self.assertEqual(round(st.uM.Jmol, 0), 8073)
        self.assertEqual(round(st.hM.Jmol, 0), 11559)
        self.assertEqual(round(st.sM.JmolK, 1), 100.8)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.16)
        self.assertEqual(round(st.cpM.JmolK, 2), 29.57)
        self.assertEqual(round(st.w, 0), 1592)

        st = pH2(T=90, P=2e7, eq="younglove")
        self.assertEqual(round(st.rho, 2), 43.51)
        self.assertEqual(round(st.rhoM, 2), 21.58)
        self.assertEqual(round(st.uM.Jmol, 1), 718.5)
        self.assertEqual(round(st.hM.Jmol, 0), 1645)
        self.assertEqual(round(st.sM.JmolK, 2), 45.65)
        self.assertEqual(round(st.cvM.JmolK, 2), 18.04)
        self.assertEqual(round(st.cpM.JmolK, 2), 31.21)
        self.assertEqual(round(st.w, 0), 1135)

        st = pH2(T=400, P=4e7, eq="younglove")
        self.assertEqual(round(st.rho, 2), 20.32)
        self.assertEqual(round(st.rhoM, 2), 10.08)
        self.assertEqual(round(st.uM.Jmol, 0), 8026)
        self.assertEqual(round(st.hM.Jmol, 0), 11994)
        self.assertEqual(round(st.sM.JmolK, 2), 89.13)
        self.assertEqual(round(st.cvM.JmolK, 2), 21.51)
        self.assertEqual(round(st.cpM.JmolK, 2), 30.04)
        self.assertEqual(round(st.w, 0), 1812)

        st = pH2(T=40, P=1e8, eq="younglove", rho0=100)
        self.assertEqual(round(st.rho, 1), 100.8)
        self.assertEqual(round(st.rhoM, 2), 49.99)
        self.assertEqual(round(st.uM.Jmol, 1), -252.1)
        self.assertEqual(round(st.hM.Jmol, 0), 1748)
        self.assertEqual(round(st.sM.JmolK, 2), 15.64)
        self.assertEqual(round(st.cvM.JmolK, 2), 14.25)
        self.assertEqual(round(st.cpM.JmolK, 2), 16.68)
        self.assertEqual(round(st.w, 0), 2394)

    def test_Assael(self):
        # Table 7, Pag 11
        self.assertEqual(round(pH2(T=298.15, rho=0).k.mWmK, 2), 192.38)
        self.assertEqual(round(pH2(T=298.15, rho=0.80844).k.mWmK, 2), 192.80)
        self.assertEqual(round(pH2(T=298.15, rho=14.4813).k.mWmK, 2), 207.85)
        self.assertEqual(round(pH2(T=35, rho=0).k.mWmK, 3), 27.222)
        self.assertEqual(round(pH2(T=35, rho=30).k.mWmK, 3), 70.334)
        self.assertEqual(round(pH2(T=18, rho=0).k.mWmK, 3), 13.643)
        self.assertEqual(round(pH2(T=18, rho=75).k.mWmK, 2), 100.52)
