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


class SF6(MEoS):
    """Multiparameter equation of state for sulfur hexafluoride"""
    name = "sulfur hexafluoride"
    CASNumber = "2551-62-4"
    formula = "SF6"
    synonym = ""
    _refPropName = "SF6"
    _coolPropName = "SulfurHexafluoride"
    rhoc = unidades.Density(742.3)
    Tc = unidades.Temperature(318.7232)
    Pc = unidades.Pressure(3754.983, "kPa")
    M = 146.0554192  # g/mol
    Tt = unidades.Temperature(223.555)
    Tb = unidades.Temperature(204.9)
    f_acent = 0.21
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 953

    _Tr = unidades.Temperature(304.013497)
    _rhor = unidades.Density(747.815849)
    _w = 0.181815238

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [11.638611086, -6.392241811],
           "ao_exp": [3.66118232, 7.87885103, 3.45981679],
           "titao": [1.617282065, 2.747115139, 4.232907175]}

    CP1 = {"ao": 3.9837756784,
           "ao_exp": [2.2181851010, -1.0921337374e1, 3.3102497939,
                      17.5189671483, 2.8903523803],
           "exp": [1114.38, 925.64, 499.26, 884.9, 1363.93]}

    CP2 = {"ao": -0.376915e-1/8.3143*146.05,
           "an": [0.305814e-2/8.3143*146.05, -0.237654e-5/8.3143*146.05],
           "pow": [1, 2]}

    guder = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur hexafluoride of "
                    "Guder and Wagner (2009)",
        "__doi__": {"autor": "Guder, C., Wagner, W.",
                    "title": "A Reference Equation of State for the "
                             "Thermodynamic Properties of Sulfur Hexafluoride "
                             "(SF6) for Temperatures from the Melting Line to "
                             "625 K and Pressures up to 150 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 38(1) (2009) 33-94",
                    "doi": "10.1063/1.3037344"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 625.0, "Pmax": 150000.0, "rhomax": 14.5,

        "nr1": [.54958259132835, -.87905033269396, -.84656969731452,
                .27692381593529, -.49864958372345e01, .48879127058055e01,
                .36917081634281e-1, .37030130305087e-3, .39389132911585e-1,
                .42477413690006e-3],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 4, 6],
        "t1": [0.125, 1.25, 1.875, 0.125, 1.5, 1.625, 1.5, 5.625, 0.625, 0.25],

        "nr2": [-.24150013863890e-1, .59447650642255e-1, -.38302880142267,
                .32606800951983, -.29955940562031e-1, -.86579186671173e-1,
                .41600684707562e01, -.41398128855814e01, -.55842159922714,
                .56531382776891, .82612463415545e-2, -.10200995338080e-1],
        "d2": [1, 2, 2, 2, 3, 6, 2, 2, 4, 4, 2, 2],
        "t2": [6, 0.25, 4.75, 5.375, 5.875, 2, 5.875, 6, 5.625, 5.75, 0, 0.5],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*12,

        "nr3": [-.21662523861406e-1, .34650943893908e-1, -.28694281385812e-1,
                .84007238998053e-2, -.26969359922498, .90415215646344e01,
                -.37233103557977e01, -.27524670823704e04, .57711861697319e04,
                -.30234003119748e04, .22252778435360e07, -.23056065559032e07,
                .63918852944475e07, -.60792091415592e07],
        "d3": [1, 3, 4, 1, 1, 4, 3, 4, 4, 4, 1, 1, 3, 3],
        "t3": [4, 1, 3, 2, 4, 3, 4, 1, 2, 3, 3, 4, 3, 4],
        "alfa3": [10, 10, 10, 10, 11, 25, 30, 30, 30, 30, 30, 30, 30, 30],
        "beta3": [150, 150, 150, 150, 225, 300, 350, 350, 350, 350, 400, 400,
                  400, 400],
        "gamma3": [1.13, 1.13, 1.13, 1.16, 1.19, 1.19, 1.16, 1.16, 1.16, 1.16,
                   1.22, 1.22, 1.22, 1.22],
        "epsilon3": [0.85, 0.85, 0.85, 0.85, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]}

    reuck = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur hexafluoride of "
                    "de Reuck (1991)",
        "__doi__": {"autor": "de Reuck, K.M., Craven, R.J.B., Cole, W.A.",
                    "title": "Report on the Development of an Equation of "
                             "State for Sulphur Hexafluoride",
                    "ref": "IUPAC Thermodynamic Tables Project Centre, 1991.",
                    "doi": ""},

        "R": 8.31448,
        "cp": CP1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 55000.0, "rhomax": 12.7,

        "nr1": [0.26945570453, -0.554046585076, -0.929624636454, .505661081063,
                -0.683495847809, 0.579161832426, -0.122636218956,
                -0.260339227668e-1, 0.222201648687e-1, -0.118992341472e-2,
                0.292000609763e-2, -0.243315775571e-2, 0.689778297550e-3],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4, 5, 10, 10, 10, 10],
        "t1": [0, 1.5, 2, 0, 1, 2, 0, 2, 0, 0.5, 1, 1.5, 2],

        "nr2": [-0.147585329235e1, 0.275952303526e1, -0.142721418498e1,
                0.598794196648e-1, 0.219991168025e-2, 0.746554473361e-2,
                0.345233637389e-2, -0.253226231963e-1, 0.433906886402e-1,
                -0.249349699078e-1, 0.338560952242e-2, 0.539985899700e-3],
        "d2": [2, 2, 2, 3, 7, 7, 9, 4, 4, 4, 6, 4],
        "t2": [3, 4, 5, 5, 1, 5, 1, 9, 14, 24, 24, 9],
        "c2": [2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 6],
        "gamma2": [1]*12}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for sulfur "
                    "hexafluoride of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": CP1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 12.65,

        "nr1": [0.12279403e1, -0.33035623e1, 0.12094019e1, -0.12316,
                0.11044657, 0.32952153e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.27017629, -0.62910351e-1, -0.3182889, -0.99557419e-1,
                -0.36909694e-1, 0.19136427e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*7}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for sulfur hexafluoride of "
                    "Polt (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von "
                             "Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},

        "R": 8.3143,
        "cp": CP2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 523.0, "Pmax": 40000.0, "rhomax": 13.133,

        "nr1": [0.131111896375, -0.792338803106, 0.580899809209,
                0.153233600406e1, -0.485096079094e1, 0.482411603806e1,
                -0.311285647219e1, 0.442141211276, 0.206313183222,
                -0.372305169645, 0.443536383059, -0.476354850910e-1,
                0.116313319336, 0.570240883234e-1, -0.152963195118,
                0.259842094503e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.131111896375, 0.792338803106, -0.580899809209,
                -0.744763581796, 0.204368923925e1, -0.129335324120e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.32678063]*6}

    eq = guder, reuck, shortSpan, polt

    # EoS in unorthodox format
    # Scalabrin, G, Bettio, L., Marchi, P. Stringari, P.
    # A Fundamental Equation of State for Sulfur Hexafluoride (SF&) in Extended
    # Equation of State Format
    # J. Phys. Chem. Ref. Data 36(2) (2007) 617-662
    # doi: 10.1063/1.2716004

    _surface = {"sigma": [0.0538, -4.064e-5], "exp": [1.271, 0.2116]}
    _dielectric = {
        "__doi__": {
            "autor": "Harvey, A.H., Mountain, R.D.",
            "title": "Correlations for the Dielectric Constants of H2S, SO2 "
                     "and SF6",
            "ref": "Int. J. Thermophys. 38 (2017) 147",
            "doi": "10.1007/s10765-017-2279-6"},
        "eq": 1,
        "a": [16.5, 0], "b": [69.2, 38.8], "c": [-177.9, -171.9],
        "Au": 0, "D": 1.2}

    _melting = {
        "eq": 2,
        "__doi__": {
            "autor": "Harvey, A.H.",
            "title": "On the Melting Curve of Sulfur Hexafluoride",
            "ref": "J. Phys. Chem. Ref. Data 46(4) (2017) 043102",
            "doi": "10.1063/1.5005537"},

        "Tmin": Tt, "Tmax": 625,
        "Tref": Tt, "Pref": 231429,
        "a2": [223.7e6], "exp2": [1.555]}
    _sublimation = {
        "eq": 2,
        "__doi__": guder["__doi__"],

        "Tref": Tt, "Pref": 231429,
        "Tmin": Tt, "Tmax": Tt,
        "a2": [-11.6942141], "exp2": [-1.07]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.09634642, 1.676662, -2.3921599, 5.86078302, -9.02978735],
        "t": [1.0, 1.5, 2.5, 4.0, 4.5]}
    _liquid_Density = {
        "eq": 3,
        "n": [2.31174688, -1.12912486, -1.439347, 0.282489982],
        "t": [0.355, 3/6, 8/6, 10/6]}
    _vapor_Density = {
        "eq": 3,
        "n": [23.68063442, 0.513062232, -24.4706238, -4.6715244, -1.7536843,
              -6.65585369],
        "t": [0.348, 1/6, 2/6, 4/6, 16/6, 34/6]}

    visco0 = {"__name__": "Quiñones-Cisneros (2012)",
              "__doi__": {
                  "autor": "Quiñones-Cisneros, S.E., Huber, M.L., Deiters, "
                           "U.K.",
                  "title": "Correlation for the Viscosity of Sulfur "
                           "Hexafluoride (SF6) from the Triple Point to 1000 K"
                           " and Pressures to 50 MPa",
                  "ref": "J. Phys. Chem. Ref. Data 41(2) (2012) 023102",
                  "doi": "10.1063/1.3702441"},

              "eq": 4, "omega": 0,

              "Toref": 318.7232,
              "no": [118.561, -378.103, 416.428, -165.295, 24.5381],
              "to": [0, 0.25, 0.5, 0.75, 1],

              "a": [-6.87811e-4, 8.22661e-4, -3.54867e-4],
              "b": [1.72737e-4, -2.02448e-4, 1.95952e-4],
              "c": [5.38783e-5, 1.63805e-6, -2.08160e-5],
              "A": [9.99563e-8, -9.64167e-9, -7.54196e-9],
              "B": [-8.98256e-8, -8.49428e-8, 0],
              "C": [-8.53432e-6, 1.14404e-5, -5.65762e-6],
              "D": [0, 0, 2.2798e-11],
              "E": [0, -5.69402e-11, 2.92190e-11]}

    _viscosity = (visco0, )

    thermo0 = {"__name__": "Assael (2012)",
               "__doi__": {
                   "autor": "Assael, M.J., Koini, I.A., Anoniadis, K.D., "
                            "Huber, M.L., Abdulagatov, I.M., Perkins, R.A.",
                   "title": "Reference Correlation of the Thermal "
                            "Conductivity of Sulfur Hexafluoride from the "
                            "Triple Point to 1000 K and up to 150 MPa",
                   "ref": "J. Phys. Chem. Ref. Data 41(2) (2012) 023104",
                   "doi": "10.1063/1.4708620"},

               "eq": 1,

               "Toref": 1, "koref": 1e-3,
               "no_num": [1461860, -18539.4, 77.7891, 0.0241059],
               "to_num": [0, 1, 2, 3],
               "no_den": [29661.7, 505.67, 1],
               "to_den": [0, 1, 2],

               "Tref_res": 318.7232, "rhoref_res": 742.297, "kref_res": 1,
               "nr": [-2.83746e-2, 2.07472e-2, -5.57180e-3, 5.32890e-3,
                      -1.61688e-3, 3.52768e-2, -4.33053e-2, 5.12084e-2,
                      -2.90262e-2, 5.98438e-3],
               "tr": [0, 0, 0, 0, 0, -1, -1, -1, -1, -1],
               "dr": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.19e-9, "gam0": 0.052, "qd": 0.35e-9,

               # Using the value in Erratum paper
               # Assael, M.J., Koini, I.A., Anoniadis, K.D., Huber, M.L.,
               # Abdulagatov, I.M., Perkins, R.A.
               # Reference Correlation of the Thermal Conductivity of Sulfur
               # Hexafluoride from the Triple Point to 1000 K and up to 150 MPa
               # J. Phys. Chem. Ref. Data 43(3) (2014) 039901
               # 10.1063/1.4885454
               "Tcref": 478.08}

    _thermal = (thermo0, )


class Test(TestCase):
    """Testing"""

    def test_guder(self):
        """Table 27, Pag 57"""
        st = SF6(T=350, P=5e6)
        self.assertEqual(round(st.rho, 7), 436.9770888)

        tau = st.Tc/350
        delta = st.rho/st.rhoc
        phi0 = st._phi0(st._constants["cp"], tau, delta)
        self.assertEqual(round(phi0["fio"], 8), 0.330559888e1)
        self.assertEqual(round(phi0["fiod"], 8), 0.169871606e1)
        self.assertEqual(round(phi0["fiodd"], 8), -0.288563626e1)
        self.assertEqual(round(phi0["fiot"], 9), 0.912770720)
        self.assertEqual(round(phi0["fiott"], 7), -0.144662979e2)
        self.assertEqual(round(phi0["fiodt"], 7), 0)

        phir = st._Helmholtz(tau, delta)
        self.assertEqual(round(phir["fir"], 9), -0.496581463)
        self.assertEqual(round(phir["fird"], 9), -0.723171558)
        self.assertEqual(round(phir["firdd"], 9), 0.405086373)
        self.assertEqual(round(phir["firt"], 8), -0.137926327e1)
        self.assertEqual(round(phir["firtt"], 8), -0.137917096e1)
        self.assertEqual(round(phir["firdt"], 8), -0.209574715e1)

        # Selected point from Table 28, Pag 71, Saturation state
        st = SF6(T=223.555, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.231424)
        self.assertEqual(round(st.Liquido.rho, 2), 1845.03)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -158.14)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.72226)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.52753)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 5), 0.83712)
        self.assertEqual(round(st.Liquido.w, 2), 552.26)
        self.assertEqual(round(st.Gas.rho, 3), 19.560)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -47.521)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.22745)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.48331)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.56309)
        self.assertEqual(round(st.Gas.w, 2), 112.84)

        st = SF6(T=250, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 0.621964)
        self.assertEqual(round(st.Liquido.rho, 2), 1705.31)
        self.assertEqual(round(st.Liquido.h.kJkg, 2), -134.81)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.62463)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.57964)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 5), 0.92638)
        self.assertEqual(round(st.Liquido.w, 2), 428.57)
        self.assertEqual(round(st.Gas.rho, 3), 50.867)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -36.589)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.23177)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.55056)
        self.assertEqual(round(st.Gas.cp.kJkgK, 5), 0.66356)
        self.assertEqual(round(st.Gas.w, 2), 111.04)

        st = SF6(T=300, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 2.468172)
        self.assertEqual(round(st.Liquido.rho, 2), 1319.00)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), -81.864)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.43688)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.68633)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 4), 1.3801)
        self.assertEqual(round(st.Liquido.w, 2), 189.80)
        self.assertEqual(round(st.Gas.rho, 2), 238.43)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -22.596)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.23932)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.70290)
        self.assertEqual(round(st.Gas.cp.kJkgK, 4), 1.2722)
        self.assertEqual(round(st.Gas.w, 3), 90.145)

        st = SF6(T=318, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 3.695565)
        self.assertEqual(round(st.Liquido.rho, 2), 917.87)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), -50.870)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.34041)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 5), 0.88030)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 3), 16.370)
        self.assertEqual(round(st.Liquido.w, 3), 70.232)
        self.assertEqual(round(st.Gas.rho, 2), 569.56)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -33.648)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.28626)
        self.assertEqual(round(st.Gas.cv.kJkgK, 5), 0.93693)
        self.assertEqual(round(st.Gas.cp.kJkgK, 3), 26.348)
        self.assertEqual(round(st.Gas.w, 3), 69.747)

        st = SF6(T=318.7, x=0.5)
        self.assertEqual(round(st.P.MPa, 6), 3.753053)
        self.assertEqual(round(st.Liquido.rho, 2), 787.75)
        self.assertEqual(round(st.Liquido.h.kJkg, 3), -44.789)
        self.assertEqual(round(st.Liquido.s.kJkgK, 5), -0.32152)
        self.assertEqual(round(st.Liquido.cv.kJkgK, 4), 1.0458)
        self.assertEqual(round(st.Liquido.cp.kJkgK, 2), 1027.68)
        self.assertEqual(round(st.Liquido.w, 3), 60.578)
        self.assertEqual(round(st.Gas.rho, 2), 696.67)
        self.assertEqual(round(st.Gas.h.kJkg, 3), -40.391)
        self.assertEqual(round(st.Gas.s.kJkgK, 5), -0.30773)
        self.assertEqual(round(st.Gas.cv.kJkgK, 4), 1.0721)
        self.assertEqual(round(st.Gas.cp.kJkgK, 2), 1273.31)
        self.assertEqual(round(st.Gas.w, 3), 62.816)

        # Selected points from Table 29, Pag 73, single phase region
        st = SF6(T=235, P=1e5)
        self.assertEqual(round(st.rho, 4), 7.6624)
        self.assertEqual(round(st.u.kJkg, 3), -52.541)
        self.assertEqual(round(st.h.kJkg, 3), -39.490)
        self.assertEqual(round(st.s.kJkgK, 5), -0.14668)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.49557)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.55917)
        self.assertEqual(round(st.w, 2), 119.82)

        st = SF6(T=245, P=5e5)
        self.assertEqual(round(st.rho, 3), 40.652)
        self.assertEqual(round(st.u.kJkg, 3), -50.563)
        self.assertEqual(round(st.h.kJkg, 3), -38.264)
        self.assertEqual(round(st.s.kJkgK, 5), -0.22772)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.53594)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.63650)
        self.assertEqual(round(st.w, 2), 112.47)

        st = SF6(T=625, P=1e6)
        self.assertEqual(round(st.rho, 3), 28.218)
        self.assertEqual(round(st.u.kJkg, 2), 236.51)
        self.assertEqual(round(st.h.kJkg, 2), 271.95)
        self.assertEqual(round(st.s.kJkgK, 5), 0.47409)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.88631)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.94703)
        self.assertEqual(round(st.w, 2), 194.24)

        st = SF6(T=290, P=2e6)
        self.assertEqual(round(st.rho, 2), 1424.28)
        self.assertEqual(round(st.u.kJkg, 3), -95.379)
        self.assertEqual(round(st.h.kJkg, 3), -93.975)
        self.assertEqual(round(st.s.kJkgK, 5), -0.47677)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.65860)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.1854)
        self.assertEqual(round(st.w, 2), 243.05)

        st = SF6(T=230, P=5e6)
        self.assertEqual(round(st.rho, 2), 1838.34)
        self.assertEqual(round(st.u.kJkg, 2), -154.40)
        self.assertEqual(round(st.h.kJkg, 2), -151.68)
        self.assertEqual(round(st.s.kJkgK, 5), -0.70521)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.54153)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.84036)
        self.assertEqual(round(st.w, 2), 554.97)

        st = SF6(T=380, P=1e7)
        self.assertEqual(round(st.rho, 2), 825.25)
        self.assertEqual(round(st.u.kJkg, 4), -3.5308)
        self.assertEqual(round(st.h.kJkg, 4), 8.5868)
        self.assertEqual(round(st.s.kJkgK, 5), -0.19035)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.76919)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.2719)
        self.assertEqual(round(st.w, 2), 142.33)

        st = SF6(T=400, P=2e7)
        self.assertEqual(round(st.rho, 2), 1123.76)
        self.assertEqual(round(st.u.kJkg, 4), -1.5810)
        self.assertEqual(round(st.h.kJkg, 3), 16.216)
        self.assertEqual(round(st.s.kJkgK, 5), -0.19659)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.77581)
        self.assertEqual(round(st.cp.kJkgK, 4), 1.0792)
        self.assertEqual(round(st.w, 2), 239.83)

        st = SF6(T=255, P=5e7)
        self.assertEqual(round(st.rho, 2), 1912.92)
        self.assertEqual(round(st.u.kJkg, 2), -144.76)
        self.assertEqual(round(st.h.kJkg, 2), -118.62)
        self.assertEqual(round(st.s.kJkgK, 5), -0.66770)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.59557)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.82513)
        self.assertEqual(round(st.w, 2), 688.44)

        st = SF6(T=300, P=1e8)
        self.assertEqual(round(st.rho, 2), 1920.29)
        self.assertEqual(round(st.u.kJkg, 2), -116.67)
        self.assertEqual(round(st.h.kJkg, 3), -64.594)
        self.assertEqual(round(st.s.kJkgK, 5), -0.56696)
        self.assertEqual(round(st.cv.kJkgK, 5), 0.66951)
        self.assertEqual(round(st.cp.kJkgK, 5), 0.86708)
        self.assertEqual(round(st.w, 2), 755.29)

    def test_shortSpan(self):
        """Table III, Pag 46"""
        st = SF6(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 0.9671)
        self.assertEqual(round(st.P.MPa, 3), 8.094)
        self.assertEqual(round(st.cp.kJkgK, 4), 0.9958)

        st2 = SF6(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 52.80)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.10913)

    def test_Quinones(self):
        """Table 8, pag 10"""
        # The dilute-gas limit differ in 6th figure
        # 0.0006 μPa·s at 300 and 400K
        self.assertEqual(round(SF6(T=300, rho=0).mu.muPas, 4), 15.2881)
        self.assertEqual(round(SF6(T=300, rho=5.92).mu.muPas, 4), 15.3037)
        self.assertEqual(round(SF6(T=300, rho=1345.1).mu.muPas, 3), 117.419)
        self.assertEqual(round(SF6(T=400, rho=0).mu.muPas, 4), 19.6790)
        self.assertEqual(round(SF6(T=400, rho=278.47).mu.muPas, 4), 24.4266)
        self.assertEqual(round(SF6(T=400, rho=1123.8).mu.muPas, 4), 84.7838)

    def test_Assael(self):
        """Table 5, from Erratum article"""
        self.assertEqual(round(SF6(T=298.15, rho=0).k.mWmK, 2), 12.95)
        self.assertEqual(round(SF6(T=298.15, rho=100).k.mWmK, 2), 14.13)
        self.assertEqual(round(SF6(T=298.15, rho=1600).k.mWmK, 2), 69.73)
        self.assertEqual(round(SF6(T=310, rho=0).k.mWmK, 2), 13.83)
        self.assertEqual(round(SF6(T=310, rho=1200).k.mWmK, 2), 48.70)
        self.assertEqual(round(SF6(T=480, rho=100).k.mWmK, 2), 28.85)

    def test_dielectric(self):
        self.assertEqual(round(SF6(T=300, rhom=10).epsilon, 5), 1.59326)
