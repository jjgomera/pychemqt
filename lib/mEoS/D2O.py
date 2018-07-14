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

from iapws._iapws import _D2O_Viscosity, _D2O_ThCond, _D2O_Tension

from lib.meos import MEoS
from lib import unidades


class D2O(MEoS):
    """Multiparameter equation of state for heavy water"""
    name = "heavy water"
    CASNumber = "7789-20-0"
    formula = "D2O"
    synonym = "deuterium oxide"
    _refPropName = "D2O"
    _coolPropName = "HeavyWater"
    Tc = unidades.Temperature(643.847)
    rhoc = unidades.Density(358)
    Pc = unidades.Pressure(21671.0, "kPa")
    M = 20.027508  # g/mol
    Tt = unidades.Temperature(276.97)
    Tb = unidades.Temperature(374.563)
    f_acent = 0.364
    momentoDipolar = unidades.DipoleMoment(1.9, "Debye")

    Fi0 = {"ao_log": [1, 3],
           "pow": [0, 1, 2, 3, 4, 5],
           "ao_pow": [-8.6739710041, 6.9611755531],
           "ao_exp": [0.00863, 0.97454, 2.0646, 0.23528, 0.29555],
           "titao": [0.4255669437, 2.6093155672, 6.0185106089, 11.3380974051,
                     29.5101165339],
           "ao_hyp": [], "hyp": []}

    Fi1 = {"ao_log": [0.5399322597e-2, 0],
           "pow": [0, 1, 2, 3, 4, 5],
           "ao_pow": [0.3087155964e2, -.3827264031e2, 0.4424799189,
                      -.1256336874e1, 0.2843343470, -.2401555088e-1],
           "tau*logtau": -.1288399716e2,
           "tau*logdelta": 0.4415884023e1,
           "ao_exp": [], "titao": [],
           "ao_hyp": [], "hyp": []}

    CP1 = {"ao": 0.39176485e1,
           "an": [-0.31123915e-3, 0.41173363e-5, -0.28943955e-8,
                  0.63278791e-12, 0.78728740],
           "pow": [1.00, 2.00, 3.00, 4.00, -0.99],
           "ao_exp": [],
           "exp": [],
           "ao_hyp": [], "hyp": []}

    herrig = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heavy water of Herrig"
        " (2017).",
        "__doi__": {
            "autor": "Herrig",
            "title": "Preliminary helmholtz equation of state for Heavy Water",
            "ref": "",
            "doi": ""},

        "R": 8.3144621, "rhoref": 17.77555*M, "Tref": 643.847,
        "cp": Fi0,
        "ref": {"Tref": 276.95, "Pref": 0.660096, "ho": 0.598, "so": 0},

        "Tmin": Tt, "Tmax": 800.0, "Pmax": 100000.0, "rhomax": 65.,
        "Pmin": 0.66103, "rhomin": 55.198,

        "nr1": [0.0105835, 0.99127253, -1.224122, 1.710643, -2.189443,
                0.1145315],
        "d1": [4.0, 1.0, 1.0, 2.0, 2.0, 3.0],
        "t1": [1.0, 0.463, 1.29, 1.307, 1.2165, 0.587],

        "nr2": [-0.89875532, -1.597051, -2.804509, 0.33016885,
                -3.396526, -0.001881],
        "c2": [1.0, 2.0, 2.0, 1.0, 2.0, 1.0],
        "d2": [1.0, 1.0, 3.0, 2.0, 2.0, 8.0],
        "t2": [2.95, 1.713, 1.929, 0.94, 3.033, 0.765],
        "gamma2": [1]*6,

        "nr3": [-0.70355957, -0.20345481, -0.70691398, 2.094255, 3.042546,
                0.8010728, 0.213384, 0.32335789, -0.0245055, 0.7380677,
                -0.21484089],
        "t3": [1.504, 2.85, 1.96, 0.969, 2.576, 2.79, 3.581, 3.67, 1.7, 1.0,
               4.1],
        "d3": [1.0, 2.0, 3.0, 1.0, 3.0, 1.0, 1.0, 2.0, 2.0, 2.0, 1.0],
        "beta3": [0.907, 0.48, 1.223, 2.61, 4.283, 1.4, 0.735, 0.24, 1067.0,
                  13.27, 1.48],
        "alfa3": [0.982, 1.34, 1.658, 1.6235, 1.4, 2.206, 0.84, 1.535, 11.33,
                  3.86, 7.56],
        "epsilon3": [2.272, 1.375, 0.648, 0.8925, 0.145, 0.291, 2.01, 1.08,
                     0.96, 0.181, 0.529],
        "gamma3": [2.263, 2.343, 0.929, 1.0, 1.383, 0.968, 1.695, 2.23, 1.07,
                   1.297, 2.41]}

    hill = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for heavy water of Hill "
        "et al. (1982).",
        "__doi__": {"autor": "Hill, P.G., MacMillan, R.D.C., and Lee, V.",
                    "title": "A Fundamental Equation of State for Heavy Water",
                    "ref": "J. Phys. Chem. Ref. Data 11, 1 (1982)",
                    "doi": "10.1063/1.555661"},

        "R": 8.3143565, "rhoref": 17.875414*M,
        "cp": Fi1,
        "ref": {"Tref": 276.95, "Pref": 0.660096, "ho": 0.598, "so": 0},

        "Tmin": Tt, "Tmax": 800.0, "Pmax": 100000.0, "rhomax": 65.,
        "Pmin": 0.66103, "rhomin": 55.198,

        "nr1": [-0.384820628204e3, 0.108213047259e4, -0.110768260635e4,
                0.164668954246e4, -0.137959852228e4, 0.598964185629e3,
                -0.100451752702e3, 0.419192736351e3, -0.107279987867e4,
                0.653852283544e3, -0.984305985655e3, 0.845444459339e3,
                -0.376799930490e3, 0.644512590492e2, -0.214911115714e3,
                0.531113962967e3, -0.135454224420e3, 0.202814416558e3,
                -0.178293865031e3, 0.818739394970e2, -0.143312594493e2,
                0.651202383207e2, -0.171227351208e3, 0.100859921516e2,
                -0.144684680657e2, 0.128871134847e2, -0.610605957134e1,
                0.109663804408e1, -0.115734899702e2, 0.374970075409e2,
                0.897967147669, -0.527005883203e1, 0.438084681795e-1,
                0.406772082680, -0.965258571044e-2, -0.119044600379e-1],
        "d1": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
               4, 4, 4, 4, 4, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8],
        "t1": [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6,
               0, 1, 2, 3, 4, 5, 6, 0, 1, 0, 1, 0, 1, 0, 1],

        "nr2": [0.382589102341e3, -0.106406466204e4, 0.105544952919e4,
                -0.157579942855e4, 0.132703387531e4, -0.579348879870e3,
                0.974163902526e2, 0.286799294226e3, -0.127543020847e4,
                0.275802674911e4, -0.381284331492e4, 0.293755152012e4,
                -0.117858249946e4, 0.186261198012e3],
        "c2": [1]*14,
        "d2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2],
        "t2": [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6],
        "gamma2": [1.5394]*14}

    eq = herrig, hill

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.80236e1, 0.23957e1, -0.42639e2, 0.99569e2, -0.62135e2],
        "t": [1.0, 1.5, 2.75, 3.0, 3.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.26406e1, 0.97090e1, -0.18058e2, 0.87202e1, -0.74487e1],
        "t": [0.3678, 1.9, 2.2, 2.63, 7.3]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.37651e1, -0.38673e2, 0.73024e2, -0.13251e3, 0.75235e2,
              -0.70412e2],
        "t": [0.409, 1.766, 2.24, 3.04, 3.42, 6.9]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "IAPWS (2007)",
              "__doi__": {
                  "autor": "Kestin, J., Sengers, J.V., Kamgar-Parsi, B., "
                           "Levelt Sengers, J.M.H.",
                  "title": "Thermophysical Properties of Fluid D2O",
                  "ref": "J. Phys. Chem. Ref. Data 13(2) (1984) 601-609",
                  "doi": "10.1063/1.555714"}}

    def _visco0(self, rho, T, fase):
        mu = _D2O_Viscosity(rho, T)
        return unidades.Viscosity(mu)

    _viscosity = visco0,

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "IAPWS (1994)",
               "__doi__": {
                  "autor": "Kestin, J., Sengers, J.V., Kamgar-Parsi, B., "
                           "Levelt Sengers, J.M.H.",
                  "title": "Thermophysical Properties of Fluid D2O",
                  "ref": "J. Phys. Chem. Ref. Data 13(2) (1984) 601-609",
                  "doi": "10.1063/1.555714"}}

    def _thermo0(self, rho, T, fase):
        k = _D2O_ThCond(rho, T)
        return unidades.ThermalConductivity(k)

    _thermal = thermo0,

    def _Surface(self, T):
        """Equation for the surface tension"""
        try:
            s = _D2O_Tension(T)
        except NotImplementedError:
            s = None
        return unidades.Tension(s)


class Test(TestCase):

    def test_D2O(self):
        # Pag 17 of IAPWS 2007 update paper
        """Table 5 pag 11"""
        fluid = D2O(eq="hill")
        Tr = 643.847
        rhor = 358
        ar = 21.671*1000/358
        sr = 21.671*1000/358./643.89
        pr = 21.671*1e6

        rho = 0.0002*rhor
        T = 0.5*Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -2.644979)
        self.assertEqual(round(state["P"]/pr, 7), 0.0004402)
        # self.assertEqual(round(state["cv"]/sr, 4), 14.2768)

        rho = 3.18*rhor
        T = 0.5*Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -0.217388)
        # self.assertEqual(round(state["P"]/pr, 7), 4.3549719)
        # self.assertEqual(round(state["cv"]/sr, 4), 41.4463)

        rho = 0.0295*rhor
        T = 0.75*Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -7.272543)
        # self.assertEqual(round(state["P"]/pr, 7), 0.0870308)
        # self.assertEqual(round(state["cv"]/sr, 4), 20.1586)

        rho = 2.83*rhor
        T = 0.75*Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -4.292707)
        # self.assertEqual(round(state["P"]/pr, 7), 4.4752958)
        # self.assertEqual(round(state["cv"]/sr, 4), 33.4367)

        rho = 0.3*rhor
        T = Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -15.163326)
        # self.assertEqual(round(state["P"]/pr, 7), 0.8014044)
        # self.assertEqual(round(state["cv"]/sr, 4), 30.8587)

        rho = 1.55*rhor
        T = Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -12.643811)
        # self.assertEqual(round(state["P"]/pr, 7), 1.0976283)
        # self.assertEqual(round(state["cv"]/sr, 4), 33.0103)

        rho = 0.4*rhor
        T = 1.2*Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -25.471535)
        # self.assertEqual(round(state["P"]/pr, 7), 1.4990994)
        # self.assertEqual(round(state["cv"]/sr, 4), 23.6594)

        rho = 1.61*rhor
        T = 1.2*Tr
        state = fluid._Helmholtz(rho, T)
        # self.assertEqual(round((state["h"]-state["P"]*1000/rho-T*state["s"])/ar, 6), -21.278164)
        # self.assertEqual(round(state["P"]/pr, 7), 4.5643798)
        # self.assertEqual(round(state["cv"]/sr, 4), 25.4800)

    def test_D2O_Viscosity(self):
        # Selected values from TAble A4, saturation state in IAPWS R4-84
        st = D2O(T=283.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 6), 0.001028)
        self.assertEqual(round(st.Liquido.mu.muPas, 0), 1678)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 9.76)

        st = D2O(T=373.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 5), 0.09634)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 328.7)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 12.62)

        st = D2O(T=473.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 3), 1.547)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 151.9)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 15.99)

        st = D2O(T=573.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 3), 8.694)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 93.6)
        self.assertEqual(round(st.Gas.mu.muPas, 2), 19.68)

        st = D2O(T=623.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 2), 16.82)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 69.0)
        self.assertEqual(round(st.Gas.mu.muPas, 1), 23.7)

        st = D2O(T=643.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 2), 21.47)
        self.assertEqual(round(st.Liquido.mu.muPas, 1), 46.0)
        self.assertEqual(round(st.Gas.mu.muPas, 1), 32.7)

        # Selected values from Table A5, pag 10
        mur = 55.2651e-6
        Tr = 643.847
        rhor = 358
        self.assertEqual(round(
            _D2O_Viscosity(3.09*rhor, 0.431*Tr)/mur, 10), 36.9123166244)
        self.assertEqual(round(
            _D2O_Viscosity(3.18*rhor, 0.5*Tr)/mur, 10), 12.4679405772)
        self.assertEqual(round(
            _D2O_Viscosity(0.0295*rhor, 0.75*Tr)/mur, 10), 0.2951479769)
        self.assertEqual(round(
            _D2O_Viscosity(0.163*rhor, 0.9*Tr)/mur, 10), 0.3619649145)
        self.assertEqual(round(
            _D2O_Viscosity(0.7*rhor, Tr)/mur, 10), 0.5528693914)
        self.assertEqual(round(
            _D2O_Viscosity(0.98*rhor, 1.1*Tr)/mur, 10), 0.7816387903)
        self.assertEqual(round(
            _D2O_Viscosity(0.8*rhor, 1.2*Tr)/mur, 10), 0.7651099154)
        self.assertEqual(round(
            _D2O_Viscosity(1.61*rhor, 1.2*Tr)/mur, 10), 1.2711900131)

    def test_D2O_ThCond(self):
        # Selected values from TAble B3, saturation state in IAPWS R4-84
        st = D2O(T=283.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 6), 0.001028)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 574.6)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 16.99)

        st = D2O(T=373.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 5), 0.09634)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 636.0)
        self.assertEqual(round(st.Gas.k.mWmK, 2), 24.82)

        st = D2O(T=473.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 3), 1.547)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 592.0)
        self.assertEqual(round(st.Gas.k.mWmK, 1), 39.00)

        st = D2O(T=573.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 3), 8.694)
        self.assertEqual(round(st.Liquido.k.mWmK, 1), 473.2)
        self.assertEqual(round(st.Gas.k.mWmK, 1), 75.2)

        st = D2O(T=623.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 2), 16.82)
        self.assertEqual(round(st.Liquido.k.mWmK, 0), 391)
        self.assertEqual(round(st.Gas.k.mWmK, 0), 143)

        st = D2O(T=643.15, x=0.5, eq="hill")
        self.assertEqual(round(st.P.MPa, 2), 21.47)
        self.assertEqual(round(st.Liquido.k.mWmK, 0), 548)
        self.assertEqual(round(st.Gas.k.mWmK, 0), 467)

        # Selected values from Table B4, pag 17
        lr = 0.742128e-3
        Tr = 643.847
        rhor = 358
        self.assertEqual(
                round(_D2O_ThCond(3.09*rhor, 0.431*Tr)/lr, 9), 762.915707396)
        self.assertEqual(
                round(_D2O_ThCond(2.95*rhor, 0.6*Tr)/lr, 9), 861.240794445)
        self.assertEqual(
                round(_D2O_ThCond(1.55*rhor, Tr)/lr, 9), 502.846952426)
        self.assertEqual(
                round(_D2O_ThCond(1.61*rhor, 1.2*Tr)/lr, 9), 471.747729424)
        self.assertEqual(
                round(_D2O_ThCond(1.37*rhor, 1.27*Tr)/lr, 9), 409.359675394)

    def test_D2O_Tension(self):
        # Selected values from table 1 in IAPWS R5-85"""
        st = D2O(T=283.15, x=0)
        self.assertEqual(round(st.sigma.mNm, 2), 74.06)

        st = D2O(T=373.15, x=0)
        self.assertEqual(round(st.sigma.mNm, 2), 58.93)

        st = D2O(T=643.15, x=0)
        self.assertEqual(round(st.sigma.mNm, 2), 0.05)
