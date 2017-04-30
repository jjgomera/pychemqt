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


from iapws._iapws import _D2O_Viscosity, _D2O_ThCond

from lib.meos import MEoS
from lib import unidades


class D2O(MEoS):
    """Multiparameter equation of state for heavy water"""
    name = "heavy water"
    CASNumber = "7789-20-0"
    formula = "D2O"
    synonym = "deuterium oxide"
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

    helmholtz1 = {
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

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for heavy water of Hill "
        "et al. (1982).",
        "__doi__": {"autor": "Hill, P.G., MacMillan, R.D.C., and Lee, V.",
                    "title": "A Fundamental Equation of State for Heavy Water",
                    "ref": "J. Phys. Chem. Ref. Data 11, 1 (1982)",
                    "doi": "10.1063/1.555661"},
        "__test__":
            # Pag 17 of IAPWS 2007 update paper
            """
            >>> st=D2O(T=0.5*D2O.Tc, rho=0.0002*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -2.644979 0.0004402 14.2768
            >>> st=D2O(T=0.5*D2O.Tc, rho=3.18*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -0.217388 4.3549719 41.4463
            >>> st=D2O(T=0.75*D2O.Tc, rho=0.0295*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -7.272543 0.0870308 20.1586
            >>> st=D2O(T=0.75*D2O.Tc, rho=2.83*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -4.292707 4.4752958 33.4367
            >>> st=D2O(T=D2O.Tc, rho=0.3*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -15.163326 0.8014044 30.8587
            >>> st=D2O(T=D2O.Tc, rho=1.55*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -12.643811 1.0976283 33.0103
            >>> st=D2O(T=1.2*D2O.Tc, rho=0.4*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -25.471535 1.4990994 23.6594
            >>> st=D2O(T=1.2*D2O.Tc, rho=1.61*D2O.rhoc)
            >>> print("%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa))
            -21.278164 4.5643798 25.4800
            """,

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

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.238, -0.152082], "exp": [1.25, 2.25]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.80236e1, 0.23957e1, -0.42639e2, 0.99569e2, -0.62135e2],
        "exp": [1.0, 1.5, 2.75, 3.0, 3.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.26406e1, 0.97090e1, -0.18058e2, 0.87202e1, -0.74487e1],
        "exp": [0.3678, 1.9, 2.2, 2.63, 7.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.37651e1, -0.38673e2, 0.73024e2, -0.13251e3, 0.75235e2,
               -0.70412e2],
        "exp": [0.409, 1.766, 2.24, 3.04, 3.42, 6.9]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "IAPWS (2007)",
              "__doi__": {"autor": "J. Kestin, J. V. Sengers, B. Kamgar‐Parsi"
                          "and J. M. H. Levelt Sengers",
                          "title": "Thermophysical Properties of Fluid D2O",
                          "ref": "J. Phys. Chem. Ref. Data 13, 601 (1984)",
                          "doi": "10.1063/1.555714"},
              "__test__":
                  # Pag 17 of IAPWS 2007 update paper
                  """
                  >>> st=D2O(T=0.431*D2O.Tc, rho=3.09*D2O.rhoc)
                  >>> print("%0.10f" % (st.mu.muPas/55.2651))
                  36.9123166244
                  >>> st=D2O(T=0.431*D2O.Tc, rho=3.23*D2O.rhoc)
                  >>> print("%0.10f" % (st.mu.muPas/55.2651))
                  34.1531546602
                  >>> st=D2O(T=0.6*D2O.Tc, rho=2.95*D2O.rhoc)
                  >>> print("%0.10f" % (st.mu.muPas/55.2651))
                  5.2437249935
                  >>> st=D2O(T=D2O.Tc, rho=0.7*D2O.rhoc)
                  >>> print("%0.10f" % (st.mu.muPas/55.2651))
                  0.5528693914
                  >>> st=D2O(T=0.9*D2O.Tc, rho=0.08*D2O.rhoc)
                  >>> print("%0.10f" % (st.mu.muPas/55.2651))
                  0.3685472578
                  >>> st=D2O(T=1.1*D2O.Tc, rho=0.98*D2O.rhoc)
                  >>> print("%0.10f" % (st.mu.muPas/55.2651))
                  0.7816387903
                  >>> st=D2O(T=1.2*D2O.Tc, rho=0.8*D2O.rhoc)
                  >>> print("%0.10f" % (st.mu.muPas/55.2651))
                  0.7651099154
                  """, }

    def _visco0(self, rho, T, fase):
        mu = _D2O_Viscosity(rho, T)
        return unidades.Viscosity(mu)

    _viscosity = visco0,

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "IAPWS (1994)",
               "__doi__": {
                   "autor": "J. Kestin, J. V. Sengers, B. Kamgar‐Parsi and J."
                   " M. H. Levelt Sengers",
                   "title": "Thermophysical Properties of Fluid D2O",
                   "ref": "J. Phys. Chem. Ref. Data 13, 601 (1984)",
                   "doi": "10.1063/1.555714"},
               "__test__":
                   # Pag 17 of IAPWS 2007 update paper
                   """
                   >>> st=D2O(T=0.431*D2O.Tc, rho=3.09*D2O.rhoc)
                   >>> print("%0.9f" % (st.k.mWmK/0.742128))
                   762.915707396
                   >>> st=D2O(T=0.431*D2O.Tc, rho=3.23*D2O.rhoc)
                   >>> print("%0.9f" % (st.k.mWmK/0.742128))
                   833.912049618
                   >>> st=D2O(T=0.6*D2O.Tc, rho=2.95*D2O.rhoc)
                   >>> print("%0.9f" % (st.k.mWmK/0.742128))
                   861.240794445
                   >>> st=D2O(T=D2O.Tc, rho=0.7*D2O.rhoc)
                   >>> print("%0.9f" % (st.k.mWmK/0.742128))
                   469.015122112
                   >>> st=D2O(T=0.9*D2O.Tc, rho=0.08*D2O.rhoc)
                   >>> print("%0.9f" % (st.k.mWmK/0.742128))
                   74.522283066
                   >>> st=D2O(T=1.1*D2O.Tc, rho=0.98*D2O.rhoc)
                   >>> print("%0.9f" % (st.k.mWmK/0.742128))
                   326.652382218
                   >>> st=D2O(T=1.2*D2O.Tc, rho=0.8*D2O.rhoc)
                   >>> print("%0.9f" % (st.k.mWmK/0.742128))
                   259.605241187
                   """, }

    def _thermo0(self, rho, T, fase):
        k = _D2O_ThCond(rho, T)
        return unidades.ThermalConductivity(k)

    _thermal = thermo0,
