#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

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


from scipy import exp

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
        "__name__": u"Helmholtz equation of state for heavy water of Hill et al. (1982).",
        "__doi__": {"autor": "Hill, P.G., MacMillan, R.D.C., and Lee, V.",
                    "title": "A Fundamental Equation of State for Heavy Water",
                    "ref": "J. Phys. Chem. Ref. Data 11, 1 (1982)",
                    "doi": "10.1063/1.555661"},
        "__test__":
            # Pag 17 of IAPWS 2007 update paper
            """
            >>> st=D2O(T=0.5*D2O.Tc, rho=0.0002*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
            -2.644979 0.0004402 14.2768
            >>> st=D2O(T=0.5*D2O.Tc, rho=3.18*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
            -0.217388 4.3549719 41.4463
            >>> st=D2O(T=0.75*D2O.Tc, rho=0.0295*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
            -7.272543 0.0870308 20.1586
            >>> st=D2O(T=0.75*D2O.Tc, rho=2.83*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
            -4.292707 4.4752958 33.4367
            >>> st=D2O(T=D2O.Tc, rho=0.3*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
            -15.163326 0.8014044 30.8587
            >>> st=D2O(T=D2O.Tc, rho=1.55*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
            -12.643811 1.0976283 33.0103
            >>> st=D2O(T=1.2*D2O.Tc, rho=0.4*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
            -25.471535 1.4990994 23.6594
            >>> st=D2O(T=1.2*D2O.Tc, rho=1.61*D2O.rhoc)
            >>> print "%0.6f %0.7f %0.4f" % (st.a*D2O.rhoc/D2O.Pc, st.Pr, st.cv.kJkgK*D2O.rhoc*D2O.Tc/D2O.Pc.kPa)
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

    eq = helmholtz1,

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
        "ao": [-0.37651e1, -0.38673e2, 0.73024e2, -0.13251e3, 0.75235e2, -0.70412e2],
        "exp": [0.409, 1.766, 2.24, 3.04, 3.42, 6.9]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "IAPWS (2007)",
              "__doi__": {"autor": "J. Kestin, J. V. Sengers, B. Kamgar‐Parsi and J. M. H. Levelt Sengers",
                          "title": "Thermophysical Properties of Fluid D2O",
                          "ref": "J. Phys. Chem. Ref. Data 13, 601 (1984)",
                          "doi": "10.1063/1.555714"},
              "__test__":
                  # Pag 17 of IAPWS 2007 update paper
                  """
                  >>> st=D2O(T=0.431*D2O.Tc, rho=3.09*D2O.rhoc)
                  >>> print "%0.10f" % (st.mu.muPas/55.2651)
                  36.9123166244
                  >>> st=D2O(T=0.431*D2O.Tc, rho=3.23*D2O.rhoc)
                  >>> print "%0.10f" % (st.mu.muPas/55.2651)
                  34.1531546602
                  >>> st=D2O(T=0.6*D2O.Tc, rho=2.95*D2O.rhoc)
                  >>> print "%0.10f" % (st.mu.muPas/55.2651)
                  5.2437249935
                  >>> st=D2O(T=D2O.Tc, rho=0.7*D2O.rhoc)
                  >>> print "%0.10f" % (st.mu.muPas/55.2651)
                  0.5528693914
                  >>> st=D2O(T=0.9*D2O.Tc, rho=0.08*D2O.rhoc)
                  >>> print "%0.10f" % (st.mu.muPas/55.2651)
                  0.3685472578
                  >>> st=D2O(T=1.1*D2O.Tc, rho=0.98*D2O.rhoc)
                  >>> print "%0.10f" % (st.mu.muPas/55.2651)
                  0.7816387903
                  >>> st=D2O(T=1.2*D2O.Tc, rho=0.8*D2O.rhoc)
                  >>> print "%0.10f" % (st.mu.muPas/55.2651)
                  0.7651099154
                  """, }

    def _visco0(self, rho, T, fase=None):
        Tr = T/643.847
        rhor = rho/358.0

        no = [1.0, 0.940695, 0.578377, -0.202044]
        fi0 = Tr**0.5/sum([n/Tr**i for i, n in enumerate(no)])

        Li = [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 5, 0, 1, 2, 3, 0, 1, 3,
              5, 0, 1, 5, 3]
        Lj = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,
              4, 5, 5, 5, 6]
        Lij = [0.4864192, -0.2448372, -0.8702035, 0.8716056, -1.051126,
               0.3458395, 0.3509007, 1.315436, 1.297752, 1.353448, -0.2847572,
               -1.037026, -1.287846, -0.02148229, 0.07013759, 0.4660127,
               0.2292075, -0.4857462, 0.01641220, -0.02884911, 0.1607171,
               -0.009603846, -0.01163815, -0.008239587, 0.004559914, -0.003886659]

        array = [lij*(1./Tr-1)**i*(rhor-1)**j for i, j, lij in zip(Li, Lj, Lij)]
        fi1 = exp(rhor*sum(array))

        return unidades.Viscosity(55.2651*fi0*fi1, "muPas")

    _viscosity = visco0,

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "IAPWS (1994)",
               "__doi__": {"autor": "J. Kestin, J. V. Sengers, B. Kamgar‐Parsi and J. M. H. Levelt Sengers",
                           "title": "Thermophysical Properties of Fluid D2O",
                           "ref": "J. Phys. Chem. Ref. Data 13, 601 (1984)",
                           "doi": "10.1063/1.555714"},
               "__test__":
                   # Pag 17 of IAPWS 2007 update paper
                   """
                   >>> st=D2O(T=0.431*D2O.Tc, rho=3.09*D2O.rhoc)
                   >>> print "%0.9f" % (st.k/0.742128)
                   762.915707396
                   >>> st=D2O(T=0.431*D2O.Tc, rho=3.23*D2O.rhoc)
                   >>> print "%0.9f" % (st.k/0.742128)
                   833.912049618
                   >>> st=D2O(T=0.6*D2O.Tc, rho=2.95*D2O.rhoc)
                   >>> print "%0.9f" % (st.k/0.742128)
                   861.240794445
                   >>> st=D2O(T=D2O.Tc, rho=0.7*D2O.rhoc)
                   >>> print "%0.9f" % (st.k/0.742128)
                   469.015122112
                   >>> st=D2O(T=0.9*D2O.Tc, rho=0.08*D2O.rhoc)
                   >>> print "%0.9f" % (st.k/0.742128)
                   74.522283066
                   >>> st=D2O(T=1.1*D2O.Tc, rho=0.98*D2O.rhoc)
                   >>> print "%0.9f" % (st.k/0.742128)
                   326.652382218
                   >>> st=D2O(T=1.2*D2O.Tc, rho=0.8*D2O.rhoc)
                   >>> print "%0.9f" % (st.k/0.742128)
                   259.605241187
                   """, }

    def _thermo0(self, rho, T, fase=None):
        rhor = rho/358
        Tr = T/643.847
        tau = Tr/(abs(Tr-1.1)+1.1)

        no = [1.0, 37.3223, 22.5485, 13.0465, 0.0, -2.60735]
        Lo = sum([Li*Tr**i for i, Li in enumerate(no)])

        nr = [483.656, -191.039, 73.0358, -7.57467]
        Lr = -167.31*(1-exp(-2.506*rhor))+sum([Li*rhor**(i+1) for i, Li in enumerate(nr)])

        f1 = exp(0.144847*Tr-5.64493*Tr**2)
        f2 = exp(-2.8*(rhor-1)**2)-0.080738543*exp(-17.943*(rhor-0.125698)**2)
        f3 = 1+exp(60*(tau-1)+20)
        f4 = 1+exp(100*(tau-1)+15)
        Lc = 35429.6*f1*f2*(1+f2**2*(5e9*f1**4/f3+3.5*f2/f4))

        Ll = -741.112*f1**1.2*(1-exp(-(rhor/2.5)**10))

        return unidades.ThermalConductivity(0.742128e-3*(Lo+Lr+Lc+Ll))

    _thermal = thermo0,
