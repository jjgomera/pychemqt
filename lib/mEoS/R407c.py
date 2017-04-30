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


from lib.meos import MEoSBlend
from lib import unidades


class R407c(MEoSBlend):
    """Multiparamenter equation of state for R407C
    (23% R32, 25% R125, 52% R134a)"""

    name = "R407C"
    CASNumber = ""
    formula = "R32+R125+R134a"
    synonym = "R407C"
    rhoc = unidades.Density(453.430936)
    Tc = unidades.Temperature(359.345)
    Pc = unidades.Pressure(4631.7, "kPa")
    M = 86.2036  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(229.52)
    f_acent = 0.363
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 62

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.4],
           "ao_pow": [2.13194, 8.05008, -14.3914],
           "ao_exp": [1.4245, 3.9419, 3.1209],
           "titao": [864/Tc, 1887/Tc, 4802/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-407C of Lemmon (2003)",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "Pseudo-Pure Fluid Equations of State for the Refrigerant Blends R-410A, R-404A, R-507A, and R-407C",
                    "ref": "Int. J. Thermophys., 24(4):991-1006, 2003.",
                    "doi": "10.1023/A:1025048800563"},
        "__test__": """
            >>> st=R407c(T=300, rhom=0)
            >>> print "%0.3g %0.1f %0.1f %0.3f %0.3f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 0.0 0.0 62.631 70.945 181.04
            >>> st=R407c(T=300, P=R407c._bubbleP(300))
            >>> print "%0.3g %0.4f %0.5f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 1.2507 13.10230 78.624 133.31 458.46
            >>> st=R407c(T=300, P=R407c._dewP(300))
            >>> print "%0.3g %0.4f %0.5f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 1.0757 0.53670 74.027 99.203 154.41
            >>> st=R407c(T=250, rhom=16)
            >>> print "%0.3g %0.3f %0.1f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            250 25.372 16.0 74.065 110.74 851.38
            """, # Table V, Pag 998

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 17.04,
        "Pmin": 19.2, "rhomin": 17.04,

        "Tj": 359.345, "Pj": 4.6317,
        "dew": {"i": [0.4*2, 0.965*2, 3.1*2, 5.0*2],
                "n": [-0.086077, -6.6364, -2.4648, -3.4776]},
        "bubble": {"i": [0.54*2, 0.925*2, 2.7*2, 4.7*2],
                   "n": [0.48722, -6.6959, -1.4165, -2.5109]},

        "nr1": [0.105880e1, -0.112018e1, 0.629064, -0.351953, 0.455978e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.241, 0.69, 2.58, 1.15, 0.248],

        "nr2": [-0.175725e1, -0.112009e1, 0.277353e-1, 0.898881, -0.117591e1,
                0.818591e-1, -0.794097e-1, -0.104047e-4, 0.233779, -0.291790,
                0.154776e-1, -0.314579e-1, -0.442552e-2, -0.101254e-1,
                0.915953e-2, -0.361575e-2],
        "d2": [1, 2, 2, 3, 3, 5, 5, 5, 1, 1, 4, 4, 2, 4, 5, 6],
        "t2": [2.15, 2.43, 5.3, 0.76, 1.48, 0.24, 2.86, 8., 3.3, 4.7, 0.45,
               8.4, 16.2, 26, 16, 8.7],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
        "gamma2": [1]*16}

    eq = helmholtz1,

    _surface = {"sigma": [0.064017], "exp": [1.2557]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.086077, -6.6364, -2.4648, -3.4776],
        "exp": [0.4, 0.965, 3.1, 5.]}
    _liquid_Pressure = {
        "eq": 5,
        "ao": [0.48722, -6.6959, -1.4165, -2.5109],
        "exp": [0.54, 0.925, 2.7, 4.7]}
