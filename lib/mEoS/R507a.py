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


class R507a(MEoSBlend):
    """Multiparameter equation of state for R507A (50% R125, 50% R143a)"""
    name = "R507A"
    CASNumber = ""
    formula = "R125+R143a"
    synonym = "R507A"
    rhoc = unidades.Density(490.7370688)
    Tc = unidades.Temperature(343.765)
    Pc = unidades.Pressure(3704.9, "kPa")
    M = 98.8592  # g/mol
    Tt = unidades.Temperature(200.0)
    Tb = unidades.Temperature(226.41)
    f_acent = 0.286
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 62
    # id = None

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.25],
           "ao_pow": [9.93541, 7.9985, -21.6054],
           "ao_exp": [0.95006, 4.1887, 5.5184],
           "titao": [364/Tc, 815/Tc, 1768/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-404A of Lemmon (2003)",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "Pseudo-Pure Fluid Equations of State for the Refrigerant Blends R-410A, R-404A, R-507A, and R-407C",
                    "ref": "Int. J. Thermophys., 24(4):991-1006, 2003.",
                    "doi": "10.1023/A:1025048800563"},
        "__test__": """
            >>> st=R507a(T=300, rhom=0)
            >>> print "%0.3g %0.1f %0.1f %0.3f %0.3f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 0.0 0.0 76.838 85.152 167.22
            >>> st=R507a(T=300, P=R507a._bubbleP(300))
            >>> print "%0.3g %0.4f %0.5f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 1.3462 10.50670 91.290 153.79 356.72
            >>> st=R507a(T=300, P=R507a._dewP(300))
            >>> print "%0.3g %0.4f %0.5f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            300 1.3450 0.73552 88.742 124.15 130.95
            >>> st=R507a(T=250, rhom=13)
            >>> print "%0.3g %0.3f %0.1f %0.3f %0.2f %0.2f" % (st.T, st.P.MPa, st.rhoM, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            250 12.916 13.0 83.541 123.21 697.89
            """, # Table V, Pag 998

        "R": 8.314472,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 50000.0, "rhomax": 14.13,
        "Pmin": 23.23, "rhomin": 14.13,

        "Tj": 343.765, "Pj": 3.7049,
        "dew": {"i": [1*2, 1.5*2, 2.1*2, 4.7*2],
                "n": [-7.5459, 2.338, -2.237, -4.1535]},
        "bubble": {"i": [1*2, 1.5*2, 2.2*2, 4.6*2],
                   "n": [-7.4853, 2.0115, -2.0141, -3.7763]},

        "nr1": [0.624982e1, -0.807855e1, 0.264843e-1, 0.286215, -0.507076e-2,
                0.109552e-1, 0.116124e-2],
        "d1": [1, 1, 1, 2, 2, 4, 6],
        "t1": [0.692, 0.943, 5.8, 0.77, 5.84, 0.24, 0.69],

        "nr2": [0.138469e1, -0.922473, -0.503562e-1, 0.822098, -0.277727,
                0.358172, -0.126426e-1, -0.607010e-2, -.815653e-1, -.233323e-1,
                .352952e-1, .159566e-1, .755927e-1, -.542007e-1, .170451e-1],
        "d2": [1, 1, 1, 2, 2, 3, 4, 7, 2, 3, 4, 4, 2, 3, 5],
        "t2": [2, 3, 7, 2.2, 4.3, 2.7, 1.2, 1.23, 12, 6, 8.5, 11.5, 13, 17, 16.2],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*15}

    eq = helmholtz1,

    _surface = {"sigma": [0.06701, -0.04297], "exp": [1.3066, 2.3145]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.5459, 2.338, -2.237, -4.1535],
        "exp": [1, 1.5, 2.1, 4.7]}
    _liquid_Pressure = {
        "eq": 5,
        "ao": [-7.4853, 2.0115, -2.0141, -3.7763],
        "exp": [1, 1.5, 2., 4.6]}
