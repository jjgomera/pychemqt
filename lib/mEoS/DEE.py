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


from lib.meos import MEoS
from lib import unidades


class DEE(MEoS):
    """Multiparameter equation of state for diethyl ether"""
    name = "diethyl ether "
    CASNumber = "60-29-7"
    formula = "C4H10O"
    synonym = ""
    rhoc = unidades.Density(264)
    Tc = unidades.Temperature(466.7)
    Pc = unidades.Pressure(3644., "kPa")
    M = 74.1216  # g/mol
    Tt = unidades.Temperature(156.92)
    Tb = unidades.Temperature(307.604)
    f_acent = 0.281
    momentoDipolar = unidades.DipoleMoment(1.151, "Debye")
    id = 162

    Fi1 = {"ao_log": [1, 3.36281],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [17.099494, -6.160844, -8.943822, 5.4621e-1, -1.6604e-2],
           "ao_exp": [], "titao": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for diethyl ether  of Thol et al. (2013).",
        "__doi__": {"autor": "Thol, M., Piazza, L., and Span, R.",
                    "title": "A New Functional Form for Equations of State for Some Weakly Associating Fluids",
                    "ref": "Int. J. Thermophys., 35(5):783-811, 2014.",
                    "doi": "10.1007/s10765-014-1633-1"},
        "__test__": """
            >>> st=DEE(T=280, rho=0.1)
            >>> print "%0.0f %0.1f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            280 0.1 0.003134775 -29.231529953 0.288943745 1.455670027 1.569255417 183.651907457
            >>> st=DEE(T=280, rho=750)
            >>> print "%0.0f %0.0f %0.9f %0.9f %0.9f %0.9f %0.9f %0.8f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            280 750 20.774279163 -396.071208713 -1.386035322 1.732367580 2.226381171 1190.89216523
            >>> st=DEE(T=400, rho=0.1)
            >>> print "%0.0f %0.1f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 0.1 0.004483768 182.466204569 0.873988712 1.841364392 1.953941717 218.048963727
            >>> st=DEE(T=400, rho=650)
            >>> print "%0.0f %0.0f %0.9f %0.9f %0.9f %0.9f %0.9f %0.8f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 650 33.847217264 -101.418571257 -0.568533506 1.949560509 2.520106238 919.575708532
            """, # Table 9, Pag 26

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 270.0, "Tmax": 500.0, "Pmax": 40000.0, "rhomax": 10.6851,
        "Pmin": 0.000001, "rhomin": 10.6851,

        "nr1": [0.376700499, -.116630334, -.73801498, -.2725701, -.4979231e-1,
                0.172267029, 0.441618910e-2],
        "d1": [1, 1, 1, 2, 3, 3, 5],
        "t1": [-0.75, -0.25, 1.25, 0.75, -1.0, -0.375, 1.25],

        "nr2": [-.153951612e1, 0.115606052e1, -.184504019e-1, -.101800599,
                -.403598704, 0.213055571e-2, -.154741976, 0.120950552e-1,
                -.143106371e-1],
        "d2": [1, 1, 2, 5, 1, 3, 4, 5, 2],
        "t2": [2.375, 3.0, 2.625, 1.875, 4.5, 5.78, 5.375, 2.75, 14.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3],
        "gamma2": [1]*9}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.3059, 1.1734, 0.7142, -4.3219],
        "exp": [1.0, 1.5, 2.2, 3.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.3275, 3.1842, -2.1407, 1.4376],
        "exp": [0.12, 0.55, 1.0, 1.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.35858, -16.843, 32.476, -33.444, -48.036],
        "exp": [0.06, 0.87, 1.3, 1.7, 5.3]}
