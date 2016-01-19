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


from lib.meos import MEoS
from lib import unidades


class R40(MEoS):
    """Multiparameter equation of state for R40"""
    name = "methyl chloride"
    CASNumber = "74-87-3"
    formula = "CH3Cl"
    synonym = "R40"
    rhoc = unidades.Density(363.219)
    Tc = unidades.Temperature(416.3)
    Pc = unidades.Pressure(6677.3, "kPa")
    M = 50.48752  # g/mol
    Tt = unidades.Temperature(175.0)
    Tb = unidades.Temperature(249.173)
    f_acent = 0.243
    momentoDipolar = unidades.DipoleMoment(1.871, "Debye")
    id = 115

    Fi1 = {"ao_log": [1, 2.92518],
           "pow": [0, 1, -1, -2, -3],
           "ao_pow": [7.499423, -2.997533, -6.0842e-2, -1.1525e-1, 1.0843e-2],
           "ao_exp": [3.764997],
           "titao": [3.7101]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R40 of Thol et al. (2013).",
        "__doi__": {"autor": "Thol, M., Piazza, L., and Span, R.",
                    "title": "A New Functional Form for Equations of State for Some Weakly Associating Fluids",
                    "ref": "Int. J. Thermophys., 35(5):783-811, 2014.",
                    "doi": "10.1007/s10765-014-1633-1"},
        "__test__": """
            >>> st=R40(T=240, rho=0.1)
            >>> print "%0.0f %0.1f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            240 0.1 0.003946471 -44.893076233 0.367418153 0.570253307 0.736404894 225.580390587
            >>> st=R40(T=240, rho=1050)
            >>> print "%0.0f %0.0f %0.9f %0.9f %0.9f %0.9f %0.9f %0.8f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            240 1050 27.686694086 -469.385073839 -1.974262988 1.057989833 1.538497999 1218.61504241
            >>> st=R40(T=400, rho=0.1)
            >>> print "%0.0f %0.1f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 0.1 0.006584665 89.522556776 0.70750466 0.790284451 0.955280932 282.066952997
            >>> st=R40(T=400, rho=900)
            >>> print "%0.0f %0.0f %0.9f %0.9f %0.9f %0.9f %0.9f %0.8f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 900 85.792012509 -200.464775222 -1.303719709 0.923813791 1.463531785 1077.73764207
            """, # Table 9, Pag 26

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 230.0, "Tmax": 630.0, "Pmax": 100000.0, "rhomax": 21.78756,
        "Pmin": 0.9, "rhomin": 21.78756,

        "nr1": [.274572058, .104235924e-1, -.125727710e1, 2.25609199e-3,
                -3.31830421e-2, .918440878e-1, .261059608e-2],
        "d1": [1, 1, 1, 2, 3, 3, 5],
        "t1": [-0.75, -0.25, 1.25, 0.75, -1.0, -0.375, 1.25],

        "nr2": [-.948880966e-1, -.843634836e-1, .226263660, -.470765940e-1,
                -.196610405, -.204318929e-1, -.692145009e-1, .148974844e-1,
                -.642544485e-2],
        "d2": [1, 1, 2, 5, 1, 3, 4, 5, 2],
        "t2": [2.375, 3.0, 2.625, 1.875, 4.5, 5.75, 5.375, 2.75, 14.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3],
        "gamma2": [1]*9}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.5074, 0.7520, -9.4148, 19.654, -20.190],
        "exp": [1.0, 1.5, 4.5, 5.8, 7.1]}
    _liquid_Density = {
        "eq": 1,
        "ao": [2.1809, 0.9228, -2.4615, 7.9722, -13.023, 9.2840],
        "exp": [0.37, 1.16, 2.0, 2.9, 3.9, 5.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.9433, -6.8001, -82.752, 202.14, -264.16, 99.135],
        "exp": [0.18, 0.9, 3.7, 4.6, 5.6, 6.7]}
