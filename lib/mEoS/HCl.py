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

from lib.meos import MEoS
from lib import unidades


class HCl(MEoS):
    """Multiparameter equation of state for hydrogen chloride"""
    name = "hydrogen chloride"
    CASNumber = "7647-01-0 "
    formula = "HCl"
    synonym = ""
    rhoc = unidades.Density(410.97)
    Tc = unidades.Temperature(324.55)
    Pc = unidades.Pressure(8263.00, "kPa")
    M = 36.460939  # g/mol
    Tt = unidades.Temperature(131.4)
    Tb = unidades.Temperature(188.199)
    f_acent = 0.12875
    momentoDipolar = unidades.DipoleMoment(1.079, "Debye")
    id = 104

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1, -1, -2],
           "ao_pow": [7.913048, -3.217575, -4.149937e-3, 8.019202e-4],
           "ao_exp": [1.054392],
           "titao": [1.241138e1]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen chloride of Thol et al. (2013).",
        "__doi__": {"autor": "Thol, M., Piazza, L., and Span, R.",
                    "title": "A New Functional Form for Equations of State for Some Weakly Associating Fluids",
                    "ref": "Int. J. Thermophys., 35(5):783-811, 2014.",
                    "doi": "10.1007/s10765-014-1633-1"},
        "__test__": """
            >>> st=HCl(T=170, rho=0.01)
            >>> print "%0.0f %0.1f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            170 0.01 0.000387586 -102.414301628 0.820319743 0.571157223 0.799483403 232.898528911
            >>> st=HCl(T=170, rho=1230)
            >>> print "%0.0f %0.0f %0.9f %0.9f %0.9f %0.9f %0.9f %0.8f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            170 1230 1.229128198 -561.339398154 -2.888318217 1.149719990 1.553944836 999.438819195
            >>> st=HCl(T=280, rho=0.1)
            >>> print "%0.0f %0.1f %0.9f %0.9f %0.9f %0.9f %0.9f %0.9f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            280 0.1 0.006381824 -14.914381738 0.580027609 0.571466689 0.800079245 298.836360637
            >>> st=HCl(T=280, rho=900)
            >>> print "%0.0f %0.0f %0.9f %0.9f %0.9f %0.9f %0.9f %0.8f" % ( \
                st.T, st.rho, st.P.MPa, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            280 900 3.421902004 -371.035989751 -2.044095830 0.961502216 2.150825813 577.782761523
            """, # Table 9, Pag 26

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 155.0, "Tmax": 330.0, "Pmax": 20000.0, "rhomax": 33.8145,
        "Pmin": 0.7, "rhomin": 33.8145,

        "nr1": [-.40937325, 0.943994574, -.178830477e1, 0.128619044,
                0.439018427e-2, 0.130480908e-1, 0.169387782e-2],
        "d1": [1, 1, 1, 2, 3, 3, 5],
        "t1": [-0.75, -0.25, 1.25, 0.75, -1.0, -0.375, 1.25],

        "nr2": [0.751559060, -.800007427, 0.430935939, 0.454319457e-2,
                -.152172259, -.436174059e-1, -.970625964e-2, 0.101144098e-1,
                0.376991644e-2],
        "d2": [1, 1, 2, 5, 1, 3, 4, 5, 2],
        "t2": [2.375, 3.0, 2.625, 1.875, 4.5, 5.78, 5.375, 2.75, 14.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3],
        "gamma2": [1]*9}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 6,
        "ao": [-0.01065138, -6.15979914, 1.55860976, -8.42734117],
        "exp": [1.0, 2.0, 6.0, 11.0]}
    _liquid_Density = {
        "eq": 2,
        "ao": [1.89232034, 0.83621066, -0.22094602, 4.70971253, -5.34396174],
        "exp": [1.0, 2.0, 4.0, 11.0, 13.0]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-2.95523223, -8.10448179, -14.78392979, -87.13352586],
        "exp": [1.29, 4.2, 11.1, 24.0]}
