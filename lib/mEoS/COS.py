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


class COS(MEoS):
    """Multiparameter equation of state for carbonyl sulfide"""
    name = "carbonyl sulfide"
    CASNumber = "463-58-1"
    formula = "COS"
    synonym = ""
    rhoc = unidades.Density(445.1565)
    Tc = unidades.Temperature(378.77)
    Pc = unidades.Pressure(6370.0, "kPa")
    M = 60.0751  # g/mol
    Tt = unidades.Temperature(134.3)
    Tb = unidades.Temperature(222.99)
    f_acent = 0.0978
    momentoDipolar = unidades.DipoleMoment(0.7152, "Debye")
    id = 219

    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [-3.6587449805, 3.7349245016],
           "ao_exp": [2.1651, 0.93456, 1.0623, 0.34269],
           "titao": [768/Tc, 1363/Tc, 3175/Tc, 12829/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for carbonyl sulfide of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids",
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"},
        "__test__": """
            >>> st=COS(T=380, rho=7*60.0751)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            380 7 6498.429 16511.877 51.563 55.861 4139.577 161.717
            """, # Table 10, Pag 842

        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 650., "Pmax": 50000.0, "rhomax": 22.52,
        "Pmin": 0.064, "rhomin": 22.5,

        "nr1": [0.94374, -2.5348, 0.59058, -0.021488, 0.082083, 0.00024689],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.21226, -0.041251, -0.22333, -0.050828, -0.028333, 0.016983],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.07246], "exp": [1.407]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.67055e1, 0.34248e1, -0.26677e1, -0.24717e1],
        "exp": [1., 1.5, 1.78, 4.8]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.76592e1, -0.19226e2, 0.27883e2, -0.23637e2, 0.99803e1],
        "exp": [0.515, 0.767, 1.034, 1.4, 1.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.32494e1, -0.71460e1, 0.35026e2, -0.34039e2, -0.64206e2, -0.15225e3],
        "exp": [0.423, 1.464, 5.3, 4.1, 7.0, 17.0]}
