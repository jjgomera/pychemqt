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
from lib.mEoS import C3


class nC15(MEoS):
    """Multiparameter equation of state for n-pentadecane"""
    name = "n-Pentadecane"
    CASNumber = "629-62-9"
    formula = "C15H32"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(226.5873288)
    Tc = unidades.Temperature(706.882)
    Pc = unidades.Pressure(1481, "kPa")
    M = 212.415  # g/mol
    Tt = unidades.Temperature(283.1)
    Tb = unidades.Temperature(543.8)
    f_acent = 0.6944
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 19

    CP1 = {"ao": 0,
           "an": [43.3493], "pow": [0.31566],
           "ao_exp": [26.192, 386.293, 97.851],
           "exp": [197.899, 1552.27, 2654.94]}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for  (2011)",
        "__doi__": {
            "autor": "Huber, M.L., Bruno, T.J., Chirico, R.D., Diky, V., "
                     "Kazakov, A.F., Lemmon, E.W., Muzny, C.D., Frenkel, M.",
            "title": "Equations of State on Demand: Application for Surrogate "
                     "Fuel Development",
            "ref": "Int. J. Thermophys. 32(3) (2011) 596-613",
            "doi": "10.1007/s10765-010-0909-3"},

        "R": 8.3144598,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 721.9, "Pmax": 654600.0, "rhomax": 893.2/M,

        "nr1": [1.58458, -3.08023, 0.302558, -0.186372, 0.0971768, .000300649],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.32, 1.23, 1.5, 1.4, 0.07, 0.8],

        "nr2": [0.961799, 0.0400807, -0.4855, -0.122765, -0.0346945, .0199921],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [2.16, 1.1, 4.1, 5.6, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (huber, )
