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


class nC14(MEoS):
    """Multiparameter equation of state for n-tetradecane"""
    name = "n-Tetradecane"
    CASNumber = "629-59-4"
    formula = "C14H30"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(227.98153796)
    Tc = unidades.Temperature(692.547)
    Pc = unidades.Pressure(1532, "kPa")
    M = 198.388  # g/mol
    Tt = unidades.Temperature(279)
    Tb = unidades.Temperature(526.7)
    f_acent = 0.6585
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 18

    CP1 = {"ao": 0,
           "an": [30.175], "pow": [0.348173],
           "ao_exp": [61.844, 419.85],
           "exp": [284.606, 1696.05]}

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
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 707.5, "Pmax": 150000.0, "rhomax": 796.1/M,

        "nr1": [1.52213, -3.02849, 0.294593, -0.183821, 0.100939, 0.000308913],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.32, 1.23, 1.5, 1.4, 0.07, 0.8],

        "nr2": [1.00568, 0.036571, -0.469823, -0.112874, -0.0293384, .0128658],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [2.16, 1.1, 4.1, 5.6, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (huber, )
