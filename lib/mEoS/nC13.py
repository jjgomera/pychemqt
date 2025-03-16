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


from lib import unidades
from lib.meos import MEoS


class nC13(MEoS):
    """Multiparameter equation of state for n-tridecane"""
    name = "n-Tridecane"
    CASNumber = "629-50-5"
    formula = "C13H28"
    synonym = ""
    _refPropName = ""
    _coolPropName = ""
    rhoc = unidades.Density(224.22906625)
    Tc = unidades.Temperature(675.634)
    Pc = unidades.Pressure(1691, "kPa")
    M = 184.361  # g/mol
    Tt = unidades.Temperature(267.8)
    Tb = unidades.Temperature(508.6)
    f_acent = 0.6156
    momentoDipolar = unidades.DipoleMoment(0, "Debye")
    id = 17

    CP1 = {"ao": 0,
           "an": [50.4642], "pow": [0.278912],
           "ao_exp": [338.579, 103.942],
           "exp": [1534.11, 2780.82]}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for n-tridecane (2011)",
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

        "Tmin": Tt, "Tmax": 690.6, "Pmax": 500000.0, "rhomax": 897.7/M,

        "nr1": [1.39164, -2.87332, 0.305134, -0.179527, 0.093016, 0.0002761],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.32, 1.23, 1.5, 1.4, 0.07, 0.8],

        "nr2": [0.986319, 0.0360788, -0.489974, -0.135264, -0.0302586,
                0.0147871],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [2.16, 1.1, 4.1, 5.6, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = (huber, )
