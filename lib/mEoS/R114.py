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


class R114(MEoS):
    """Multiparameter equation of state for R114"""
    name = "1,2-dichloro-1,1,2,2-tetrafluoroethane"
    CASNumber = "76-14-2"
    formula = "CClF2CClF2"
    synonym = "R114"
    _refPropName = "R114"
    _coolPropName = "R114"
    rhoc = unidades.Density(579.969)
    Tc = unidades.Temperature(418.83)
    Pc = unidades.Pressure(3257.0, "kPa")
    M = 170.921  # g/mol
    Tt = unidades.Temperature(180.63)
    Tb = unidades.Temperature(276.741)
    f_acent = 0.25253
    momentoDipolar = unidades.DipoleMoment(0.658, "Debye")
    id = 231

    CP1 = {"ao": 0.97651380e-1,
           "an": [0.3240861e-2, -0.5895364e-5, 0.6737929e-8, -0.3546364e-11],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    platzer = {
        "__type__": "Helmholtz",
        "__name__": "Bender equation of state for R-114 of Platzer (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., Maurer, G.",
                    "title": "Thermophysical Properties of Refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 273.15, "Tmax": 507.0, "Pmax": 21000.0, "rhomax": 8.942,
        "Pmin": 0.2, "rhomin": 10.4,

        "nr1": [-0.340776521414, 0.323001398420, -0.424950537596e-1,
                0.107938879710e1, -0.199243619673e1, -0.155135133506,
                -0.121465790553, -0.165038582393e-1, -0.186915808643,
                0.308074612567, 0.115861416115, 0.276358316589e-1,
                0.108043243088, 0.460683793064e-1, -0.174821616881,
                0.317530854287e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.340776521414, -0.323001398420, 0.424950537596e-1,
                -0.166940100976e1, 0.408693082002e1, -0.241738963889e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.21103865]*6}

    eq = platzer,

    _surface = {"sigma": [0.05239], "exp": [1.258]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.72195e1, 0.16357e1, -0.14576e1, -0.69580e1, 0.57181e1],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.43023, 0.22722e2, -0.27118e2, 0.13247e2, -0.90529e1],
        "t": [0.095, 0.93, 1.1, 2.0, 3.0]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.46609, -6.8355, -167.15, 1.5805e4, -3.1859e4, 2.1548e4],
        "t": [0.09, 0.76, 4.0, 6.5, 7.0, 8.0]}
