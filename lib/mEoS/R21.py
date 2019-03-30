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


from lib import unidades
from lib.meos import MEoS


class R21(MEoS):
    """Multiparameter equation of state for R21"""
    name = "dichlorofluoromethane"
    CASNumber = "75-43-4"
    formula = "CHCl2F"
    synonym = "R21"
    _refPropName = "R21"
    _coolPropName = "R21"
    rhoc = unidades.Density(526.0138)
    Tc = unidades.Temperature(451.48)
    Pc = unidades.Pressure(5181.2, "kPa")
    M = 102.9227  # g/mol
    Tt = unidades.Temperature(142.8)
    Tb = unidades.Temperature(282.01)
    f_acent = 0.2061
    momentoDipolar = unidades.DipoleMoment(1.37, "Debye")
    id = 642

    CP1 = {"ao": 0.2376576/8.31451*102.92,
           "an": [0.12714330e-2/8.31451*102.92, 0.32413520e-6/8.31451*102.92,
                  -2.4924280e-9/8.31451*102.92, 1.7172080e-12/8.31451*102.92],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    platzer = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-21 of Platzer (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., Maurer, G.",
                    "title": "Thermophysical Properties of Refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 200.0, "Tmax": 473.19, "Pmax": 137900.0, "rhomax": 15.36,
        "Pmin": 0.6828e-4, "rhomin": 16.519,

        "nr1": [-.44386484873e2, 9.26505600935, -.551709104376, .504676623431,
                -.732431415692, -.868403860387, .146234705555, -.280576335053,
                0.864743656093, -0.270767233732e1, 0.330476390706e1,
                -0.210878239171, 0.449531449589, 0.120779813143,
                -0.277297953777, 0.305441291172e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.443864848730e2, -0.926505600935e1, 0.551709104376,
                0.121128809552e1, 0.167119476587, -0.504876793028e-1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.07470252]*6}

    eq = platzer,

    _surface = {"sigma": [0.06924], "exp": [1.259]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.70336e1, 0.15672e1, -0.33932e1, 0.17582e1, -0.86765e1],
        "t": [1.0, 1.5, 3.0, 7.0, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.33546, 0.18208e2, -0.26400e2, 0.10586e2],
        "t": [0.09, 0.78, 0.92, 1.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.38213, -0.55559e1, -0.15886e2, -0.44766e2, -0.27606e3],
        "t": [0.09, 0.667, 2.5, 6.0, 15.0]}
