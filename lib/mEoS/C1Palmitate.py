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


class C1Palmitate(MEoS):
    """Multiparameter equation of state for methyl palmitate"""
    name = "methyl palmitate"
    CASNumber = "112-39-0"
    formula = "C17H34O2"
    synonym = ""
    _refPropName = "MPALMITA"
    _coolPropName = "MethylPalmitate"
    rhoc = unidades.Density(242.59424202)
    Tc = unidades.Temperature(755.0)
    Pc = unidades.Pressure(1350.0, "kPa")
    M = 270.45066  # g/mol
    Tt = unidades.Temperature(302.71)
    Tb = unidades.Temperature(602.3)
    f_acent = 0.91
    momentoDipolar = unidades.DipoleMoment(1.54, "Debye")

    f = 8.314472
    CP1 = {"an": [120.529/f],
           "pow": [0.0801627],
           "ao_exp": [345.62/f, 289.038/f, 301.639/f],
           "exp": [2952.37, 734.653, 1593.55]}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl linoleate of Huber"
                    "et al. (2009).",
        "__doi__": {"autor": "Huber, M.L., Lemmon, E.W., Kazakov, A., Ott, "
                             "L.S., Bruno, T.J.",
                    "title": "Model for the Thermodynamic Properties of a "
                             "Biodiesel Fuel",
                    "ref": "Energy Fuels, 23 (7) (2009) 3790–3797",
                    "doi": "10.1021/ef900159g"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 302.71, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 3.36,

        "nr1": [0.4282821e-1, 2.443162, -3.75754, -0.1588526, 0.4055990e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.36, 1.22, 1.45, 0.7],

        "nr2": [-1.52409, -0.7686167, 1.79995, -1.590967, -0.1267681e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [3.0, 3.9, 2.2, 2.9, 1.25],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.2198347e1, -0.7737211, -0.4314520],
        "d3": [1, 1, 3],
        "t3": [2.6, 3.0, 3.2],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76],
        "nr4": []}

    eq = huber,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.13378e2, 0.12258e2, -0.12205e2, -0.58118e1, -0.25451e1],
        "t": [1.0, 1.5, 2.04, 4.3, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.54203, 0.13191e2, -0.20107e2, 0.11328e2, -0.60761],
        "t": [0.18, 0.5, 0.7, 0.9, 1.5]}
    _vapor_Density = {
        "eq": 2,
        "n": [-11.612, 1.63e2, -4.7913e2, 7.2986e2, -4.8202e2, -1.8198e2],
        "t": [0.65, 1.78, 2.15, 2.7, 3.1, 9.8]}
