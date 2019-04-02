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
from lib.mEoS import C3


class RC318(MEoS):
    """Multiparameter equation of state for RC318"""
    name = "octafluorocyclobutane"
    CASNumber = "406-58-6"
    formula = "cyclo-C4F8"
    synonym = "RC318"
    _refPropName = "RC318"
    _coolPropName = "RC318"
    rhoc = unidades.Density(619.973)
    Tc = unidades.Temperature(388.38)
    Pc = unidades.Pressure(2777.5, "kPa")
    M = 200.0312  # g/mol
    Tt = unidades.Temperature(233.35)
    Tb = unidades.Temperature(267.175)
    f_acent = 0.3553
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 692

    # Factor neccessary because Cp is not reduced by R
    f = 200.04/8.31451
    CP1 = {"ao": 0.121*f,
           "an": [0.2903e-2*f, -0.25327e-5*f, 0.77191e-9*f], "pow": [1, 2, 3],
           "ao_exp": [], "exp": []}

    platzer = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-C318 of Platzer (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., Maurer, G.",
                    "title": "Thermophysical Properties of Refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 623.0, "Pmax": 60000.0, "rhomax": 8.6452,
        "Pmin": 19.461, "rhomin": 8.6452,

        "nr1": [-0.104729119796e1, 0.138034128822e1, -0.333625769594,
                0.109415569278e1, -0.268265237178e1, 0.173403063905e1,
                -0.163611246876e1, 0.304834499143, 0.102771551787,
                -0.232367895587e-1, 0.166151957803, -0.250103914479e-1,
                0.935680977639e-1, 0.431929127445e-1, -0.133439845861,
                0.255416632165e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.104729119796e1, -0.138034128822e1, 0.333625769594,
                -0.510485781618, 0.181840728111e1, -0.138530893970e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.99943992]*6}

    eq = platzer,

    _surface = {"sigma": [0.0507], "exp": [1.25]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.78467e1, 0.24555e1, -0.30824e1, -0.58263e1, 0.35483e1],
        "t": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "n": [-0.30181, 0.29345e1, -0.13741e1, 0.14650e1, 0.16963],
        "t": [0.11, 0.32, 0.57, 0.84, 2.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.24491e2, 0.53255e2, -0.38863e2, -0.24938e2, -0.90092e2],
        "t": [0.61, 0.77, 0.92, 3.3, 7.5]}

    trnECS = {"__name__": "Huber (2003)",

              "__doi__": {
                  "autor": "Huber, M.L., Laesecke, A., Perkins, R.A.",
                  "title": "Model for the Viscosity and Thermal Conductivity "
                           "of Refrigerants, Including a New Correlation for "
                           "the Viscosity of R134a",
                  "ref": "Ind. Eng. Chem. Res., 42(13) (2003) 3163-3178",
                  "doi": "10.1021/ie0300880"},

              "eq": "ecs",

              "ref": C3,
              "visco": "visco1",
              "thermo": "thermo0",

              "ek": 299.76, "sigma": 0.5947, "omega": 5,

              "psi": [1.21141, -3.37573e-2], "psi_d": [0, 1],
              "fint": [1.35697e-3, -1.11635e-7], "fint_t": [0, 1],
              "chi": [1.5249, -0.147564], "chi_d": [0, 1],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03, "Xio": 0.194e-9,
              "gam0": 0.0496, "qd": 3.56085e-10, "Tcref": 1.5*Tc}

    _viscosity = trnECS,
    _thermal = trnECS,
