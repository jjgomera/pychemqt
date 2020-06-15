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


class MD4M(MEoS):
    """Multiparameter equation of state for tetradecamethylhexasiloxane"""
    name = "tetradecamethylhexasiloxane"
    CASNumber = "107-52-8"
    formula = "C14H42O5Si6"
    synonym = "MD4M"
    _refPropName = "MD4M"
    _coolPropName = "MD4M"
    rhoc = unidades.Density(285.6576532213632)
    Tc = unidades.Temperature(653.2)
    Pc = unidades.Pressure(877.47, "kPa")
    M = 458.99328  # g/mol
    Tt = unidades.Temperature(214.15)
    Tb = unidades.Temperature(532.723)
    f_acent = 0.825
    momentoDipolar = unidades.DipoleMoment(1.308, "Debye")

    f = 8.314472
    CP1 = {"ao": -20.071/f,
           "an": [2228.5e-3/f, -1311.4e-6/f, 286.2e-9/f],
           "pow": [1, 2, 3]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD4M of Colonna  (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A., "
                             "Lemmon, E.W.",
                    "title": "Multiparameter Equations of State for Selected "
                             "Siloxanes",
                    "ref": "Fluid Phase Equilibria, 244:193-211, 2006.",
                    "doi":  "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.09,

        "nr1": [1.18492421, -1.87465636, -0.65713510e-1, -0.61812689,
                0.19535804, 0.50678740e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [1.23544082, 0.49462708e-1, -0.73685283, -0.19991438,
                -0.55118673e-1, 0.28325885e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = colonna,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.10532e2, 0.33939e1, -0.89744e1, -0.56150e1],
        "t": [1.0, 1.5, 2.75, 5.1]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.10453e1, 0.55476, 0.44536e1, -0.76179e1, 0.46237e1],
        "t": [0.235, 0.6, 0.95, 1.35, 1.7]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.10890e1, -0.84374e1, -0.35615e2, -0.73478e3, 0.19915e4,
              -0.16317e4],
        "t": [0.231, 0.8, 2.9, 7.7, 9.0, 10.0]}
