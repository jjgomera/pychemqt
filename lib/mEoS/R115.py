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


class R115(MEoS):
    """Multiparameter equation of state for R115"""
    name = "chloropentafluoroethane"
    CASNumber = "76-15-3"
    formula = "CClF2CF3"
    synonym = "R115"
    rhoc = unidades.Density(614.77633568)
    Tc = unidades.Temperature(353.1)
    Pc = unidades.Pressure(3129.0, "kPa")
    M = 154.466416  # g/mol
    Tt = unidades.Temperature(173.75)
    Tb = unidades.Temperature(233.9)
    f_acent = 0.248
    momentoDipolar = unidades.DipoleMoment(0.52, "Debye")
    id = 229

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [7.142, 10.61],
           "exp": [289, 1301],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 2.4409547,
           "an": [0.053544743, -0.81861429e-4, 0.10410538e-6, -0.71645701e-10],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-115 of McLinden and Lemmon (2013)",
        "__doi__": {"autor": "McLinden, M.O. and Lemmon, E.W.",
                    "title": "Thermodynamic Properties of R-227ea, R-365mfc, R-115, and R-13I1",
                    "ref": "to be submitted to J. Chem. Eng. Data, 2013.",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 60000.0, "rhomax": 11.3,
        "Pmin": 2.2, "rhomin": 11.3,

        "nr1": [1.20873, -3.54460, 0.745302, 0.114128, 0.000436572],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.988385, 1.13878, -0.0215633, -0.630230, 0.0167901,
                -0.149412, -0.0271153],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Bender equation of state for R-115 of Platzer et al. (1990).",
        "__doi__": {"autor": "Platzer, B., Polt, A., and Maurer, G.",
                    "title": "Thermophysical properties of refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},

        "R": 8.31451,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 200.0 , "Tmax": 450.0, "Pmax": 7000.0, "rhomax": 10.7,
        "Pmin": 6.213, "rhomin": 10.743,

        "nr1": [-0.377294477051, -0.695891789165e-1, 0.206972205161,
                0.266609543946, -0.117158857583e1, 0.817521154071,
                -0.978729789251, -0.174482448760, 0.143598704796e1,
                -0.265460417723e1, 0.165212655822e1, -0.588257570097,
                0.738774518022, 0.296779702685, -0.534330750773, 0.659766160237e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.377294477051, 0.695891789165e-1, -0.206972205161,
                -0.350603135603, 0.108682541098e1, -0.619304197853],
        "d2": [2]*6,
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [0, 0, 0, 2, 2, 2],
        "gamma2": [1.50553819]*6}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.04771], "exp": [1.246]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.77016e1, 0.43462e1, -0.40020e1, -0.65510e1, 0.39278e1],
        "exp": [1.0, 1.5, 1.9, 5.2, 6.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18245e2, -0.57373e2, 0.78511e2, -0.50979e2, 0.14361e2],
        "exp": [0.556, 0.75, 0.95, 1.2, 1.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.35696e1, -0.83593e1, -0.34007e3, 0.40114e3, 0.84442e2, -0.22137e3],
        "exp": [0.421, 1.5, 4.7, 5.0, 5.4, 6.0]}
