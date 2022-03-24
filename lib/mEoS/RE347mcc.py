#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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


class RE347mcc(MEoS):
    """Multiparameter equation of state for RE347mcc"""
    name = "methyl-heptafluoropropyl-ether"
    CASNumber = "375-03-1"
    formula = "CF3CF2CF2OCH3"
    synonym = "HFE-7000"
    _refPropName = "RE347MCC"
    _coolPropName = ""
    rhoc = unidades.Density(524.143687088)
    Tc = unidades.Temperature(437.7)
    Pc = unidades.Pressure(2476.2, "kPa")
    M = 200.0548424  # g/mol
    Tt = unidades.Temperature(250)
    Tb = unidades.Temperature(307.349)
    f_acent = 0.411
    momentoDipolar = unidades.DipoleMoment(3.13, "Debye")

    CP1 = {"ao": 13.09,
           "ao_exp": [13.78, 14.21], "exp": [2045, 850]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for RE347mcc of Zhou (2012)",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., Mahmoud, A.M.",
                    "title": "Equations of state for RE245cb2, RE347mcc, "
                             "RE245fa2 and R1216",
                    "ref": "Preliminary equation",
                    "doi":  ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 20000.0, "rhomax": 7.662,

        "nr1": [0.0330627, 2.606165, -4.902937, 2.228012, 1.494115, -2.420459,
                0.160067],
        "d1": [4, 1, 1, 1, 2, 2, 3],
        "t1": [1, 0.34, 0.77, 1.02, 0.79, 1.017, 0.634],

        "nr2": [1.383893, -2.092005, -0.5904708],
        "d2": [2, 1, 2],
        "t2": [1.35, 2.25, 2.5],
        "c2": [1, 2, 2],
        "gamma2": [1]*6,

        "nr3": [-0.701794, 2.765425, 0.6860982, -2.208170, 0.1739594,
                -0.9028007, -0.0213123],
        "d3": [1, 1, 2, 2, 3, 3, 1],
        "t3": [2, 1.66, 1.33, 2.0, 1.87, 1.75, 1.05],
        "alfa3": [0.593, 1.36, 1.73, 1.483, 0.617, 1.596, 9.64],
        "beta3": [0.0872, 1.176, 1.53, 0.78, 0.088, 1.04, 263.0],
        "gamma3": [1.06, 1.22, 0.92, 1.08, 1.21, 0.85, 1.12],
        "epsilon3": [1.12, 0.79, 1.055, 0.5, 0.84, 0.85, 0.91]}

    eq = zhou,

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.05031], "exp": [1.232]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-8.0456, 2.6285, -2.7498, -5.4277, -4.3693],
        "t": [1.0, 1.5, 2., 4.25, 12.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.5144, 2.3745, -2.6363, 2.0830, 0.50537],
        "t": [0.29, 0.85, 1.5, 2.2, 9.]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.0640, -6.4226, -18.982, -58.689, -117.64, -253.93],
        "t": [0.321, 0.96, 2.75, 5.9, 12., 22.]}
