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


class EthylBenzene(MEoS):
    """Multiparameter equation of state for ethylbenzene"""
    name = "ethylbenzene"
    CASNumber = "100-41-4"
    formula = "C8H10"
    synonym = ""
    _refPropName = "EBENZENE"
    _coolPropName = "EthylBenzene"
    rhoc = unidades.Density(291.)
    Tc = unidades.Temperature(617.12)
    Pc = unidades.Pressure(3622.4, "kPa")
    M = 106.165  # g/mol
    Tt = unidades.Temperature(178.2)
    Tb = unidades.Temperature(409.314)
    f_acent = 0.305
    momentoDipolar = unidades.DipoleMoment(0.6, "Debye")
    id = 45

    Fi1 = {"ao_log": [1, 4.2557889],
           "pow": [0, 1],
           "ao_pow": [5.70409, -0.52414353],
           "ao_exp": [9.7329909, 11.201832, 25.440749],
           "titao": [585/Tc, 4420/Tc, 1673/Tc]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylbenzene of Zhou et "
                    "al. (2012).",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., and Wu, J.",
                    "title": "Thermodynamic Properties of o-Xylene, m-Xylene, "
                             "p-Xylene, and Ethylbenzene",
                    "ref": "J. Phys. Chem. Ref. Data 41, 023103 (2012).",
                    "doi": "10.1063/1.3703506"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 60000.0, "rhomax": 9.124,
        "Pmin": 0.000004002, "rhomin": 9.123,

        "nr1": [0.0018109418, -0.076824284, 0.041823789, 1.5059649, -2.4122441,
                -0.47788846, 0.18814732],
        "d1": [5, 1, 4, 1, 1, 2, 3],
        "t1": [1, 1, 0.92, 0.27, 0.962, 1.033, 0.513],

        "nr2": [-1.0657412, -0.20797007, 1.1222031, -0.99300799, -0.027300984],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.31, 3.21, 1.26, 2.29, 1.0],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.3757894, -0.44477155, -0.07769742, -2.16719],
        "d3": [1, 1, 3, 3],
        "t3": [0.6, 3.6, 2.1, 0.5],
        "alfa3": [1.178, 1.07, 1.775, 15.45],
        "beta3": [2.437, -1.488, -4, -418.6],
        "gamma3": [1.2667, 0.4237, 0.8573, 1.15],
        "epsilon3": [0.5494, 0.7235, 0.493, 0.8566]}

    eq = zhou,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.8411, 2.5921, -3.502, -2.7613],
        "exp": [1.0, 1.5, 2.5, 5.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [3.5146, -3.7537, 5.476, -3.4724, 1.2141],
        "exp": [0.43, 0.83, 1.3, 1.9, 3.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.2877, -3.6071, -15.878, -53.363, -128.57],
        "exp": [0.42, 0.98, 2.48, 5.9, 13.4]}
