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


class mXylene(MEoS):
    """Multiparameter equation of state for m-xylene"""
    name = "m-xylene"
    CASNumber = "108-38-3"
    formula = "C8H10"
    synonym = "1,3-dimethylbenzene"
    rhoc = unidades.Density(282.929725)
    Tc = unidades.Temperature(616.89)
    Pc = unidades.Pressure(3534.6, "kPa")
    M = 106.165  # g/mol
    Tt = unidades.Temperature(225.3)
    Tb = unidades.Temperature(412.214)
    f_acent = 0.326
    momentoDipolar = unidades.DipoleMoment(0.3, "Debye")
    id = 43

    Fi1 = {"ao_log": [1, 1.169909],
           "pow": [0, 1],
           "ao_pow": [12.652887, -0.45975624],
           "ao_exp": [4.44312, 2.862794, 24.83298, 16.26077],
           "titao": [160/Tc, 190/Tc, 1333/Tc, 3496/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for m-xylene of Zhou et al. (2012).",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., and Wu, J.",
                    "title": "Thermodynamic Properties of o-Xylene, m-Xylene, p-Xylene, and Ethylbenzene",
                    "ref": "J. Phys. Chem. Ref. Data 41, 023103 (2012).",
                    "doi": "10.1063/1.3703506"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 200000.0, "rhomax": 8.677,
        "Pmin": 0.003123, "rhomin": 8.677,

        "nr1": [0.000012791017, 0.041063111, 1.505996, -2.3095875, -0.46969,
                0.171031],
        "d1": [8, 4, 1, 1, 2, 3],
        "t1": [1.0, 0.91, 0.231, 0.772, 1.205, 0.323],

        "nr2": [-1.001728, -0.3945766, 0.6970578, -0.3002876, -0.024311],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.7, 3.11, 0.768, 4.1, 0.818],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.815488, -0.330647, -0.123393, -0.54661],
        "d3": [1, 1, 3, 3],
        "t3": [2.0, 2.9, 3.83, 0.5],
        "alfa3": [1.0244, 1.3788, 0.9806, 6.3563],
        "beta3": [1.66, 1.9354, 1.0323, 78],
        "gamma3": [1.1013, 0.6515, 0.4975, 1.26],
        "epsilon3": [0.713, 0.9169, 0.6897, 0.7245]}

    eq = helmholtz1,

    _surface = {"sigma": [0.0661], "exp": [1.29]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.5635, 1.2857, -3.2346, -1.9018],
        "exp": [1.0, 1.5, 3.1, 5.6]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.43346, 3.8716, -3.0144, 1.619],
        "exp": [0.16, 0.6, 1.0, 1.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-1.1597, -6.0358, -16.712, -45.482, -98.418],
        "exp": [0.26, 0.78, 2.6, 5.7, 11.7]}
