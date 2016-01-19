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


class oXylene(MEoS):
    """Multiparameter equation of state for o-xylene """
    name = "o-xylene "
    CASNumber = "95-47-6"
    formula = "C8H10"
    synonym = "1,2-dimethylbenzene"
    rhoc = unidades.Density(285)
    Tc = unidades.Temperature(630.259)
    Pc = unidades.Pressure(3737.5, "kPa")
    M = 106.165  # g/mol
    Tt = unidades.Temperature(247.985)
    Tb = unidades.Temperature(417.521)
    f_acent = 0.312
    momentoDipolar = unidades.DipoleMoment(0.63, "Debye")
    id = 42

    Fi1 = {"ao_log": [1, 2.748798],
           "pow": [0, 1],
           "ao_pow": [10.137376, -0.91282993],
           "ao_exp": [4.754892, 6.915052, 15.84813, 10.93886],
           "titao": [225/Tc, 627/Tc, 1726/Tc, 4941/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for o-xylene of Zhou et al. (2012).",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., and Wu, J.",
                    "title": "Thermodynamic Properties of o-Xylene, m-Xylene, p-Xylene, and Ethylbenzene",
                    "ref": "J. Phys. Chem. Ref. Data 41, 023103 (2012).",
                    "doi": "10.1063/1.3703506"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 70000.0, "rhomax": 8.648,
        "Pmin": 0.0228, "rhomin": 8.647,

        "nr1": [0.0036765156, -0.13918171, 0.014104203, 1.5398899, -2.3600925,
                -0.44359159, 0.19596977],
        "d1": [5, 1, 4, 1, 1, 2, 3],
        "t1": [1.0, 0.6, 0.91, 0.3, 0.895, 1.167, 0.435],

        "nr2": [-1.0909408, -0.21890801, 1.1179223, -0.93563815, -0.018102996],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.766, 3.8, 1.31, 3.0, 0.77],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.4172368, -0.57134695, -0.081944041, -40.682878],
        "d3": [1, 1, 3, 3],
        "t3": [1.41, 4.8, 1.856, 2.0],
        "alfa3": [1.1723, 1.095, 1.6166, 20.4],
        "beta3": [2.442, 1.342, 3.0, 450.0],
        "gamma3": [1.2655, 0.3959, 0.7789, 1.162],
        "epsilon3": [0.552, 0.728, 0.498, 0.894]}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.2834, -1.5813, 7.6516, -7.9953, -2.2277],
        "exp": [1.0, 1.5, 1.9, 2.4, 6.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.9743, 16.511, -52.934, 87.962, -71.719, 22.569],
        "exp": [0.3, 0.96, 1.4, 1.9, 2.4, 3.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-1.29038, -33.3428, 142.046, -292.211, 293.950, -159.504, -88.2170],
        "exp": [0.32, 1.14, 1.7, 2.2, 2.8, 3.5, 9.8]}
