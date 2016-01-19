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


class pXylene(MEoS):
    """Multiparameter equation of state for p-xylene"""
    name = "p-xylene"
    CASNumber = "106-42-3"
    formula = "C8H10"
    synonym = "1,4-dimethylbenzene"
    rhoc = unidades.Density(286)
    Tc = unidades.Temperature(616.168)
    Pc = unidades.Pressure(3531.5, "kPa")
    M = 106.165  # g/mol
    Tt = unidades.Temperature(286.4)
    Tb = unidades.Temperature(411.47)
    f_acent = 0.324
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 44

    Fi1 = {"ao_log": [1, 4.2430504],
           "pow": [0, 1],
           "ao_pow": [5.9815241, -0.52477835],
           "ao_exp": [5.2291378, 19.549862, 16.656178, 5.9390291],
           "titao": [414/Tc, 1256/Tc, 2649/Tc, 6681/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for p-xylene of Zhou et al. (2012).",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W., and Wu, J.",
                    "title": "Thermodynamic Properties of o-Xylene, m-Xylene, p-Xylene, and Ethylbenzene",
                    "ref": "J. Phys. Chem. Ref. Data 41, 023103 (2012).",
                    "doi": "10.1063/1.3703506"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 200000.0, "rhomax": 8.166,
        "Pmin": 0.580, "rhomin": 8.165,

        "nr1": [0.0010786811, -0.103161822, 0.0421544125, 1.47865376, -2.4266,
                -0.46575193, 0.190290995],
        "d1": [5, 1, 4, 1, 1, 2, 3],
        "t1": [1.0, 0.83, 0.83, 0.281, 0.932, 1.1, 0.443],

        "nr2": [-1.06376565, -0.209934069, 1.25159879, -0.951328356,
                -0.0269980032],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.62, 2.5, 1.2, 3.0, 0.778],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.37103180, -0.494160616, -0.0724317468, -3.69464746],
        "d3": [1, 1, 3, 3],
        "t3": [1.13, 4.5, 2.2, 2.0],
        "alfa3": [1.179, 1.065, 1.764, 13.675],
        "beta3": [2.445, 1.483, 4.971, 413.0],
        "gamma3": [1.267, 0.4242, 0.864, 1.1465],
        "epsilon3": [0.54944, 0.7234, 0.4926, 0.8459]}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.7221, 1.5789, -13.035, 18.453, -11.345],
        "exp": [1.0, 1.5, 3.8, 4.6, 5.5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.1783, 3.4488, -2.3906, 1.5933],
        "exp": [0.15, 0.5, 0.9, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-6.17784, -0.38825, -19.0575, -541.124, 1251.55, -920.22],
        "exp": [0.653, 0.17, 2.6, 7.8, 8.9, 10.]}
