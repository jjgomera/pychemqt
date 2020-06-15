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


class C1Linolenate(MEoS):
    """Multiparameter equation of state for methyl linolenate"""
    name = "methyl linolenate"
    CASNumber = "301-00-8"
    formula = "C19H32O2"
    synonym = ""
    _refPropName = "MLINOLEN"
    _coolPropName = "MethylLinolenate"
    rhoc = unidades.Density(247.798121314)
    Tc = unidades.Temperature(772.0)
    Pc = unidades.Pressure(1369.0, "kPa")
    M = 292.45618  # g/mol
    Tt = unidades.Temperature(218.65)
    Tb = unidades.Temperature(629.13)
    f_acent = 1.14
    momentoDipolar = unidades.DipoleMoment(1.54, "Debye")

    f = 8.314472
    CP1 = {"an": [79.5913/f], "pow": [0.214648],
           "ao_exp": [290.379/f, 81.4323/f, 474.881/f],
           "exp": [1213.24, 578.752, 2799.79]}

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

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 3.29,

        "nr1": [0.4070829e-1, 2.412375, -3.756194, -0.1526466, 0.4682918e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.15, 1.24, 1.6, 1.28],

        "nr2": [-1.470958, -0.76455, 1.908964, -1.629366, -0.1242073e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.9, 3.15, 2.16, 2.8, 1.4],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [2.180707, -0.7537264, -0.4347781],
        "d3": [1, 1, 3],
        "t3": [2.5, 3.0, 3.1],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76]}

    eq = huber,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.14239e2, 0.85361e1, -0.29678e2, 0.29106e2, -0.16922e2],
        "t": [1.0, 1.5, 2.86, 3.5, 4.6]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.16258e2, -0.19465e3, 0.14231e4, -0.22420e4, 0.10012e4],
        "t": [0.681, 1.26, 1.58, 1.7, 1.8]}
    _vapor_Density = {
        "eq": 2,
        "n": [-11.463, 45.192, -65.779, -0.18386e4, 0.40689e4, -0.25124e4],
        "t": [0.65, 1.55, 1.8, 6.6, 7.2, 7.8]}
