#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


class R1216(MEoS):
    """Multiparameter equation of state for R1216"""
    name = "hexafluoropropene"
    CASNumber = "116-15-4"
    formula = "C3F6"
    synonym = "R1216"
    _refPropName = "R1216"
    _coolPropName = ""
    rhoc = unidades.Density(583.40757266496)
    Tc = unidades.Temperature(358.9)
    Pc = unidades.Pressure(3149.528, "kPa")
    M = 150.0225192  # g/mol
    Tt = unidades.Temperature(117.654)
    Tb = unidades.Temperature(242.81)
    f_acent = 0.333
    momentoDipolar = unidades.DipoleMoment(1.088, "Debye")
    id = 669

    CP1 = {"ao": 5.878676,
           "ao_exp": [9.351559, 9.192089, 7.983222],
           "exp": [561, 1486, 7595]}

    zhou = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R1216 of Zhou (2010)",
        "__doi__": {"autor": "Zhou, Y., Lemmon, E.W.",
                    "title": "Preliminary equation, 2010.",
                    "ref": "",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 12000.0, "rhomax": 15,

        "nr1": [.37582356e-1, .14558246e1, -.2701615e1, -.3357347, .1885495],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.3, 1.0, 1.35, 0.4],

        "nr2": [-0.1689206, 0.1122147e1, -0.6405048, -0.25931535e-1,
                0.42940852, -0.10163408e1, -0.43691328e-1],
        "d2": [3, 2, 2, 7, 1, 1, 1],
        "t2": [1., 1.68, 2.36, 0.615, 1.32, 2.12, 3.],
        "c2": [2, 1, 2, 1, 1, 2, 3],
        "gamma2": [1]*7,

        "nr3": [0.12530663e1, -0.54254994, -0.15327764, -0.92102535e-2],
        "d3": [1, 1, 3, 3],
        "t3": [0.82, 2.85, 2.83, 1.67],
        "alfa3": [0.9665, 1.503, 0.97, 5.87],
        "beta3": [1.24, 0.776, 0.86, 478],
        "gamma3": [1.284, 0.42, 0.434, 1.074],
        "epsilon3": [0.67, 0.925, 0.75, 0.73],
        "nr4": []}

    eq = (zhou, )

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.053876, 0.038318], "exp": [1.0944, 2.3239]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.9011, 3.1506, -3.0852, -4.2112, -15.438],
        "t": [1.0, 1.5, 2., 4.5, 19.]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.7159, 2.3953, -5.8035, 10.749, -10.537, 4.7535],
        "t": [0.31, 0.97, 1.7, 2.4, 3.2, 4.1]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.4969, -5.8935, -16.846, -55.082, -140.43],
        "t": [0.353, 1.05, 2.74, 6., 13.3]}
