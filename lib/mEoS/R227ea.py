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


class R227ea(MEoS):
    """Multiparameter equation of state for R227ea"""
    name = "1,1,1,2,3,3,3-heptafluoropropane"
    CASNumber = "431-89-0"
    formula = "CF3CHFCF3"
    synonym = "R227ea"
    rhoc = unidades.Density(594.25)
    Tc = unidades.Temperature(374.9)
    Pc = unidades.Pressure(2925.0, "kPa")
    M = 170.02886  # g/mol
    Tt = unidades.Temperature(146.35)
    Tb = unidades.Temperature(256.81)
    f_acent = 0.357
    momentoDipolar = unidades.DipoleMoment(1.456, "Debye")
    id = 671
    # id = 1872

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [11.43, 12.83], "exp": [403, 1428],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-227ea of McLinden and Lemmon (2013).",
        "__doi__": {"autor": "McLinden, M.O. and Lemmon, E.W.",
                    "title": "Thermodynamic Properties of R-227ea, R-365mfc, R-115, and R-13I1",
                    "ref": "to be submitted to J. Chem. Eng. Data, 2013.",
                    "doi": ""},
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 475.0, "Pmax": 60000.0, "rhomax": 11.05,
        "Pmin": 0.0073, "rhomin": 11.05,

        "nr1": [2.024341, -2.605930, 0.4957216, -0.8240820, 0.06543703],
        "d1": [1, 1, 2, 2, 4],
        "t1": [0.34, 0.77, 0.36, 0.9, 1],

        "nr2": [-1.02461, .6247065, .2997521, -.353917, -1.232043, -.8824483],
        "d2": [1, 3, 6, 6, 2, 3],
        "t2": [2.82, 2.1, 0.9, 1.13, 3.8, 2.75],
        "c2": [1, 1, 1, 1, 2, 2],
        "gamma2": [1]*6,

        "nr3": [0.1349661, -0.2662928, 0.1764733, 0.01536163, -0.004667185,
                -11.70854, 0.9114512],
        "d3": [1, 2, 1, 1, 4, 2, 1],
        "t3": [1.5, 1.5, 2.5, 5.4, 4, 1, 3.5],
        "alfa3": [0.83, 2.19, 2.44, 3.65, 8.88, 8.23, 2.01],
        "beta3": [1.72, 5.2, 2.31, 1.02, 5.63, 50.9, 1.56],
        "gamma3": [0.414, 1.051, 1.226, 1.7, 0.904, 1.42, 0.926],
        "epsilon3": [1.13, 0.71, 1.2, 1.7, 0.546, 0.896, 0.747]}

    eq = helmholtz1,

    _surface = {"sigma": [0.06127, -0.009516, -0.00192],
                "exp": [1.192, 0.9795, 1.421]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.77961e1, 0.21366e1, -0.26023e1, -0.57444e1, 0.23982e1],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.20032e1, 0.49235, 0.13738, 0.21057, -0.12834],
        "exp": [0.345, 0.74, 1.2, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.2135e1, -.68425e1, -.21447e2, -.20457e3, .51795e3, -.45908e3],
        "exp": [0.324, 1.03, 3.0, 7.4, 9.0, 10.0]}
