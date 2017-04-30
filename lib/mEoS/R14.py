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


class R14(MEoS):
    """Multiparameter equation of state for R14"""
    name = "tetrafluoromethane"
    CASNumber = "75-73-0"
    formula = "CF4"
    synonym = "R14"
    rhoc = unidades.Density(625.66)
    Tc = unidades.Temperature(227.51)
    Pc = unidades.Pressure(3750.0, "kPa")
    M = 88.0046  # g/mol
    Tt = unidades.Temperature(89.54)
    Tb = unidades.Temperature(145.10)
    f_acent = 0.1785
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 218

    CP1 = {"ao": 3.9465247,
           "an": [-.88586725e-2, 0.13939626e-3, -0.30056204e-6, 0.20504001e-9],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-14 of Platzer et al. (1990)",
        "__doi__": {"autor": "Platzer, B., Polt, A., and Maurer, G.",
                    "title": "Thermophysical properties of refrigerants",
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""},
        "R": 8.31451,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 120.0, "Tmax": 623.0, "Pmax": 51000.0, "rhomax": 20.764,
        "Pmin": 0.64144, "rhomin": 20.764,

        "nr1": [-.334698748966, .586690904687, -.147068929692, .103999039623e1,
                -.245792025288e1, .799614557889, -.749498954929, .152177772502,
                -.293408331764, .717794502866, -.426467444199e-1, .226562749365,
                -0.391091694003, -0.257394804936e-1, 0.554844884782e-1,
                0.610988261204e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [.334698748966, -.586690904687, .147068929692, -.190315426142,
                .716157133959, -.703161904626],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.99832625]*6}

    eq = helmholtz1,

    _surface = {"sigma": [0.0423], "exp": [1.24]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.61905e1, -0.91398e1, 0.12192e2, -0.47215e1, -0.20439e1],
        "exp": [1.0, 1.5, 1.64, 2.5, 7.3]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.10612e1, 0.44343e1, -0.38753e1, 0.29825e1, 0.30746],
        "exp": [0.1, 0.24, 0.4, 0.6, 3.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.55804e2, .10868e3, -.64257e2, -.11954e4, .36688e4, -.25956e4],
        "exp": [0.713, 0.84, 1.0, 5.8, 6.3, 6.6]}
