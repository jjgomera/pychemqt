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


class D5(MEoS):
    """Multiparameter equation of state for decamethylcyclopentasiloxane"""
    name = "decamethylcyclopentasiloxane"
    CASNumber = "541-02-6"
    formula = "C10H30O5Si5"
    synonym = "D5"
    _refPropName = "D5"
    _coolPropName = "D5"
    rhoc = unidades.Density(292.570762680819)
    Tc = unidades.Temperature(619.23462341)
    Pc = unidades.Pressure(1161.46, "kPa")
    M = 370.7697  # g/mol
    Tt = unidades.Temperature(226.0)
    Tb = unidades.Temperature(484.05)
    f_acent = 0.658
    momentoDipolar = unidades.DipoleMoment(1.349, "Debye")
    # id=1671

    f = 8.314472
    CP1 = {"ao": -34.898/f,
           "an": [1861.5e-3/f, -1403.4e-6/f, 500.0e-9/f],
           "pow": [1, 2, 3]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of "
                    "Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A., "
                             "Lemmon, E.W.",
                    "title": "Multiparameter Equations of State for Selected "
                             "Siloxanes",
                    "ref": "Fluid Phase Equilibria, 244:193-211, 2006.",
                    "doi":  "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 300, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.83,

        "nr1": [1.40844725, -2.29248044, 0.42851607, -0.73506382, 0.16103808,
                0.29643278e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.82412481, 0.15214274, -0.68495890, -0.55703624e-1,
                0.13055391e-1, -0.31853761e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = colonna,

    _surface = {
        "__doi__": {
            "autor": "Mulero, A., Cachadiña, I.",
            "title": "Recommended Correlations for the Surface Tension of "
                     "Several Fluids Included in the REFPROP Program",
            "ref": "J. Phys. Chem. Ref. Data 43(2) (2014) 023104",
            "doi": "10.1063/1.4878755"},
        "sigma": [0.04408], "exp": [1.357]}

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.99967e1, 0.70091e1, -0.72265e1, -0.62938e1],
        "t": [1.0, 1.5, 1.87, 3.8]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.303988e3, -0.110342e4, 0.134359e4, -0.705243e3, 0.164540e3],
        "t": [0.57, 0.65, 0.73, 0.84, 0.96]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.37577e1, -0.47669e1, -0.24233e2, -0.29872e3, 0.34441e3,
              -0.32498e3],
        "t": [0.459, 1.02, 2.6, 6.7, 7.7, 11.0]}
