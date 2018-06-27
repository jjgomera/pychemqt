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


class D6(MEoS):
    """Multiparameter equation of state for dodecamethylcyclohexasilosane"""
    name = "dodecamethylcyclohexasiloxane"
    CASNumber = "540-97-6"
    formula = "C12H36Si6O6"
    synonym = "D6"
    _refPropName = "D6"
    rhoc = unidades.Density(279.09572983533354)
    Tc = unidades.Temperature(645.78)
    Pc = unidades.Pressure(961.0, "kPa")
    M = 444.924  # g/mol
    Tt = unidades.Temperature(270.2)
    Tb = unidades.Temperature(518.11)
    f_acent = 0.736
    momentoDipolar = unidades.DipoleMoment(1.559, "Debye")
    # id=1674

    CP1 = {"ao": 468.7,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [425104546.6, 0, 3151243909.0, 0],
           "hyp": [786.8, 0, 1792.1, 0]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hexamethyldisiloxane of "
                    "Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., and Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes:"
                             " [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and "
                             "[O-Si-(CH3)2]6",
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi":  "10.1016/j.fluid.2007.10.001"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.246,
        "Pmin": 0.00016, "rhomin": 2.245,

        "nr1": [1.69156186, -3.37962568, 0.38609039, 0.64598995e-1,
                0.10589012, 0.45456825e-4],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.74169279, -0.88102648e-1, -0.17373336, -0.10951368,
                -0.62695695e-1, 0.37459986e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = colonna,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.96557e1, 0.62155, 0.17863e1, -0.10496e2, -0.84102e1],
        "exp": [1.0, 1.5, 1.72, 3.18, 11.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.42563e2, -0.15707e3, 0.29502e3, -0.24191e3, 0.65145e2],
        "exp": [0.537, 0.68, 0.85, 1.0, 1.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.20930e1, -0.94442e1, -0.44731e2, -0.57898e2, -0.35144e2,
               -0.29661e3],
        "exp": [0.338, 1.02, 3.46, 7.1, 7.4, 15.0]}
