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


class MD3M(MEoS):
    """Multiparameter equation of state for dodecamethylpentasiloxane"""
    name = "dodecamethylpentasiloxane"
    CASNumber = "141-63-9"
    formula = "C12H36Si5O4"
    synonym = "MD3M"
    _refPropName = "MD3M"
    _coolPropName = "MD3M"
    rhoc = unidades.Density(263.9218791237794)
    Tc = unidades.Temperature(628.36)
    Pc = unidades.Pressure(945.0, "kPa")
    M = 384.839  # g/mol
    Tt = unidades.Temperature(192.0)
    Tb = unidades.Temperature(503.03)
    f_acent = 0.722
    momentoDipolar = unidades.DipoleMoment(1.223, "Debye")

    CP1 = {"ao": 463.2,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [609372332.2, 0, 4290277999.0, 0],
           "hyp": [908.5, 0, 2117.1, 0]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD3M of Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., and Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes:"
                             " [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and "
                             "[O-Si-(CH3)2]6",
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi":  "10.1016/j.fluid.2007.10.001"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 2.54,
        "Pmin": 0.4e-12, "rhomin": 2.54,

        "nr1": [1.20540386, -2.42914797, 0.69016432, -0.69268041, 0.18506046,
                0.31161436e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.99862519, 0.74229034e-1, -0.80259136, -0.20865337,
                -0.36461791e-1, 0.19174051e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = colonna,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.92608e1, 0.15861e1, -0.32859e1, -0.75194e1, -0.34883e1],
        "t": [1.0, 1.5, 2.46, 3.7, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.74156, 0.21723e1, 0.66412e2, -0.17125e3, 0.10848e3],
        "t": [0.22, 0.51, 5.5, 6.0, 6.4]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.19054e1, -0.74526e1, -0.10520e3, 0.24548e3, -0.23783e3,
              -0.21226e3],
        "t": [0.332, 0.88, 3.25, 4.0, 4.6, 12.0]}
