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


class MD2M(MEoS):
    """Multiparameter equation of state for decamethyltetrasiloxane"""
    name = "decamethyltetrasiloxane"
    CASNumber = "141-62-8"
    formula = "C10H30Si4O3"
    synonym = "MD2M"
    _refPropName = "MD2M"
    _coolPropName = "MD2M"
    rhoc = unidades.Density(284.1716396703609)
    Tc = unidades.Temperature(599.40)
    Pc = unidades.Pressure(1227.0, "kPa")
    M = 310.685  # g/mol
    Tt = unidades.Temperature(205.2)
    Tb = unidades.Temperature(467.51)
    f_acent = 0.668
    momentoDipolar = unidades.DipoleMoment(1.12, "Debye")
    # id = 1837

    CP1 = {"ao": 331.9,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [329620742.8, 0, 2556558319.0, 0],
           "hyp": [795.1, 0, 1813.8, 0]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MD2M of Colonna (2006)",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., and Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes:"
                             " [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and "
                             "[O-Si-(CH3)2]6",
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi":  "10.1016/j.fluid.2007.10.001"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.033,
        "Pmin": 0.0000005, "rhomin": 3.032,

        "nr1": [1.33840331, -2.62939393, 0.4398383, -0.53496715, 0.1818844,
                0.40774609e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [1.13444506, 0.5774631e-1, -0.5917498, -0.11020225,
                -0.34942635e-1, 0.7646298e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = colonna,

    _vapor_Pressure = {
        "eq": 3,
        "n": [-0.10029e2, 0.44434e1, -0.36753e1, -0.68925e1, -0.32211e1],
        "t": [1.0, 1.5, 2.06, 3.5, 10.0]}
    _liquid_Density = {
        "eq": 1,
        "n": [0.12608e2, -0.32120e2, 0.33559e2, -0.11695e2, 0.76192],
        "t": [0.48, 0.64, 0.8, 1.0, 2.6]}
    _vapor_Density = {
        "eq": 2,
        "n": [-0.24684e1, -0.71262e1, -0.27601e2, -0.49458e2, -0.24106e2,
               -0.19370e3],
        "t": [0.376, 0.94, 2.9, 5.9, 6.2, 13.0]}
