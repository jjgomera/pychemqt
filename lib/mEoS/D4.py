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


class D4(MEoS):
    """Multiparameter equation of state for octamethylcyclotetrasiloxane"""
    name = "octamethylcyclotetrasiloxane"
    CASNumber = "556-67-2"
    formula = "C8H24O4Si4"
    synonym = "D4"
    rhoc = unidades.Density(307.0335906736056)
    Tc = unidades.Temperature(586.49127187)
    Pc = unidades.Pressure(1332.0, "kPa")
    M = 296.61576  # g/mol
    Tt = unidades.Temperature(290.25)
    Tb = unidades.Temperature(448.504)
    f_acent = 0.592
    momentoDipolar = unidades.DipoleMoment(1.090, "Debye")
    # id=1430

    CP1 = {"ao": -18.256,
           "an": [1427.2e-3, -990.20e-6, 300.0e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for D4 of Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., Guardone, A., "
                             "Lemmon, E.W.",
                    "title": "Multiparameter Equations of State for Selected "
                             "Siloxanes",
                    "ref": "Fluid Phase Equilibria, 244:193-211, 2006.",
                    "doi":  "10.1016/j.fluid.2006.04.015"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": 300.0, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.21,
        "Pmin": 0.0696, "rhomin": 3.2,

        "nr1": [1.05392408, -2.22981918, 0.77573923, -0.6937405, 0.18721557,
                0.42193330e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.70301835, 0.47851888e-1, -0.8025348, -0.18968872,
                -0.22211781e-1, 0.60103354e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = colonna,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.87935e1, 0.27204e1, -0.48174e1, -0.69086e1],
        "exp": [1.0, 1.5, 2.2, 4.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.14563e1, -0.94215, 0.45065e1, -0.27688e1, 0.8745],
        "exp": [0.24, 0.5, 0.75, 1.0, 2.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.16204e1, -0.57888e1, -0.24291e2, 0.53567e2, -0.12135e3,
               -0.10976e4],
        "exp": [0.31, 0.78, 2.5, 4.4, 5.0, 15.0]}

    visco0 = {"eq": 5, "omega": 3,
              "__doi__": {"autor": "T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E",
                          "title": "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties",
                          "ref": "Ind. Eng. Chem. Res., 1988, 27 (4), pp 671–679",
                          "doi": "10.1021/ie00076a024"},
              "__name__": "Chung (1988)",
              "w": 0.592, "mur": 0.0, "k": 0.0}

    _viscosity = visco0,
#    _thermal=visco0,
