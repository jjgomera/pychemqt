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


class RE245cb2(MEoS):
    """Multiparameter equation of state for RE245cb2"""
    name = "methyl-pentafluoroethyl-ether"
    CASNumber = "22410-44-2"
    formula = "CF3CF2OCH3"
    synonym = "HFE-245cb2"
    rhoc = unidades.Density(499.507581544)
    Tc = unidades.Temperature(406.813)
    Pc = unidades.Pressure(2886.4, "kPa")
    M = 150.047336  # g/mol
    Tt = unidades.Temperature(250)
    Tb = unidades.Temperature(278.76)
    f_acent = 0.354
    momentoDipolar = unidades.DipoleMoment(2.785, "Debye")
    id = 671
    # id = 1817

    CP1 = {"ao": 10.196438,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [10.214789, 10.503071, 0.98682562, 0],
           "hyp": [814, 2031, 3040, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for RE245cb2 of Zhou et al. (2010)",
        "__doi__": {"autor": "Zhou, Y. and Lemmon, E.W.",
                    "title": "preliminary equation, 2010.",
                    "ref": "",
                    "doi":  ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 400000.0, "rhomax": 10.02,
        "Pmin": 8.272, "rhomin": 10.,

        "nr1": [0.41453162e-1, 0.15010352e1, -0.23142144e1, -0.471412, 0.17182],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.25, 0.786, 1.32, 0.338],

        "nr2": [-0.98793, -0.392049, 0.6848583, -0.32413816, -0.2414796e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.82, 2., 1., 3., 0.766],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.82792487, -0.31833343, -0.11929747, -0.65010212],
        "d3": [1, 1, 3, 3],
        "t3": [1.75, 3.5, 3.86, 2.75],
        "alfa3": [1.023, 1.384, 0.998, 6.9],
        "beta3": [1.727, 1.543, 1.075, 88],
        "gamma3": [1.1, 0.64, 0.5, 1.26],
        "epsilon3": [0.713, 0.917, 0.69, 0.743]}

    eq = helmholtz1,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.8026, 1.8804, -2.8375, -4.3077],
        "exp": [1, 1.5, 2.5, 5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [1.8378, 2.5311, -7.084, 18.678, -30.228, 22.985],
        "exp": [0.32, 1.08, 1.9, 2.8, 3.8, 4.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-1.5224, -5.7245, -15.972, -50.473, -6.8916],
        "exp": [0.286, 0.82, 2.5, 5.6, 7.3]}

    visco0 = {"eq": 5, "omega": 3,
              "__doi__": {"autor": "T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E",
                          "title": "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties",
                          "ref": "Ind. Eng. Chem. Res., 1988, 27 (4), pp 671–679",
                          "doi": "10.1021/ie00076a024"},
              "__name__": "Chung (1988)",
              "w": 0.354, "mur": 0.0, "k": 0.0}

    _viscosity = visco0,
