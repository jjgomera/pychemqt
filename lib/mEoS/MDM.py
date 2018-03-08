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


class MDM(MEoS):
    """Multiparamenter equation of state for octamethyltrisiloxane"""
    name = "octamethyltrisiloxane"
    CASNumber = "107-51-7"
    formula = "C8H24O2Si3"
    synonym = "MDM"
    rhoc = unidades.Density(256.73940949935815)
    Tc = unidades.Temperature(564.09)
    Pc = unidades.Pressure(1415.0, "kPa")
    M = 236.531  # g/mol
    Tt = unidades.Temperature(187.2)
    Tb = unidades.Temperature(425.66)
    f_acent = 0.529
    momentoDipolar = unidades.DipoleMoment(1.079, "Debye")
    # id = 1893

    CP1 = {"ao": 275.1,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [266040871.9, 0, 2051643622.0, 0],
           "hyp": [802.6, 0, 1829.6, 0]}

    colonna = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for MDM of Colonna (2006).",
        "__doi__": {"autor": "Colonna, P., Nannan, N.R., and Guardone, A.",
                    "title": "Multiparameter equations of state for siloxanes:"
                             " [(CH3)3-Si-O1/2]2-[O-Si-(CH3)2]i=1,…,3, and "
                             "[O-Si-(CH3)2]6",
                    "ref": "Fluid Phase Equilibria 263:115-130, 2008",
                    "doi":  "10.1016/j.fluid.2007.10.001"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 673.0, "Pmax": 30000.0, "rhomax": 3.94,
        "Pmin": 0.0000008, "rhomin": 3.93,

        "nr1": [1.19735372, -2.40380622, 0.3256564, -0.19971259, 0.11206277,
                0.15893999e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.51234323, -0.20660361e-1, -0.38978114, -0.1186931,
                -0.37203537e-1, 0.18359984e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.0],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = colonna,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.85589e1, 0.20278e1, -0.28501e1, -0.64397e1, -0.85460e1],
        "exp": [1.0, 1.5, 2.3, 4.0, 13.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.54145, -0.27650e-1, 0.41558e1, -0.19104e1, 0.67606],
        "exp": [0.12, 0.36, 0.6, 0.8, 2.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.16483e1, -0.71410e1, -0.23088e2, -0.70554e2, 0.19938e1,
               -0.20193e3],
        "exp": [0.296, 0.905, 2.8, 5.9, 12.0, 13.0]}

    visco0 = {"eq": 5, "omega": 3,
              "__doi__": {"autor": "T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E",
                          "title": "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties",
                          "ref": "Ind. Eng. Chem. Res., 1988, 27 (4), pp 671–679",
                          "doi": "10.1021/ie00076a024"},
              "__name__": "Chung (1988)",
              "w": 0.531, "mur": 0.0, "k": 0.0}

    _viscosity = visco0,
