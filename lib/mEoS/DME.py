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


class DME(MEoS):
    """Multiparameter equation of state for dimethylether"""
    name = "dimethylether"
    CASNumber = "115-10-6"
    formula = "CH3-O-CH3"
    synonym = "R-170"
    rhoc = unidades.Density(273.6465336)
    Tc = unidades.Temperature(400.378)
    Pc = unidades.Pressure(5336.8, "kPa")
    M = 46.06844  # g/mol
    Tt = unidades.Temperature(131.66)
    Tb = unidades.Temperature(248.368)
    f_acent = 0.196
    momentoDipolar = unidades.DipoleMoment(1.301, "Debye")
    id = 133

    Fi1 = {"ao_log": [1, 3.039],
           "pow": [0, 1],
           "ao_pow": [-1.928925, 3.150284],
           "ao_exp": [2.641, 2.123, 8.992, 6.191],
           "titao": [361/Tc, 974/Tc, 1916/Tc, 4150/Tc],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for DME of Wu et al. (2011).",
        "__doi__": {"autor": "Wu, J., Zhou, Y., and Lemmon, E.W.",
                    "title": "An Equation of State for the Thermodynamic Properties of Dimethyl Ether",
                    "ref": "J. Phys. Chem. Ref. Data 40, 023104 (2011)",
                    "doi":  "10.1063/1.3582533"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1.0, "ho": 23242, "so": 131.3883},

        "Tmin": Tt, "Tmax": 525.0, "Pmax": 40000.0, "rhomax": 19.15,
        "Pmin": 0.0022, "rhomin": 19.15,

        "nr1":  [0.29814139e-1, 0.14351700e1, -0.26496400e1, -0.29515532,
                 0.17035607],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1.0, 0.4366, 1.011, 1.137, 0.45],

        "nr2": [-0.94642918, -0.99250514e-1, 0.11264071e1, -0.76936548,
                -0.20717696e-1, 0.24527037],
        "d2": [1, 3, 2, 2, 7, 1],
        "t2": [2.83, 1.5, 1.235, 2.675, 0.7272, 1.816],
        "c2": [2, 2, 1, 2, 1, 1],
        "gamma2": [1]*6,

        "nr3": [0.11863438e1, -0.49398368, -0.16388716, -0.27583584e-1],
        "d3": [1, 1, 3, 3],
        "t3": [1.783, 3.779, 3.282, 1.059],
        "alfa3": [0.965336, 1.508580, 0.963855, 9.726430],
        "beta3": [1.287190, 0.806235, 0.777942, 197.681000],
        "gamma3": [1.277720, 0.430750, 0.429607, 1.138490],
        "epsilon3": [0.672698, 0.924246, 0.750815, 0.800022],
        "nr4": []}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for DME of Ihmels and Lemmon (2007)",
        "__doi__": {"autor": "Ihmels, E.C. and Lemmon, E.W.",
                    "title": "Experimental densities, vapor pressures, and critical point, and a fundamental equation of state for dimethyl ether",
                    "ref": "Fluid Phase Equilibria. 10/2007; 260(1):36-48",
                    "doi":  "10.1016/j.fluid.2006.09.016"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1.0, "ho": 23242, "so": 131.3883},

        "Tmin": Tt, "Tmax": 1350.0, "Pmax": 1000000.0, "rhomax": 73.96,
        "Pmin": 0.61166, "rhomin": 55.497,

        "nr1":  [1.22690, -2.47245, 0.119889, 0.0000354],
        "d1": [1, 1, 3, 8],
        "t1": [0.21, 1.0, 0.5, 1.0],

        "nr2": [0.567139, 0.166649, -0.078412, -0.289066, -0.031272, -0.065607],
        "d2": [2, 1, 5, 1, 4, 3],
        "t2": [1.4, 3.1, 1.5, 5.0, 5.9, 3.7],
        "c2": [1, 1, 1, 2, 2, 2],
        "gamma2": [1]*6}

    eq = helmholtz1, helmholtz2

    _surface = {"sigma": [0.063157], "exp": [1.2595]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.112782, 1.971239, -2.276083, -2.215774],
        "exp": [1, 1.5, 2.5, 5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [7.884834, -10.516328, 5.39142, 0.40489],
        "exp": [0.54, 0.74, 0.95, 11.43]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-4.136444, -4.302025, -12.03214, -39.527936, -89.4768],
        "exp": [1.467, 4.2, 8.0, 17.0, 36.0]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Meng (2012)",
              "__doi__": {"autor": "Meng, X., Zhang, J., Wu, J., Liu, Z.",
                          "title": "Experimental Measurement and Modeling of the Viscosity of Dimethyl Ether",
                          "ref": "J. Chem. Eng. Data, 2012, 57 (3), pp 988–993",
                          "doi":  "10.1021/je201297j"},

              "ek": 317.937, "sigma": 0.446704,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.14508011,
              "n_ideal": [0.294261, -0.377826, -0.491673],
              "t_ideal": [0, 1, 2],

              "Tref_res": 400.378, "rhoref_res": 5.94*M, "etaref_res": 1,
              "n_packed": [],
              "t_packed": [],
              "n_poly": [-2.70002, 4.44583, -104.998, 78.27474, 41.3751,
                         -175.055, 62.81975, 0.21302, 112.3219, 6.50681],
              "t_poly": [-5.92, -4.36, -2.93, -1.64, -7.86, -4.25, -4.79,
                         -5.87, -3.11, -0.45],
              "d_poly": [3, 3, 3, 4, 5, 2, 2, 5, 2, 1],
              "g_poly": [0]*10,
              "c_poly": [0, 0, 1, 1, 2, 1, 1, 0, 2, 0]}

    _viscosity = visco0,
