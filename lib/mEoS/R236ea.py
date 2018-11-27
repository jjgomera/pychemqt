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
from lib.mEoS.R134a import R134a


class R236ea(MEoS):
    """Multiparameter equation of state for R236ea"""
    name = "1,1,1,2,3,3-hexafluoropropane"
    CASNumber = "431-63-0"
    formula = "CF3CHFCHF2"
    synonym = "R236ea"
    _refPropName = "R236EA"
    _coolPropName = "R236EA"
    rhoc = unidades.Density(565.)
    Tc = unidades.Temperature(412.44)
    Pc = unidades.Pressure(3420.0, "kPa")
    M = 152.0384  # g/mol
    Tt = unidades.Temperature(170.0)
    Tb = unidades.Temperature(279.322)
    f_acent = 0.369
    momentoDipolar = unidades.DipoleMoment(1.129, "Debye")
    # id = 1873

    Fi1 = {"ao_log": [1, 2.762],
           "pow": [0, 1],
           "ao_pow": [-14.121424135, 10.2355589225],
           "ao_exp": [0.7762, 10.41, 12.18, 3.332],
           "titao": [144/Tc, 385/Tc, 1536/Tc, 7121/Tc]}

    CP1 = {"ao": 5.30694,
           "an": [0.03973, -1.859e-5], "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    rui = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R236ea of Rui (2013)",
        "__doi__": {"autor": "Rui, X., Pan, J., Wang, Y.",
                    "title": "An Equation of State for Thermodynamic "
                             "Properties of 1,1,1,2,3,3-Hexafluoropropane "
                             "(R236ea)",
                    "ref": "Fluid Phase Equilibria 341 (2013) 78-85",
                    "doi": "10.1016/j.fluid.2012.12.026"},

        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1.,
                "ho": 56317.4970978844, "so": 282.8465334259},

        "Tmin": 240.0, "Tmax": 412.0, "Pmax": 6000.0, "rhomax": 10.5,
        "Pmin": 0.02, "rhomin": 11.7,

        "nr1": [0.051074, 2.5584, -2.9180, -0.71485, 0.15534],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1., 0.264, 0.5638, 1.306, 0.2062],

        "nr2": [-1.5894, -0.784, 0.85767, -0.67235, -0.017953],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.207, 2.283, 1.373, 2.33, 0.6376],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.3165, -0.42023, -0.28053, -1.4134, -0.0000062617],
        "d3": [1, 1, 3, 3, 2],
        "t3": [1.08, 1.67, 3.502, 4.357, 0.6945],
        "alfa3": [1.019, 1.341, 1.034, 5.264, 24.44],
        "beta3": [1.3, 2.479, 1.068, 79.85, 49.06],
        "gamma3": [1.13, 0.6691, 0.465, 1.28, 0.8781],
        "epsilon3": [0.7119, 0.9102, 0.678, 0.7091, 1.727]}

    ecs = {"__type__": "ECS",
           "__name__": "Thermodynamic Extended Corresponding States model w/ T- and rho-dependent shape factors.",
           "__doc__":  u"""Huber, M.L. and Ely, J.F., "A predictive extended corresponding states model for pure and mixed refrigerants including an equation of state for R134a," Int. J. Refrigeration, 17:18-31, 1994.""",
           "cp": CP1,
           "ref": R134a,
           "eq": "helmholtz1",
           # "eq": "MBWR",
           "R": 8.314471,

            "Tmin": 242.0, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 10.465,
#            "Pmin": aaaaaaa, "rhomin": aaaaaaa,

           "ft": [-0.67786992, -0.52182651],
           "ft_add": [], "ft_add_exp": [],
           "fd": [0.113833347e-1], "fd_exp": [1],
           "ht": [0.142369159e1, 0.870214752e-1],
           "ht_add": [0.195298641e-1], "ht_add_exp": [1],
           "hd": [], "hd_exp": []}

    # eq = rui, ecs
    eq = rui,

    _surface = {"sigma": [0.306974, -0.247277], "exp": [1.12614, 1.09899]}
    _vapor_Pressure = {
        "eq": 3,
        "n": [-7.9095, 2.3374, -2.6453, -5.7058],
        "t": [1, 1.5, 2.15, 4.75]}
    _liquid_Density = {
        "eq": 1,
        "n": [1.6074, 1.5021, -1.106, 0.91146],
        "t": [0.31, 0.75, 1.3, 1.9]}
    _vapor_Density = {
        "eq": 2,
        "n": [-2.7426, -6.2268, -15.109, -49.524, -114.11],
        "t": [0.376, 1.1, 2.7, 5.5, 11]}

    trnECS = {"eq": "ecs",
              "__name__": "Extended Corresponding States model",
              "__doc__": """Huber, M.L., Laesecke, A., and Perkins, R.A., Model for the viscosity and thermal conductivity of refrigerants, including a new correlation for the viscosity of R134a, Ind.Eng.Chem.Res. 42: 3163-3178 (2003).""",

              "ref": R134a,
              "ref_eq": "helmholtz1",
              "eq_visco": "visco0",
              "eq_thermo": "thermo0",

              "f_int": [1.32e-3],
              "psi": [1.0],
              "phi": [1.0],

              "critical": 3,
              "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
              "Xio": 0.194e-9, "gam0": 0.0496, "qd": 1.5e-9, "Tcref": 579.49}

#    _viscosity=trnECS,
#    _thermal=trnECS,
