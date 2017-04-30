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
from lib.mEoS.R113 import R113


class C4F10(MEoS):
    """Multiparameter equation of state for perfluorobutane"""
    name = "perfluorobutane"
    CASNumber = "355-25-9"
    formula = "C4F10"
    synonym = ""
    rhoc = unidades.Density(599.8356)
    Tc = unidades.Temperature(386.326)
    Pc = unidades.Pressure(2323.4, "kPa")
    M = 238.03  # g/mol
    Tt = unidades.Temperature(145.0)
    Tb = unidades.Temperature(271.061)
    f_acent = 0.371
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 693

    CP1 = {"ao": 2.0150709,
           "R": 8.31451,
           "an": [0.96863193e-1, -0.99797537e-4, 0.37348060e-7],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    ecs = {
        "__type__": "ECS",
        "__name__": "Thermodynamic Extended Corresponding States model w/ T- and rho-dependent shape factors.",
        "__doc__":  u"""Huber, M.L. and Ely, J.F., "A predictive extended corresponding states model for pure and mixed refrigerants including an equation of state for R134a," Int. J. Refrigeration, 17:18-31, 1994.""",
        "cp": CP1,
        "ref": R113,
        "eq": "helmholtz1",
        "R": 8.314471,

        "Tmin": 189., "Tmax": 500., "Pmax": 30000.0, "rhomax": 7.64,
#        "Pmin": 0.61166, "rhomin": 55.497,

        "ft": [0.776042865e-2, -0.641975631],
        "ft_add": [], "ft_add_exp": [],
        "fd": [], "fd_exp": [],
        "ht": [0.278313281e-2, -0.593657910],
        "ht_add": [], "ht_add_exp": [],
        "hd": [-0.236093735e-2], "hd_exp": [1]}

    eq = ecs,

    _surface = {"sigma": [0.04429], "exp": [1.242]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.72217e1, -0.18886e2, 0.47288e2, -0.29794e2, -0.50457e1],
        "exp": [1.0, 1.5, 1.65, 1.8, 4.8]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.36787e1, -0.20581e1, 0.98945, 0.60478e-1],
        "exp": [0.4, 0.6, 0.8, 2.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.42967e1, -0.10715e2, -0.33457e2, -0.72206e2],
        "exp": [0.486, 1.7, 4.2, 8.0]}

    trnECS = {"eq": "ecs",
              "__name__": "Extended Corresponding States model",
              "__doc__": """Huber, M.L., Laesecke, A., and Perkins, R.A., Model for the viscosity and thermal conductivity of refrigerants, including a new correlation for the viscosity of R134a, Ind.Eng.Chem.Res. 42: 3163-3178 (2003).""",

              "ref": R113,
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

#if __name__ == "__main__":
#    import doctest
#    doctest.testmod()
#
#    cyc5=C4F10(T=300., P=0.1)
#    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
#
