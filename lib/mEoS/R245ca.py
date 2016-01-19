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
from .R134a import R134a


class R245ca(MEoS):
    """Multiparameter equation of state for R245ca"""
    name = "1,1,2,2,3-pentafluoropropane"
    CASNumber = "679-86-7"
    formula = "CHF2CF2CH2F"
    synonym = "R245ca"
    rhoc = unidades.Density(525.4679248)
    Tc = unidades.Temperature(447.57)
    Pc = unidades.Pressure(3940.7, "kPa")
    M = 134.04794  # g/mol
    Tt = unidades.Temperature(191.5)
    Tb = unidades.Temperature(298.412)
    f_acent = 0.355
    momentoDipolar = unidades.DipoleMoment(1.740, "Debye")
    id = 693

    CP1 = {"ao": 8.888,
           "an": [], "pow": [],
           "ao_exp": [0.8843, 5.331, 14.46],
           "exp": [865, 2830, 1122],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": -3.8444,
           "an": [5.24008e-1, -3.74976e-4], "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-245ca of Zhou and Lemmon (2013).",
        "__doi__": {"autor": "Zhou, Y. and Lemmon, E.W.",
                    "title": "unpublished equation, 2013.",
                    "ref": "",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 10000.0, "rhomax": 12.21,
        "Pmin": 0.0708, "rhomin": 12.21,

        "nr1": [0.04489247, 1.526476, -2.408320, -0.5288088, 0.18222346],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.26, 1., 1.2, 0.67],

        "nr2": [-1.063228, -0.223149, 1.18738, -0.9772383, -0.02296938],
        "d2": [1, 3, 2, 2, 7],
        "t2": [1.92, 2., 1.5, 1.93, 1.06],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [1.364444, -0.5080666, -0.06649496, -1.128359],
        "d3": [1, 1, 3, 3],
        "t3": [0.17, 3.9, 1., 1.],
        "alfa3": [1.16, 1.1, 1.64, 13.8],
        "beta3": [2.4, 1.5, 4.2, 379],
        "gamma3": [1.265, 0.42, 0.864, 1.15],
        "epsilon3": [0.55, 0.724, 0.524, 0.857]}

    ecs = {"__type__": "ECS",
           "__name__": "Thermodynamic Extended Corresponding States model w/ T- and rho-dependent shape factors.",
           "__doc__":  u"""Huber, M.L. and Ely, J.F., "A predictive extended corresponding states model for pure and mixed refrigerants including an equation of state for R134a," Int. J. Refrigeration, 17:18-31, 1994.""",
           "cp": CP2,
           "ref": R134a,
           "eq": "helmholtz1",
           "R": 8.314471,

            "Tmin": 200.0, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 11.995,
#            "Pmin": aaaaaaa, "rhomin": aaaaaaa,

           "ft": [-0.241011472, -0.788477331],
           "ft_add": [], "ft_add_exp": [],
           "fd": [], "fd_exp": [],
           "ht": [0.160567866e1, -0.727455038],
           "ht_add": [], "ht_add_exp": [],
           "hd": [], "hd_exp": []}

    eq = helmholtz1,

    _surface = {"sigma": [0.069297, -0.022419], "exp": [1.2795, 3.1368]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.8807, 2.1026, -3.0768, -4.9894],
        "exp": [1.0, 1.5, 2.5, 4.95]}
    _liquid_Density = {
        "eq": 1,
        "ao": [4.0176, -4.7916, 7.8662, -7.1049, 3.1949],
        "exp": [0.48, 1.0, 1.62, 2.3, 3.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-4.65885, -1.03328, -13.5047, -48.4225, -104.097],
        "exp": [0.5, 1.09, 2.1, 5.1, 10.4]}

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

#if __name__ == "__main__":
#    import doctest
#    doctest.testmod()
#
#    cyc5=R245ca(T=300., P=0.1)
#    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
