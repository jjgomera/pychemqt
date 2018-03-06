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


class C1Oleate(MEoS):
    """Multiparameter equation of state for methyl oleate"""
    name = "methyl oleate"
    CASNumber = "112-62-9"
    formula = "C19H36O2"
    synonym = ""
    rhoc = unidades.Density(241.000222029)
    Tc = unidades.Temperature(782.0)
    Pc = unidades.Pressure(1246.0, "kPa")
    M = 296.48794   # g/mol
    Tt = unidades.Temperature(253.47)
    Tb = unidades.Temperature(627.18)
    f_acent = 0.91
    momentoDipolar = unidades.DipoleMoment(1.63, "Debye")
    # id = 919

    CP1 = {"ao": 0.0,
           "an": [90.2385], "pow": [0.146118],
           "ao_exp": [234.797, 335.768, 431.66],
           "exp": [613.529, 1405.31, 2867.76],
           "ao_hyp": [], "hyp": []}

    huber = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methyl oleate of Huber "
                    "et al. (2009).",
        "__doi__": {"autor": "Huber, M.L., Lemmon, E.W., Kazakov, A., Ott, "
                             "L.S., and Bruno, T.J.",
                    "title": "Model for the Thermodynamic Properties of a "
                             "Biodiesel Fuel",
                    "ref": "Energy Fuels, 23 (7) (2009) 3790–3797",
                    "doi": "10.1021/ef900159g"},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 50000.0, "rhomax": 3.05,
        "Pmin": 0.0000000004, "rhomin": 3.05,

        "nr1": [0.4596121e-1, 2.2954, -3.554366, -0.2291674, 0.6854534e-1],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.34, 1.14, 1.4, 0.6],

        "nr2": [-1.535778, -0.7334697, 1.712700, -1.471394, -0.1724678e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [3.3, 4.1, 1.9, 3.8, 1.3],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1]*5,

        "nr3": [0.2115470e1, -0.7555374, -0.4134269],
        "d3": [1, 1, 3],
        "t3": [3.4, 3.8, 4.0],
        "alfa3": [1.1, 1.6, 1.1],
        "beta3": [0.9, 0.65, 0.75],
        "gamma3": [1.14, 0.65, 0.77],
        "epsilon3": [0.79, 0.9, 0.76],
        "nr4": []}

    eq = huber,

    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.13900e2, 0.16246e2, -0.15568e2, -0.73568e1, -0.48739e1],
        "exp": [1.0, 1.5, 1.93, 4.2, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.19920e2, 0.12230e3, -0.23582e3, 0.21009e3, -0.73435e2],
        "exp": [0.461, 0.6, 0.75, 0.91, 1.05]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-13.426, 1.8069e2, -1.1108e3, 1.3265e3, -4.6421e2, -2.1070e2],
        "exp": [0.667, 1.71, 2.2, 2.46, 3.0, 9.7]}

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2010)",
               "__doi__": {"autor": "Perkins, R.A. and Huber, M.L.",
                           "title": "Measurement and Correlation of the Thermal Conductivities of Biodiesel Constituent Fluids: Methyl Oleate and Methyl Linoleate",
                           "ref": "Energy Fuels, 2011, 25 (5), pp 2383–2388",
                           "doi": "10.1021/ef200417x"},
               "__test__": """
                    >>> st=C1Oleate(T=450, P=1e2)
                    >>> print "%0.0f %0.4f %0.6g %0.6g" % (st.T, st.P.MPa, st.rho, st.k.WmK)
                    450 0.0001 0.00792667 0.0110996
                    >>> st=C1Oleate(T=450, P=1e6)
                    >>> print "%0.0f %0.0f %0.6g %0.6g" % (st.T, st.P.MPa, st.rho, st.k.WmK)
                    450 1 764.716 0.123794
                    >>> st=C1Oleate(T=450, P=2e7)
                    >>> print "%0.0f %0.0f %0.6g %0.6g" % (st.T, st.P.MPa, st.rho, st.k.WmK)
                    450 20 787.08 0.133856
                    """, # Table 3, Pag 2386

               "Tref": 782.0, "kref": 1,
               "no": [-0.27125000e-3, 0.25936500e-2, 0.35024100e-1, -0.90227300e-2],
               "co": [0, 1, 2, 3],

               "Trefb": 782.0, "rhorefb": 0.81285*M, "krefb": 1.,
               "nb": [-0.41010600e-1, 0.32844300e-1, -0.41850600e-2, 0.0, 0.0,
                      0.60665700e-1, -0.49840700e-1, 0.12175200e-1, 0.0, 0.0],
               "tb": [0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
               "db": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 8.75e-10, "Tcref": 1173.}

    thermo1 = {"eq": 5, "omega": 3,
               "__name__": "Chung (1988)",
               "__doi__": {"autor": "T-H. Chung, Ajlan, M., Lee, L.L. and Starling, K.E.",
                           "title": "Generalized Multiparameter Correlation for Nonpolar and Polar Fluid Transport Properties",
                           "ref": "Ind. Eng. Chem. Res., 1988, 27 (4), pp 671–679",
                           "doi": "10.1021/ie00076a024"},
               "w": 0.91, "mur": 0.0, "k": 0.0}

    _thermal = thermo0, thermo1
