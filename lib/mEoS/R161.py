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


from unittest import TestCase

from lib.meos import MEoS
from lib import unidades


class R161(MEoS):
    """Multiparameter equation of state for R161"""
    name = "fluoroethane"
    CASNumber = "353-36-6"
    formula = "C2H5F"
    synonym = "R161"
    _refPropName = "R161"
    _coolPropName = "R161"
    rhoc = unidades.Density(302.0010921)
    Tc = unidades.Temperature(375.25)
    Pc = unidades.Pressure(5046.0, "kPa")
    M = 48.0595  # g/mol
    Tt = unidades.Temperature(130.0)
    Tb = unidades.Temperature(235.6)
    f_acent = 0.216
    momentoDipolar = unidades.DipoleMoment(1.9397, "Debye")
    id = 247

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-6.9170707460, 5.4837900434],
           "ao_exp": [1.08888, 1.80842, 8.72417, 5.67715],
           "titao": [329/Tc, 742/Tc, 1644/Tc, 3922/Tc]}

    Fi2 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-6.9187, 5.4788],
           "ao_exp": [2.059, 9.253, 6.088],
           "titao": [420/Tc, 1548/Tc, 3882/Tc]}

    CP1 = {"ao": 3.985,
           "an": [], "pow": [],
           "ao_exp": [2.077, 9.265, 6.054], "exp": [420, 1548, 3882],
           "ao_hyp": [], "hyp": []}

    qi = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-161 of Qi (2012)",
        "__doi__": {"autor": "Qi, H., Fang, D., Gao, K., Meng, X., Wu, J.",
                    "title": "Compressed Liquid Densities and Helmholtz Energy"
                             " Equation of State for Fluoroethane (R161)",
                    "ref": "Int. J. Thermophys. 37(3) (2016) 55",
                    "doi": "10.1007/s10765-016-2061-1"},

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "IIR",

        # "rhoc": 6.2839,

        "Tmin": Tt, "Tmax": 420.0, "Pmax": 100000.0, "rhomax": 18,
        # "Pmin": 0.005512, "rhomin": 19.91,

        "nr1": [0.005145283, -0.001882274, 1.884722, -3.1819965, -0.24432415,
                0.27792467],
        "d1": [5, 4, 1, 1, 2, 3],
        "t1": [1, 0.68, 0.32, 0.92, 1.23, 0.846],

        "nr2": [-0.4414064, -0.402065, 0.24171113, -0.16603585, -0.03440867,
                -0.000099185],
        "d2": [1, 3, 2, 2, 7, 5],
        "t2": [4.208, 3.06, 1.85, 4.28, 1.003, 1.12],
        "c2": [2, 2, 1, 2, 1, 1],
        "gamma2": [1]*6,

        "nr3": [1.0146668, -0.03542609, -0.006038245, -0.025437558,
                -0.00515678, 0.006396804],
        "d3": [1, 1, 3, 3, 2, 2],
        "t3": [1.055, 0.8, 4.08, 1.6, 3.85, 0.57],
        "alfa3": [0.96212, 3.2147, 2.6288, 0.8657, 2.3839, 1.7814],
        "beta3": [0.62848, 4.5968, 4.9696, 0.239, 0.788, 7.0874],
        "gamma3": [1.9363, 1.5054, 1.3691, 2.3594, 0.5581, 0.6326],
        "epsilon3": [0.70192, 1.23824, 0.73324, 0.6258, 1.564, 1.4861]}

    wu = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-161 of Wu and "
                    "Zhou (2012)",
        "__doi__": {"autor": "Wu, J., Zhou, Y.",
                    "title": "An Equation of State for Fluoroethane (R161)",
                    "ref": "Int. J. Thermophys. 33(2) (2012) 220-234",
                    "doi": "10.1007/s10765-011-1151-3"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 28559.6, "so": 167.205},

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 5000.0, "rhomax": 20.0,
        "Pmin": 0.005512, "rhomin": 19.91,

        "nr1": [1.511, -2.3, -0.457, 0.1683, 0.04133],
        "d1": [1, 1, 2, 3, 4],
        "t1": [0.37, 0.97, 1.14, 0.744, 1.],

        "nr2": [0.62187, -0.0265, -1.03, -0.285, -0.476],
        "d2": [2, 7, 1, 2, 3],
        "t2": [1.26, 1., 1.8, 3., 2.25],
        "c2": [1, 1, 2, 2, 2],
        "gamma2": [1]*5,

        "nr3": [0.82, -0.3532, -0.116, -0.0220583, -1.63148],
        "d3": [11, 1, 3, 3, 3],
        "t3": [1, 1.2, 5.3, 1, 4],
        "alfa3": [0.96, 1.35, 1.26, 1.23, 16.8],
        "beta3": [2.7, 5.2, 3.9, 4.7, 413],
        "gamma3": [0.9, 0.69, 0.67, 0.67, 1.15],
        "epsilon3": [0.683, 0.892, 0.785, 1.33, 0.86]}

    refprop = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-161 of Lemmon "
                    "(2005).",
        "__doi__": {"autor": "Lemmon, E.W.",
                    "title": "preliminary equation, 2005.",
                    "ref": "",
                    "doi": ""},

        "R": 8.314472,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 50000.0, "rhomax": 20.0,
        "Pmin": 0.006, "rhomin": 19.95,

        "nr1": [0.75688, -1.4110, -0.63922, 0.055685, 0.00028395],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [.73357, .67596, .011369, -.56406, -.094362, -.1678, .00034215],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 2],
        "gamma2": [1]*7}

    eq = wu, refprop

    _surface = {"sigma": [0.05385], "exp": [1.111]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.75224e1, 0.29140e1, -0.30129e1, -0.44497e1, 0.24207e1],
        "exp": [1.0, 1.5, 2.3, 6.0, 7.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-.22587e2, .13424e3, -.2671e3, .3389e3, -.31059e3, .13009e3],
        "exp": [0.56, 0.7, 0.9, 1.2, 1.5, 1.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.62548e1, 0.10499e2, -0.20353e2, -0.36709e2, -0.86781e2],
        "exp": [0.56, 1.3, 1.7, 5.0, 11.0]}
