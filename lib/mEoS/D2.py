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


class D2(MEoS):
    """Multiparameter equation of state for deuterium"""
    name = "deuterium"
    CASNumber = "7782-39-0"
    formula = "D2"
    synonym = ""
    rhoc = unidades.Density(69.405886)
    Tc = unidades.Temperature(38.34)
    Pc = unidades.Pressure(1679.6, "kPa")
    M = 4.0282  # g/mol
    Tt = unidades.Temperature(18.724)
    Tb = unidades.Temperature(23.661)
    f_acent = -0.136
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [-2.0677351753, 2.4237151502],
           "ao_exp": [-3.54145, 3.0326, -3.52422, -1.73421, -3.57135, 2.14858,
                      6.23107, -3.30425, 6.23098, -3.57137, 3.32901, 0.97782],
           "titao": [7174.1/Tc, 8635/Tc, 902.7/Tc, 181.1/Tc, 438.5/Tc,
                     5034.2/Tc, 269.9/Tc, 229.9/Tc, 666.4/Tc, 452.8/Tc,
                     192.0/Tc, 1187.6/Tc]}

    CP1 = {"ao": 2.4512991,
           "an": [0.43563077e-02, -0.53169470e-03, 0.17067184e-04,
                  -0.53819932e-08, 0.89310438e-12],
           "pow": [1, 1.5, 2, 3, 4],
           "ao_exp": [0.18403263e2, -0.21257617e2, 0.41091635e1],
           "exp": [319, 361, 518],
           "ao_hyp": [], "hyp": []}

    richardson = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for deuterium of Richardson et al. (2014).",
        "__doi__": {"autor": "Richardson, I.A., Leachman, J.W., and Lemmon, E.W.",
                    "title": "Fundamental Equation of State for Deuterium ",
                    "ref": "J. Phys. Chem. Ref. Data 43, 013103 (2014)",
                    "doi": "10.1063/1.4864752"},
        "__test__":
            # Table 7, Pag 12
            """
            >>> st=D2(T=18.724, x=0.5)
            >>> print "%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, \
                st.Gas.cv.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            18.724 17.189 174.630 0.455 -30.450 286.160 -1.415 15.494 3.355 3.143 5.627 5.364 1085.50 250.95
            >>> st=D2(T=20, x=0.5)
            >>> print "%0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, \
                st.Gas.cv.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            20 29.425 171.660 0.737 -23.077 291.490 -1.038 14.690 3.371 3.167 5.852 5.464 1058.60 257.98
            >>> st=D2(T=25, x=0.5)
            >>> print "%0.0f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, \
                st.Gas.cv.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            25 146.40 158.870 3.144 9.282 307.750 0.370 12.309 3.442 3.318 6.996 6.197 937.27 278.61
            >>> st=D2(T=30, x=0.5)
            >>> print "%0.0f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, \
                st.Gas.cv.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            30 445.75 143.160 9.055 49.449 313.220 1.757 10.550 3.576 3.568 9.069 8.013 777.76 288.11
            >>> st=D2(T=35, x=0.5)
            >>> print "%0.0f %0.1f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, \
                st.Gas.cv.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            35 1036.7 120.480 23.053 105.180 297.410 3.329 8.821 3.846 3.954 16.103 15.430 561.44 286.83
            >>> st=D2(T=38, x=0.5)
            >>> print "%0.0f %0.1f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, \
                st.Gas.h.kJkg, st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, \
                st.Gas.cv.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            38 1600.6 88.550 50.659 172.270 244.900 5.009 6.920 4.187 4.300 132.060 137.920 374.12 297.24
            """,

        "R": 8.3144621,
        "cp": Fi1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 2000000.0, "rhomax": 43.351,
        "Pmin": 17.189, "rhomin": 43.351,

        "nr1": [0.006267958, 10.53609, -10.14149, 0.356061, 0.1824472,
                -1.129638, -0.0549812, -0.6791329],
        "d1": [4, 1, 1, 2, 3, 1, 3, 2],
        "t1": [1, 0.462, 0.5584, 0.627, 1.201, 0.309, 1.314, 1.1166],

        "nr2": [1.347918, -0.8657582, 1.719146, -1.917977, 0.1233365,
                -0.07936891],
        "d2": [2, 2, 1, 1, 3, 2],
        "t2": [1.25, 1.25, 1.395, 1.627, 1.0, 2.5],
        "c2": [1, 1, 2, 2, 2, 2],
        "gamma2": [1]*6,

        "nr3": [1.686617, -4.240326, 1.857114, -0.5903705, 1.520171,
                2.361373, -2.297315],
        "d3": [1, 1, 2, 3, 3, 1, 3],
        "t3": [0.635, 0.664, 0.7082, 2.25, 1.524, 0.67, 0.709],
        "alfa3": [0.868, 0.636, 0.668, 0.65, 0.745, 0.782, 0.693],
        "beta3": [0.613, 0.584, 0.57, 1.056, 1.01, 1.025, 1.029],
        "gamma3": [0.6306, 0.711, 0.6446, 0.8226, 0.992, 1.2184, 1.203],
        "epsilon3": [1.46, 1.7864, 1.647, 0.541, 0.969, 1.892, 1.076],
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for deuterium of McCarty (1989)",
        "__doi__": {"autor": "McCarty, R.D.",
                    "title": "Correlations for the Thermophysical Properties of Deuterium",
                    "ref": "National Institute of Standards and Technology, Boulder, CO, 1989",
                    "doi": ""},

        "R": 8.31434,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 423.0, "Pmax": 320000.0, "rhomax": 43.38,
        "Pmin": 19.462, "rhomin": 43.365,

        "b": [None, 0.4894244053982e-4, 0.5600164604601e-1, -0.6301493491211,
              0.2538329946038e1, 0.1723475985309e3, 0.2956238369436e-4,
              -0.3926317169317e-2, 0.1195764193293e-1, 0.1136916678824e5,
              -0.1916378195727e-6, 0.3153535946452e-3, 0.2122937335070e-1,
              -0.1057999371607e-5, -0.6722062598854e-4, -0.3030166828627,
              0.1980817195099e-5, -0.1453922641871e-7, 0.1780919116891e-3,
              -0.1823145348424e-5, -0.1135358616578e5, -0.1943542941899e4,
              -0.3632847669580e2, 0.1087745118380e3, -0.4078276062687e-1,
              0.6460021864005e-2, -0.4480242189217e-4, -0.2475011206216e-3,
              -0.8834384656760e-8, -0.1081622159862e-8, -0.1478159334303e-10,
              0.7926922356112e-11, 0.5721547329378e-11]}

    eq = richardson, MBWR

    _surface = {"sigma": [0.009376], "exp": [1.258]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-5.5706, 1.7631, -0.5458, 1.2154, -1.1556],
        "exp": [1., 1.5, 2.83, 4.06, 5.4]}
    _liquid_Density = {
        "eq": 1,
        "ao": [3.3769, -5.3693, 11.943, -17.361, 15.170, -6.3079],
        "exp": [0.512, 1.12, 1.8, 2.55, 3.4, 4.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.8111, -7.3624, 2.2294, -21.443, 12.796, -31.334],
        "exp": [0.528, 2.03, 3.6, 5.0, 6.5, 9.0]}
