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


class nC4(MEoS):
    """Multiparameter equation of state for n-butane"""
    name = "n-butane"
    CASNumber = "106-97-8"
    formula = "CH3-(CH2)2-CH3"
    synonym = "R-600"
    rhoc = unidades.Density(228.)
    Tc = unidades.Temperature(425.125)
    Pc = unidades.Pressure(3796.0, "kPa")
    M = 58.1222  # g/mol
    Tt = unidades.Temperature(134.895)
    Tb = unidades.Temperature(272.660)
    f_acent = 0.201
    momentoDipolar = unidades.DipoleMoment(0.05, "Debye")
    id = 6
    _Tr = unidades.Temperature(406.785141)
    _rhor = unidades.Density(230.384826)
    _w = 0.194240287

    Fi1 = {"ao_log": [1, 3.24680487],
           "pow": [0, 1],
           "ao_pow": [12.54882924, -5.46976878],
           "ao_exp": [5.54913289, 11.4648996, 7.59987584, 9.66033239],
           "titao": [0.7748404445, 3.3406025522, 4.9705130961, 9.9755537783],
           "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, 3.33944],
           "pow": [0, 1],
           "ao_pow": [20.884143364, -91.638478026],
           "ao_exp": [], "titao": [],
           "ao_hyp": [9.44893, 6.89406, 24.4618, 14.7824],
           "hyp": [1.101487798, 0.43195766, 4.502440459, 2.124516319]}

    Fi3 = {"ao_log": [1, 3.240207],
           "pow": [0, 1],
           "ao_pow": [-5.404217, 4.91136],
           "ao_exp": [5.513671, 7.388450, 10.250630, 11.061010],
           "titao": [327.55988/Tc, 1319.06935/Tc, 4138.63184/Tc, 1864.36783/Tc],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": -1.3491511376e1,
           "an": [3.8802310194e5, -1.5444296890e5, 2.8455082239e3,
                  6.6142595353e-2, -2.4307965028e-5, 1.5044248429e-10],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-8.3933423467], "exp": [3000],
           "ao_hyp": [], "hyp": []}

    CP6 = {"ao": 0.801601/8.3143*58.124,
           "an": [0.655936e-3/8.3143*58.124, 0.12277e-4/8.3143*58.124,
                  -0.165626e-7/8.3143*58.124, 0.67736e-11/8.3143*58.124],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Buecker and Wagner (2006)",
        "__doi__": {"autor": "Bücker, D., Wagner, W.",
                    "title": "Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Butane and Isobutane",
                    "ref": "J. Phys. Chem. Ref. Data 35, 929 (2006)",
                    "doi": "10.1063/1.1901687"},
        "__test__":
            # Table 44, Pag 974
            """
            >>> st=nC4(T=134.895, x=0.5)
            >>> print "%0.6g %0.8f %0.4f %0.6f %0.5g %0.5g %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            134.895 0.00000067 734.9588 0.000034 -721.64 -225.72 -3.036 0.640 1.441 0.963 1.973 1.106 1826.82 148.87
            >>> st=nC4(T=156, x=0.5)
            >>> print "%0.6g %0.6f %0.3f %0.5f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            156 0.000020 715.291 0.00091 -679.81 -201.60 -2.748 0.317 1.443 1.034 1.992 1.178 1693.20 159.37
            >>> st=nC4(T=170, x=0.5)
            >>> print "%0.6g %0.6f %0.3f %0.5f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            170 0.000116 702.237 0.00478 -651.81 -184.83 -2.576 0.171 1.447 1.078 2.009 1.221 1610.11 165.94
            >>> st=nC4(T=180, x=0.5)
            >>> print "%0.6g %0.6f %0.3f %0.5f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            180 0.000335 692.876 0.01301 -631.64 -172.50 -2.461 0.090 1.453 1.109 2.024 1.253 1552.13 170.43
            >>> st=nC4(T=300, x=0.5)
            >>> print "%0.6g %0.5f %0.2f %0.4f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            300 0.25760 570.68 6.5164 -367.83 -8.24 -1.349 -0.150 1.729 1.602 2.451 1.811 890.88 202.15
            >>> st=nC4(T=400, x=0.5)
            >>> print "%0.6g %0.4f %0.2f %0.3f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            400 2.4954 408.48 73.077 -80.60 113.39 -0.542 -0.057 2.173 2.210 3.838 3.623 318.35 154.77
            >>> st=nC4(T=410, x=0.5)
            >>> print "%0.6g %0.4f %0.2f %0.3f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            410 2.9578 377.13 95.371 -42.76 116.18 -0.452 -0.064 2.247 2.306 4.677 4.840 246.50 141.36
            >>> st=nC4(T=420, x=0.5)
            >>> print "%0.6g %0.4f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            420 3.4897 327.77 135.00 3.37 108.06 -0.344 -0.095 2.357 2.444 8.852 10.719 165.59 124.46
            >>> st=nC4(T=424, x=0.5)
            >>> print "%0.6g %0.4f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            424 3.7262 284.03 173.21 32.66 91.15 -0.277 -0.139 2.454 2.549 34.430 44.895 127.24 115.73
            >>> st=nC4(T=425, x=0.5)
            >>> print "%0.6g %0.4f %0.2f %0.2f %0.2f %0.2f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f %0.2f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            425 3.7881 250.17 205.54 50.78 73.93 -0.235 -0.180 2.534 2.589 375.35 460.13 114.85 112.63
            """
            # Table 45, Pag 980
            """
            >>> st=nC4(T=200, P=1e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.6g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            200 674.06 -590.85 -590.7 -2.2462 1.4733 2.062 1438.79
            >>> st=nC4(T=425, P=5e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            425 8.619 183.13 241.14 0.45027 2.1279 2.3062 244.71
            >>> st=nC4(T=425, P=1e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            425 18.184 175.73 230.73 0.3333 2.1469 2.3739 233.55
            >>> st=nC4(T=370, P=1.5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            370 473.15 -181.62 -178.45 -0.79043 2.0102 3.05 508.67
            >>> st=nC4(T=425, P=2e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            425 41.696 157.9 205.86 0.18997 2.195 2.6031 206.82
            >>> st=nC4(T=400, P=3e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 416.28 -90.374 -83.167 -0.55176 2.1613 3.623 348.01
            >>> st=nC4(T=425, P=4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            425 323.28 6.3288 18.702 -0.31172 2.364 7.3685 173.79
            >>> st=nC4(T=575, P=4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            575 55.37 518.95 591.19 0.87953 2.7242 3.0176 265.11
            >>> st=nC4(T=425, P=5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            425 368.63 -12.88 0.68419 -0.36083 2.2866 4.1643 268.65
            >>> st=nC4(T=425, P=1e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            425 427.92 -41.303 -17.934 -0.43378 2.2319 3.1867 454.17
            >>> st=nC4(T=500, P=6.9e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            500 515.49 90.412 224.26 -0.17641 2.533 2.986 923.31
            """,

        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 575., "Pmax": 200000.0, "rhomax": 13.86,
        "Pmin": 0.000653, "rhomin": 12.645,

        "nr1": [0.25536998241635e1, -0.44585951806696e1, 0.82425886369063,
                0.11215007011442, -0.35910933680333e-1, 0.16790508518103e-1,
                0.32734072508724e-1],
        "d1": [1, 1, 1, 2, 3, 4, 4],
        "t1": [0.50, 1.00, 1.50, 0.00, 0.50, 0.50, 0.75],
        "nr2": [0.95571232982005, -0.10003385753419e1, 0.85581548803855e-1,
                -0.25147918369616e-1, -0.15202958578918e-2, 0.47060682326420e-2,
                -0.97845414174006e-1, -0.48317904158760e-1, 0.17841271865468,
                0.18173836739334e-1, -0.11399068074953, 0.19329896666669e-1,
                0.11575877401010e-2, 0.15253808698116e-3, -0.43688558458471e-1,
                -0.82403190629989e-2],
        "d2": [1, 1, 2, 7, 8, 8, 1, 2, 3, 3, 4, 5, 5, 10, 2, 6],
        "t2": [2.00, 2.50, 2.50, 1.50, 1.00, 1.50, 4.00, 7.00, 3.00, 7.00,
               3.00, 1.00, 6.00, 0.00, 6.00, 13.00],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*16,

        "nr3": [-0.28390056949441e-1, 0.14904666224681e-2],
        "d3": [1, 2],
        "t3": [2., 0.],
        "alfa3": [10, 10],
        "beta3": [150, 200],
        "gamma3": [1.16, 1.13],
        "epsilon3": [0.85, 1.]}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for butane of Younglove and Ely (1987)",
        "__doi__": {"autor": "Younglove, B.A. and Ely, J.F.",
                    "title": "Thermophysical Properties of Fluids. II. Methane, Ethane, Propane, Isobutane, and Normal Butane ",
                    "ref": "J. Phys. Chem. Ref. Data 16, 577 (1987)",
                    "doi": "10.1063/1.555785"},

        "R": 8.31434,
        "cp": CP4,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 19275.7, "so": 309.909},

        "Tmin": 134.86, "Tmax": 500., "Pmax": 70000.0, "rhomax": 13.2,
        "Pmin": 6.736e-4, "rhomin": 12.65,

        "b": [None, 0.153740104603e-1, -0.160980034611, -0.979782459010e1,
              0.499660674504e3, -0.102115607687e7, 0.236032147756e-2,
              -0.137475757093e1, -0.907038733865e3, 0.385421748213e6,
              -0.349453710700e-4, 0.157361122714, 0.102301474068e3,
              0.182335737331e-1, -0.404114307787e1, 0.187979855783e1,
              0.362088795040, -0.738762248266e-2, -0.218618590563e1,
              0.118802729027, 0.706854198713e6, -0.219469885796e9,
              -0.182454361268e5, 0.206790377277e10, 0.111757550145e3,
              0.558779925986e5, -0.159579054026e2, -0.148034214622e7,
              -0.245206328201, 0.218305259309e3, -0.923990627338e-4,
              -0.205267776639e1, 0.387639044820e2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004",
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032-3091",
                    "doi": "10.1021/je300655b"},
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 575., "Pmax": 69000.0, "rhomax": 13.2,
        "Pmin": 0.000653, "rhomin": 12.645,

        "nr1": [0.10626277411455e1, -0.28620951828350e1, 0.88738233403777,
                -0.12570581155345, 0.10286308708106, 0.25358040602654e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.32325200233982, -0.37950761057432e-1, -0.32534802014452,
                -0.79050969051011e-1, -0.20636720547775e-1, 0.57053809334750e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Miyamoto and Watanabe (2001)",
        "__doi__": {"autor": "Miyamoto, H. and Watanabe, K.",
                    "title": "A Thermodynamic Property Model for Fluid-Phase n-Butane",
                    "ref": "Int. J. Thermophys., 22(2):459-475, 2001.",
                    "doi": "10.1023/A:1010722814682"},
        "R": 8.314472,
        "cp": Fi3,
        "ref": "IIR",

        "Tmin": 134.87, "Tmax": 589., "Pmax": 69000.0, "rhomax": 13.15,
        "Pmin": 0.000688, "rhomin": 12.652,

        "nr1": [2.952054e-1, -1.32636, -2.031317e-3, 2.240301e-1,
                -3.635425e-2, 1.905841e-3, 7.409154e-5, -1.401175e-6],
        "d1": [1, 1, 2, 2, 3, 5, 8, 8],
        "t1": [-0.25, 1.5, -0.75, 0, 1.25, 1.5, 0.5, 2.5],

        "nr2": [-2.492172, 2.386920, 1.424009e-3, -9.393388e-3, 2.616590e-3,
                -1.977323e-1, -3.809534e-2, 1.523948e-3, -2.391345e-2,
                -9.535229e-3, 3.928384e-5],
        "d2": [3, 3, 8, 5, 6, 1, 5, 7, 2, 3, 15],
        "t2": [1.5, 1.75, -0.25, 3, 3, 4, 2, -1, 2, 19, 5],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*11}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for butane of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"},
        "__test__": """
            >>> st=nC4(T=700, rho=200, eq=4)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            3.2176 18.416 3.5758
            >>> st2=nC4(T=750, rho=100, eq=4)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            213.77 0.37465
            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": 134.86, "Tmax": 750., "Pmax": 100000.0, "rhomax": 13.20,
        "Pmin": 0.00064578, "rhomin": 12.671,

        "nr1": [0.10626277e1, -0.28620952e1, 0.88738233, -0.12570581,
                0.10286309, 0.25358041e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.323252, -0.37950761e-1, -0.32534802, -0.79050969e-1,
                -0.20636721e-1, 0.57053809e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Polt et al. (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP6,
        "ref": "NBP",

        "Tmin": 140.0, "Tmax": 589., "Pmax": 30000.0, "rhomax": 12.81,
        "Pmin": 0.00161, "rhomin": 12.573,

        "nr1": [-0.504188295325, 0.541067401063, -0.760421383062e-1,
                0.846035653528, -0.191317317203e1, 0.521441860186,
                -0.783511318207, 0.689697797175e-1, 0.947825461055e-1,
                -0.141401831669, 0.382675021672, -0.423893176684e-1,
                0.677591792029e-1, 0.567943363340e-1, -0.131517698401,
                0.221136942526e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.504188295325, -0.541067401063, 0.760421383062e-1,
                -0.619109535460e-1, 0.423035373804, -0.390505508895],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.08974964]*6}

    helmholtz6 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for butane of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"},
        "R": 8.314472,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [1.18936994, 1.05407451, -3.24964532, 8.25263908e-2,
                2.76467405e-4, -8.09869214e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-9.38097492e-2, 1.46213532e-1, 4.01168502e-1, -1.28716120e-2,
                -2.75191070e-1, -1.62708971e-2, -7.04082962e-2, -2.32871995e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, MBWR, GERG, helmholtz3, helmholtz4, helmholtz5, helmholtz6

    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.0557549],  "expt0": [-1.], "expd0": [1.],
                   "a1": [20.611, 0.02], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [66.64, 24.44, -7461.2, -1983.6],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.00066566,
                "Tmin": 134.895, "Tmax": 575.0,
                "a1": [-558558235.4, 558558236.4], "exp1": [0, 2.206],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _surface = {"sigma": [0.05138], "exp": [1.209]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.71897e1, 0.26122e1, -0.21729e1, -0.27230e1],
        "exp": [1, 1.5, 2., 4.5]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.52341e1, -0.62011e1, 0.36063e1, 0.22137],
        "exp": [0.44, 0.6, 0.76, 5.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.27390e1, -0.57347e1, -0.16408e2, -0.46986e2, -0.10090e3],
        "exp": [0.39, 1.14, 3.0, 6.5, 14.0]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.17067154, -0.48879666, 0.039038856],
              "__name__": "Vogel (1999)",
              "__doi__": {"autor": "Vogel, E., Kuechenmeister, C., and Bich, E.",
                          "title": "Viscosity for n-Butane in the Fluid Region",
                          "ref": "High Temp. - High Pressures, 31(2):173-186, 1999.",
                          "doi": "10.1068/htrt154"},

              "ek": 280.51, "sigma": 0.57335,
              "Tref": 1, "rhoref": 1.*M, "etaref": 1.,
              "n_chapman": 0.1628213/M**0.5,

              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455, -0.34664158],
              "t_virial": [0.0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 280.51, "etaref_virial": 0.1135034,

              "Tref_res": 425.125, "rhoref_res": 3.92*M, "etaref_res": 1,
              "n_packed": [2.30873963359, 2.03404037254],
              "t_packed": [0, 0.5],
              "n_poly": [-54.7737770846, 58.0898623034, 0, 35.2658446259,
                         -39.6682203832, 0, -1.83729542151, 0, 0,
                         -0.833262985358, 1.93837020663, 0, -188.075903903],
              "t_poly": [0, -1, -2, 0, -1, -2, 0, -1, -2, 0, -1, -2, 0],
              "d_poly": [2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1],
              "g_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
              "n_num": [188.075903903],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [1, 0]}

    visco1 = {"eq": 2, "omega": 2,
              "__name__": "Younglove (1987)",
              "__doi__": {"autor": "Vogel, E., Kuechenmeister, C., Bich, E., and Laesecke, A.",
                          "title": "Reference Correlation of the Viscosity of Propane",
                          "ref": "J. Phys. Chem. Ref. Data 27, 947 (1998)",
                          "doi": "10.1063/1.556025"},

              "ek": 440., "sigma": 0.503103,
              "n_chapman": 0.20352457/M**0.5,
              "F": [0.1630521851e1, 0.0, 1.40, 425.16],
              "E": [-0.2724386845e2, 0.8012766611e3, 0.2503978646e2,
                    -0.1309704275e5, -0.8313305258e-1, 0.6636975027e2,
                    0.9849317662e4],
              "rhoc": 3.920}

    visco2 = {"eq": 4, "omega": 1,
              "__name__": "Quiñones-Cisneros (2006)",
              "__doi__": {"autor": "S.E.Quiñones-Cisneros and U.K. Deiters",
                          "title": "Generalization of the Friction Theory for Viscosity Modeling",
                          "ref": "J. Phys. Chem. B, 2006, 110 (25), pp 12820–12834",
                          "doi": "10.1021/jp0618577"},

              "Tref": 425.125, "muref": 1.0,
              "ek": 440., "sigma": 0.503103, "n_chapman": 0,
              "n_ideal": [18.3983, -57.1255, 49.3197],
              "t_ideal": [0, 0.25, 0.5],
              "a": [-1.34110938674421e-5, -8.56587924603951e-5, -6.45720639242339e-13],
              "b": [1.49859653515567e-4, -1.71133855507542e-4, 7.37953726544736e-13],
              "c": [3.53018109777015e-7, -1.93040375218067e-5, -1.26469933968355e-14],
              "A": [-3.63389393526204e-9, -7.73717469888952e-10, 0.0],
              "B": [3.70980259815724e-8, 2.07658634467549e-9, 0.0],
              "C": [-1.12495594619911e-7, 7.66906137372152e-8, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2002)",
               "__doi__": {"autor": "Perkins, R.A, Ramires, M.L.V., Nieto de Castro, C.A. and Cusco, L.",
                           "title": "Measurement and Correlation of the Thermal Conductivity of Butane from 135 K to 600 K at Pressures to 70 MPa",
                           "ref": "J. Chem. Eng. Data, 2002, 47 (5), pp 1263–1271",
                           "doi": "10.1021/je0101202"},

               "Tref": 425.16, "kref": 1.,
               "no": [1.62676e-3, 9.75703e-4, 2.89887e-2],
               "co": [0, 1, 2],

               "Trefb": 425.16, "rhorefb": 3.92, "krefb": 1.,
               "nb": [-3.04337e-2, 4.18357e-2, 1.65820e-1, -1.47163e-1,
                      -1.48144e-1, 1.33542e-1, 5.25500e-2, -4.85489e-2,
                      -6.29367e-3, 6.44307e-3],
               "tb": [0, 1, 0, 1, 0, 1, 0, 1, 0, 1],
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.875350e-9, "Tcref": 637.68}

    thermo1 = {"eq": 2, "omega": 2,
               "__name__": "Younglove (1987)",
               "__doi__": {"autor": "Younglove, B.A. and Ely, J.F.",
                           "title": "Thermophysical Properties of Fluids. II. Methane, Ethane, Propane, Isobutane, and Normal Butane",
                           "ref": "J. Phys. Chem. Ref. Data 16, 577 (1987)",
                           "doi": "10.1063/1.555785"},

               "visco": visco1,
               "n_chapman": 2.0352526600e-1,
               "G": [0.1530992335e1, -0.2114511021],
               "E": [0.4024170074e-2, 0.1561435847e1, -0.6004381127e3,
                     -0.7547260841e-3, -0.2069676662e-1, 0.9382534978e2,
                     -0.1711371457, 0.3647724935e2],

               "critical": 2,
               "X": [0.000769608, 13.2533, 0.485554, 1.01021],
               "Z": 9.10218e-10}

    _thermal = thermo0, thermo1
