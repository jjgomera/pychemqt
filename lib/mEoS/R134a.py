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


class R134a(MEoS):
    """Multiparameter equation of state for R134a"""
    name = "1,1,1,2-tetrafluoroethane"
    CASNumber = "811-97-2"
    formula = "CF3CH2F"
    synonym = "R134a"
    rhoc = unidades.Density(511.9)
    Tc = unidades.Temperature(374.21)
    Pc = unidades.Pressure(4059.28, "kPa")
    M = 102.032  # g/mol
    Tt = unidades.Temperature(169.85)
    Tb = unidades.Temperature(247.076)
    f_acent = 0.32684
    momentoDipolar = unidades.DipoleMoment(2.058, "Debye")
    id = 236
    # id = 1235

    Fi1 = {"ao_log": [1, -1.629789],
           "pow": [0, 1, -0.5, -0.75],
           "ao_pow": [-1.019535, 9.047135, -9.723916, -3.92717],
           "ao_exp": [], "titao": []}

    Fi2 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.25],
           "ao_pow": [10.78497786, 8.612977410, -24.37548384],
           "ao_exp": [7.451784998, -4.239239505, 6.445739825],
           "titao": [-4.103830338, -2.561528683, -2.084607363]}

    CP2 = {"ao": 19.4006,
           "an": [0.258531, -1.29665e-4], "pow": [1, 2],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    # TDOO: Add Huber-Ely meos, file in todo folder
    tillner = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Tillner-Roth & Baehr (1994).",
        "__doi__": {"autor": "Tillner-Roth, R. and Baehr, H.D.",
                    "title": "An international standard formulation of the thermodynamic properties of 1,1,1,2-tetrafluoroethane (HFC-134a) for temperatures from 170 K to 455 K at pressures up to 70 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 23, 657 (1994)",
                    "doi": "10.1063/1.555958"},
        "__test__":
            #Table 8, Pag 696
            """
            >>> st=R134a(T=169.85, x=0.5)
            >>> print "%0.2f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            169.85 0.00039 1591.1 0.02817 71.454 334.94 0.4126 1.9639 0.7922 0.5029 1.1838 0.5853 1119.9 126.79
            >>> st=R134a(T=180, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            180 0.00113 1564.2 0.07701 83.482 340.88 0.4813 1.9113 0.7912 0.5267 1.1871 0.6096 1068.3 130.05
            >>> st=R134a(T=200, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            200 0.00631 1510.4 0.38977 107.39 353.05 0.6073 1.8356 0.8015 0.5732 1.2057 0.5386 967.60 135.98
            >>> st=R134a(T=220, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            220 0.02443 1455.1 1.3850 131.77 365.65 0.7234 1.7865 0.8193 0.6203 1.2331 0.7109 869.85 141.00
            >>> st=R134a(T=250, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            250 0.11561 1367.8 5.9545 169.56 384.60 0.8841 1.7442 0.8514 0.6961 1.2864 0.8044 728.38 145.98
            >>> st=R134a(T=300, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            300 0.70282 1199.6 34.192 237.18 413.26 1.1286 1.7155 0.9144 0.8426 1.4324 1.0438 497.88 143.87
            >>> st=R134a(T=.350, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            350 2.4610 951.31 140.99 316.49 429.02 1.3674 1.6889 1.0036 1.0300 1.9613 1.8493 254.05 120.33
            >>> st=R134a(T=360, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            360 3.0404 870.11 193.58 336.05 427.07 1.4207 1.6735 1.0389 1.0853 2.4367 2.6063 196.04 111.24
            >>> st=R134a(T=370, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            370 3.7278 740.31 293.89 360.64 417.68 1.4856 1.6398 1.1145 1.1690 5.1048 6.8621 127.23 99.370
            >>> st=R134a(T=374, x=0.5)
            >>> print "%0.0f %0.5f %0.1f %0.5f %0.3f %0.2f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.2f" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            374 4.0416 587.91 434.05 380.85 399.50 1.5387 1.5885 1.2120 1.2409 101.66 137.23 92.401 91.389
            """
            #Table 9, Pag 700
            """
            >>> st=R134a(T=170, P=1e4)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            170 1590.7 71.636 0.4136 0.7921 1.1838 1119.2
            >>> st=R134a(T=210, P=2e4)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            210 1483.0 119.52 0.6664 0.8098 1.2186 918.36
            >>> st=R134a(T=305, P=3e4)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            305 1.2137 431.60 2.0221 0.7662 0.8501 165.14
            >>> st=R134a(T=440, P=4e4)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            440 1.1177 561.31 2.3484 0.9864 1.0688 196.72
            >>> st=R134a(T=260, P=6e4)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            260 2.8904 394.24 1.8338 0.6929 0.7842 151.67
            >>> st=R134a(T=370, P=1e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            370 3.3472 489.70 2.0970 0.8776 0.9631 180.23
            >>> st=R134a(T=230, P=1.6e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            230 1427.0 144.24 0.7784 0.8295 1.2490 822.77
            >>> st=R134a(T=200, P=2e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            200 1510.7 107.47 0.6070 0.8016 1.2055 968.48
            >>> st=R134a(T=275, P=3e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            275 14.679 400.05 1.7305 0.7625 0.8981 147.40
            >>> st=R134a(T=300, P=4e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 17.859 420.31 1.7796 0.7926 0.9207 153.82
            >>> st=R134a(T=235, P=5e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            235 1413.3 150.63 0.8048 0.8348 1.2568 801.11
            >>> st=R134a(T=175, P=6e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            175 1578.2 77.821 0.4473 0.7911 1.1843 1096.0
            >>> st=R134a(T=300, P=8e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 1200.2 237.19 1.1283 0.9143 1.4312 499.04
            >>> st=R134a(T=220, P=1e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            220 1457.1 132.16 0.7221 0.8195 1.2315 874.95
            >>> st=R134a(T=310, P=1.2e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            310 1161.9 251.69 1.1748 0.9285 1.4758 454.86
            >>> st=R134a(T=460, P=1.4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            460 39.636 574.72 2.0934 1.0324 1.1453 192.05
            >>> st=R134a(T=265, P=2e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            265 1327.3 189.62 0.9566 0.8689 1.3123 672.83
            >>> st=R134a(T=275, P=2.5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            275 1297.2 202.99 1.0047 0.8809 1.3342 631.65
            >>> st=R134a(T=180, P=3e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            180 1568.1 84.816 0.4781 0.7921 1.1846 1079.9
            >>> st=R134a(T=275, P=3.5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            275 1300.9 203.23 1.0027 0.8808 1.3293 639.71
            >>> st=R134a(T=370, P=4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            370 787.42 356.29 1.4729 1.0747 3.1454 157.31
            >>> st=R134a(T=270, P=5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            270 1321.8 197.03 0.9758 0.8748 1.3115 672.94
            >>> st=R134a(T=175, P=6e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            175 1584.8 80.250 0.4417 0.7928 1.1803 1115.9
            >>> st=R134a(T=200, P=1e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            200 1525.7 111.68 0.5958 0.8045 1.1953 1010.3
            >>> st=R134a(T=265, P=1.2e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            265 1357.1 192.75 0.9403 0.8696 1.2807 740.01
            >>> st=R134a(T=320, P=1.4e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            320 1202.9 266.29 1.1868 0.9343 1.3799 549.36
            >>> st=R134a(T=450, P=1.6e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            450 718.77 466.57 1.7036 1.1007 1.6893 243.70
            >>> st=R134a(T=220, P=2e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            220 1490.4 140.01 0.6992 0.8247 1.2076 962.36
            >>> st=R134a(T=355, P=3e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            355 1177.0 316.95 1.2971 0.9750 1.3680 570.01
            >>> st=R134a(T=185, P=4e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            185 1598.4 107.49 0.4744 0.8063 1.1659 1179.8
            >>> st=R134a(T=270, P=5e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            270 1425.2 213.72 0.9167 0.8842 1.2345 904.12
            >>> st=R134a(T=460, P=7e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            460 1113.9 471.45 1.5916 1.1023 1.4020 617.03
            """,

        "R": 8.314471,
        "Tref": 374.18, "rhoref": 508,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 455.0, "Pmax": 70000.0, "rhomax": 15.60,
        "Pmin": 0.3896, "rhomin": 15.5942,

        "nr1": [0.5586817e-1, 0.498223, 0.2458698e-1, 0.8570145e-3,
                0.4788584e-3, -0.1800808e1, 0.2671641, -0.4781652e-1],
        "d1": [2, 1, 3, 6, 6, 1, 1, 2],
        "t1": [-0.5, 0, 0, 0, 1.5, 1.5, 2, 2],

        "nr2": [0.1423987e-1, 0.3324062, -0.7485907e-2, 0.1017263e-3,
                -0.5184567, -0.8692288e-1, 0.2057144, -0.5000457e-2,
                0.4603262e-3, -0.3497836e-2, 0.6995038e-2, -0.1452184e-1,
                -0.1285458e-3],
        "d2": [5, 2, 2, 4, 1, 4, 1, 2, 4, 1, 5, 3, 10],
        "t2": [1, 3, 5, 1, 5, 5, 6, 10, 10, 10, 18, 22, 50],
        "c2": [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4],
        "gamma2": [1]*13}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-134a of Huber and McLinden (1992)",
        "__doi__": {"autor": "Huber, M.L. and McLinden, M.O.",
                    "title": "Thermodynamic properties of R134a (1,1,1,2-tetrafluoroethane)",
                    "ref": "International Refrigeration Conference, West Lafayette, IN, July 14-17, 453-462, 1992.",
                    "doi": "10.1016/j.fluid.2004.06.028"},
        "__test__": """
            >>> st=R134a(T=-103.3+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -103.30 0.00039 1591.2 0.028 71.89 335.08 0.4143 1.9638 1.147 0.585 1135 127
            >>> st=R134a(T=-80+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -80.00 0.00369 1526.2 0.235 99.65 349.03 0.5674 1.8585 1.211 0.637 999 134
            >>> st=R134a(T=-50+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -50.00 0.02948 1443.1 1.651 136.21 367.83 0.7432 1.7812 1.229 0.712 858 142
            >>> st=R134a(T=-25+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -25.00 0.10637 1371.2 5.505 167.43 383.57 0.8755 1.7465 1.270 0.788 742 146
            >>> st=R134a(T=0+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            0.00 0.29269 1293.7 14.420 200.00 398.68 1.0000 1.7274 1.335 0.883 626 147
            >>> st=R134a(T=25+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            25.00 0.66526 1206.3 32.318 234.47 412.44 1.1197 1.7166 1.425 1.012 508 144
            >>> st=R134a(T=50+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            50.00 1.3177 1102.0 66.164 271.59 423.63 1.2373 1.7078 1.569 1.218 387 137
            >>> st=R134a(T=75+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.3f %0.3f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            75.00 2.3639 963.3 133.31 313.14 429.26 1.3580 1.6916 1.915 1.730 262 122
            >>> st=R134a(T=100+273.15, x=0.5, eq=1)
            >>> print "%0.2f %0.5f %0.1f %0.3f %0.2f %0.2f %0.4f %0.4f %0.0f %0.0f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.w, st.Gas.w)
            100.00 3.9721 646.7 399.3 374.03 407.08 1.5207 1.6093 105 94
            """, # Table 9, Pag 459

        "R": 8.314471,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 70000.0, "rhomax": 15.60,
        "Pmin": 0.3922, "rhomin": 15.60,

        "b": [None, 0.965209362217e-1, -0.401824768889e1, 0.395239532858e2,
              0.134532868960e4, -0.139439741347e7, -0.309281355175e-2,
              0.292381512283e1, -0.165146613555e4, 0.150706003118e7,
              0.534973948313e-4, 0.543933317622, -0.211326049762e3,
              -0.268191203847e-1, -0.541067125950, -0.851731779398e3,
              0.205188253646, -0.733050188093e-2, 0.380655963862e1,
              -0.105832087589, -0.679243084424e6, -0.126998378601e9,
              -0.426234431829e5, 0.101973338234e10, -0.186699526782e3,
              -0.933426323419e5, -0.571735208963e1, -0.176762738787e6,
              -0.397282752308e-1, 0.143016844796e2, 0.803085294260e-4,
              -0.171959073552, 0.226238385661e1]}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-134a of Span and "
                    "Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. "
                             "III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1) (2003) 111-162",
                    "doi": "10.1023/A:1022362231796"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 15.6,
        "Pmin": 0.38818, "rhomin": 15.588,

        "nr1": [0.106631890000e1, -0.244959700000e1, 0.446457180000e-1,
                0.756568840000e-1, 0.206520890000e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.42006912, 0.76739111, 0.17897427e-2, -0.36219746,
                -0.6780937e-1, -0.10616419, -0.18185791e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    astina = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Astina and Sato (2004)",
        "__doi__": {"autor": "Astina, I.M. and Sato, H.",
                    "title": "A Fundamental Equation of State for 1,1,1,2-Tetrafluoroethane with an Intermolecular Potential Energy Background and Relialbe Ideal-Gas Properties",
                    "ref": "Fluid Phase Equilib., 221:103-111, 2004.",
                    "doi": "10.1016/j.fluid.2004.03.004"},

        "R": 8.314472,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 460.0, "Pmax": 70000.0, "rhomax": 15.58,
        "Pmin": 0.327, "rhomin": 15.58,

        "nr1": [1.832124209, -2.940698861, 5.156071823e-1, 2.756965911e-1,
                1.225264939, -6.486749497e-1, -9.286738053e-1, 3.920381291e-1,
                1.056692108e-1],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 4],
        "t1": [0.5, 1.125, 3.25, 0.5, 1.875, 2.75, 1.625, 2.125, 1.125],

        "nr2": [-7.586523371e-1, -1.198140136, -2.878260390e-1,
                -9.723032379e-2, 5.307113358e-2, -4.681610582e-2,
                -9.604697902e-3, 6.668035048e-3, 2.361266290e-3],
        "d2": [1, 2, 3, 2, 3, 4, 4, 5, 6],
        "t2": [3.75, 1.5, 0.75, 9, 8.5, 5.5, 32, 23, 31],
        "c2": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*9}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-134a of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.314471,
        "cp": Fi1,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [1.08605179, 1.03772416, -2.92069735, 9.15573346e-2,
                2.40541430e-4, -2.00239570e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-1.61424796e-2, -2.15499979e-1, 3.11819936e-1, 1.12867938e-3,
                -0.283454532, -4.21157950e-2, -8.08314045e-2, -1.59762784e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = tillner, MBWR, shortSpan, astina, sun
    _PR = 0.001032

    _surface = {"sigma": [0.05801], "exp": [1.241]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.77513e1, 0.29263e1, -0.26622e1, -0.39711e1],
        "exp": [1.0, 1.5, 1.9, 4.25]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.12449e2, -0.41023e2, 0.73641e2, -0.64635e2, 0.22551e2],
        "exp": [0.5, 0.7, 0.9, 1.1, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.29174e1, -0.72542e1, -0.23306e2, 0.59840e1, -0.71821e2],
        "exp": [0.383, 1.21, 3.3, 5.6, 7.0]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Huber (2003)",
              "__doi__": {"autor": "Huber, M.L., Laesecke, A., and Perkins, R.A.",
                          "title": "Model for the Viscosity and Thermal Conductivity of Refrigerants, Including a New Correlation for the Viscosity of R134a",
                          "ref": "Ind. Eng. Chem. Res., 2003, 42 (13), pp 3163–3178",
                          "doi": "10.1021/ie0300880"},

              "ek": 299.363, "sigma": 0.468932,
              "collision": [0.355404, -0.464337, 0.257353e-1],
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.215729/M**0.5,

              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 299.363, "etaref_virial": 0.0620984,

              "Tref_res": 374.21, "rhoref_res": 5.0170613*M, "etaref_res": 1000,
              "n_packed": [3.163695635587490, -0.8901733752064137e-1,
                           0.1000352946668359],
              "t_packed": [0, 1, 2],
              "n_poly": [-0.2069007192080741e-1, 0.3560295489828222e-3,
                         0.2111018162451597e-2, 0.1396014148308975e-1,
                         -0.4564350196734897e-2, -0.3515932745836890e-2,
                         -0.2147633195397038],
              "t_poly": [0, -6, -2, -0.5, 2, 0, 0],
              "d_poly": [1, 2, 2, 2, 2, 3, 0],
              "g_poly": [0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 0],
              "n_num": [0.2147633195397038],
              "t_num": [0],
              "d_num": [0],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    visco1 = {"eq": 4, "omega": 1,
              "__name__": "Quiñones-Cisneros (2006)",
              "__doi__": {"autor": "S.E.Quiñones-Cisneros and U.K. Deiters",
                          "title": "Generalization of the Friction Theory for Viscosity Modeling",
                          "ref": "J. Phys. Chem. B, 2006, 110 (25), pp 12820–12834",
                          "doi": "10.1021/jp0618577"},

              "Tref": 374.21, "muref": 1.0,
              "ek": 299.363, "sigma": 0.468932, "n_chapman": 0,
              "n_ideal": [31.2515, -89.6122, 73.0823],
              "t_ideal": [0, 0.25, 0.5],

              "a": [1.07271318464787e-4, -4.41655360682255e-5, 0.0],
              "b": [1.66457266522365e-4, -4.80292908400793e-5, 0.0],
              "c": [8.08333416284215e-5, -4.90359549823121e-5, 0.0],
              "A": [-2.12476175599662e-8, 2.81647242085073e-9, 0.0],
              "B": [1.35593527573090e-8, 0.0, 3.17549774078234e-10],
              "C": [0.0, 4.81768878752129e-7, -1.17148596093671e-7],
              "D": [0.0, 0.0, 0.0]}

    visco2 = {"eq": 1, "omega": 1,
              "__name__": "Laesecke (2003)",
              "__doi__": {"autor": "Laesecke, A.",
                          "title": "Data reassessment and full surface correlation of the viscosity of HFC-134a (1,1,1,2-tetrafluoroethane)",
                          "ref": "J. Phys. Chem. Ref. Data.",
                          "doi": ""},
              "__doc__": """, ",""",
              "ek": 288.82, "sigma": 0.50647,
              "collision": [0.355404, -0.464337, 0.257353e-1],
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.215729/M**0.5,

              "n_virial": [-0.17999496e1, 0.46692621e2, -0.53460794e3,
                           0.33604074e4, -0.13019164e5, 0.33414230e5,
                           -0.58711743e5, 0.71426686e5, -0.59834012e5,
                           0.33652741e5, -0.12027350e5, 0.24348205e4,
                           -0.20807957e3],
              "t_virial": [0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5,
                           -5, -5.5, -6],
              "Tref_virial": 288.82, "etaref_virial": 0.07823693,

              "Tref_res": 374.18, "rhoref_res": 4.9788302*M, "etaref_res": 1000,
              "n_packed": [3.07383, 0.482539055],
              "t_packed": [0, 1],
              "n_poly": [-0.331249e-1, -0.468509e-3, 0.306398],
              "t_poly": [0, 0, 0],
              "d_poly": [1, 2, 0],
              "g_poly": [0, 0, -1],
              "c_poly": [0, 0, 0],
              "n_num": [-0.306398, 0.215221],
              "t_num": [0, 0],
              "d_num": [0, 1],
              "g_num": [0, 0],
              "c_num": [0, 0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2000)",
               "__doi__": {"autor": "Perkins, R.A., Laesecke, A., Howley, J., Ramires, M.L.V., Gurova, A.N., and Cusco, L.",
                           "title": "Experimental thermal conductivity values for the IUPAC round-robin sample of 1,1,1,2-tetrafluoroethane (R134a)",
                           "ref": "NIST Interagency/Internal Report (NISTIR) - 6605",
                           "doi": ""},

               "Tref": 1., "kref": 1.,
               "no": [-1.05248e-2, 8.00982e-5],
               "co": [0, 1],

               "Trefb": 339.173, "rhorefb": 5.049886, "krefb": 2.055e-3,
               "nb": [1.836526, 5.126143, -1.436883, 6.261441e-1],
               "tb": [0]*4,
               "db": [1, 2, 3, 4],
               "cb": [0]*4,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 5.285356e-10, "Tcref": 561.411}

    _thermal = thermo0,


# class Test(TestCase):

    # def test_shortSpan(self):
        # # Table III, Pag 117
        # st = R134a(T=500, rho=500, eq="shortSpan")
        # self.assertEqual(round(st.cp0.kJkgK, 4), 1.1577)
        # self.assertEqual(round(st.P.MPa, 3), 14.656)
        # self.assertEqual(round(st.cp.kJkgK, 4), 1.6129)

        # st2 = R134a(T=600, rho=100, eq="shortSpan")
        # self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 181.97)
        # self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.41386)
