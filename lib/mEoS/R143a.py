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


class R143a(MEoS):
    """Multiparameter equation of state for R143a"""
    name = "1,1,1-trifluoroethane "
    CASNumber = "420-46-2"
    formula = "CF3CH3"
    synonym = "R143a"
    rhoc = unidades.Density(431.00006645)
    Tc = unidades.Temperature(345.857)
    Pc = unidades.Pressure(3761.0, "kPa")
    M = 84.041  # g/mol
    Tt = unidades.Temperature(161.34)
    Tb = unidades.Temperature(225.909)
    f_acent = 0.2615
    momentoDipolar = unidades.DipoleMoment(2.34, "Debye")
    id = 243

    Fi1 = {"ao_log": [1, -1],
           "pow": [0, 1, -0.33],
           "ao_pow": [5.903087, 7.307253, -16.59105],
           "ao_exp": [4.4402, 3.7515], "titao": [1791/Tc, 823/Tc]}

    CP2 = {"ao": 0.10002060,
           "an": [-0.96337511e-3, 0.31822397e3, 0.46917620e-1],
           "pow": [1.5, -1.25, 1],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP3 = {"ao": 1.838736,
           "an": [3.01994e-2, -1.78455e-5, 4.42442e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, -0.8999794],
           "pow": [0, 1, -1.5, 1.25, -1],
           "ao_pow": [-5.556942, 8.93748, 1.652398, -0.6827433, -8.113464],
           "ao_exp": [], "titao": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-143a of Lemmon and Jacobsen (2000).",
        "__doi__": {"autor": "Lemmon, E.W. and Jacobsen, R.T.",
                    "title": "An International Standard Formulation for the Thermodynamic Properties of 1,1,1-Trifluoroethane (HFC-143a) for Temperatures From 161 to 450 K and Pressures to 50 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 29, 521 (2000)",
                    "doi": "10.1063/1.1318909"},

        "__test__":
            #Table, Pag 541
            """
            >>> st=R143a(T=161.34, x=0.5)
            >>> print "%0.3f %0.5f %0.5g %0.5f %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -111.810 0.00107 1330.5 0.06754 52.52 319.59 0.31417 1.9695 0.8138 0.5283 1.211 0.6299 1058.1 137.6
            >>> st=R143a(T=-75+273.15, x=0.5)
            >>> print "%i %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -75 0.02189 1239.5 1.1379 97.831 342.6 0.56688 1.8021 0.8288 0.6304 1.26 0.7421 882.8 149.1
            >>> st=R143a(T=-50+273.15, x=0.5)
            >>> print "%i %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -50 0.08874 1173.9 4.2098 130.05 358.58 0.71968 1.7438 0.8608 0.7046 1.318 0.8331 763.9 154.2
            >>> st=R143a(T=-25+273.15, x=0.5)
            >>> print "%i %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -25 0.26144 1103.3 11.716 163.93 373.98 0.86289 1.7093 0.8991 0.7863 1.392 0.9492 644.7 155.8
            >>> st=R143a(T=0+273.15, x=0.5)
            >>> print "%i %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            0 0.61967 1024.3 27.306 200 387.81 1 1.6876 0.9408 0.8756 1.495 1.109 524.3 153.1
            >>> st=R143a(T=25+273.15, x=0.5)
            >>> print "%i %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            25 1.26157 930.22 57.653 239.19 398.54 1.1349 1.6693 0.9873 0.9737 1.669 1.369 399.5 144.5
            >>> st=R143a(T=50+273.15, x=0.5)
            >>> print "%i %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            50 2.30735 802.97 120.31 283.9 402.43 1.2748 1.6416 1.051 1.093 2.118 2.07 262.7 127.9
            >>> st=R143a(T=70+273.15, x=0.5)
            >>> print "%i %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            70 3.55268 600.85 270.1 333.2 385.42 1.4172 1.5694 1.198 1.272 7.72 11.5 122.4 104.2
            >>> st=R143a(T=345.857, x=0.5)
            >>> print "%0.3f %0.5f %0.5g %0.5g %0.5g" % (\
                st.T.C, st.P.MPa, st.rho, st.h.kJkg, st.s.kJkgK)
            72.707 3.76100 431 358.91 1.4906
            """
            #Table , Pag 544
            """
            >>> st=R143a(T=0+273.15, P=1e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0 3.788 375.11 401.51 1.9057 0.7918 0.9036 171.5
            >>> st=R143a(T=-100+273.15, P=2e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.5g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            -100 1302.2 66.816 66.97 0.3997 0.8117 1.22 1002.5
            >>> st=R143a(T=300+273.15, P=5e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 8.8753 690.02 746.36 2.5782 1.274 1.378 246.1
            >>> st=R143a(T=0+273.15, P=1e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0 1026.5 199.05 200.03 0.99874 0.9406 1.49 529
            >>> st=R143a(T=-30+273.15, P=1.5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            -30 1122.3 156.09 157.43 0.83188 0.8913 1.368 679.4
            >>> st=R143a(T=40+273.15, P=2e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            40 863.21 262.5 264.82 1.216 1.02 1.849 324.6
            >>> st=R143a(T=50+273.15, P=2.5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            50 808.38 280.31 283.4 1.2725 1.048 2.07 270.5
            >>> st=R143a(T=60+273.15, P=3e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            60 736.18 300.25 304.32 1.3343 1.089 2.582 208.1
            >>> st=R143a(T=0+273.15, P=4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0 1042.2 196.56 200.39 0.98946 0.9398 1.456 563.2
            >>> st=R143a(T=-100+273.15, P=5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.5g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            -100 1309.2 65.652 69.471 0.39291 0.815 1.214 1020.7
            >>> st=R143a(T=300+273.15, P=6e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 113.58 672.29 725.12 2.302 1.292 1.469 238.1
            >>> st=R143a(T=100+273.15, P=8e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            100 616.43 359.58 372.56 1.5068 1.117 2.3 201.4
            >>> st=R143a(T=0+273.15, P=1e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0 1068.3 192.37 201.73 0.97356 0.9401 1.411 620.1
            >>> st=R143a(T=-60+273.15, P=1.5e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            -60 1231.4 111.77 123.95 0.63496 0.8523 1.262 894.5
            >>> st=R143a(T=300+273.15, P=2e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 386.5 632.94 684.68 2.1235 1.321 1.629 271.5
            >>> st=R143a(T=0+273.15, P=2.5e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            0 1116.5 184.56 206.95 0.94246 0.9436 1.355 729.5
            >>> st=R143a(T=-80+273.15, P=5e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.5g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            -80 1324.2 79.79 117.55 0.46889 0.8378 1.201 1111.3
            >>> st=R143a(T=300+273.15, P=1e8)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g" % (\
                st.T.C, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 890.5 569.53 681.83 1.9149 1.368 1.595 671.2
            """,

        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1., "ho": 33936.4, "so": 198.961},

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 100000.0, "rhomax": 15.85,
        "Pmin": 1.0749, "rhomin": 15.832,

        "nr1": [.77736443e1, -.870185e1, -.27779799, .1460922, .89581616e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.67, 0.833, 1.7, 1.82, 0.35],

        "nr2": [-0.20552116, 0.10653258, 0.23270816e-1, -0.13247542e-1,
                -0.42793870e-1, 0.36221685, -0.25671899, -0.92326113e-1,
                0.83774837e-1, 0.17128445e-1, -0.17256110e-1, 0.49080492e-2],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [3.9, 0.95, 0, 1.19, 7.2, 5.9, 7.65, 7.5, 7.45, 15.5, 22, 19],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-143a of Outcalt and McLinden (1996)",
        "__doi__": {"autor": "Outcalt, S.L. and McLinden, M.O.",
                    "title": "An equation of state for the thermodynamic properties of R143a (1,1,1-trifluoroethane)",
                    "ref": "Int. J. Thermophys., 18(6):1445-1463, 1997.",
                    "doi": "10.1007/BF02575344"},

        "R": 8.314471,
        "cp": CP3,

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 15.84,
        "Pmin": 1.069, "rhomin": 15.8328,

        "b": [None, -0.240561786316e-1, 0.262345913719e1, -0.650858041394e2,
              0.995952053681e4, -0.147536464961e7, 0.135498153308e-2,
              -0.281726617426e1, 0.134371062574e4, 0.850286316514e6,
              -0.180516636446e-3, 0.618889066246, -0.223083798271e3,
              -0.119095922349e-1, -0.173933336877e1, -0.420847601180e3,
              0.213502079796, -0.565708555185e-2, 0.185442296800e1,
              -0.520377059921e-1, -0.846735696108e6, -0.207964483848e8,
              -0.349977290513e5, 0.576427827667e9, -0.389131863941e3,
              0.103074054089e5, -0.447627052215e1, -0.106673161101e6,
              -0.219511369081e-1, 0.642186519493e1, -0.938317030843e-4,
              -0.478594713528e-1, -0.206555883874e1]}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-143a of Li et al. (1999).",
        "__doi__": {"autor": "Li, J., Tillner-Roth, R., Sato, H., and Watanabe, K.",
                    "title": "An Equation of State for 1,1,1-Trifluoroethane (R-143a)",
                    "ref": "Int. J. Thermophys., 20(6):1639-1651, 1999.",
                    "doi": "10.1023/A:1022645626800"},
        "R": 8.31451,
        "cp": Fi2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 50000.0, "rhomax": 15.84,
        "Pmin": 1.0808, "rhomin": 15.819,

        "nr1": [.1606645e-1, .4163515e1, -.5031058e1, -.1920208e-1, .1470093e-2],
        "d1": [5, 1, 1, 2, 4],
        "t1": [0, 0.5, 0.75, 2.5, 2.5],

        "nr2": [0.1775429, -0.7316069e-2, -0.9555916e-1, -0.5822518,
                -0.4211022e-3, -0.2059847e-1, 0.3711325e-1, 0.1799723e-3,
                -0.4145922e-1, 0.7682566e-4, -0.2089695e-2, 0.1958633e-2,
                -0.3198325e-5, -0.5376561e-2],
        "d2": [3, 8, 3, 1, 10, 1, 4, 8, 2, 12, 8, 2, 5, 3],
        "t2": [0.25, 0.25, 2, 3, 3, 8, 8, 8, 10, 8, 17, 20, 35, 27],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4],
        "gamma2": [1]*14}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-143a of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. III. Results for Polar Fluids",
                    "ref": "Int. J. Thermophys., 24(1):111-162, 2003.",
                    "doi": "10.1023/A:1022362231796"},
        "__test__": """
            >>> st=R143a(T=700, rho=200, eq=3)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            1.2785 20.152 1.6702
            >>> st2=R143a(T=750, rho=100, eq=3)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            201.13 0.47846
            """, # Table III, Pag 117

        "R": 8.31451,
        "cp": CP2,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 15.82,
        "Pmin": 1.072, "rhomin": 15.816,

        "nr1": [.10306886e1, -.29497307e1, .6943523, .71552102e-1, .19155982e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.79764936e-1, 0.56859424, -0.90946566e-2, -0.24199452,
                -0.70610813e-1, -0.75041709e-1, -0.16411241e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1, MBWR, helmholtz2, helmholtz3

    _surface = {"sigma": [0.05416], "exp": [1.255]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73938e1, 0.19948e1, -0.18487e1, -0.41927e1, 0.14862e1],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.21135e1, 0.10200e2, -0.30836e2, 0.39909e2, -0.18557e2],
        "exp": [0.348, 1.6, 2.0, 2.4, 2.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.28673e1, -.63818e1, -.16314e2, -.45947e2, -.13956e1, -.24671e3],
        "exp": [0.384, 1.17, 3.0, 6.2, 7.0, 15.0]}
