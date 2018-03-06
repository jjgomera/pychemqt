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

from scipy import exp

from lib.meos import MEoS
from lib import unidades


class Ethylene(MEoS):
    """Multiparameter equation of state for ehylene"""
    name = "ethylene"
    CASNumber = "74-85-1"
    formula = "CH2=CH2"
    synonym = "R-1150"
    rhoc = unidades.Density(214.23999998)
    Tc = unidades.Temperature(282.35)
    Pc = unidades.Pressure(5041.8, "kPa")
    M = 28.05376  # g/mol
    Tt = unidades.Temperature(103.986)
    Tb = unidades.Temperature(169.379)
    f_acent = 0.0866
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 22
    _Tr = unidades.Temperature(273.316763)
    _rhor = unidades.Density(216.108926)
    _w = 0.085703183

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [8.68815523, -4.47960564],
           "ao_exp": [2.49395851, 3.0027152, 2.5126584, 3.99064217],
           "titao": [4.43266896, 5.74840149, 7.8027825, 15.5851154]}

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [1]*12,
           "exp": [4353.90715, 2335.22515, 1930.91322, 1471.92565, 4464.69725,
                   1778.39697, 1365.45204, 1356.81905, 4469.01375, 1188.47564,
                   4300.67034, 2077.67413],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 3.554495281,
           "an": [0.5603615762e6, -0.2141069802e5, 0.2532008897e3,
                  -0.9951927478e-2, 0.5108931070e-4, -0.1928667482e-7],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-0.2061703241e2],
           "exp": [3000],
           "ao_hyp": [], "hyp": []}

    smukala = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Smukala et al. (2000)",
        "__doi__": {"autor": "Smukala, J., Span, R., Wagner, W.",
                    "title": "New equation of state for ethylene covering the fluid region from the melting line to 450 K at pressures up to 300 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 29, 1053 (2000)",
                    "doi": "10.1063/1.1329318"},
        "__test__":
            # Table 32, Pag 1093
            """
            >>> st=Ethylene(T=103.989, x=0.5)
            >>> print "%0.3f %0.6f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            103.989 0.000122 654.6 0.00396 -819.13 -251.60 -4.8014 0.6561 1.622 0.89014 2.4295 1.1868 1766.6 202.68
            >>> st=Ethylene(T=120, x=0.5)
            >>> print "%0.0f %0.6f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            120 0.001368 634.17 0.03852 -780.21 -232.74 -4.4533 0.10889 1.553 0.89369 2.4271 1.192 1660.4 217.52
            >>> st=Ethylene(T=140, x=0.5)
            >>> print "%0.0f %0.6f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            140 0.01185 608.02 0.28758 -731.85 -209.71 -4.0807 -0.35111 1.4655 0.90645 2.4081 1.2125 1520.9 233.95
            >>> st=Ethylene(T=160, x=0.5)
            >>> print "%0.0f %0.6f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            160 0.056236 580.87 1.2123 -683.70 -188.03 -3.7597 -0.66179 1.3948 0.93356 2.4067 1.2599 1376.7 247.39
            >>> st=Ethylene(T=180, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            180 0.18185 552.20 3.5889 -635.17 -168.72 -3.4752 -0.88387 1.3465 0.97832 2.4413 1.3468 1227 256.99
            >>> st=Ethylene(T=200, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            200 0.45549 521.22 8.4936 -585.35 -152.89 -3.2155 -1.0533 1.3214 1.0431 2.5287 1.492 1070 261.95
            >>> st=Ethylene(T=220, x=0.5)
            >>> print "%0.0f %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            220 0.95664 486.67 17.452 -533 -141.96 -2.9709 -1.1935 1.3196 1.132 2.7003 1.7383 903.5 241.54
            >>> st=Ethylene(T=240, x=0.5)
            >>> print "%0.0f %0.4f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            240 1.7731 446.11 33.066 -476.16 -138.29 -2.7314 -1.3236 1.3435 1.254 3.0431 2.2112 724.4 254.88
            >>> st=Ethylene(T=260, x=0.5)
            >>> print "%0.0f %0.4f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            260 3.0036 393.47 61.542 -410.63 -147.57 -2.4812 -1.4695 1.4069 1.4388 3.9458 3.5116 524.13 240.47
            >>> st=Ethylene(T=280, x=0.5)
            >>> print "%0.0f %0.4f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            280 4.7836 290.7 140.7 -313.28 -203.54 -2.1411 -1.7492 1.7785 1.9809 19.563 29.261 246.68 208.88
            >>> st=Ethylene(T=282, x=0.5)
            >>> print "%0.0f %0.4f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            282 5.0023 253.12 175.8 -287.09 -232.15 -2.0508 -1.856 2.2089 2.3981 146.97 225.24 188.90 191.32
            """
            # Table 33, Pag 1097
            """
            >>> st=Ethylene(T=120, P=1e5)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            120 627.76 -768.13 -767.97 -4.3546 1.5303 2.422 1626.6
            >>> st=Ethylene(T=450, P=1e5)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            450 0.75081 138.75 271.94 0.73428 1.7662 2.065 394.35
            >>> st=Ethylene(T=200, P=5e5)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            200 521.3 -586.28 -585.32 -3.2158 1.3214 2.528 1070.5
            >>> st=Ethylene(T=300, P=1e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 11.932 -97.036 -13.23 -0.70588 1.264 1.6432 320.09
            >>> st=Ethylene(T=120, P=1.5e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            120 635.02 -780.78 -778.42 -4.4581 1.554 2.424 1669
            >>> st=Ethylene(T=400, P=2e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            400 17.627 41.788 155.25 -0.41873 1.6102 1.9804 365.48
            >>> st=Ethylene(T=255, P=3e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            255 411.09 -436.46 -429.17 -2.5531 1.381 3.515 591.37
            >>> st=Ethylene(T=120, P=4e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            105 655.24 -817.9 -811.79 -4.7896 1.6188 2.4244 1781.9
            >>> st=Ethylene(T=280, P=5e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            280 311.34 -340.71 -324.65 -2.1843 1.5996 9.1274 306.79
            >>> st=Ethylene(T=105, P=6e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            105 656.17 -818.49 -809.35 -4.7954 1.6194 2.4215 1792.4
            >>> st=Ethylene(T=450, P=7e6)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            450 57.569 98.708 220.3 -0.61249 1.8068 2.3187 379.44
            >>> st=Ethylene(T=200, P=1e7)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            200 535.31 -596.46 -577.78 -3.268 1.3312 2.4163 1174.5
            >>> st=Ethylene(T=450, P=1.5e7)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            450 132 50.02 163.66 -0.93729 1.8444 2.6487 392.55
            >>> st=Ethylene(T=110, P=2e7)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            110 656.51 -810.64 -780.17 -4.7224 1.609 2.4043 1827.8
            >>> st=Ethylene(T=300, P=2.5e7)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 417 -377.14 -317.19 -2.3407 1.4528 2.6835 756.22
            >>> st=Ethylene(T=130, P=5e7)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            130 647.55 -773.34 -696.13 -4.4046 1.5638 2.3411 1826.4
            >>> st=Ethylene(T=150, P=1e8)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 649.39 -743.25 -589.25 -4.1916 1.5535 2.2574 1891.7
            >>> st=Ethylene(T=450, P=2e8)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            450 512.76 -174.18 215.86 -1.9089 2.0216 2.5659 1444.4
            >>> st=Ethylene(T=310, P=3e8)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            310 621.79 -479.7 2.7729 -2.9385 1.7006 2.2342 1964.4
            >>> st=Ethylene(T=450, P=3e8)
            >>> print "%0.0f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            450 562.26 -195.31 338.25 -2.0494 2.0671 2.5681 1740.6
            """,

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 300000.0, "rhomax": 27.03,
        "Pmin": 0.12196, "rhomin": 23.334,

        "nr1": [0.18617429100670e1, -0.30913708460844e1, -0.17384817095516,
                0.80370985692840e-1, 0.23682707317354, 0.21922786610247e-1],
        "d1": [1, 1, 1, 2, 2, 4],
        "t1": [0.5, 1., 2.5, 0.0, 2.0, 0.5],

        "nr2": [0.11827885813193, -0.21736384396776e-1, 0.44007990661139e-1,
                0.12554058863881, -0.13167945577241, -0.52116984575897e-2,
                0.15236081265419e-3, -0.24505335342756e-4, 0.28970524924022,
                -0.18075836674288, 0.15057272878461, -0.14093151754458,
                0.22755109070253e-1, 0.14026070529061e-1, 0.61697454296214e-2,
                -0.41286083451333e-3, .12885388714785e-1, -.69128692157093e-1,
                0.10936225568483, -0.81818875271794e-2, -0.56418472117170e-1,
                0.16517867750633e-2, .95904006517001e-2, -.26236572984886e-2],
        "d2": [1, 1, 3, 4, 5, 7, 10, 11, 1, 1, 2, 2, 4, 4, 6, 7, 4, 5, 6, 6,
               7, 8, 9, 10],
        "t2": [1., 4., 1.25, 2.75, 2.25, 1., 0.75, 0.5, 2.5, 3.5, 4., 6., 1.5,
               5., 4.5, 15., 20., 23., 22., 29., 19., 15., 13., 10.],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4,
               4, 4, 4],
        "gamma2": [1]*24,

        "nr3": [-0.50242414011355e2, 0.74846420119299e4, -0.68734299232625e4,
                -0.93577982814338e3, 0.94133024786113e3],
        "d3": [2, 2, 2, 3, 3],
        "t3": [1., 0., 1., 2., 3.],
        "alfa3": [25.]*5,
        "beta3": [325, 300, 300, 300, 300],
        "gamma3": [1.16, 1.19, 1.19, 1.19, 1.19],
        "epsilon3": [1.]*5}

    jahangiri = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Jahangiri et al. (1986)",
        "__doi__": {"autor": "Jahangiri, M., Jacobsen, R.T, Stewart, R.B., and McCarty, R.D.",
                    "title": "Thermodynamic properties of ethylene from the freezing line to 450 K at pressures to 260 MPa",
                    "ref": "J. Phys. Chem. Ref. Data 15, 593 (1986)",
                    "doi": "10.1063/1.555753"},
        "__test__":
            # Table 18, Pag 613
            """
            >>> st=Ethylene(T=200, x=0.5, eq=2)
            >>> print "%0.1f %0.5f %0.3f %0.3f" % (st.T, st.cpM0.kJkmolK, st.hM0.kJkmol, st.sM0.kJkmolK)
            200.0 35.352 -3801.885 -15.367
            >>> st=Ethylene(T=250, x=0.5, eq=2)
            >>> print "%0.1f %0.5f %0.3f %0.3f" % (st.T, st.cpM0.kJkmolK, st.hM0.kJkmol, st.sM0.kJkmolK)
            250.0 38.624 -1958.438 -7.153
            >>> st=Ethylene(T=300, x=0.5, eq=2)
            >>> print "%0.1f %0.5f %0.3f %0.3f" % (st.T, st.cpM0.kJkmolK, st.hM0.kJkmol, st.sM0.kJkmolK)
            300.0   43.027  79.437  0.266
            >>> st=Ethylene(T=350, x=0.5, eq=2)
            >>> print "%0.1f %0.5f %0.3f %0.3f" % (st.T, st.cpM0.kJkmolK, st.hM0.kJkmol, st.sM0.kJkmolK)
            350.0 47.964 2353.064 7.265
            >>> st=Ethylene(T=400, x=0.5, eq=2)
            >>> print "%0.1f %0.5f %0.3f %0.3f" % (st.T, st.cpM0.kJkmolK, st.hM0.kJkmol, st.sM0.kJkmolK)
            400.0 52.989 4877.189 13.999
            >>> st=Ethylene(T=450, x=0.5, eq=2)
            >>> print "%0.1f %0.5f %0.3f %0.3f" % (st.T, st.cpM0.kJkmolK, st.hM0.kJkmol, st.sM0.kJkmolK)
            450.0 57.844 7649.04 20.523
            >>> st=Ethylene(T=500, x=0.5, eq=2)
            >>> print "%0.1f %0.5f %0.3f %0.3f" % (st.T, st.cpM0.kJkmolK, st.hM0.kJkmol, st.sM0.kJkmolK)
            500.0 62.411 10656.721 26.856
            """
            # Table 24, Pag 635
            """
            >>> st=Ethylene(T=200.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f" % (st.T, st.virialB.ccg*st.M)
            200.15 -310.248
            >>> st=Ethylene(T=250.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f" % (st.T, st.virialB.ccg*st.M)
            250.15 -199.921
            >>> st=Ethylene(T=300.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f" % (st.T, st.virialB.ccg*st.M)
            300.15 -138.087
            >>> st=Ethylene(T=350.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f" % (st.T, st.virialB.ccg*st.M)
            350.15 -98.356
            >>> st=Ethylene(T=400.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f" % (st.T, st.virialB.ccg*st.M)
            400.15 -66.019
            >>> st=Ethylene(T=450.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f" % (st.T, st.virialB.ccg*st.M)
            450.15 -50.099
            """
            # Table 25, Pag 637
            """
            >>> st=Ethylene(T=238.18, x=0.5, eq=2)
            >>> print "%0.2f %0.3f %0.1f %0.1f" % (st.T, st,P.MPa, st.Liquido.hM.Jmol, st.Gas.hM.Jmol)
            238.18 1.681 16108.4 25739.7
            >>> st=Ethylene(T=273.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f %0.1f %0.1f" % (st.T, st,P.MPa, st.Liquido.hM.Jmol, st.Gas.hM.Jmol)
            273.15 4.099 19639.9 24839.8
            >>> st=Ethylene(T=282.15, x=0.5, eq=2)
            >>> print "%0.2f %0.3f %0.1f %0.1f" % (st.T, st,P.MPa, st.Liquido.hM.Jmol, st.Gas.hM.Jmol)
            282.15 5.018 24644.4 23021.6
            """,

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 29610, "so": 219.225},

        "Tmin": Tt, "Tmax": 450.0, "Pmax": 260000.0, "rhomax": 26.67,
        "Pmin": 0.1225, "rhomin": 23.348,

        "nr1": [0.324893703388e1, -0.101727886161e2, 0.738660405252e1,
                -0.156891635862e1, -0.888451428662e-1, 0.602106814262e-1,
                0.107832458846, -0.200402521069e-1, 0.195049141244e-2,
                0.671800640346e-1, -0.420045146918e-1, -0.162050762577e-2,
                0.555515679497e-3, 0.758367114630e-3, -0.287854402074e-3],
        "d1": [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 6, 6, 6],
        "t1": [0.5, 1, 1.25, 1.75, 4, 2, 4, 5, 6, 0.25, 3, 0.25, 0.5, 2.5, 3],

        "nr2": [0.6258987063e-1, -0.641843116e-1, -0.1368693752, 0.517920766,
                -0.3026331319, 0.7757213872, -0.2639890864e1, 0.2927563554e1,
                -0.1066267599e1, -0.538047154e-1, 0.127792108, -0.745015231e-1,
                -0.1624304356e-1, 0.1476032429, -0.2003910489, 0.2926905618,
                -0.1389040901, 0.5913513541e1, -0.380037013e2, 0.969194057e2,
                -0.1226256839e3, 0.7702379476e2, -0.1922684672e2,
                -0.3800045701e-2, 0.1118003813e-1, 0.2945841426e-2],
        "d2": [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4,
               4, 4, 8, 8, 8],
        "t2": [0.5, 1, 0.5, 2, 4, 3, 4, 5, 6, 2, 3, 4, 1.5, 0.5, 1.5, 4, 5, 1,
               2, 3, 4, 5, 6, 0.5, 1, 5],
        "c2": [3, 3, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 3, 2, 2, 2, 2, 4, 4, 4, 4,
               4, 4, 2, 2, 2],
        "gamma2": [1]*26}

    shortSpan = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for ethylene of Span "
                    "and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. "
                             "II. Results for nonpolar fluids.",
                    "ref": "Int. J. Thermophys. 24 (1) (2003) 41-109",
                    "doi": "10.1023/A:1022310214958"},

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",

        "Tmin": Tt, "Tmax": 750.0, "Pmax": 100000.0, "rhomax": 27.03,
        "Pmin": 0.12123, "rhomin": 23.34,

        "nr1": [0.9096223, -0.24641015e1, 0.56175311, -0.19688013e-1,
                0.78831145e-1, 0.21478776e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.23151337, -0.37804454e-1, -0.20122739, -0.44960157e-1,
                -0.2834296e-1, 0.12652824e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethylene of McCarty and Jacobsen (1981)",
        "__doi__": {"autor": "McCarty, R.D., Jacobsen, R.T.",
                    "title": "An Equation of State for Fluid Ethylene",
                    "ref": "Natl. Bur. Stand., Tech. Note 1045, 1981.",
                    "doi": ""},
        "__test__":
            # Table, Pag 138
            """
            >>> st=Ethylene(T=110, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 23.2516 110 7369.71 87.438 43.54 72.05 1706.12
            >>> st=Ethylene(T=150, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 21.4687 150 10056.28 108.312 40.64 65.72 1507.50
            >>> st=Ethylene(T=200, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 19.0802 200 13387.64 127.463 37.33 68.11 1172.01
            >>> st=Ethylene(T=250, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 16.1946 250 16970.99 143.418 38.08 77.31 811.65
            >>> st=Ethylene(T=300, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 11.5343 300 21586.21 160.141 42.53 118.76 418.90
            >>> st=Ethylene(T=350, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 5.3640 350 27794.36 179.364 44.3 92.15 309.36
            >>> st=Ethylene(T=400, P=1e7, eq=1)
            >>> print "%0.0f %0.4f %0.0f %0.2f %0.3f %0.2f %0.2f %0.2f" % (\
                st.P.MPa, st.rho, st.T, st.hM.kJkmol, st.s.kJkmolK, st.cv.kJkmolK, st.cp.kJkmolK, st.Liquido.w, st.Gas.w)
            10 3.7626 400 31676.15 189.761 46.49 70.37 347.16
            """,
        "R": 8.31434,
        "cp": CP2,
        "ref": "IIR",

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 40000.0, "rhomax": 23.343,
        "Pmin": 0.1213, "rhomin": 23.343,

        "b": [None, -0.2146684366683e-1, 0.1791433722534e1, -0.3675315603930e2,
              0.3707178934669e4, -0.3198282566709e6, 0.5809379774732e-3,
              -0.7895570824899, 0.1148620375835e3, 0.2713774629193e6,
              -0.8647124319107e-4, 0.1617727266385, -0.2731527496271e2,
              -0.2672283641459e-2, -0.4752381331990e-1, -0.6255637346217e2,
              0.4576234964434e-2, -0.7534839269320e-4, 0.1638171982209,
              -0.3563090740740e-2, -0.1833000783170e6, -0.1805074209985e8,
              -0.4794587918874e4, 0.3531948274957e8, -0.2562571039155e2,
              0.1044308253292e4, -0.1695303363659, -0.1710334224958e4,
              -0.2054114462372e-3, 0.6727558766661e-1, -0.1557168403328e-5,
              -0.1229814736077e-3, 0.4234325938573e-3]}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.31451,
        "cp": CP1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 29610, "so": 219.225},

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [8.42278605e-1, 8.65139678e-1, -2.79801027, 6.74520156e-2,
                2.42445468e-4, -2.74767618e-3],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-1.48602227e-2, 1.29307481e-1, 3.74759088e-1, -1.25336440e-2,
                -2.33507187e-1, 1.38862785e-2, -4.88033330e-2, -2.38141707e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = smukala, MBWR, jahangiri, shortSpan, sun

    _surface = {"sigma": [0.0477], "exp": [1.17]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [10.725], "expt1": [0], "expd1": [1],
                   "a2": [55.19, 49.5, -2045, -1154.],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.9, 2.9]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 1000,
                "Tmin": Tt, "Tmax": 450.0,
                "a1": [0.1225e-3, 0.357924e3, -0.357924e3],
                "exp1": [0, 0.20645e1, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.3905741, 1.4060338, -1.6589923, 1.0278028, -2.5071716],
        "exp": [1.0, 1.5, 2.5, 3.0, 4.5]}
    _liquid_Density = {
        "eq": 4,
        "ao": [1.8673079, -0.61533892, -0.058973772, 0.10744720],
        "exp": [1.029, 1.5, 4.0, 6.0]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-1.9034556, -0.75837929, -3.7717969, -8.7478586, -23.885296,
               -54.197979],
        "exp": [1.047, 2.0, 3.0, 7.0, 14.5, 28.0]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "Holland (1983)",
              "__doi__": {"autor": "Holland, P.M., Eaton, B.E., and Hanley, H.J.M.",
                          "title": "A Correlation of the Viscosity and Thermal Conductivity Data of Gaseous and Liquid Ethylene",
                          "ref": "J. Phys. Chem. Ref. Data 12, 917 (1983)",
                          "doi": "10.1063/1.555701"},
              "__test__":
                  # Table 5, pag 924
                  """
                  >>> st=Ethylene(T=110, P=1e5)
                  >>> print "%0.1f" % st.mu.muPas*10
                  5660.5
                  >>> st=Ethylene(T=140, P=1e6)
                  >>> print "%0.1f" % st.mu.muPas*10
                  2769.8
                  >>> st=Ethylene(T=200, P=5e6)
                  >>> print "%0.1f" % st.mu.muPas*10
                  1223.7
                  >>> st=Ethylene(T=300, P=1e5)
                  >>> print "%0.1f" % st.mu.muPas*10
                  103.8
                  >>> st=Ethylene(T=130, P=1e7)
                  >>> print "%0.1f" % st.mu.muPas*10
                  3278.5
                  >>> st=Ethylene(T=300, P=5e7)
                  >>> print "%0.1f" % st.mu.muPas*10
                  759.0
                  >>> st=Ethylene(T=500, P=1e5)
                  >>> print "%0.1f" % st.mu.muPas*10
                  165.1
                  >>> st=Ethylene(T=500, P=5e7)
                  >>> print "%0.1f" % st.mu.muPas*10
                  394.1
                  """}

    def _visco0(self, rho, T, fase):
        GV = [-3.5098225018e6, 2.5008406184e6, -5.8365540744e5, 4.5549146583e3,
              2.2881683403e4, -4.7318682077e3, 4.5022249258e2, -2.1490688088e1,
              4.1649263233e-1]
        muo = 0
        for i in range(-3, 6):
            muo += GV[i+3]*T**(i/3.)

        mu1 = 0
        tita = (rho-self.rhoc)/self.rhoc
        j = [0, -4.8544486732, 1.3033585236e1, 2.7808928908e4, -1.8241971308e3,
             1.5913024509, -2.0513573927e2, -3.9478454708e4]
        deltamu = exp(j[1]+j[4]/T)*(exp(rho.gcc**0.1*(j[2]+j[3]/T**1.5)+tita*rho.gcc**0.5*(j[5]+j[6]/T+j[7]/T**2))-1.)

        return unidades.Viscosity((muo+mu1+deltamu)*1e-7, "Pas")

    visco1 = {"eq": 2, "omega": 2,
              "__name__": "NIST",
              "__doi__": {"autor": "",
                          "title": "Coefficients are taken from NIST14, Version 9.08",
                          "ref": "",
                          "doi": ""},

              "ek": 224.7, "sigma": 0.4163,
              "n_chapman": 0.141374566253583/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-8.03553028329404, -439.8962514, 8.69536237617, 5773.08496161,
                    0.267589139152, -34.39391627, 66.4795135739],
              "rhoc": 7.63299886259}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "Holland (1983)",
               "__doi__": {"autor": "Holland, P.M., Eaton, B.E., and Hanley, H.J.M.",
                           "title": "A Correlation of the Viscosity and Thermal Conductivity Data of Gaseous and Liquid Ethylene",
                           "ref": "J. Phys. Chem. Ref. Data 12, 917 (1983)",
                           "doi": "10.1063/1.555701"},
               "__test__":
                   # Table 6, pag 927
                   """
                   >>> st=Ethylene(T=110, P=1e5)
                   >>> print "%0.2f" % st.k.mWmK
                   261.77
                   >>> st=Ethylene(T=140, P=1e6)
                   >>> print "%0.2f" % st.k.mWmK
                   223.14
                   >>> st=Ethylene(T=200, P=5e6)
                   >>> print "%0.2f" % st.k.mWmK
                   158.50
                   >>> st=Ethylene(T=300, P=1e5)
                   >>> print "%0.2f" % st.k.mWmK
                   20.56
                   >>> st=Ethylene(T=130, P=1e7)
                   >>> print "%0.2f" % st.k.mWmK
                   244.97
                   >>> st=Ethylene(T=300, P=5e7)
                   >>> print "%0.2f" % st.k.mWmK
                   129.32
                   >>> st=Ethylene(T=500, P=1e5)
                   >>> print "%0.2f" % st.k.mWmK
                   49.95
                   >>> st=Ethylene(T=500, P=5e7)
                   >>> print "%0.2f" % st.k.mWmK
                   93.57
                   """}

    def _thermo0(self, rho, T, fase):
        GT = [-2.903423528e5, 4.680624952e5, -1.8954783215e5, -4.8262235392e3,
              2.243409372e4, -6.6206354818e3, 8.9937717078e2, -6.0559143718e1,
              1.6370306422]
        lo = 0
        for i in range(-3, 6):
            lo += GT[i+3]*T**(i/3.)

        tita = (rho.gcc-0.221)/0.221
        j = [0, -1.304503323e1, 1.8214616599e1, 9.903022496e3, 7.420521631e2,
             -3.0083271933e-1, 9.6456068829e1, 1.350256962e4]
        l1 = exp(j[1]+j[4]/T)*(exp(rho.gcc**0.1*(j[2]+j[3]/self.T**1.5)+tita*rho.gcc**0.5*(j[5]+j[6]/T+j[7]/T**2))-1.)

        lc = 0
        # FIXME: no sale
#        deltarho=(self.rho/self.M-0.221)/0.221
#        deltaT=(self.T-282.34)/282.34
#        xkt=(1.0/self.rho/self.M/self.derivative("P", "rho", "T")*1e3)**0.5
#        b=abs(deltarho)/abs(deltaT)**1.19
#        xts=(self.rho/self.M)**2*xkt*5.039/.221**2
#        g=xts*abs(deltaT)**1.19
#        xi=0.69/(b**2*5.039/g/Boltzmann/282.34)**0.5
#        f=exp(-18.66*deltaT**2-4.25*deltarho**4)
#        c=(self.M/self.rho.gcc/Avogadro/Boltzmann/self.T)**0.5
#        lc=c*Boltzmann*self.T**2/6.0/pi/self.mu.muPas/xi*self.dpdT**2*self.kappa**0.5*f
#        print lo, l1
        return unidades.ThermalConductivity(lo+l1+lc, "mWmK")


    thermo1 = {"eq": 1, "critical": 0,
               "__name__": "NIST14",
               "__doi__": {"autor": "",
                           "title": "Coefficients are taken from NIST14, Version 9.08",
                           "ref": "",
                           "doi": ""},

               "Tref": 224.7, "kref": 1e-3,
               "no": [1.35558587, -0.14207565869509, 1],
               "co": [0, -1, -96],

               "Trefb": 282.350007277, "rhorefb": 7.63299886259, "krefb": 1e-3,
               "nb": [15.3064493136, 25.0280721432, -15.4526955192,
                      0.8590418672, 3.32700049633, -0.333048907849],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6}

    _thermal = thermo0, thermo1

# TODO: Add MBWR equation of Younglove


class Test(TestCase):

    def test_shortSpan(self):
        # Table III, Pag 46
        st = Ethylene(T=700, rho=200, eq="shortSpan")
        self.assertEqual(round(st.cp0.kJkgK, 4), 2.7682)
        self.assertEqual(round(st.P.MPa, 3), 48.416)
        self.assertEqual(round(st.cp.kJkgK, 4), 3.0651)

        st2 = Ethylene(T=750, rho=100, eq="shortSpan")
        self.assertEqual(round(st2.h.kJkg-st.h.kJkg, 2), 174.10)
        self.assertEqual(round(st2.s.kJkgK-st.s.kJkgK, 5), 0.47681)
