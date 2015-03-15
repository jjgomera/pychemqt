#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class CH4(MEoS):
    """Multiparameter equation of state for methane"""
    name = "methane"
    CASNumber = "74-82-8"
    formula = "CH4"
    synonym = "R-50"
    rhoc = unidades.Density(162.66)
    Tc = unidades.Temperature(190.564)
    Pc = unidades.Pressure(4599.2, "kPa")
    M = 16.0428  # g/mol
    Tt = unidades.Temperature(90.694)
    Tb = unidades.Temperature(111.667)
    f_acent = 0.01142
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 2
    _Tr = unidades.Temperature(186.659809)
    _rhor = unidades.Density(163.413536)
    _w = 0.010528102

    CP1 = {"ao": 4.00160,
           "an": [], "pow": [],
           "ao_exp": [0.008449, 4.6942, 3.4865, 1.6572, 1.4115],
           "exp": [648, 1957, 3895, 5705, 15080],
           "ao_hyp": [], "hyp": []}
           
    Fi1 = {"ao_log": [1, 3.00160],
           "pow": [0, 1],
           "ao_pow": [9.91243972, -6.33270087],
           "ao_exp": [0.008449, 4.6942, 3.4865, 1.6572, 1.4115],
           "titao": [648/Tc, 1957/Tc, 3895/Tc, 5705/Tc, 15080/Tc], 
           "ao_hyp": [], "hyp": []}

    Fi2 = {"R": 8.314510, 
           "ao_log": [1, 3.00088],
           "pow": [0, 1],
           "ao_pow": [19.597508817, -83.959667892],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [0.76315, 0.0046, 8.74432, -4.46921],
           "hyp": [4.306474465, 0.936220902, 5.577233895, 5.722644361]}

    CP2 = {"ao": 4.00088,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.76315, 0.0046, 8.74432, -4.46921],
           "hyp": [4.306474465*Tc, 0.936220902*Tc, 5.577233895*Tc, 5.722644361*Tc]}

    CP3 = {"ao": 3.5998324,
           "an": [2.614717613495e-1, -5.671028952515e-2, 4.105505612671e-3],
           "pow": [1./3, 2./3, 1],
           "ao_exp": [4.7206715],
           "exp": [2009.15202],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 0.15438149595e2,
           "an": [-0.18044750507e7, 0.77426666393e5, -0.13241658754e4,
                  -0.51479005257e-1, 0.10809172196e-3, -0.65501783437e-7],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-0.67490056171e1], "exp": [3000],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methane of Setzmann and Wagner (1991)",
        "__doi__": {"autor": "Setzmann, U., Wagner, W.",
                    "title": "A New Equation of State and Tables of Thermodynamic Properties for Methane Covering the Range from the Melting Line to 625 K at Pressures up to 1000 MPa", 
                    "ref": "J. Phys. Chem. Ref. Data 20, 1061 (1991)",
                    "doi": "10.1063/1.555898"}, 
        "__test__": 
            # Table 39, Pag 1117
            """
            >>> st=CH4(T=90.694, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            90.694 0.011696 451.48 0.25074 -982.76 -438.5 -7.3868 -1.3857 2.1677 1.5735 3.3678 2.11 1538.6 249.13
            
            >>> st=CH4(T=100, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            100 0.034376 438.89 0.67457 -951.21 -420.73 -7.0562 -1.7514 2.1136 1.5887 3.4084 2.1458 1452 260.09

            >>> st=CH4(T=126, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            126 0.28667 400.52 4.7434 -859.94 -378.42 -6.2511 -2.4295 1.995 1.6591 3.6102 2.3638 1190.8 281.28

            >>> st=CH4(T=160, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            160 1.5921 336.31 25.382 -726.14 -353.87 -5.3391 -3.0124 1.9037 1.8473 4.4354 3.4189 795.43 283.01

            >>> st=CH4(T=162, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            162 1.7235 331.57 27.66 -717.27 -354.28 -5.2864 -3.0457 1.9028 1.8655 4.5448 3.5665 769.03 282.03

            >>> st=CH4(T=164, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            164 1.8626 326.64 30.132 -708.21 -355.01 -5.2334 -3.0797 1.9027 1.8852 4.6702 3.7374 742.1 280.91

            >>> st=CH4(T=166, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            166 2.0096 321.5 32.822 -698.95 -356.07 -5.18 -3.1145 1.9037 1.9066 4.8154 3.9373 714.59 279.65

            >>> st=CH4(T=168, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            168 2.1647 316.14 35.758 -689.46 -357.52 -5.1261 -3.1503 1.9059 1.93 4.9854 4.174 686.42 278.23
            
            >>> st=CH4(T=170, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            170 2.3283 310.5 38.974 -679.7 -359.4 -5.0715 -3.1874 1.9095 1.9556 5.1872 4.4585 657.52 276.66

            >>> st=CH4(T=172, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            172 2.5007 304.56 42.514 -669.64 -361.77 -5.0159 -3.226 1.9149 1.984 5.4311 4.8066 627.77 274.93

            >>> st=CH4(T=174, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            174 2.6822 298.25 46.434 -659.64 -361.77 -5.0159 -3.226 1.9225 2.0157 5.7318 5.2416 597.05 273.02

            >>> st=CH4(T=176, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            176 2.8732 291.5 50.808 -648.37 -368.31 -4.9009 -3.3096 1.933 2.0515 6.1124 5.8 565.18 270.92

            >>> st=CH4(T=178, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            178 3.074 284.21 55.74 -637.01 -372.71 -4.8406 -3.3558 1.9473 2.0925 6.6105 6.5421 531.94 268.6

            >>> st=CH4(T=180, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            180 3.2852 273.23 61.375 -625 -378.11 -4.7778 -3.4062 1.9669 2.1404 7.2923 7.574 497.01 266.04

            >>> st=CH4(T=182, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            182 3.5071 267.33 67.94 -612.14 -384.81 -4.7112 -3.4621 1.9944 2.1977 8.2863 9.103 459.94 263.17

            >>> st=CH4(T=184, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            184 3.7405 257.14 75.81 -598.08 -393.29 -4.6393 -3.5262 2.0346 2.2688 9.8808 11.592 420 259.89

            >>> st=CH4(T=186, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            186 3.986 244.93 85.704 -582.19 -404.47 -4.5586 -3.6032 2.0978 2.362 12.883 16.333 375.88 255.97
            
            >>> st=CH4(T=188, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            188 4.2448 228.93 99.377 -562.88 -420.58 -4.4612 -3.7044 2.213 2.5001 20.738 28.774 324.57 250.72
            
            >>> st=CH4(T=190, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            190 4.5186 200.78 125.18 -532.67 -451.91 -4.3082 -3.8831 2.6022 2.8546 94.012 140.81 250.31 238.55

            >>> st=CH4(T=190.564, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Gas.rho, st.Gas.h.kJkg, st.Gas.s.kJkgK)
            190.564 4.5992 162.66 -495.35 -4.1144
            """
            # Table 39, Pag 1117
            """
            >>> st=CH4(T=115, P=25000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            115 0.42278 -447.44 -388.31 -1.2861 1.5667 2.1012 280.5

            >>> st=CH4(T=620, P=25000)
            >>> print "%0.6g %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            620 0.07780 565.93 887.27 2.693 2.8267 3.3451 616.69

            >>> st=CH4(T=90.704, P=50000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            90.704 451.49 -982.78 -982.67 -7.3867 2.1677 3.3676 1538.8
            
            >>> st=CH4(T=230, P=50000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            230 0.42032 -267.02 -148.06 -0.19648 1.5947 2.1175 397.03
            
            >>> st=CH4(T=620, P=50000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            620 0.15559 565.86 887.22 2.3337 2.8267 3.3454 616.74

            >>> st=CH4(T=190, P=100000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            190 1.023 -331.02 -233.27 -0.96123 1.5707 2.104 360.51

            >>> st=CH4(T=190, P=2000000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            190 24.354 -362.79 -280.67 -2.6867 1.6906 2.7325 328.4

            >>> st=CH4(T=180, P=4000000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            180 285.08 -644.73 -630.7 -4.8236 1.9295 6.1588 552.74

            >>> st=CH4(T=190, P=4000000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            190 70.942 -423.54 -367.16 -3.4054 2.0278 6.887 279.41

            >>> st=CH4(T=150, P=4500000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 365.67 -778.92 -766.62 -5.6533 1.9241 3.8443 989.75

            >>> st=CH4(T=185, P=4500000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            185 266.01 -618.18 -601.26 -4.6722 1.9687 7.5262 477.08

            >>> st=CH4(T=190, P=4500000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            190 115.69 -477.41 -438.51 -3.8118 2.6453 55.55 246.05
            """, 
            
        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO",
        
        "Tmin": Tt, "Tmax": 625.0, "Pmax": 1000000.0, "rhomax": 40.072, 
        "Pmin": 11.696, "rhomin": 28.142, 

        "nr1": [0.43679010280e-1, 0.67092361990, -0.17655778590e01,
                0.85823302410, -0.12065130520e01, 0.51204672200,
                -0.40000107910e-3, -0.12478424230e-1, 0.31002697010e-1,
                0.17547485220e-2, -0.31719216050e-5, -0.22403468400e-5,
                0.29470561560e-6],
        "d1": [1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 8, 9, 10],
        "t1": [-0.5, 0.5, 1., 0.5, 1., 1.5, 4.5, 0., 1., 3., 1., 3., 3.],

        "nr2": [0.18304879090, 0.15118836790, -0.42893638770, 0.68940024460e-1,
                -0.14083139960e-1, -0.30630548300e-1, -0.29699067080e-1,
                -0.19320408310e-1, -0.11057399590, 0.99525489950e-1,
                0.85484378250e-2, -0.61505556620e-1, -0.42917924230e-1,
                -0.18132072900e-1, 0.34459047600e-1, -0.23859194500e-2,
                -0.11590949390e-1, 0.66416936020e-1, -0.23715495900e-1,
                -0.39616249050e-1, -0.13872920440e-1, 0.33894895990e-1,
                -0.29273787530e-2],
        "d2": [1, 1, 1, 2, 4, 5, 6, 1, 2, 3, 4, 4, 3, 5, 5, 8, 2, 3, 4, 4, 4, 5, 6],
        "t2": [0., 1., 2., 0., 0., 2., 2., 5., 5., 5., 2., 4., 12., 8., 10.,
               10., 10., 14., 12., 18., 22., 18., 14.],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*23,

        "nr3": [0.93247999460e-4, -0.62871715180e01,  0.12710694670e02,
                -0.64239534660e01],
        "d3": [2, 0, 0, 0],
        "t3": [2., 0., 1., 2.],
        "alfa3": [10, 40, 40, 40],
        "beta3": [200, 250, 250, 250],
        "gamma3": [1.07, 1.11, 1.11, 1.11],
        "epsilon3": [1]*4}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethane of Younglove and Ely (1987)",
        "__doc__":  u"""Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
        
        "Tmin": 90.68, "Tmax": 600.0, "Pmax": 200000.0, "rhomax": 36.2029, 
        "Pmin": 11.744, "rhomin": 28.147, 

        "R": 8.31434,
        "cp": CP4,
        "b": [None, 0.9898937956e-4, 0.2199608275, -0.5322788000e1,
              0.2021657962e3, -0.2234398926e5, 0.106794028e-3, 0.1457922469e-2,
              -0.9265816666e1, 0.2915364732e4, 0.2313546209e-5, 0.1387214274e-2,
              0.4780467451e-1, 0.1176103833e-3, -0.198209673e-2, -0.2512887756,
              0.9748899826e-4, -0.1202192137e-5, 0.4128353939e-3, -0.7215842918e-5,
              0.5081738255e4, -0.9198903192e6, -0.2732264677e2, 0.7499024351e6,
              0.1114060908e-1, 0.1083955159e2, -0.4490960312e-3, -0.1380337847e2,
              -0.2371902232e-6, 0.3761652197e-3, -0.2375166954e-8,
              -0.1237640790e-6, 0.6766926453e-5]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methane of Kunz and Wagner (2004).",
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": Fi2,

        "Tmin": 90.6941, "Tmax": 625.0, "Pmax": 1000000.0, "rhomax": 40.072, 
#        "Pmin": 73.476, "rhomin": 29.249, 

        "nr1":  [0.57335704239162, -0.16760687523730e1, 0.23405291834916,
                 -0.21947376343441, 0.16369201404128e-1, 0.15004406389280e-1],
        "d1": [1, 1, 2, 2, 4, 4],
        "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

        "nr2": [0.98990489492918e-1, 0.58382770929055, -0.7478686756039,
                0.30033302857974, 0.20985543806568, -0.18590151133061e-1,
                -0.15782558339049, 0.12716735220791, -0.32019743894346e-1,
                -0.68049729364536e-1, 0.24291412853736e-1, 0.51440451639444e-2,
                -0.19084949733532e-1, 0.55229677241291e-2, -0.44197392976085e-2,
                0.40061416708429e-1, -0.33752085907575e-1, -0.25127658213357e-2],
        "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
        "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
               14, 11.5, 26, 28, 30, 16],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
        "gamma2": [1]*18}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methane of Friend et al. (1989)",
        "__doc__":  u"""Friend, D.G., Ely, J.F., and Ingham, H., "Thermophysical Properties of Methane," J. Phys. Chem. Ref. Data, 18(2):583-638, 1989.""",
        "R": 8.31451,
        "cp": CP3,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 100000.0, "rhomax": 29.714, 
        "Pmin": 11.694, "rhomin": 28.145, 

        "nr1": [0.384436099659, -0.179692598800e1, 0.329444947369,
                0.226312728442e-1, 0.759236768798e-1, 0.693758447259e-1,
                0.241163263947e-1, 0.107009920854e-1, -0.380933275164e-1,
                0.471537561143e-3, 0.556607678805e-3, 0.548759346533e-6,
                -0.999632699967e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-0.128087979280, 0.380198873377e-1, 0.139226650551,
                -0.874996348859e-1, -0.334894165760e-2, -0.517576297122e-1,
                0.252835179116e-1, 0.518703205950e-3, -0.166770594525e-2,
                -0.607401927389e-3, -0.972915359991e-4, -0.298844010462e-4,
                -0.130940111124e-1, 0.198175833798e-1, 0.208465762327e-1,
                -0.358025052631e-1, -0.203486851741, 0.215964755088,
                -0.429340628249e-2],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*19}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for methane of Span and Wagner (2003)",
        "__doc__":  u"""Span, R., Wagner, W. Equations of state for technical applications. II. Results for nonpolar fluids. Int. J. Thermophys. 24 (2003), 41 – 109.""",
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 40.072, 
        "Pmin": 11.661, "rhomin": 28.167, 

        "nr1": [0.89269676, -0.25438282e1, 0.64980978, 0.20793471e-1,
                0.70189104e-1, 0.23700378e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.16653334, -0.43855669e-1, -0.1572678, -0.35311675e-1,
                -0.29570024e-1, 0.14019842e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1, MBWR, GERG, helmholtz3, helmholtz4

    _surface = {"sigma": [0.0308936, 0.0249105, -0.0068276],
                "exp": [1.25, 2.25, 3.25]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [6.5443, 0.0133], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [8.4578, 3.7196, -352.97, -100.65],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 11.696,
                "Tmin": Tt, "Tmax": 625.0,
                "a1": [1, 0.247568e5, -0.736602e4, -0.247568e5, 0.736602e4],
                "exp1": [0, 1.85, 2.1, 0, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 11.696,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-12.84], "exp2": [1],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 6,
        "ao": [-6.036219, 1.409353, -0.4945199, -1.443048],
        "exp": [2, 3, 4, 9]}
    _liquid_Density = {
        "eq": 3,
        "ao": [1.9906389, -0.78756197, 0.036976723],
        "exp": [0.354, 0.5, 2.5]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-1.880284, -2.8526531, -3.000648, -5.251169, -13.191859, -37.553961],
        "exp": [1.062, 2.5, 4.5, 7.5, 12.5, 23.5]}

    visco0 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quinones-Cisneros, M.L. Huber and U.K. Deiters, "Reference Correlation for the Viscosity of Methane," in preparation, for submission to J. Phys. Chem. Reference Data, 2007.""",
              "__name__": "Quiñones-Cisneros (2007)",
              "Tref": 190.564, "muref": 1.0,
              "ek": 174., "sigma": 0.36652, "n_chapman": 0,
              "n_ideal": [0.028790445329809258e3, -0.08883896490106571e3,
                          0.0854278871311819e3, -0.018038099301409677e3],
              "t_ideal": [0, 0.25, 0.5, 0.75],

              "a": [-3.388499774849180e-5, 1.357751056436950e-5, 0.0],
              "b": [2.757578135745610e-5, -3.437257168370160e-5, 0.0],
              "c": [2.891185964062900e-5, -9.980820268031560e-6, -2.073683069612580e-7],
              "A": [1.599928656708460e-8, -1.914111948640950e-9, 0.0],
              "B": [-2.591543952601510e-9, 3.260205780076830e-9, 0.0],
              "C": [-3.328762135404730e-8, 1.806469586530360e-7, 4.380464338668540e-10],
              "D": [3.151938769973220e-12, 0.0, 0.0]}

    visco1 = {"eq": 1, "omega": 1,
              "collision": [0.215309028, -0.46256942, 0.051313823,
                            0.030320660, -0.0070047029],
              "__name__": "Vogel (2000)",
              "__doc__": """Vogel, E., Wilhelm, J., Kuechenmeister, C., and Jaesche, M., "High-precision viscosity measurements on methane," High Temp. - High Pressures, 32(1):73-81, 2000.""",
              "ek": 160.78, "sigma": 0.37333,
              "n_chapman": 0.0855422/M**0.5,

              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 159.7, "etaref_virial": 0.0306525,

              "Tref_res": 190.564, "rhoref_res": 10.139*M, "etaref_res": 1,
              "n_packed": [3.10860501398],
              "t_packed": [0],
              "n_poly": [-3.02256904347, 17.6965130175, 3.11150846518,
                         -21.5685107769, 0.672852409238, 10.2387524315,
                         -1.09330775541, -1.20030749419, -21.1009923406],
              "t_poly": [0, -1, 0, -1, 0, -1, 0, -1, 0],
              "d_poly": [2, 2, 3, 3, 4, 4, 5, 5, 1],
              "g_poly": [0, 0, 0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "n_num": [21.1009923406],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    visco2 = {"eq": 2, "omega": 2,
              "__name__": "Younglove (1987)",
              "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
              "ek": 168., "sigma": 0.368,
              "n_chapman": 0.1069188/M**0.5,
              "F": [0.16969859271, -0.13337234608e-1, 0.140e1, 0.168e3],
              "E": [-0.1620427429e2, 0.4270589027e3, 0.1402596278e2,
                    -0.3916837745e4, -0.347709909e-1, 0.2136542674e2,
                    0.1436802482e4],
              "rhoc": 10.15}

    visco3 = {"eq": 1, "omega": 2,
              "__name__": "Friend (1989)",
              "__doc__": """Friend, D.G., Ely, J.F., and Ingham, H., "Tables for the Thermophysical Properties of Methane," NIST Technical Note 1325, 1989.""",
              "Tref": 174., "etaref": 10.0,
              "ek": 174., "sigma": 0.36652,
              "n_chapman": 0.14105376/M**0.5,

              "Tref_res": 190.551, "rhoref_res": 10.139*M, "etaref_res": 12.149,
              "n_num": [0.41250137, -0.14390912, 0.10366993, 0.40287464,
                        -0.24903524, -0.12953131, 0.06575776, 0.02566628,
                        -0.03716526],
              "t_num": [0, -1, 0, -1, -1.5, 0, -2, 0, -1],
              "d_num": [1, 1, 2, 2, 2, 3, 3, 4, 4],
              "g_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "c_num": [0, 0, 0, 0, 0, 0, 0, 0, 0],
              "n_den": [1.0, -0.38798341, 0.03533815],
              "t_den": [0, 0, -1.0],
              "d_den": [0, 1, 1],
              "g_den": [0, 0, 0],
              "c_den": [0, 0, 0]}

    _viscosity = visco0, visco1, visco2, visco3

    thermo0 = {"eq": 1, "critical": "CH4",
               "__name__": "Friend (1989)",
               "__doc__": """Friend, D.G., Ely, J.F., and Ingham, H., "Tables for the Thermophysical Properties of Methane," NIST Technical Note 1325, 1989.""",

               "Tref": 174., "kref": 1e-3,
               "no": [1.45885, -0.4377162, 0],
               "co": [0, -1, -96],

               "Trefb": 190.551, "rhorefb": 10.139, "krefb": 6.29638e-3,
               "nb": [1.5554612, 1., 2.4149207, 0.55166331, -0.52837734,
                      0.073809553, 0.24465507, -0.047613626],
               "tb": [0, 0, 0, 0, 0, -1, 0, -1],
               "db": [2, 0, 1, 3, 4, 4, 5, 5],
               "cb": [0, -99, 0, 0, 0, 0, 0, 0]}

    thermo1 = {"eq": 2, "omega": 2,
               "__name__": "Younglove (1987)",
               "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
               "visco": visco2,
               "n_chapman": 0.1069188,
               "G": [0.1346953698e1, -0.3254677753],
               "E": [0.2325800819e-2, -0.2477927999, 0.3880593713e2,
                     -0.1579519146e-6, 0.3717991328e-2, -0.9616989434,
                     -0.3017352774e-1, 0.4298153386],

               "critical": 2,
               "X": [37.42368, 3.16714, 0.78035, 0.60103],
               "Z": 6.512707e-10}

    _thermal = thermo0, thermo1
