#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Ar(MEoS):
    """Multiparamter equation of state for argon"""
    name = "argon"
    CASNumber = "7440-37-1"
    formula = "Ar"
    synonym = "R-740"
    rhoc = unidades.Density(535.6)
    Tc = unidades.Temperature(150.687)
    Pc = unidades.Pressure(4863, "kPa")
    M = 39.948  # g/mol
    Tt = unidades.Temperature(83.8058)
    Tb = unidades.Temperature(87.302)
    f_acent = -0.00219
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 98
    _Tr = unidades.Temperature(147.707801)
    _rhor = unidades.Density(540.014968)
    _w = 0.000305675

    Fi1 = {"ao_log": [1, 1.5],
           "pow": [0, 1],
           "ao_pow": [8.31666243, -4.94651164],
           "ao_exp": [], "titao": []}

    CP1 = {"ao": 2.5,
           "an": [], "pow": [], "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}
           
    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "FEQ Helmholtz equation of state for argon of Tegeler et al. (1999).",
        "__doi__": {"autor": "Tegeler, Ch., Span, R., Wagner, W.",
                    "title": "A New Equation of State for Argon Covering the Fluid Region for Temperatures From the Melting Line to 700 K at Pressures up to 1000 MPa", 
                    "ref": "J. Phys. Chem. Ref. Data 28, 779 (1999)",
                    "doi": "10.1063/1.556037"}, 
        "__test__": 
            #Table 33, Pag 828
            """
            >>> st=Ar(T=83.8058, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            83.8058 0.068891 1416.77 4.0546 -276.56 -112.85 -2.544 -0.59044 0.5496 0.32471 1.1157 0.55503 862.43 168.12

            >>> st=Ar(T=90, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            90 0.13351 1378.63 7.4362 -269.61 -110.55 -2.4645 -0.69718 0.52677 0.33094 1.1212 0.57569 819.45 172.83

            >>> st=Ar(T=120, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            120 1.213 1162.82 60.144 -233.48 -106.71 -2.1274 -1.071 0.45763 0.38934 1.3324 0.86265 584.19 185.09

            >>> st=Ar(T=130, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            130 2.0255 1068.13 103.56 -219.29 -109.83 -2.0197 -1.1777 0.4492 0.42745 1.5638 1.1717 487.88 184.85

            >>> st=Ar(T=140, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            140 3.1682 943.71 178.86 -202.29 -117.65 -1.9023 -1.2978 0.45984 0.49404 2.2247 2.1036 371.63 181.5

            >>> st=Ar(T=142, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            142 3.4435 911.61 201.37 -198.23 -120.25 -1.8756 -1.3265 0.46729 0.51706 2.5349 2.5648 344.14 180

            >>> st=Ar(T=144, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            144 3.7363 874.98 228.48 -193.76 -123.46 -1.8467 -1.3584 0.47972 0.54719 3.0262 3.3149 313.8 177.93

            >>> st=Ar(T=146, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            146 4.0479 831.38 262.63 -188.68 -127.57 -1.8142 -1.3956 0.50257 0.58923 3.9312 4.7346 278.88 174.89

            >>> st=Ar(T=148, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            148 4.3797 775.03 309.6 -182.49 -133.29 -1.7749 -1.4424 0.55094 0.6568 6.2097 8.383 236.08 169.81
            
            >>> st=Ar(T=150, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            150 4.7346 680.43 394.5 -173.01 -143.6 -1.7145 -1.5185 0.70603 0.82182 23.582 35.468 174.74 157.01
            """
            #Table 33, Pag 828
            """
            >>> st=Ar(T=83.814, P=1e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            83.814 1416.8 -276.61 -276.54 -2.544 0.54961 1.1156 862.52
            >>> st=Ar(T=700, P=1e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            700 0.68619 63.355 209.09 0.44677 0.31223 0.5205 492.95
            >>> st=Ar(T=150, P=5e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 16.605 -110.45 -80.334 -0.70404 0.32098 0.55987 224.97
            >>> st=Ar(T=150, P=5e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 16.605 -110.45 -80.334 -0.70404 0.32098 0.55987 224.97
            >>> st=Ar(T=170, P=1e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            170 29.723 -105.62 -71.972 -0.78987 0.32356 0.57801 238.88
            >>> st=Ar(T=125, P=2e6)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            125 1122.34 -228.45 -226.66 -2.0773 0.45179 1.4048 544.65
            >>> st=Ar(T=135, P=3e6)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            135 1020.52 -214.63 -211.69 -1.9694 0.44845 1.7159 445.83
            >>> st=Ar(T=150, P=3e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 128.93 -124.64 -101.37 -1.1772 0.39203 1.0311 205.67
            >>> st=Ar(T=140, P=4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            140 968.76 -207.81 -203.68 -1.9185 0.45035 1.9268 403.8
            >>> st=Ar(T=145, P=4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            145 862.44 -196.53 -191.89 -1.8358 0.48302 3.1513 306.38
            >>> st=Ar(T=150, P=4e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 209.45 -134.47 -115.38 -1.3116 0.46106 1.8982 193.39
            >>> st=Ar(T=125, P=5e6)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            125 1150.27 -231.07 -226.72 -2.0989 0.45236 1.2969 586.37
            >>> st=Ar(T=150, P=5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 765.37 -186.33 -179.79 -1.7622 0.52622 5.1511 248.19
            >>> st=Ar(T=150, P=1e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            150 964.88 -203.07 -192.7 -1.8855 0.43203 1.5594 445.1
            """, 
            
        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 2000., "Pmax": 1000000.0, "rhomax": 50.65, 
        "Pmin": 68.891, "rhomin": 35.465, 

        "nr1": [0.887223049900e-1, 0.705148051673, -0.168201156541e1,
                 -0.149090144315, -0.120248046009, -0.121649787986,
                 0.400359336268, -0.271360626991, 0.242119245796,
                 0.578895831856e-2, -0.410973356153e-1, 0.247107615416e-1],
        "d1": [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4],
        "t1": [0., 0.25, 1., 2.75, 4.0, 0., 0.25, 0.75, 2.75, 0.0, 2.0, 0.75],

        "nr2": [-0.321813917507, 0.332300176958, 0.310199862873e-1,
                -0.307770860024e-1, 0.938911374196e-1, -0.906432106820e-1,
                -0.457783492767e-3, -0.826597290252e-4, 0.130134156031e-3,
                -0.113978400020e-1, -0.244551699605e-1, -0.643240671760e-1,
                0.588894710937e-1, -0.649335521130e-3, -0.138898621584e-1,
                0.404898392969, -0.386125195947, -0.188171423322,
                0.159776475965, 0.539855185139e-1, -0.289534179580e-1,
                -0.130254133814e-1, 0.289486967758e-2, -0.226471343048e-2,
                0.176164561964e-2],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
               3, 3, 4, 4],
        "d2": [1, 1, 3, 4, 4, 5, 7, 10, 10, 2, 2, 4, 4, 8, 3, 5, 5, 6, 6, 7, 7,
               8, 9, 5, 6],
        "t2": [3., 3.5, 1., 2., 4., 3., 0., 0.5, 1., 1., 7., 5., 6., 6., 10.,
               13., 14., 11., 14., 8., 14., 6., 7., 24., 22.],
        "gamma2": [1]*25,

        "nr3": [0.585524544828e-2, -0.6925190827, 0.153154900305e1,
                -0.273804474498e-2],
        "d3": [2, 1, 2, 3],
        "t3": [3, 1, 0, 0],
        "alfa3": [20]*4,
        "beta3": [250, 375, 300, 225],
        "gamma3": [1.11, 1.14, 1.17, 1.11],
        "epsilon3": [1, 1, 1, 1],
        "nr4": []}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for argon of Kunz and Wagner (2004).",
        "__doc__": u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 700., "Pmax": 1000000.0, "rhomax": 50.65, 
        "Pmin": 68.891, "rhomin": 35.465, 

        "nr1": [0.85095714803969, -0.24003222943480e1, 0.54127841476466,
                0.16919770692538e-1, 0.68825965019035e-1, 0.21428032815338e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.17429895321992, -0.33654495604194e-1, -0.13526799857691,
                -0.16387350791552e-1, -0.24987666851475e-1, 0.88769204815709e-2],
        "c2": [1, 1, 2, 2, 3, 3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for argon of Stewart and Jacobsen (1989).",
        "__doi__": {"autor": "Stewart, R.B. and Jacobsen, R.T.",
                    "title": "Thermodynamic Properties of Argon from the Triple Point to 1200 K at Pressures to 1000 MPa", 
                    "ref": "J. Phys. Chem. Ref. Data, 18(2):639-798, 1989",
                    "doi": "10.1063/1.555829"}, 
        "__test__": 
            #Table 14, Pag 379
            """
            >>> st=Ar(T=83.804, x=0.5, eq=3)
            >>> print "%0.6g %0.4g %0.5g %0.5g %0.5g %0.5g %0.4g %0.5g %0.4g %0.4g %0.3g %0.3g" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, \
                st.Liquido.cpM.JmolK, st.Liquido.w, st.Gas.w)
            83.804 0.06895 35.475 0.10152 -4835.9 1701.4 53.29 131.3 21.34 42.61 853 208
            >>> st=Ar(T=90, x=0.5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.5g %0.4g %0.4g %0.3g %0.3g" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, \
                st.Liquido.cpM.JmolK, st.Liquido.w, st.Gas.w)
            90 0.13362 34.538 0.18651 -4568.1 1777.5 56.35 126.86 20.59 43.49 811 186
            >>> st=Ar(T=100, x=0.5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.5g %0.4g %0.4g %0.3g %0.3g" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, \
                st.Liquido.cpM.JmolK, st.Liquido.w, st.Gas.w)
            100 0.32401 32.918 0.42327 -4120.9 1876.3 61 120.98 19.55 45.48 742 180
            >>> st=Ar(T=110, x=0.5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.4g %0.5g %0.4g %0.4g %0.3g %0.3g" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, \
                st.Liquido.cpM.JmolK, st.Liquido.w, st.Gas.w)
            110 0.66575 31.133 0.83561 -3648.1 1931.3 65.41 116.13 18.75 48.45 668 181
            >>> st=Ar(T=120, x=0.5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.0f %0.0f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, st.Gas.cvM.JmolK,\
                st.Liquido.cpM.JmolK, st.Gas.cpM.JmolK, st.Liquido.w, st.Gas.w)
            120 1.2139 29.123 1.5090 -3138.4 1917.7 69.68 111.81 18.16 16.75 53.27 36.15 586 182
            >>> st=Ar(T=130, x=0.5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.0f %0.0f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, st.Gas.cvM.JmolK,\
                st.Liquido.cpM.JmolK, st.Gas.cpM.JmolK, st.Liquido.w, st.Gas.w)
            130 2.027 26.748 2.597 -2570 1798.1 73.99 107.59 17.88 18.05 62.79 48.3 490 182
            >>> st=Ar(T=140, x=0.5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.0f %0.0f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, st.Gas.cvM.JmolK,\
                st.Liquido.cpM.JmolK, st.Gas.cpM.JmolK, st.Liquido.w, st.Gas.w)
            140 3.1704 23.59 4.4877 -1883.5 1489.2 78.74 102.83 18.44 20.15 90.97 84.59 368 179
            >>> st=Ar(T=150, x=0.5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.1f %0.1f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.0f %0.0f" % (\
                st.T, st.P.MPa, st.Liquido.rhoM, st.Gas.rhoM, st.Liquido.hM.Jmol, st.Gas.hM.Jmol, \
                st.Liquido.sM.JmolK, st.Gas.sM.JmolK, st.Liquido.cvM.JmolK, st.Gas.cvM.JmolK,\
                st.Liquido.cpM.JmolK, st.Gas.cpM.JmolK, st.Liquido.w, st.Gas.w)
            150 4.7363 16.973 9.7709 -707.3 487.7 86.28 94.25 23.74 26 762.2 1098.27 198 171
            """
            #Table 15, Pag 684
            """
            >>> st=Ar(T=84, P=8e4, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            84 35.447 -4829.6 -4827.3 53.39 21.31 42.64 852
            >>> st=Ar(T=300, P=8e4, eq=3)
            >>> print "%0.6g %0.5f %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            300 0.03209 3736.4 6229.5 156.81 12.48 20.82 323
            >>> st=Ar(T=1200, P=8e4, eq=3)
            >>> print "%0.6g %0.5f %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            1200 0.00802 14965 24944 185.64 12.47 20.79 645
            >>> st=Ar(T=84, P=1e5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            84 35.448 -4829.8 -4827 53.39 21.31 42.63 852
            >>> st=Ar(T=150, P=1.5e5, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            150 0.12154 1845.3 3079.4 137.02 12.58 21.23 227
            >>> st=Ar(T=116, P=1e6, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            116 29.967 -3380.9 -3347.5 67.97 18.37 50.98 621
            >>> st=Ar(T=150, P=2e6, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            150 1.8988 1482.8 2536.2 112.99 14.26 30.07 214
            >>> st=Ar(T=138, P=3e6, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            138 24.393 -2164.9 -2041.9 77.65 18.16 80.50 400
            >>> st=Ar(T=150, P=3e6, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            150 3.2266 1217.3 2147.1 107.70 15.78 41.36 205
            >>> st=Ar(T=144, P=4e6, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.4g %0.4g %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            144 22.202 -1757.7 -1577.5 80.64 18.97 110.7 326
            >>> st=Ar(T=150, P=4e6, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.4g %0.4g %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            150 5.2481 822.44 1584.6 102.31 18.61 76.39 193
            >>> st=Ar(T=150, P=5e6, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.4g %0.4g %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            150 18.975 -1229.6 -966.11 84.46 21.06 209.4 245
            >>> st=Ar(T=100, P=1e7, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.4g %0.4g %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            100 33.825 -4265.2 -3969.5 59.62 19.92 42.91 802
            >>> st=Ar(T=1200, P=1e7, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.4g %0.4g %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            1200 0.98207 14899 25082 145.44 12.53 20.95 660
            >>> st=Ar(T=150, P=1e8, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.4g %0.4g %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            150 34.219 -3331.3 -408.95 67.05 19.55 36.23 956
            >>> st=Ar(T=450, P=1e9, eq=3)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.4g %0.4g %0.0f" % (\
                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
            450 41.651 2711.6 26720 80.63 20.67 28.69 1815
            """, 

        "R": 8.31434,
        "cp": CP1,
        "ref": {"Tref": 298.15, "Pref": 101.325, "ho": 6197, "so": 154.732}, 
        "Tc": 150.6633, "Pc": 4860, "rhoc": 13.29, "Tt": 83.804, 
        
        "Tmin": 83.804, "Tmax": 1200., "Pmax": 1000000.0, "rhomax": 45.814, 
        "Pmin": 68.961, "rhomin": 35.475, 

        "nr1": [0.7918675715, -0.1633346151e1, -0.439530293, 0.1033899999,
                0.2061801664, -0.2888681776, 0.439801055, -0.8429550391e-1,
                -0.2155658654, 0.4786509099, -0.3525884593, 0.3015073692e-1,
                0.2987679059e-1, -0.1522568583e-1, 0.7435785786e-3],
        "d1": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 6],
        "t1": [0.25, 1, 3, 4, 0.25, 1, 2.5, 3.5, 0.75, 1, 1.5, 2.5, 1, 2, 2],

        "nr2": [0.7099541624e-1, -0.2904237185e-1, -0.6223078525e-1,
                0.1410895187e-3, -0.1481241783e-2, 0.3023342784e-1,
                -0.6126784685e-1, 0.270996709e-1, 0.9411034405e-1,
                -0.7291645114e-2, -0.1586314976e-2, 0.9510948813e-3,
                0.7786181844e-3],
        "c2": [3, 3, 2, 4, 6, 3, 3, 3, 2, 2, 4, 2, 2],
        "d2": [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 8, 8],
        "t2": [5, 7, 5, 22, 16, 10, 14, 16, 4, 8, 10, 5, 6],
        "gamma2": [1]*13,

        "nr3": [],
        "nr4": []}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for argon of Span and Wagner (2003).",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.", 
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"}, 
        "__test__": """
            >>> st=Ar(T=700, rho=200, eq=4)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            0.5203 31.922 0.5630
            >>> st2=Ar(T=750, rho=100, eq=4)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            25.97 0.18479
            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 750., "Pmax": 100000.0, "rhomax": 50.65, 
        "Pmin": 69.026, "rhomin": 35.498, 

        "nr1": [0.85095715, -0.24003223e1, 0.54127841, 0.16919771e-1,
                0.68825965e-1, 0.21428033e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.17429895, -0.33654496e-1, -0.135268, -0.16387351e-1,
                -0.24987667e-1, 0.88769205e-2],
        "c2": [1, 1, 2, 2, 3, 3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "BWR  MBWR equation of state for argon of Younglove (1982).",
        "__doc__": """Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
        "R": 8.31434,
        "cp": CP1,

        "Tmin": 83.80, "Tmax": 400., "Pmax": 101000.0, "rhomax": 50.65, 
        "Pmin": 68.906, "rhomin": 35.4, 

        "b": [None, -0.6569731294e-3, 0.1822957801, -0.3649470141e1,
              0.1232012107e3, -0.8613578274e4, 0.7978579691e-4, -0.2911489110e-1,
              0.7581821758e1, 0.8780488169e4, 0.1423145989e-6, 0.1674146131e-2,
              -0.3200447909, 0.2561766372e-4, -0.5475934941e-3, -0.4505032058,
              0.2013254653e-4, -0.1678941273e-6, 0.4207329271e-3, -0.5444212996e-5,
              -0.8004855011e4, -0.1319304201e6, -0.4954923930e2, 0.8092132177e5,
              -0.9870104061e-1, 0.2020441562e1, -0.1637417205e-3, -0.7038944136,
              -0.1154324539e-6, 0.1555990117e-4, -0.1492178536e-9,
              -0.1001356071e-7, 0.2933963216e-6]}

    eq = helmholtz1, MBWR, GERG, helmholtz3, helmholtz4
    _PR = -0.0034

    _dielectric = {
        "eq": 3, "Tref": 273.16, "rhoref": 1000.,
        "a0": [],  "expt0": [], "expd0": [],
        "a1": [4.1414], "expt1": [0], "expd1": [1],
        "a2": [1.597, 0.262, -117.9], "expt2": [0, 1, 0], "expd2": [2, 2, 3.1]}
    _melting = {
        "eq": 1, "Tref": Tt, "Pref": 68.891,
        "Tmin": 83.8058, "Tmax": 700.0, 
        "a1": [1, -7476.26651, 9959.06125, 7476.26651, -9959.06125],
        "exp1": [0, 1.05, 1.275, 0, 0],
        "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {
        "eq": 3, "Tref": Tt, "Pref": 68.891,
        "Tmin": 83.8058, "Tmax": 83.8058, 
        "a1": [], "exp1": [],
        "a2": [-11.1307], "exp2": [1],
        "a3": [], "exp3": []}
    _surface = {"sigma": [0.037898063], "exp": [1.278]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-5.9409785, 1.3553888, -0.4649761, -1.5399043],
        "exp": [1., 1.5, 2., 4.5]}
    _liquid_Density = {
        "eq": 3,
        "ao": [1.5004264, -0.3138129, 0.086461622, -0.041477525],
        "exp": [0.334, 2./3, 7./3, 4]}
    _vapor_Density = {
        "eq": 5,
        "ao": [-0.29182e1, 0.97930e-1, -0.13721e1, -0.22898e1],
        "exp": [0.72, 1.25, 0.32, 4.34]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Lemmon (2004)",
               "__doi__": {"autor": "Lemmon, E.W. and Jacobsen, R.T.",
                            "title": "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air", 
                            "ref": "Int. J. Thermophys., 25:21-69, 2004.",
                            "doi": "10.1023/B:IJOT.0000022327.04529.f3"}, 
               "__test__": """
                    >>> st=Ar(T=100, rhom=0)
                    >>> print "%0.5f" % st.mu.muPas
                    8.18940
                    >>> st=Ar(T=300, rhom=0)
                    >>> print "%0.4f" % st.mu.muPas
                    22.7241
                    >>> st=Ar(T=100, rhom=33)
                    >>> print "%0.3f" % st.mu.muPas
                    184.232
                    >>> st=Ar(T=200, rhom=10)
                    >>> print "%0.4f" % st.mu.muPas
                    25.5662
                    >>> st=Ar(T=300, rhom=5)
                    >>> print "%0.4f" % st.mu.muPas
                    26.3706
                    >>> st=Ar(T=150.69, rhom=13.4)
                    >>> print "%0.4f" % st.mu.muPas
                    27.6101
                    """, # Table V, Pag 28

              "ek": 143.2, "sigma": 0.335,
              "n_poly": [12.19, 13.99, 0.005027, -18.93, -6.698, -3.827],
              "t_poly": [0.42, 0.0, 0.95, 0.5, 0.9, 0.8],
              "d_poly": [1, 2, 10, 5, 1, 2],
              "g_poly": [0, 0, 0, 1, 1, 1],
              "c_poly": [0, 0, 0, 2, 4, 4]}

    visco1 = {"eq": 3,
              "__doc__": """Younglove, B.A. and Hanley, H.J.M., The Viscosity and Thermal Conductivity Coefficients of Gaseous and Liquid Argon," J. Phys. Chem. Ref. Data, 15(4):1323-1337, 1986.""",
              "__name__": "Younglove (1986)",
              "Tref": 1, "muref": 1.0,
              "n_poly": [-0.8973188257e5, 0.8259113473e5, -0.2766475915e5,
                         0.3068539784e4, 0.4553103615e3, -0.1793443839e3,
                         0.2272225106e2, -0.1350672796e1, 0.3183693230e-1],
              "t_poly": [-1., -2./3, -1./3, 0, 1./3, 2./3, 1., 4./3, 5./3],
              "n_num": [0.5927733783, -0.4251221169e2, -0.2698477165e-1,
                        0.3727762288e2, -0.3958508720e4, 0.3636730841e-2,
                        -0.2633471347e1, 0.2936563322e3, -0.3811869019e-4,
                        0.4451947464e-1, -0.5385874487e1, 1, -0.1115054926e-1, 
                        -0.1328893444e1],
              "t_num": [0, -1, 0, -1, -2, 0, -1, -2, 0, -1, -2, 0, 0, -1],
              "d_num": [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 0, 1, 1],
              "n_den": [1.0, -0.1115054926e-1, -0.1328893444e1],
              "t_den": [0, 0, -1],
              "d_den": [0, 1, 1]}
              
    visco2 = {"eq": 2, "omega": 2,
              "collision": [25.7830291943396, -234.320222858983, 814.636688705024,
                            -1452.04353466585, 1467.17535558104, -870.164951237067,
                            313.024934147423, -61.2072628957372, 5.07700488990665],
              "__name__": "Younglove (1982)",
              "__doc__": """Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
              "ek": 152.8, "sigma": 0.3297,
              "n_chapman": 0.16871158559818,
              "t_chapman": 0,
              "F": [5.85384107393e-3, -3.09546765250e-3, 1.4, 152.8],
              "E": [-12.313579086, 40.136071933, 11.6160872385243,
                    -413.04094973717, 4.13624595833e-2, 7.96883967907912,
                    234.196850483958],
              "rhoc": 13.4424752177831}

    visco3 = {"eq": 1, "omega": 1,
              "__doc__": """Lemmon, E.W. and Jacobsen, R.T, unpublished equation, 2001""",
              "__name__": "Lemmon (2001)",
              "collision": [0.5136, -0.5218, 0.8852e-1, 0.3445e-2, -0.2289e-2],
              "ek": 104.9, "sigma": 0.3478,
              "n_poly": [0.648491572951e1, 0.921829714883e1, 0.251873627628e1,
                         0.710198542697e-1, 0.233411121182, -0.365386072189e1],
              "t_poly": [0, 0, 1.1303, 0.2879, 4.4346, 2.2411],
              "d_poly": [1, 2, 4, 7, 4, 3],
              "g_poly": [0, 0, 0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0, 0, 0]}

    _viscosity = visco0, visco1, visco2, visco3

    thermo0 = {"eq": 1,
               "__name__": "Lemmon (2004)",
               "__doi__": {"autor": "Lemmon, E.W. and Jacobsen, R.T.",
                            "title": "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air", 
                            "ref": "Int. J. Thermophys., 25:21-69, 2004.",
                            "doi": "10.1023/B:IJOT.0000022327.04529.f3"}, 
               "__test__": """
                    >>> st=Ar(T=100, rhom=0)
                    >>> print "%0.5f" % st.k.mWmK
                    6.36587
                    >>> st=Ar(T=300, rhom=0)
                    >>> print "%0.4f" % st.k.mWmK
                    17.8042
                    >>> st=Ar(T=100, rhom=33)
                    >>> print "%0.3f" % st.k.mWmK
                    111.266
                    >>> st=Ar(T=200, rhom=10)
                    >>> print "%0.4f" % st.k.mWmK
                    26.1377
                    >>> st=Ar(T=300, rhom=5)
                    >>> print "%0.4f" % st.k.mWmK
                    23.2302
                    >>> st=Ar(T=150.69, rhom=13.4)
                    >>> print "%0.4f" % st.k.mWmK
                    856.793
                    """, # Table V, Pag 28

               "Tref": 150.687, "kref": 1e-3,
               "no": [0.8158, -0.432],
               "co": [-97, -0.77],

               "Trefb": 150.687, "rhorefb": 13.40742965, "krefb": 1e-3,
               "nb": [13.73, 10.07, 0.7375, -33.96, 20.47, -2.274, -3.973],
               "tb": [0.0, 0.0, 0.0, 0.8, 1.2, 0.8, 0.5],
               "db": [1, 2, 4, 5, 6, 9, 1],
               "cb": [0, 0, 0, 2, 2, 2, 4],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.13e-9, "gam0": 0.55e-1, "qd": 0.32e-9, "Tcref": 301.374}

    thermo1 = {"eq": 3,
               "__name__": "Younglove (1982)",
               "__doc__": """Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",

               "ek": 152.8, "sigma": 0.3297,
               "Nchapman": 0.16871158559818,
               "tchapman": 0,
               "b": [2.64712871543e-2, -.216629583011974, 0.709700888884514,
                     -1.21908891344223, 1.20168985706305, -.700084760049098,
                     0.24816605762696, -4.79479287295e-2, 3.93679190444e-3],
               "F": [9.64428741429e-4, 3.02391316601e-4, 1, 152.8],
               "E": [-33.327027332, -355.59415848, 22.2441164817987,
                     1663.62775376509, 0, 0, 0],
               "rhoc": 25.0325423049965,
               "ff": 1.7124,
               "rm": 0.00000003669}

    thermo2 = {"eq": 1,
               "__name__": "Perkins (1991)",
               "__doc__": """Perkins, R.A., Friend, D.G., Roder, H.M., and Nieto de Castro, C.A., "Thermal Conductivity Surface of Argon:  A Fresh Analysis," Int. J. Thermophys., 12(6):965-984, 1991""",

               "Tref": 1.0, "kref": 1e-3,
               "no": [.1225067272e5, -.9096222831e4, .2744958263e4,
                      -.4170419051e3, .2527591169e2, .1604421067e1,
                      -.2618841031, .1381696924e-1, -.2463115922e-3],
               "co": [-1, -2./3, -1./3, 0, 1./3, 2./3, 1., 4./3, 5./3],

               "Trefb": 1., "rhorefb": 1., "krefb": 1.,
               "nb": [0.757894e-3, 0.612624e-4, -0.205353e-5,  0.745621e-7],
               "tb": [0, 0, 0, 0],
               "db": [1, 2, 3, 4],
               "cb": [0, 0, 0, 0],

               "critical": 4,
               "Tcref": 150.86, "Pcref": 4905.8, "rhocref": 13.41, "kcref": 1e-3,
               "gamma": 1.02,
               "expo": 0.46807, "alfa": 39.8, "beta": 5.45, "Xio": 6.0795e-1}

    _thermal = thermo0, thermo1, thermo2


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    for eq in (0, 2, 3, 4):
        st=Ar(T=300, P=1e6, eq=eq)
        print "%0.6g %0.5g %0.1f %0.3f %0.3f %0.3f %0.3f %0.2f" % (\
            st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)

#            >>> st=Ar(T=84, P=8e4, eq=3)
#            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.2f %0.2f %0.0f" % (\
#                st.T, st.rhoM, st.uM.Jmol, st.hM.Jmol, st.sM.JmolK, st.cvM.JmolK, st.cpM.JmolK, st.w)
#            84 35.447 -4829.6 -4827.3 53.39 21.31 42.64 852
