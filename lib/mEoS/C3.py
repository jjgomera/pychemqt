#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class C3(MEoS):
    """Multiparameter equation of state for propane"""
    name = "propane"
    CASNumber = "74-98-6"
    formula = "CH3CH2CH3"
    synonym = "R-290"
    rhoc = unidades.Density(220.4781)
    Tc = unidades.Temperature(369.89)
    Pc = unidades.Pressure(4251.2, "kPa")
    M = 44.09562  # g/mol
    Tt = unidades.Temperature(85.525)
    Tb = unidades.Temperature(231.036)
    f_acent = 0.1521
    momentoDipolar = unidades.DipoleMoment(0.084, "Debye")
    id = 4
    _Tr = unidades.Temperature(354.964211)
    _rhor = unidades.Density(221.906745)
    _w = 0.149041513

    Fi1 = {"ao_log": [1, 3],
           "pow": [0, 1],
           "ao_pow": [-4.970583, 4.29352],
           "ao_exp": [3.043, 5.874, 9.337, 7.922],
           "titao": [393/Tc, 1237/Tc, 1984/Tc, 4351/Tc], 
           "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, 3.02939],
           "pow": [0, 1],
           "ao_pow": [31.602908195, -84.463284382],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [6.60569, 3.197, 19.1921, -8.37267],
           "hyp": [1.297521801, 0.543210978, 2.583146083, 2.777773271]}

    CP2 = {"ao": 4.02256195,
           "an": [], "pow": [],
           "ao_exp": [2.90591124, 4.68495401, 10.2971154, 8.08977905],
           "exp": [388.87291, 1145.03868, 1880.40472, 4228.18881],
           "ao_hyp": [], "hyp": []}

    CP3 = {"ao": 4.02939,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [6.60569, 3.197, 19.1921, -8.37267],
           "hyp": [1.297521801*Tc, 0.543210978*Tc, 2.583146083*Tc, 2.777773271*Tc]}

    CP4 = {"ao": 4.021394,
           "an": [], "pow": [],
           "ao_exp": [2.889980, 4.474243, 1.048251e1, 8.139803],
           "exp": [387.69088, 1129.1386, 1864.95906, 4224.43701],
           "ao_hyp": [], "hyp": []}

    CP5 = {"ao": -5.4041204338,
           "an": [3.1252450099e6, -1.1415253638e5, 1.4971650720e3,
                  3.9215452897e-2, -2.1738913926e-5, 4.8274541303e-9],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [3.1907016349], "exp": [1500],
           "ao_hyp": [], "hyp": []}

    CP6 = {"ao": 4.02939,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.1521038e7, 0.1290245e6, 0.1751511e8, -0.8835886e7],
           "hyp": [0.479856e3, 0.2008930e3, 0.9553120e4, 0.1027290e4]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Lemmon et al. (2009)",
        "__doi__": {"autor": "Lemmon, E.W., McLinden, M.O., Wagner, W.",
                    "title": "Thermodynamic Properties of Propane.  III.  A Reference Equation of State for Temperatures from the Melting Line to 650 K and Pressures up to 1000 MPa", 
                    "ref": "J. Chem. Eng. Data, 2009, 54 (12), pp 3141–3180",
                    "doi": "10.1021/je900217v"}, 
        "__test__":
            # Table 5, Pag AH
            """
            >>> st=C3(T=200, rhom=14)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % (\
                st.T, st.rhoM, st.P.MPa, st.cvM.JmolK, st.cpM.JmolK, st.w)
            200.0 14.0 2.3795138 61.078424 93.475362 1381.9552
            >>> st=C3(T=300, rhom=12)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % (\
                st.T, st.rhoM, st.P.MPa, st.cvM.JmolK, st.cpM.JmolK, st.w)
            300.0 12.0 19.053797 73.972542 108.61529 958.4052
            >>> st=C3(T=300, rhom=0.4)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % (\
                st.T, st.rhoM, st.P.MPa, st.cvM.JmolK, st.cpM.JmolK, st.w)
            300.0 0.4 0.84694991 69.021875 85.753997 221.88959
            >>> st=C3(T=400, rhom=5)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % (\
                st.T, st.rhoM, st.P.MPa, st.cvM.JmolK, st.cpM.JmolK, st.w)
            400.0 5.0 6.646284 97.017439 271.07044 194.65847
            >>> st=C3(T=369.9, rhom=5)
            >>> print "%0.1f %0.1f %0.8g %0.8g %0.8g %0.8g" % (\
                st.T, st.rhoM, st.P.MPa, st.cvM.JmolK, st.cpM.JmolK, st.w)
            369.9 5.0 4.2519399 117.71621 753625 130.898
            """
            # Table A1, Pag A1
            """
            >>> st=C3(T=C3.Tt, x=0.5)
            >>> print "%0.6g %0.4e %0.2f %0.2e %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -187.625 1.7203e-10 733.13 1.07e-08 -196.64 366.26 -1.396 5.186 1.355 0.6907 1.916 0.8792 2136.4 143.3
            >>> st=C3(T=-110+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5f %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -110 0.0011644 654.05 0.03790 -43.988 445.38 -0.1298 2.87 1.343 0.9614 2.032 1.151 1611.8 191.7
            >>> st=C3(T=-50+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.1f %0.4g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            -50 0.070569 589.9 1.727 82.753 516.48 0.5298 2.473 1.428 1.182 2.212 1.397 1212.5 216.5
            >>> st=C3(T=273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            0 0.47446 528.59 10.351 200 574.87 1 2.372 1.572 1.427 2.493 1.739 885.5 221.3
            >>> st=C3(T=20+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            20 0.83646 500.06 18.082 251.64 595.95 1.18 2.354 1.647 1.544 2.666 1.949 752.9 217.2
            >>> st=C3(T=40+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            40 1.3694 467.46 30.165 307.15 614.21 1.359 2.34 1.732 1.678 2.913 2.263 617.0 208.6
            >>> st=C3(T=60+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            60 2.1168 427.97 49.493 368.14 627.36 1.543 2.321 1.832 1.836 3.337 2.841 474.2 194.1
            >>> st=C3(T=80+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            80 3.1319 373.29 84.406 438.93 628.73 1.742 2.279 1.969 2.057 4.545 4.707 314.9 171.6
            >>> st=C3(T=85+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            85 3.4361 353.96 98.818 459.81 624.75 1.798 2.259 2.023 2.144 5.433 6.182 269.1 164.1
            >>> st=C3(T=90+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            90 3.7641 328.83 119 483.71 616.47 1.862 2.227 2.107 2.26 7.623 9.888 218.3 155.5
            >>> st=C3(T=95+273.15, x=0.5)
            >>> print "%0.6g %0.5g %0.2f %0.5g %0.5g %0.5g %0.4g %0.4g %0.4g %0.4g %0.4g %0.4g %0.1f %0.1f" % (\
                st.T.C, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            95 4.1195 286.51 156.31 516.33 595.81 1.948 2.164 2.302 2.467 23.59 36.07 158.1 144.1
            """
            , 
        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1, "ho": 26148.48, "so": 157.9105}, 

        "Tmin": Tt, "Tmax": 650.0, "Pmax": 1000000.0, "rhomax": 20.6, 
        "Pmin": 0.00000017, "rhomin": 16.63, 

        "nr1":  [0.42910051e-1, 0.17313671e1, -0.24516524e1, 0.34157466,
                 -0.46047898],
        "d1": [4, 1, 1, 2, 2],
        "t1": [1, 0.33, 0.8, 0.43, 0.9],

        "nr2": [-0.66847295, 0.20889705, 0.19421381, -0.22917851, -0.60405866,
                0.66680654e-1],
        "d2": [1, 3, 6, 6, 2, 3],
        "t2": [2.46, 2.09, 0.88, 1.09, 3.25, 4.62],
        "c2": [1, 1, 1, 1, 2, 2],
        "gamma2": [1]*6,

        "nr3": [0.17534618e-1, 0.33874242, 0.22228777, -0.23219062,
                -0.92206940e-1, -0.47575718, -0.17486824e-1],
        "d3": [1, 1, 1, 2, 2, 4, 1],
        "t3": [0.76, 2.50, 2.75, 3.05, 2.55, 8.40, 6.75],
        "alfa3": [0.963, 1.977, 1.917, 2.307, 2.546, 3.28, 14.6],
        "beta3": [2.33, 3.47, 3.15, 3.19, 0.92, 18.8, 547.8],
        "gamma3": [0.684, 0.829, 1.419, 0.817, 1.500, 1.426, 1.093],
        "epsilon3": [1.283, 0.6936, 0.788, 0.473, 0.8577, 0.271, 0.948]}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Buecker and Wagner (2006)",
        "__doi__": {"autor": "Bücker, D., Wagner, W.",
                    "title": "Reference Equations of State for the Thermodynamic Properties of Fluid Phase n-Butane and Isobutane", 
                    "ref": "J. Phys. Chem. Ref. Data 35, 929 (2006)",
                    "doi": "10.1063/1.1901687"}, 
                    
        "R": 8.314472,
        "cp": CP2,

        "Tmin": Tt, "Tmax": 500.0, "Pmax": 100000.0, "rhomax": 17.41, 
        "Pmin": 0.00000017, "rhomin": 16.62, 

        "nr1": [.21933784906951e1, -.38432884604893e1, .56820219711755,
                .11235233289697, -.13246623110619e-1, .14587076590314e-1,
                .19654925217128e-1],
        "d1": [1, 1, 1, 2, 3, 4, 4],
        "t1": [0.50, 1.00, 1.50, 0.00, 0.50, 0.50, 0.75],

        "nr2": [.73811022854042, -.85976999811290, .14331675665712,
                -.23280677426427e-1, -.98713669399783e-4, .45708225999895e-2,
                -.27766802597861e-1, -.10523131087952, .97082793466314e-1,
                .20710537703751e-1, -.54720320371501e-1, .64918009057295e-3,
                .74471355056336e-2, -.27504616979066e-3, -.77693374632348e-2,
                -.17367624932157e-2],
        "d2": [1, 1, 2, 7, 8, 8, 1, 2, 3, 3, 4, 5, 5, 10, 2, 6],
        "t2": [2.00, 2.50, 2.50, 1.50, 1.00, 1.50, 4.00, 7.00, 3.00, 7.00,
               3.00, 1.00, 6.00, 0.00, 6.00, 13.00],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*16,

        "nr3": [-.38248057095416e-1, -.68797254435490e-2],
        "d3": [1, 2],
        "t3": [2., 0.],
        "alfa3": [10, 10],
        "beta3": [150, 200],
        "gamma3": [1.16, 1.13],
        "epsilon3": [0.85, 1.]}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for propane of Younglove and Ely (1987)",
        "__doc__":  u"""Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
        "R": 8.31434,
        "cp": CP5,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 17.36, 
        "Pmin": 1.685e-7, "rhomin": 16.617, 

        "b": [None, -0.2804337729e-2, 0.1180666107e1, -0.3756325860e2,
              0.5624374521e4, -0.9354759605e6, -0.4557405505e-3, 0.1530044332e1,
              -0.1078107476e4, 0.2218072099e6, 0.6629473971e-4, -0.6199354447e-1,
              0.6754207966e2, 0.6472837570e-2, -0.6804325262, -0.9726162355e2,
              0.5097956459e-1, -0.1004655900e-2, 0.4363693352, -0.1249351947e-1,
              0.2644755879e6, -0.7944237270e8, -0.7299920845e4, 0.5381095003e9,
              0.3450217377e2, 0.9936666689e4, -0.2166699036e1, -0.1612103424e6,
              -0.3633126990e-2, 0.1108612343e2, -0.1330932838e-3,
              -0.3157701101e-1, 0.1423083811e1]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": 85.48, "Tmax": 500.0, "Pmax": 100000.0, "rhomax": 17.41, 
#        "Pmin": 73.476, "rhomin": 29.249, 

        "nr1": [0.10403973107358e1, -0.28318404081403e1, 0.84393809606294,
                -0.76559591850023e-1, 0.94697373057280e-1, 0.24796475497006e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.27743760422870, -0.43846000648377e-1, -0.26991064784350,
                -0.69313413089860e-1, -0.29632145981653e-1, 0.14040126751380e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Miyamoto and Watanabe (2001)",
        "__doc__":  u"""Miyamoto, H. and Watanabe, K. "A Thermodynamic Property Model for Fluid-Phase Propane," Int. J. Thermophys., 21(5):1045-1072, 2000.""",
        "R": 8.314472,
        "cp": CP4,

        "Tmin": Tt, "Tmax": 623.0, "Pmax": 103000.0, "rhomax": 17.41, 
        "Pmin": 0.00000017, "rhomin": 16.64, 

        "nr1":  [2.698378e-1, -1.339252, -2.273858e-2, 2.414973e-1,
                 -3.321461e-2, 2.203323e-3, 5.935588e-5, -1.137457e-6],
        "d1": [1, 1, 2, 2, 3, 5, 8, 8],
        "t1": [-0.25, 1.5, -0.75, 0, 1.25, 1.5, 0.5, 2.5],

        "nr2": [-2.379299, 2.337373, 1.242344e-3, -7.352787e-3, 1.965751e-3,
                -1.402666e-1, -2.093360e-2, -2.475221e-4, -1.482723e-2,
                -1.303038e-2, 3.634670e-5],
        "d2": [3, 3, 8, 5, 6, 1, 5, 7, 2, 3, 15],
        "t2": [1.5, 1.75, -0.25, 3, 3, 4, 2, -1, 2, 19, 5],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*11}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for propane of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.", 
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"}, 
#        "__test__": """
#            >>> st=C3(T=700, rho=200, eq=5)
#            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
#            3.2350 21.175 3.5658
#            >>> st2=C3(T=750, rho=100, eq=5)
#            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
#            212.66 0.41879
#            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": CP6,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 17.36, 
        "Pmin": 0.00000015304, "rhomin": 16.706, 

        "nr1":  [0.10403973e1, -0.28318404e1, 0.8439381, -0.76559592e-1,
                 0.94697373e-1, 0.24796475e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.2774376, -0.43846001e-1, -0.26991065, -0.69313413e-1,
                -0.29632146e-1, 0.14040127e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz6 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.314472,
        "cp": Fi1,
        "ref": {"Tref": 273.15, "Pref": 1, "ho": 26148.48, "so": 157.9105}, 

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [9.70439249e-1, 9.73671323e-1, -2.96661981, 7.84340496e-2,
                2.78440866e-4, -6.77622221e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-8.56371936e-2, 1.77467443e-1, 3.91636018e-1, -8.03312946e-3,
                -2.60385851e-1, -1.91104746e-2, -6.31331470e-2, -2.27769095e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, MBWR, helmholtz2, GERG, helmholtz4, helmholtz5, helmholtz6

    _surface = {"sigma": [0.05666, -0.005291], "exp": [1.265, 1.265]}
    _dielectric = {"eq": 4, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.1573071],  "expt0": [-1.], "expd0": [1.],
                   "a1": [15.85, 0.036], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [172.75, 505.67, -388.21, -2078.8],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.35, 2.35]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.00000017,
                "Tmin": Tt, "Tmax": 2000.0,
                "a1": [-4230000000000, 4230000000001], "exp1": [0, 1.283],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.7722, 1.6938, -1.3341, -3.1876, 0.94937],
        "exp": [1, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [1.82205, 0.65802, 0.21109, 0.083973],
        "exp": [0.345, 0.74, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-2.4887, -5.1069, -12.174, -30.495, -52.192, -134.89],
        "exp": [0.3785, 1.07, 2.7, 5.5, 10., 20.]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.25104574, -0.47271238, 0, 0.060836515],
              "__name__": "Vogel (1998)",
              "__doc__": """Vogel, E., Kuechenmeister, C., Bich, E., and Laesecke, A., "Reference Correlation of the Viscosity of Propane," J. Phys. Chem. Ref. Data, 27(5):947-970, 1998.""",
              "ek": 263.88, "sigma": 0.49748,
              "Tref": 1, "rhoref": 1.*M, "etaref": 1.,
              "n_chapman": 0.141824/M**0.5,

              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0.0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 263.88, "etaref_virial": 0.0741445,

              "Tref_res": 369.82, "rhoref_res": 5.*M, "etaref_res": 1,
              "n_packed": [0.250053938863e1, 0.215175430074e1],
              "t_packed": [0, 0.5],
              "n_poly": [0.359873030195e2, -0.180512188564e3, 0.877124888223e2,
                         -0.105773052525e3, 0.205319740877e3, -0.129210932610e3,
                         0.589491587759e2, -0.129740033100e3, 0.766280419971e2,
                         -0.959407868475e1, 0.210726986598e2,
                         -0.143971968187e2, -0.161688405374e4, ],
              "t_poly": [0, -1, -2, 0, -1, -2, 0, -1, -2, 0, -1, -2, 0],
              "d_poly": [2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1],
              "g_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
              "n_num": [0.161688405374e4],
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
              "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
              "ek": 358.9, "sigma": 0.47,
              "n_chapman": 0.177273976/M**0.5,
              "F": [0, 0, 1.12, 359.],
              "E": [-14.113294896, 968.22940153, 13.686545032, -12511.628378,
                    0.0168910864, 43.527109444, 7659.45434720],
              "rhoc": 5.0}

    visco2 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quiñones-Cisneros (2006)",
              "Tref": 369.825, "muref": 1.0,
              "ek": 358.9, "sigma": 0.47, "n_chapman": 0,
              "n_ideal": [12.3057, -42.5793, 40.3486],
              "t_ideal": [0, 0.25, 0.5],
              "a": [-9.34267734206329e-6, -4.93309341792654e-5, 1.46749885301233e-13],
              "b": [9.60710434008784e-5, -8.18030722274335e-5, 3.00126073333685e-12],
              "c": [7.68800436177747e-5, -4.18871321795657e-5, -7.2008794976648e-15],
              "A": [-8.49308621313605e-9, -4.91414639525551e-10, 0.0],
              "B": [2.08794813407621e-8, 9.21785453914614e-10, 0.0],
              "C": [-4.05944109221870e-7, 1.31730904193479e-7, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"eq": 1,
               "__name__": "Marsh (2002)",
               "__doc__": """Marsh, K., Perkins, R., and Ramires, M.L.V., "Measurement and Correlation of the Thermal Conductivity of Propane from 86 to 600 K at Pressures to 70 MPa," J. Chem. Eng. Data, 47(4):932-940, 2002""",

               "Tref": 369.85, "kref": 1.,
               "no": [-1.24778e-3, 8.16371e-3, 1.99374e-2],
               "co": [0, 1, 2],

               "Trefb": 369.85, "rhorefb": 5., "krefb": 1.,
               "nb": [-3.69500e-2, 4.82798e-2, 1.48658e-1, -1.35636e-1,
                      -1.19986e-1, 1.17588e-1, 4.12431e-2, -4.36911e-2,
                      -4.86905e-3, 6.16079e-3],
               "tb": [0, 1]*5,
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.716635e-9, "Tcref": 554.73}

    thermo1 = {"eq": 2, "omega": 2,
               "__name__": "Younglove (1987)",
               "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
               "visco": visco1,
               "n_chapman": 1.77273976e-1,
               "G": [0.1422605e1, -0.179749],
               "E": [0.3113890422e-2, -0.225755973, 0.5674370999e2,
                     -0.7840963643e-4, 0.2291785465e-1, -0.2527939890e1,
                     -0.6265334654e-1, 0.2518064809e1],

               "critical": 2,
               "X": [3.98, 5.45, 0.468067, 1.08],
               "Z": 8.117e-10}

    _thermal = thermo0, thermo1

if __name__ == "__main__":

    st=C3(T=85.525, x=0.5)
    print "%0.6g %0.6f %0.5f %0.5f %0.5g %0.5g %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.2f %0.2f" % (\
        st.T.K, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
        st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
        st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
