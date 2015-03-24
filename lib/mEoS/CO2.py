#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class CO2(MEoS):
    """Multiparameter equation of state for carbon dioxide"""
    name = "carbon dioxide"
    CASNumber = "124-38-9"
    formula = "CO2"
    synonym = "R-744"
    rhoc = unidades.Density(467.6)
    Tc = unidades.Temperature(304.1282)
    Pc = unidades.Pressure(7.3773, "MPa")
    M = 44.0098  # g/mol
    Tt = unidades.Temperature(216.592)
    Tb = unidades.Temperature(194.686)
    f_acent = 0.22394
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 49

    CP1 = {"ao": 3.5,
           "an": [], "pow": [],
           "ao_exp": [1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678],
           "exp": [958.49956, 1858.80115, 2061.10114, 3443.89908, 8238.20035],
           "ao_hyp": [], "hyp": []}
    
    Fi1 = {"ao_log": [1, 2.5],
           "pow": [0, 1],
           "ao_pow": [8.37304456, -3.70454304],
           "ao_exp": [1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678],
           "titao": [3.15163, 6.11190, 6.77708, 11.32384, 27.08792], 
           "ao_hyp": [], "hyp": []}

           
    Fi2 = {"ao_log": [1, 2.50002],
           "pow": [0, 1],
           "ao_pow": [11.925152758, -16.118762264],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [2.04452, -1.06044, 2.03366, 0.01393],
           "hyp": [3.022758166, -2.844425476, 1.589964364, 1.12159609]}

    CP2 = {"ao": 3.50002,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [2.04452, -1.06044, 2.03366, 0.01393],
           "hyp": [3.022758166*Tc, -2.844425476*Tc, 1.589964364*Tc, 1.12159609*Tc]}

    CP3 = {"ao": 3.5,
           "an": [], "pow": [],
           "ao_exp": [2, 1, 1], "exp": [960.11, 1932, 3380.2],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Span and Wagner (1996)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "A New Equation of State for Carbon Dioxide Covering the Fluid Region from the Triple‐Point Temperature to 1100 K at Pressures up to 800 MPa", 
                    "ref": "J. Phys. Chem. Ref. Data 25, 1509 (1996)",
                    "doi": "10.1063/1.555991"}, 
        "__test__": 
            # Table 34, Pag 1560
            """
            >>> st=CO2(T=216.592, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            216.592 0.51796 1178.46 13.761 -426.74 -76.364 -2.2177 -0.59999 0.97466 0.62921 1.9532 0.90872 975.85 222.78
            >>> st=CO2(T=230, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            230 0.89291 1128.68 23.271 -400.21 -72.178 -2.1003 -0.67406 0.95667 0.67004 1.997 1.0053 879.09 223.57
            >>> st=CO2(T=240, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            240 1.2825 1088.87 33.295 -379.94 -70.293 -2.0155 -0.72532 0.94535 0.70534 2.051 1.1033 806.38 222.96
            >>> st=CO2(T=250, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            250 1.785 1045.97 46.644 -359.07 -69.736 -1.9323 -0.77492 0.93643 0.74591 2.132 1.2366 731.78 221.22
            >>> st=CO2(T=260, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            260 2.4188 998.89 64.417 -337.34 -70.862 -1.8495 -0.82456 0.93227 0.79426 2.2554 1.4295 652.58 218.19
            >>> st=CO2(T=270, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            270 3.2033 945.83 88.374 -314.37 -74.223 -1.7658 -0.87641 0.93959 0.85168 2.4534 1.7307 565.46 213.75
            >>> st=CO2(T=280, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            280 4.1607 883.58 121.74 -289.48 -80.84 -1.6792 -0.93401 0.96046 0.92316 2.8141 2.2769 471.54 207.72
            >>> st=CO2(T=290, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            290 5.3177 804.67 171.96 -261.15 -93.025 -1.5846 -1.0049 0.99373 1.026 3.6756 3.6142 371.95 199.45
            >>> st=CO2(T=300, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            300 6.7131 679.24 268.58 -223.4 -119.7 -1.4631 -1.1175 1.1199 1.2476 8.6979 11.921 245.67 185.33
            >>> st=CO2(T=304.1282, x=0.5)
            >>> print "%0.7g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Gas.rho, st.Gas.h.kJkg, st.Gas.s.kJkgK)
            304.1282 7.3773 467.6 -174.53 -1.3054
            """
            # Table 35, Pag 1562
            """
            >>> st=CO2(T=190, P=5e4)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            190 1.4089 -121.78 -86.286 -0.22345 0.54661 0.7466 218.9
            >>> st=CO2(T=1100, P=5e4)
            >>> print "%0.6g %0.5f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            1100 0.24057 675.87 883.71 1.5137 1.0702 1.2592 494.54
            >>> st=CO2(T=280, P=1e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            280 1.9021 -68.784 -16.209 -0.05256 0.639 0.8333 261.03
            >>> st=CO2(T=300, P=101325)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 1.7966 -55.76 0.63726 0.00307 0.65935 0.85262 269.38
            >>> st=CO2(T=1100, P=2e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            1100 0.96197 675.75 883.66 1.2517 1.0702 1.2594 494.75
            >>> st=CO2(T=500, P=5e5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            500 5.3126 93.079 187.19 0.1757 0.82672 1.0207 340.25
            >>> st=CO2(T=220, P=1e6)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            220 1167.03 -420.8 -419.95 -2.1884 0.97034 1.9589 953.55
            >>> st=CO2(T=300, P=2e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 39.42 -69.155 -18.42 -0.60581 0.71116 1.0206 254.15
            >>> st=CO2(T=500, P=5e6)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            500 54.826 81.644 172.84 -0.28196 0.84035 1.084 337.43
            >>> st=CO2(T=300, P=1e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 801.62 -257.46 -244.98 -1.5496 0.94964 2.9906 414.28
            >>> st=CO2(T=1000, P=2e7)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            1000 101.27 553.1 750.6 0.24652 1.0516 1.2727 501.16
            >>> st=CO2(T=300, P=5e7)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            300 1028.94 -310.05 -261.46 -1.7458 0.92303 1.7514 827.23
            >>> st=CO2(T=240, P=1e8)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            240 1257.21 -424.46 -344.92 -2.2162 1.0037 1.6854 1265.1
            >>> st=CO2(T=600, P=2e8)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            600 876.69 8.1703 236.3 -0.96106 1.0016 1.4019 942
            >>> st=CO2(T=310, P=6e8)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            310 1444.64 -386.39 28.939 -2.1927 1.1442 1.5427 1903
            >>> st=CO2(T=1100, P=8e8)
            >>> print "%0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            1100 1092.77 545.32 1277.4 -0.43587 1.202 1.4286 1542.2
            """, 

        "R": 8.31451,
        "cp": Fi1,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 2000., "Pmax": 800000.0, "rhomax": 37.24, 
        "Pmin": 517.95, "rhomin": 26.777, 

        "nr1": [0.38856823203161, 0.29385475942740e1, -0.55867188534934e1,
                -0.76753199592477, 0.31729005580416, 0.54803315897767,
                0.12279411220335],
        "d1": [1, 1, 1, 1, 2, 2, 3],
        "t1": [0.0, 0.75, 1.0, 2.0, 0.75, 2.0, 0.75],

        "nr2": [0.21658961543220e1, 0.15841735109724e1, -0.23132705405503,
                0.58116916431436e-1, -0.55369137205382, 0.48946615909422,
                -0.24275739843501e-1, 0.62494790501678e-1, -0.12175860225246,
                -0.37055685270086, -0.16775879700426e-1, -0.11960736637987,
                -0.45619362508778e-1, 0.35612789270346e-1, -0.74427727132052e-2,
                -0.17395704902432e-2, -0.21810121289527e-1, 0.24332166559236e-1,
                -0.37440133423463e-1, 0.14338715756878, -0.13491969083286,
                -0.23151225053480e-1, 0.12363125492901e-1, 0.21058321972940e-2,
                -0.33958519026368e-3, 0.55993651771592e-2, -0.30335118055646e-3],
        "d2": [1, 2, 4, 5, 5, 5, 6, 6, 6, 1, 1, 4, 4, 4, 7, 8, 2, 3, 3, 5, 5,
               6, 7, 8, 10, 4, 8],
        "t2": [1.5, 1.5, 2.5, 0, 1.5, 2, 0, 1, 2, 3, 6, 3, 6, 8, 6, 0, 7, 12,
               16, 22, 24, 16, 24, 8, 2, 28, 14],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4,
               4, 4, 4, 4, 5, 6],
        "gamma2": [1]*27,

        "nr3": [-0.21365488688320e3, 0.26641569149272e5, -0.24027212204557e5,
                -0.28341603423999e3, 0.21247284400179e3],
        "d3": [2, 2, 2, 3, 3],
        "t3": [1., 0., 1., 3., 3.],
        "alfa3": [25, 25, 25, 15, 20],
        "beta3": [325, 300, 300, 275, 275],
        "gamma3": [1.16, 1.19, 1.19, 1.25, 1.22],
        "epsilon3": [1.]*5,

        "nr4": [-0.66642276540751, 0.72608632349897, 0.55068668612842e-1],
        "a4": [3.5, 3.5, 3.],
        "b4": [0.875, 0.925, 0.875],
        "beta4": [0.3]*3,
        "A": [0.7]*3,
        "B": [0.3, 0.3, 1.],
        "C": [10., 10., 12.5],
        "D": [275]*3}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for carbon dioxide of Ely et al. (1987)",
        "__doc__":  u"""Ely, J.F., Magee, J.W., and Haynes, W.M., "Thermophysical properties for special high CO2 content mixtures," Research Report RR-110, Gas Processors Association, Tulsa, OK, 1987.""",
        "R": 8.31434,
        "cp": CP3,

        "Tmin": 216.58, "Tmax": 440.1, "Pmax": 40000.0, "rhomax": 27.778, 
        "Pmin": 518.2, "rhomin": 26.778, 

        "b": [None, -0.981851065838e-2, 0.995062267309, -0.228380160313e2,
              0.281827634529e4, -0.347001262699e6, 0.394706709102e-3,
              -0.325550000110, 0.484320083063e1, -0.352181542995e6,
              -0.324053603343e-4, 0.468596684665e-1, -0.754547012075e1,
              -0.381894354016e-4, -0.442192933859e-1, 0.516925168095e2,
              0.212450985237e-2, -0.261009474785e-4, -0.888533388977e-1,
              0.155226179403e-2, 0.415091004940e6, -0.110173967489e8,
              0.291990583344e4, 0.143254606508e8, 0.108574207533e2,
              -0.247799657039e3, 0.199293590763e-1, 0.102749908059e3,
              0.377618865158e-4, -0.332276512346e-2, 0.179196707121e-7,
              0.945076627807e-5, -0.123400943061e-2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Kunz and Wagner (2004)",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 1100., "Pmax": 800000.0, "rhomax": 37.24, 
#        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1": [0.52646564804653, -0.14995725042592e1, 0.27329786733782,
                0.12949500022786],
        "d1": [1, 1, 2, 3],
        "t1": [0, 1.25, 1.625, 0.375],

        "nr2": [0.15404088341841, -0.58186950946814, -0.18022494838296,
                -0.95389904072812e-1, -0.80486819317679e-2, -0.3554775127309e-1,
                -0.28079014882405, -0.82435890081677e-1, 0.10832427979006e-1,
                -0.67073993161097e-2, -0.46827907600524e-2, -0.28359911832177e-1,
                0.19500174744098e-1, -0.21609137507166, 0.43772794926972,
                -0.22130790113593, 0.15190189957331e-1, -0.15380948953300e-1],
        "d2": [3, 3, 4, 5, 6, 6, 1, 4, 1, 1, 3, 3, 4, 5, 5, 5, 5, 5],
        "t2": [0.375, 1.375, 1.125, 1.375, 0.125, 1.625, 3.75, 3.5, 7.5, 8, 6,
               16, 11, 24, 26, 28, 24, 26],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 5, 5, 5, 6, 6],
        "gamma2": [1]*18}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Ely et al. (1987)",
        "__doi__": {"autor": "Ely, J.F., Magee, J.W., and Haynes, W.M.",
                    "title": "Thermophysical properties for special high CO2 content mixtures", 
                    "ref": "Research Report RR-110, Gas Processors Association, Tulsa, OK, 1987.",
                    "doi": ""}, 
        "R": 8.31434,
        "cp": CP3,

        "Tmin": 216.58, "Tmax": 1000., "Pmax": 100000.0, "rhomax": 26.776, 
        "Pmin": 518.03, "rhomin": 26.776, 

        "nr1": [0.485497428986, -0.191900462349e1, 0.451739876847,
                0.838475229022e-2, 0.310719428397, -0.183619563850,
                0.448878785519e-1, -0.362211893044e-1, -0.169827491865e-1,
                0.803504394396e-3, 0.320223641512e-3, -0.658956249553e-5,
                -0.461991678692e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-0.385989029443, 0.131878614095, 0.109639470331,
                -0.310044422115e-1, -0.989797992915e-1, -0.222934996927e-1,
                -0.225488505376e-1, -0.595661202393e-2, -0.219959964099e-1,
                0.140330955537e-1, -0.315424157971e-2, 0.443394060420e-3,
                -0.487628903103e-2, -0.311643343682e-1, 0.226083669848e-1,
                0.186651858191e-1, -0.399277963883, 0.464945130861,
                -0.817090055061e-1],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5, ],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*19}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for carbon dioxide of Span and Wagner (2003)",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. III. Results for Polar Fluids," Int. J. Thermophys., 24(1):111-162, 2003.""",
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 600., "Pmax": 100000.0, "rhomax": 37.24, 
        "Pmin": 517.86, "rhomin": 26.795, 

        "nr1": [0.89875108, -0.21281985e1, -0.6819032e-1, 0.76355306e-1,
                0.22053253e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.41541823, 0.71335657, 0.30354234e-3, -0.36643143,
                -0.14407781e-2, -0.89166707e-1, -0.23699887e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for carbon dioxide of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [-4.71122371e-1, 9.13375599e-1, -1.96793707, 6.89687161e-2,
                2.15658922e-4, 9.51876380e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-4.91366518e-3, 7.32487713e-1, 8.70918629e-1, -5.35917679e-3,
                -4.03818537e-1, -2.40820897e-2, -1.04239403e-1, -2.16335828e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, MBWR, GERG, helmholtz3, helmholtz4, helmholtz5

    _surface = {"sigma": [0.084497], "exp": [1.28]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [7.3455, 0.00335], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [83.93, 145.1, -578.8, -1012.],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.55, 2.55]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 517.95,
                "Tmin": Tt, "Tmax": 1100.0,
                "a1": [1], "exp1": [0],
                "a2": [1955.539, 2055.4593], "exp2": [1, 2],
                "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 517.950,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-14.740846, 2.4327015, -5.3061778],
                    "exp2": [1, 1.9, 2.9],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.0602087, 1.9391218, -1.6463597, -3.2995634],
        "exp": [1, 1.5, 2., 4.]}
    _liquid_Density = {
        "eq": 4,
        "ao": [1.92451080, -0.62385555, -0.32731127, 0.39245142],
        "exp": [1.02, 1.5, 5.0, 5.5]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-1.7074879, -0.8227467, -4.6008549, -10.111178, -29.742252],
        "exp": [1.02, 1.5, 3.0, 7.0, 14.0]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.235156, -0.491266, 5.211155e-2, 5.347906e-2,
                            -1.537102e-2],
              "__name__": "Fenghour (1998)",
              "__doc__": """Fenghour, A., Wakeham, W.A., Vesovic, V., "The Viscosity of Carbon Dioxide," J. Phys. Chem. Ref. Data, 27:31-44, 1998.""",
              "bmega_b": [0.235156, -0.491266, 5.211155e-2, 5.347906e-2,
                          -1.537102e-2],
              "ek": 251.196, "sigma": 1.,
              "n_chapman": 1.00697/M**0.5,
              "Tref": 1.,

              "Tref_res": 251.196, "rhoref": 0.0227222*M,
              "n_poly": [0.4071119e-2, 0.7198037e-4, 0.2411697e-16,
                         0.2971072e-22, -0.1627888e-22],
              "t_poly": [0, 0, -3, 0, -1],
              "d_poly": [1, 2, 6, 8, 8],
              "g_poly": [0, 0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0, 0]}

    visco1 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quiñones-Cisneros (2006)",
              "ek": 251.196, "sigma": 0.3751, "n_chapman": 0,
              "Tref": 304.1282, "muref": 1.0,

              "n_ideal": [69.18424, -215.8618, 210.94362, -49.0494],
              "t_ideal": [0, 0.25, 0.5, 0.75],

              "a": [1.19805e-4,  -1.25861e-4, 5.48871e-5],
              "b": [3.15921e-5, -2.60469e-5, 7.09199e-6],
              "c": [1.80689e-5, -7.41742e-6, 0.0],
              "A": [-2.31066e-9, 0.0, 5.42486e-10],
              "B": [1.04558e-8, -2.20758e-9, 0.0],
              "C": [1.03255e-6, -8.56207e-7, 3.84384e-7],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 1,
               "__name__": "Vesovic (1990)",
               "__doc__": """Vesovic, V., Wakeham, W.A., Olchowy, G.A., Sengers, J.V., Watson, J.T.R., and Millat, J., "The transport properties of carbon dioxide," J. Phys. Chem. Ref. Data, 19:763-808, 1990""",

               "Tref": 251.196, "kref": 1e-3,
               "no": [7.5378307, 4.8109652e-2],
               "co": [0.5, -99],
               "noden": [0.4226159, 0.6280115, -0.5387661, 0.6735941,
                         -0.4362677, 0.2255388],
               "toden": [0, -1, -2, -3, -6, -7],

               "Trefb": 1., "rhorefb": 2.272221e-2, "krefb": 1e-3,
               "nb": [2.447164e-2, 8.705605e-5, -6.547950e-8, 6.594919e-11],
               "tb": [0, 0, 0, 0],
               "db": [1, 2, 3, 4],
               "cb": [0, 0, 0, 0],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 1.5e-10, "gam0": 0.052, "qd": 0.4e-9, "Tcref": 450.}

    _thermal = thermo0,
