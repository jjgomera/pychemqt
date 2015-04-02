#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp

from lib import unidades
from lib.iapws import _Viscosity, _ThCond, _Dielectric
from lib.meos import MEoS


class H2O(MEoS):
    """Multiparameter equation of state for water (including IAPWS95)

#    >>> water=H2O(T=300, rho=996.5560)
#    >>> print("%0.10f %0.8f %0.5f %0.9f" % (water.P.MPa, water.cv.kJkgK, water.w, water.s.kJkgK))
#    0.0992418350 4.13018112 1501.51914 0.393062643
#
#    >>> water=H2O(T=500, rho=0.435)
#    >>> print("%0.10f %0.8f %0.5f %0.9f" % (water.P.MPa, water.cv.kJkgK, water.w, water.s.kJkgK))
#    0.0999679423 1.50817541 548.31425 7.944882714
#
#    >>> water=H2O(T=900., P=700e6)
#    >>> print("%0.4f %0.8f %0.5f %0.8f" % (water.rho, water.cv.kJkgK, water.w, water.s.kJkgK))
#    870.7690 2.66422350 2019.33608 4.17223802
#
#    >>> water=H2O(T=300., P=0.1e6)
#    >>> print("%0.2f %0.5f %0.2f %0.2f %0.5f %0.4f %0.1f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    300.00 0.10000 996.56 112.65 0.39306 4.1806 1501.5 -0.066682
#
#    >>> water=H2O(T=500., P=0.1e6)
#    >>> print("%0.2f %0.5f %0.5f %0.1f %0.4f %0.4f %0.2f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    500.00 0.10000 0.43514 2928.6 7.9447 1.9813 548.31 -0.0094137
#
#    >>> water=H2O(T=450., x=0.5)
#    >>> print("%0.2f %0.5f %0.4f %0.1f %0.4f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.virialB))
#    450.00 0.93220 9.5723 1761.8 4.3589 -0.013028
#
#    >>> water=H2O(P=1.5e6, rho=1000.)
#    >>> print("%0.2f %0.4f %0.1f %0.3f %0.5f %0.4f %0.1f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    286.44 1.5000 1000.0 57.253 0.19931 4.1855 1462.1 -0.085566
#
#    >>> water=H2O(h=3000e3, s=8e3)
#    >>> print("%0.2f %0.5f %0.5f %0.1f %0.4f %0.4f %0.2f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    536.24 0.11970 0.48547 3000.0 8.0000 1.9984 567.04 -0.0076606
#
#    >>> water=H2O(h=150e3, s=0.4e3)
#    >>> print("%0.2f %0.5f %0.2f %0.2f %0.5f %0.4f %0.1f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    301.27 35.50549 1011.48 150.00 0.40000 4.0932 1564.1 -0.065238
#
#    >>> water=H2O(T=450., rho=300)
#    >>> print("%0.2f %0.5f %0.2f %0.2f %0.4f %0.6f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    450.00 0.93220 300.00 770.82 2.1568 0.010693 -0.013028
#
#    >>> water=H2O(rho=300., P=0.1e6)
#    >>> print("%0.2f %0.5f %0.2f %0.2f %0.4f %0.7f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    372.76 0.10000 300.00 420.56 1.3110 0.0013528 -0.025144
#
#    >>> water=H2O(h=1500e3, P=0.1e6)
#    >>> print("%0.2f %0.5f %0.4f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    372.76 0.10000 1.2303 1500.0 4.2068 0.47952 -0.025144
#
#    >>> water=H2O(s=5e3, P=3.5e6)
#    >>> print("%0.2f %0.4f %0.3f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    515.71 3.5000 25.912 2222.8 5.0000 0.66921 -0.0085877
#
#    >>> water=H2O(T=500., u=900e3)
#    >>> print("%0.2f %0.2f %0.2f %0.2f %0.1f %0.4f %0.4f %0.1f %0.7f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    500.00 108.21 903.62 900.00 1019.8 2.4271 4.1751 1576.0 -0.0094137
#
#    >>> water=H2O(P=0.3e6, u=1550.e3)
#    >>> print("%0.2f %0.5f %0.4f %0.1f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    406.67 0.30000 3.3029 1550.0 1640.8 4.3260 0.49893 -0.018263
#
#    >>> water=H2O(rho=300, h=1000.e3)
#    >>> print("%0.2f %0.4f %0.2f %0.2f %0.1f %0.4f %0.6f %0.7f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    494.92 2.3991 300.00 992.00 1000.0 2.6315 0.026071 -0.0097064
#
#    >>> water=H2O(rho=30, s=8.e3)
#    >>> print("%0.2f %0.3f %0.3f %0.1f %0.1f %0.4f %0.4f %0.2f %0.9f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    1562.42 21.671 30.000 4628.5 5350.9 8.0000 2.7190 943.53 0.000047165
#
#    >>> water=H2O(rho=30, s=4.e3)
#    >>> print("%0.2f %0.4f %0.3f %0.1f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    495.00 2.4029 30.000 1597.3 1677.4 4.0000 0.39218 -0.0097015
#
#    >>> water=H2O(rho=300, u=1000.e3)
#    >>> print("%0.2f %0.4f %0.3f %0.1f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    496.44 2.4691 300.000 1000.0 1008.2 2.6476 0.02680 -0.0096173
#
#    >>> water=H2O(s=3.e3, h=1000.e3)
#    >>> print("%0.2f %0.6f %0.5f %0.2f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    345.73 0.034850 0.73526 952.60 1000.0 3.0000 0.29920 -0.034124
#
#    >>> water=H2O(u=995.e3, h=1000.e3)
#    >>> print("%0.2f %0.4f %0.2f %0.2f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    501.89 2.7329 546.58 995.00 1000.0 2.6298 0.00866 -0.009308
#
#    >>> water=H2O(u=1000.e3, s=3.e3)
#    >>> print("%0.2f %0.6f %0.5f %0.2f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    371.24 0.094712 1.99072 1000.00 1047.6 3.0000 0.28144 -0.025543

    """
    name = "water"
    CASNumber = "7732-18-5"
    formula = "H2O"
    synonym = "R-718"
    Tc = unidades.Temperature(647.096)
    rhoc = unidades.Density(322.)
    Pc = unidades.Pressure(22064.0, "kPa")
    M = 18.015268  # g/mol
    Tt = unidades.Temperature(273.16)
    Tb = unidades.Temperature(373.1243)
    f_acent = 0.3443
    momentoDipolar = unidades.DipoleMoment(1.855, "Debye")
    id = 62

    Fi1 = {"ao_log": [1, 3.00632],
           "pow": [0, 1],
           "ao_pow": [-8.3204464837497, 6.6832105275932],
           "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
           "titao": [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]}
           
    Fi2 = {"ao_log": [1, 3.00392],
           "pow": [0, 1],
           "ao_pow": [8.203520690, -11.996306443],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [0.01059, -0.98763, 3.06904, 0],
           "hyp": [0.415386589, 1.763895929, 3.874803739, 0]}

    Fi3 = {"ao_log": [1, 3.00632],
           "pow": [0, 1],
           "ao_pow": [-8.318441, 6.681816],
           "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
           "titao": [1.287202151, 3.537101709, 7.740210774, 9.243749421, 27.5056402]}

    Fi4 = {"ao_log": [1, 3.00632],
           "pow": [0, 1],
           "ao_pow": [-8.3177095, 6.6815049],
           "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
           "titao": [1.287202151, 3.537101709, 7.740210774, 9.243749421, 27.5056402]}


    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for water of Wagner and Pruß (2002).",
        "__doi__": {"autor": u"Wagner, W., Pruß, A.",
                    "title": "The IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water Substance for General and Scientific Use", 
                    "ref": "J. Phys. Chem. Ref. Data 31, 387 (2002)",
                    "doi": "10.1063/1.1461829"}, 
        "__test__": 
            # Table 6.6, Pag 436
            """
            >>> wt=H2O()
            >>> tau=wt.Tc/500
            >>> delta=838.025/wt.rhoc
            >>> print "%0.9g %0.9g %0.9g %0.9g %0.9g %0.9g" % wt._phi0(wt._constants["cp"], tau, delta)
            2.04797734 9.04611106 -1.93249185 0.384236747 -0.147637878 0
            >>> print "%0.9g %0.9g %0.9g %0.9g %0.9g %0.9g" % wt._phir(tau, delta)[:6]
            -3.42693206 -5.81403435 -2.23440737 -0.36436665 0.856063701 -1.12176915
            >>> tau=wt.Tc/647
            >>> delta=358/wt.rhoc
            >>> print "%0.9g %0.9g %0.9g %0.9g %0.9g %0.9g" % wt._phi0(wt._constants["cp"], tau, delta)
            -1.56319605 9.80343918 -3.43316334 0.899441341 -0.808994726 0
            >>> print "%0.9g %0.9g %0.9g %0.9g %0.9g %0.9g" % wt._phir(tau, delta)[:6]
            -1.21202657 -3.21722501 -9.96029507 -0.714012024 0.475730696 -1.3321472
            """
            #Table 13.1, Pag 486
            """
            >>> st=H2O(T=273.16, x=0.5)
            >>> print "%0.6g %0.3g %0.6g %0.3g %0.3f %0.6g %0.4f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            273.16 0.000612 999.793 0.00485 0.001 2500.92 -0.0000 9.1555 4.2174 1.4184 4.2199 1.8844 1402.3 409

            >>> st=H2O(T=300, x=0.5)
            >>> print "%0.6g %0.4g %0.6g %0.4g %0.3f %0.6g %0.4f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            300 0.003537 996.513 0.02559 112.565 2549.85 0.3931 8.5174 4.1305 1.4422 4.1809 1.9141 1501.4 427.89

            >>> st=H2O(T=400, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.3f %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            400 0.24577 937.486 1.3694 532.953 2715.7 1.6013 7.0581 3.6324 1.6435 4.2555 2.2183 1509.5 484.67
            
            >>> st=H2O(T=500, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.3f %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            500 2.6392 831.313 13.199 975.431 2802.48 2.581 6.2351 3.2255 2.2714 4.6635 3.4631 1239.6 504.55

            >>> st=H2O(T=600, x=0.5)
            >>> print "%0.6g %0.5g %0.6g %0.5g %0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            600 12.345 649.411 72.842 1505.36 2677.81 3.519 5.4731 3.0475 3.3271 6.9532 9.1809 749.57 457.33

            >>> st=H2O(T=646, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            646 21.775 402.96 243.46 1963.49 2238.06 4.2214 4.6465 4.5943 5.1457 204.58 385.23 297.13 331.61
                
            >>> st=H2O(T=647, x=0.5)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
                st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
                st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
            647 22.038 357.34 286.51 2029.44 2148.56 4.3224 4.5065 6.2344 6.274 3905.2 5334.1 251.19 285.32
            """
            #Table 13.2, Pag 495
            """
            >>> st=H2O(T=290, P=50000)
            >>> print "%0.6g %0.5g %0.5g %0.5g %0.4g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            290 998.78 70.725 70.775 0.2513 4.1682 4.1868 1472.2

            >>> st=H2O(T=600, P=15000000)
            >>> print "%0.6g %0.6g %0.6g %0.6g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            600 659.407 1474.91 1497.65 3.4994 3.0282 6.583 787.74

            >>> st=H2O(T=640, P=22500000)
            >>> print "%0.6g %0.6g %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            640 531.385 1744.16 1786.5 3.9443 3.1724 12.142 529.59
            
            >>> st=H2O(T=580, P=40000000)
            >>> print "%0.6g %0.6g %0.6g %0.6g %0.5g %0.5g %0.5g %0.5g" % (\
                st.T, st.rho, st.u.kJkg, st.h.kJkg, st.s.kJkgK, st.cv.kJkgK, st.cp.kJkgK, st.w)
            580 753.362 1306.43 1359.53 3.2061 2.9966 4.9863 1099.6
            """, 
            
        "R": 8.314371357587,
        "cp": Fi1,
        "ref": {"Tref": Tt, "Pref": 0.611655, "ho": 0.611872, "so": 0}, 

        "Tmin": Tt, "Tmax": 2000., "Pmax": 2000000.0, "rhomax": 73.96, 
        "Pmin": 0.61248, "rhomin": 55.49696, 

        "nr1": [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
                0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
                0.88089493102134e-2],
        "d1": [1, 1, 1, 2, 2, 3, 4],
        "t1": [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1],

        "nr2": [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
                -0.19232721156002, -0.25709043003438, 0.16074868486251,
                -0.4009282892587e-1, 0.39343422603254e-6, -0.75941377088144e-5,
                0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
                0.36582165144204e-6, -0.13251180074668e-11, -0.62639586912454e-9,
                -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
                -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
                -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
                0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
                0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
                -0.20393486513704e-1, -0.16554050063734e-2, 0.19955571979541e-2,
                0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
                0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
                -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
                0.31777497330738, -0.11841182425981],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
               2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,
               6, 6],
        "d2": [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
               4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14,
               3, 6, 6, 6],
        "t2": [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10,
               10, 3, 7, 10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,
               23, 10, 50, 44, 46, 50],
        "gamma2": [1]*44,

        "nr3": [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4],
        "d3": [3]*3,
        "t3": [0, 1, 4],
        "alfa3": [20]*3,
        "beta3": [150, 150, 250],
        "gamma3": [1.21, 1.21, 1.25],
        "epsilon3": [1.]*3,

        "nr4": [-0.14874640856724, 0.31806110878444],
        "a4": [3.5, 3.5],
        "b4": [0.85, 0.95],
        "B": [0.2, 0.2],
        "C": [28, 32],
        "D": [700, 800],
        "A": [0.32, .32],
        "beta4": [0.3, 0.3]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032-3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 1350.0, "Pmax": 1000000.0, "rhomax": 73.96, 
        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1": [0.82728408749586, -0.18602220416584e1, -0.11199009613744e1,
                0.15635753976056, 0.87375844859025, -0.36674403715731,
                0.53987893432436e-1],
        "d1": [1, 1, 1, 2, 2, 3, 4],
        "t1": [0.5, 1.25, 1.875, 0.125, 1.5, 1, 0.75],

        "nr2": [0.10957690214499e1, 0.53213037828563e-1, 0.13050533930825e-1,
                -0.41079520434476, 0.14637443344120, -0.55726838623719e-1,
                -0.11201774143800e-1, -0.66062758068099e-2, 0.46918522004538e-2],
        "c2": [1, 1, 1, 2, 2, 2, 3, 5, 5],
        "d2": [1, 5, 5, 1, 2, 4, 4, 1, 1],
        "t2": [1.5, 0.625, 2.625, 5, 4, 4.5, 3, 4, 6],
        "gamma2": [1]*9}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Saul and Wagner-58 coeff (1989).",
        "__doi__": {"autor": "Saul, A. and Wagner, W.",
                    "title": "A Fundamental Equation for Water Covering the Range from the Melting Line to 1273 K at Pressures up to 25000 MPa", 
                    "ref": "J. Phys. Chem. Ref. Data 18, 1537 (1989)",
                    "doi": "10.1063/1.555836"}, 
        "R": 8.31434,
        "cp": Fi3,
        "ref": {"Tref": Tt, "Pref": 611.655, "ho": 0.611872, "so": 0}, 

        "Tmin": Tt, "Tmax": 1273., "Pmax": 400000.0, "rhomax": 55.49, 
        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1": [0.8216377478, -0.2543894379, -0.08830868648, -0.8903097248e-6, 
                -0.1241333357e-5, 0.2895590286e-8, 0.1403610309e-10, 
                0.8183943371e-12, -0.2397905287e-12],
        "d1": [1, 1, 2, 5, 8, 11, 11, 13, 13],
        "t1": [0, 2, 0, 9, 0, 0, 12, 7, 13],

        "nr2": [-0.7519743341, -0.4151278588, -0.103051374e1, -0.1648036888e1, 
                -0.4686350251, 0.3560258142, -0.6364658294, 0.2227482363, 
                -0.8954849939e-1, 0.1557686788e-2, 0.1347719088e-2, 
                -0.1301353385e-2, 0.9987368673e-6, 0.2263629476e-3, 
                0.289330495e-5, 0.1995437169, -0.2707767662e-1, 
                0.1849068216e-1, -0.4402394357e-2, -0.8546876737e-1, 
                0.1220538576, -0.2562237041, 0.2555034636, -0.6323203907e-1, 
                0.3351397575e-4, -0.6152834985e-1, -0.3533048208e-3,
                0.3146309259e-1, -0.2261795983e-2, 0.18689702e-3, 
                -0.1384614556e-2, 0.2713160073e-2, -0.4866118539e-2, 
                0.3751789129e-2, -0.5692669373e-3, -0.5876414555, 0.5687838346,
                -0.1642158198, 0.5878635885, -0.2844301931, -0.2049198337, 
                -0.4039233716e-2, 0.5459049594e-1, -0.8914260146e-2, 
                0.4974411254e-2],
        "c2": [1]*15+[2]*20+[3]*10,
        "d2": [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 11, 1, 1, 1, 1, 2, 2,
               4, 5, 6, 6, 7, 7, 8, 10, 10, 11, 11, 11, 11, 11, 2, 2, 2, 3, 3,
               4, 4, 5, 5, 5],
        "t2": [0, 1, 3, 1, 5, 5, 2, 3, 5, 6, 4, 1, 8, 0, 1, 0, 9, 10, 11, 0, 
               8, 5, 4, 2, 12, 3, 10, 3, 2, 8, 0, 1, 3, 4, 6, 13, 14, 15, 14,
               16, 13, 26, 15, 23, 25],
        "gamma2": [1]*45, 
        
        "nr5": [-0.709318338e-2, 0.1718796342e-1, -0.1482653038e-1, 
                0.4517292884e-2], 
        "d5": [1, 2, 3, 4], 
        "t5": [50, 40, 32, 26]
        }

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Saul and Wagner-38 coeff (1989).",
        "__doi__": {"autor": "Saul, A. and Wagner, W.",
                    "title": "A Fundamental Equation for Water Covering the Range from the Melting Line to 1273 K at Pressures up to 25000 MPa", 
                    "ref": "J. Phys. Chem. Ref. Data 18, 1537 (1989)",
                    "doi":  "10.1063/1.555836"}, 
        "R": 8.31434,
        "cp": Fi4,
        "ref": {"Tref": Tt, "Pref": 611.655, "ho": 0.611872, "so": 0}, 

        "Tmin": Tt, "Tmax": 1273., "Pmax": 400000.0, "rhomax": 55.49, 
        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1": [0.2330009013, -0.1402091128e1, 0.1172248041, -0.1850749499,
                0.1770110422, 0.5525151794e-1, -0.341325738e-3, 0.8557274367e-3,
                0.3716900685e-3, -0.1308871233e-3, 0.3216895199e-4,
                0.2785881034e-6],
        "d1": [1, 1, 2, 2, 2, 2, 3, 5, 5, 6, 7, 8],
        "t1": [0, 2, 0, 1, 2, 3, 5, 0, 1, 3, 2, 5],

        "nr2": [-0.352151113, 0.7881914536e-1, -0.151966661e-1, -0.1068458586,
                -0.2055046288, 0.9146198012, 0.3213343569e-3, -0.1133591391e1,
                -0.3107520749, 0.1217901527e1, -0.4481710831, 0.5494218772e-1,
                -0.8665222096e-4, 0.3844084088e-1, 0.9853044884e-2,
                -0.1767598472e-1, 0.1488549222e-2, -0.3070719069e-2,
                0.388080328e-2, -0.2627505215e-2, 0.5258371388e-3, -0.1716396901,
                0.7188823624e-1, 0.5881268357e-1, -0.145593888e-1, -0.12161394e-1],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
               2, 3, 3, 3, 3, 3],
        "d2": [1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 11, 11, 11, 11,
               11, 2, 2, 3, 3, 5],
        "t2": [5, 7, 9, 5, 4, 6, 13, 5, 2, 3, 2, 0, 11, 1, 4, 0, 0, 3, 5, 6,
               7, 13, 14, 15, 24, 15],
        "gamma2": [1]*26}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for water of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.314371357587,
        "cp": Fi1,
        "ref": {"name": "CUSTOM",
            "Tref": Tt, "Pref": 611.655, "ho": 0.611872, "so": 0}, 

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [3.46821920e-1, 5.03423025e-1, -3.51059570e-1, 5.07004866e-2,
                1.99939129e-4, -5.69888763e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [-1.96198912e-1, -2.02509554, -1.09353609, 7.25785202e-2,
                2.16072642e-1, -1.01542630e-1, 7.46926106e-2, 2.18830463e-3],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, GERG, helmholtz2, helmholtz3, helmholtz4
    _PR = 0.0043451

    _surface = {"sigma": [0.2358, -0.147375], "exp": [1.256, 2.256]}
    _melting = {"Tmin": 251.165, "Tmax": 370.0}
    _sublimation = {"Tmin": 50.0, "Tmax": Tt}
#    _sublimation={"eq": 2, "Tref": 1, "Pref": 0.133332237, "a1": [-0.212144006e2, 0.273203819e2, -0.61059813e1], "exp1": [-0.9933333333, 0.206667, 0.703333], "a2": [], "exp2": [], "a3": [], "exp3": []}

    _vapor_Pressure = {
        "eq": 6,
        "ao": [-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719,
               1.80122502],
        "exp": [2, 3, 6, 7, 8, 15]}
    _liquid_Density = {
        "eq": 2,
        "ao": [1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352,
               -6.74694450e5],
        "exp": [1, 2, 5, 16, 43, 110]}
    _vapor_Density = {
        "eq": 4,
        "ao": [-2.0315024, -2.6830294, -5.38626492, -17.2991605, -44.7586581,
               -63.9201063],
        "exp": [1, 2, 4, 9, 18.5, 35.5]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "IAPWS (1997)",
              "__code__": (_Viscosity, )}

    visco1 = {"eq": 4, "omega": 1,
              "__name__": u"Quiñones-Cisneros (2006)",
              "__doi__": {"autor": "S.E.Quiñones-Cisneros and U.K. Deiters",
                          "title": "Generalization of the Friction Theory for Viscosity Modeling", 
                          "ref": "J. Phys. Chem. B, 2006, 110 (25), pp 12820–12834",
                          "doi": "10.1021/jp0618577"}, 

              "Tref": 647.096, "muref": 1.0,
              "ek": 809.1, "sigma": 0.2641, "n_chapman": 0,
              "n_ideal": [151.138, -444.318, 398.262, -81.7008],
              "t_ideal": [0, 0.25, 0.5, 0.75],

              "a": [-1.17407105202836e-5, -3.78854818708520e-7, 3.56742875797909e-8],
              "b": [1.62216397984014e-6, -8.36595322447571e-6, 9.10862531286788e-8],
              "c": [1.92706925578893e-5, -1.28679815491711e-5, 0.0],
              "A": [-3.30144899918610e-10, 0.0, 1.02931444103415e-11],
              "B": [5.03139997945133e-10, 1.82304182380560e-10, 0.0],
              "C": [8.01449084635477e-10, 5.65613687804585e-9, 1.10163426018591e-10],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1

    def _visco0(self, rho, T, fase):
        """IAPWS, Release on the IAPWS Formulation 2008 for the Viscosity of Ordinary Water Substance (International Association for the Properties of Water and Steam, 2008)"""
        ref = H2O()
        ref._ref("OTO")
        estado = ref._Helmholtz(rho, 1.5*647.096)
        drho = 1/estado["dpdrho"]*1e3
        return _Viscosity(rho, T, fase=fase, drho=drho)

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "IAPWS (1997)", 
               "__code__": (_ThCond, )}

    _thermal = thermo0,

    def _thermo0(self, rho, T, fase):
        """IAPWS, Release on the IAPWS Formulation 2011 for the Thermal Conductivity of Ordinary Water Substance"""
        ref = H2O()
        ref._ref("OTO")
        estado = ref._Helmholtz(rho, 1.5*647.096)
        drho = 1/estado["dpdrho"]*1e3
        return _ThCond(rho, T, fase, drho)
        
    def _Dielectric(self, rho, T):
        return unidades.Dimensionless(_Dielectric(rho, T))

    @classmethod
    def _Melting_Pressure(cls, T):
#        if 251.165 <= T <= 273.16:
#            Tita = T/cls.Tt
#            a = [0.119539337e7, 0.808183159e5, 0.33382686e4]
#            expo = [3., 0.2575e2, 0.10375e3]
#            suma = 1
#            for ai, expi in zip(a, expo):
#                suma += ai*(1-Tita**expi)
#            P1 = suma*611.657
#        else:
#            P1 = None
        
        if 251.165 <= T <= 256.164:
            Tref = 251.165
            Pref = 208566.
            Tita = T/Tref
            P2 = Pref*(1-0.299948*(1-Tita**60.))
        elif 256.164 < T <= 273.31:
            Tref = 256.164
            Pref = 350100.
            Tita = T/Tref
            P2 = Pref*(1-1.18721*(1-Tita**8.))
        elif 273.31 < T <= 355:
            Tref = 273.31
            Pref = 632400.
            Tita = T/Tref
            P2 = Pref*(1-1.07476*(1-Tita**4.6))
        elif 355. < T:
            Tref = 355
            Pref = 2216000.
            Tita = T/Tref
            P2 = Pref*exp(1.73683*(1-1./Tita)-0.544606e-1*(1-Tita**5)+0.806106e-7*(1-Tita**22))
        return unidades.Pressure(P2, "kPa")

    @classmethod
    def _Sublimation_Pressure(cls, T):
        Pref = 611.657
        Tita = T/cls.Tt
        a = [-0.212144006e2, 0.273203819e2, -0.61059813e1]
        expo = [0.333333333e-2, 1.20666667, 1.70333333]
        suma = 0
        for ai, expi in zip(a, expo):
            suma += ai*Tita**expi
        return unidades.Pressure(exp(suma/Tita)*Pref)


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()
    
#    st=H2O(T=300, x=0.5)
#    print "%0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (
#        st.T, st.P.kPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg,
#        st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, 
#        st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)


#    water=H2O(T=273.16, P=612)
#    print water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.cv.kJkgK, water.cp.kJkgK, water.w
#    water=H2O(T=273.16, P=612, eq="GERG")
#    print water.T, water.P.MPa, water.rho, water.h.kJkg, water.s, water.cv, water.cp, water.w
    

#    water=H2O(P=1e6, x=1)
#    print("%0.2f %0.4f %0.3f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))

#    water=H2O(s=5e3, P=3.5e6)
#    print("%0.2f %0.4f %0.3f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    515.71 3.5000 25.912 2222.8 5.0000 0.66921 -0.0085877

#    water=H2O(rho=300., P=0.1e6)
#    print("%0.2f %0.5f %0.2f %0.2f %0.4f %0.7f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    372.76 0.10000 300.00 420.56 1.3110 0.0013528 -0.025144

#    water=H2O(h=1500e3, P=0.1e6)
#    print("%0.2f %0.5f %0.4f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    372.76 0.10000 1.2303 1500.0 4.2068 0.47952 -0.025144

#    cyc5=H2O(T=500., rho=838.025, recursion=False)
#    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
#    print cyc5.cp.kJkgK, cyc5.cp0.kJkgK

#    water = H2O(T=300, x=0.5)
#    print water.v0, water.rho0
#    print water.h0, water.u0, water.s0
#    print water.a0, water.g0, water.gamma0
#    print water.cp0.kJkgK, water.cv0.kJkgK, water.cp0_cv


#    water=H2O(T=300, rho=1000)
#    print water.joule.Katm, water.virialB, water.virialC, water.Hvap.kJkg
#    0.0992418350 4.13018112 1501.51914 0.393062643
#
#    >>> water=H2O(T=500, rho=0.435)
#    >>> print("%0.10f %0.8f %0.5f %0.9f" % (water.P, water.cv, water.w, water.s))
#    0.0999679423 1.50817541 548.31425 7.944882714
#
#    water=H2O(T=647.096, rho=100.0)
#    print water.P.MPa, water.rho, water.h, water.s
#    870.7690 2.66422350 2019.33608 4.17223802
#
#    >>> water=H2O(T=300., P=0.1)
#    >>> print("%0.2f %0.5f %0.2f %0.2f %0.5f %0.4f %0.1f %0.6f" % (water.T, water.P, water.rho, water.h, water.s, water.cp, water.w, water.virialB))
#    300.00 0.10000 996.56 112.65 0.39306 4.1806 1501.5 -0.066682
#
#    >>> water=H2O(T=500., P=0.1)
#    >>> print("%0.2f %0.5f %0.5f %0.1f %0.4f %0.4f %0.2f %0.7f" % (water.T, water.P, water.rho, water.h, water.s, water.cp, water.w, water.virialB))
#    500.00 0.10000 0.43514 2928.6 7.9447 1.9813 548.31 -0.0094137
#
#    >>> water=H2O(T=450., x=0.5)
#    >>> print("%0.2f %0.5f %0.4f %0.1f %0.4f %0.6f" % (water.T, water.P, water.rho, water.h, water.s, water.virialB))
#    450.00 0.93220 9.5723 1761.8 4.3589 -0.013028

#    print H2O._Melting_Pressure(252)

#{'P': 10513642.277577657, 'T': 600.0}
#{'P': 12996040.325499836, 'T': 600.0}
#T 600.0 1
#{'P': 16064562.56384358, 'T': 600.0}
#/usr/lib/python2.7/dist-packages/scipy/optimize/minpack.py:152: RuntimeWarning: The iteration is not making good progress, as measured by the 
#  improvement from the last ten iterations.
#  warnings.warn(msg, RuntimeWarning)
#{'P': 19857600.000000007, 'T': 600.0}
#T 600.0 1
#{'P': 19898125.714285713, 'T': 600.0}

#{'P': 3338.5258541981202, 'T': 600.0}
#{'P': 20000000.0, 'T': 639.43649632653057}
#    water=H2O(P=2384230.3837643899, T=700.0)
#    print water.P.MPa, water.T, water.rho, water.x

#    water=H2O(T=600, x=0)
#    print water.P, water.T, water.x    
#    water=H2O(P=12996040.325499836, x=0.5)
#    print water.P, water.T, water.x

#    water=H2O(T=500., u=900e3)
#    print water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB
#    500.00 108.21 903.62 900.00 1019.8 2.4271 4.1751 1576.0 -0.0094137
#
#    >>> water=H2O(P=0.3, u=1550.)
#    >>> print("%0.2f %0.5f %0.4f %0.1f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    406.67 0.30000 3.3029 1550.0 1640.8 4.3260 0.49893 -0.018263
#
#    >>> water=H2O(rho=300, h=1000.)
#    >>> print("%0.2f %0.4f %0.2f %0.2f %0.1f %0.4f %0.6f %0.7f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    494.92 2.3991 300.00 992.00 1000.0 2.6315 0.026071 -0.0097064
#
#    >>> water=H2O(rho=30, s=8.)
#    >>> print("%0.2f %0.3f %0.3f %0.1f %0.1f %0.4f %0.4f %0.2f %0.9f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.cp.kJkgK, water.w, water.virialB))
#    1562.42 21.671 30.000 4628.5 5350.9 8.0000 2.7190 943.53 0.000047165
#
#    >>> water=H2O(rho=30, s=4.)
#    >>> print("%0.2f %0.4f %0.3f %0.1f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    495.00 2.4029 30.000 1597.3 1677.4 4.0000 0.39218 -0.0097015
#
#    >>> water=H2O(rho=300, u=1000.)
#    >>> print("%0.2f %0.4f %0.3f %0.1f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    496.44 2.4691 300.000 1000.0 1008.2 2.6476 0.02680 -0.0096173
#{'s': 11000.0, 'T': 411.99634285714285}
#{'P': 22064000.0, 'T': 600.0}
#79.6790793628

#    water=H2O(P=1e5, h=5000000)
#    print water.P.bar, water.T, water.rho, water.x, water.s.kJkgK, water.h.kJkg
#    water2=H2O(P=water.P, T=water.T)
#    print water2.P.bar, water2.T, water2.rho, water2.x, water2.s.kJkgK, water2.h.kJkg
#    water2=H2O(s=water.s, T=water.T)
#    print water2.P.bar, water2.T, water2.rho, water2.x, water2.s.kJkgK, water2.h.kJkg

#    water=H2O(P=1e5, x=0.5)
#    print water.P.bar, water.T, water.rho, water.x, water.s.kJkgK, water.h.kJkg
#    water=H2O(P=water.P, h=water.h)
#    print water.P.bar, water.T, water.rho, water.x, water.s.kJkgK, water.h.kJkg
#    water2=H2O(s=water.s, T=water.T)
#    print water2.P.bar, water2.T, water2.rho, water2.x, water2.s.kJkgK, water2.h.kJkg

#    water=H2O(T=280, x=0)
#    print water.P, water.T, water.rho, water.x, water.h
    
#    water=H2O(P=612.48, T=650)
#    print("%0.2f %0.6f %0.5f %0.2f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    345.73 0.034850 0.73526 952.60 1000.0 3.0000 0.29920 -0.034124
#
#    >>> water=H2O(u=995., h=1000.)
#    >>> print("%0.2f %0.4f %0.2f %0.2f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    501.89 2.7329 546.58 995.00 1000.0 2.6298 0.00866 -0.009308
#
#    >>> water=H2O(u=1000., s=3.)
#    >>> print("%0.2f %0.6f %0.5f %0.2f %0.1f %0.4f %0.5f %0.6f" % (water.T, water.P.MPa, water.rho, water.u.kJkg, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
#    371.24 0.094712 1.99072 1000.00 1047.6 3.0000 0.28144 -0.025543

#    wt=H2O()
#    tau=wt.Tc/500
#    delta=838.025/wt.rhoc
#    print wt._phi0(wt._constants["cp"], tau, delta)

    for eq in (0, 1, 3, 4):
        st=H2O(T=500, x=0.5,  eq=eq)
        print "%0.6g %0.5g %0.5g %0.5g %0.2f %0.6g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g %0.5g" % (\
            st.T, st.P.MPa, st.Liquido.rho, st.Gas.rho, st.Liquido.h.kJkg, st.Gas.h.kJkg, \
            st.Liquido.s.kJkgK, st.Gas.s.kJkgK, st.Liquido.cv.kJkgK, st.Gas.cv.kJkgK, \
            st.Liquido.cp.kJkgK, st.Gas.cp.kJkgK, st.Liquido.w, st.Gas.w)
