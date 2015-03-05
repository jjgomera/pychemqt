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
    >>> water=H2O(s=5e3, P=3.5e6)
    >>> print("%0.2f %0.4f %0.3f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))
    515.71 3.5000 25.912 2222.8 5.0000 0.66921 -0.0085877

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

    CP1 = {"ao": 4.00632,
           "an": [], "pow": [],
           "ao_exp": [0.012436, 0.97315, 1.27950, 0.96956, 0.24873],
           "exp": [833, 2289, 5009, 5982, 17800],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 4.00392,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.01059, -0.98763, 3.06904, 0],
           "hyp": [0.415386589*Tc, 1.763895929*Tc, 3.874803739*Tc, 0]}
           
    Fi0 = {"ao_log": [1, 3.00632],
           "pow": [0, 1],
           "ao_pow": [-8.3204464837497, 6.6832105275932],
           "ao_exp": [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
           "titao": [1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": u"Helmholtz equation of state for water of Wagner and Pruß (2002).",
        "__doc__":  u"""Wagner, W., Pruß, A. The IAPWS formulation 1995 for the thermodyamic properties of ordinary water substance for general and scientific use. J. Phys. Chem. Ref. Data 31 (2002), 387 – 535.""",
        "R": 8.314371357587,
#        "cp": CP1,
        "cp": Fi0,

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
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP2,

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
        "__name__": "Helmholtz equation of state for water of Saul and Wagner (1989).",
        "__doc__":  u"""Saul, A. and Wagner, W., "A Fundamental Equation for Water Covering the Range From the Melting Line to 1273 K at Pressures up to 25000 MPa," J. Phys. Chem. Ref. Data, 18(4):1537-1564, 1989.""",
        "R": 8.31434,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 1273., "Pmax": 400000.0, "rhomax": 55.49, 
        "Pmin": 0.61166, "rhomin": 55.497, 

        "nr1": [0.2330009013, -0.1402091128e-1, 0.1172248041, -0.1850749499,
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
               7, 13, 14, 13, 24, 15],
        "gamma2": [1]*26}

    eq = helmholtz1, GERG, helmholtz2
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
              "__doc__": u"""S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": u"Quiñones-Cisneros (2006)",
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
    
#    water=H2O(T=1500, s=13000.)
#    print("%0.2f %0.4f %0.3f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))

    water=H2O(P=1e6, x=1)
    print("%0.2f %0.4f %0.3f %0.1f %0.4f %0.5f %0.7f" % (water.T, water.P.MPa, water.rho, water.h.kJkg, water.s.kJkgK, water.x, water.virialB))

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
