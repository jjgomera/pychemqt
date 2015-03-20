#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class nC5(MEoS):
    """Multiparameter equation of state for n-pentane

    >>> pentano=nC5(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (pentano.T, pentano.rho, pentano.u.kJkg, pentano.h.kJkg, pentano.s.kJkgK, pentano.cv.kJkgK, pentano.cp.kJkgK, pentano.w)
    300.0 618.99 -22.72 -22.56 0.44047 1.7397 2.3237 1002.8
    """
    name = "pentane"
    CASNumber = "109-66-0"
    formula = "CH3-(CH2)3-CH3"
    synonym = "R-601"
    rhoc = unidades.Density(232.)
    Tc = unidades.Temperature(469.7)
    Pc = unidades.Pressure(3370.0, "kPa")
    M = 72.14878  # g/mol
    Tt = unidades.Temperature(143.47)
    Tb = unidades.Temperature(309.21)
    f_acent = 0.251
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 8
    _Tr = unidades.Temperature(449.271155)
    _rhor = unidades.Density(233.873368)
    _w = 0.247058753

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.2857243e6, 0.1542722e8, 0.1051520e9, 0],
           "hyp": [0.17867e3, 0.840538e3, 0.177425e4, 0]}
           
    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [8.95043, 21.836, 33.4032, 0],
           "hyp": [0.380391739, 1.789520971, 3.777411113, 0]}

    CP2 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [8.95043, 21.836, 33.4032, 0],
           "hyp": [0.380391739*Tc, 1.789520971*Tc, 3.777411113*Tc, 0]}

    CP3 = {"ao": 10.288132,
           "an": [-0.2695377e-1, 0.20951065e-3, -0.27910773e-6, 0.12266269e-9],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 22.5012/8.3159524*4.184,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [2.057417e8/8.3159524*4.184, 2.972927e7/8.3159524*4.184, 0, 0],
           "hyp": [1.71958e3, 8.02069e2, 0, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for pentane of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.", 
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"}, 
        "__test__": """
            >>> st=nC5(T=700, rho=200)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            3.2053 13.454 3.6052
            >>> st2=nC5(T=750, rho=100)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            213.42 0.34915
            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 11.2, 
        "Pmin": 0.76322e-4, "rhomin": 10.566, 

        "nr1": [0.10968643e1, -0.29988888e1, 0.99516887, -0.16170709,
                0.11334460, 0.26760595e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40979882, -0.40876423e-1, -0.38169482, -0.10931957,
                -0.32073223e-1, 0.16877016e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for \
                    Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": 143.47, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 10.57, 
#        "Pmin": 73.476, "rhomin": 29.249, 

        "nr1": [0.10968643098001e1, -0.29988888298061e1, 0.99516886799212,
                -0.16170708558539, 0.11334460072775, 0.26760595150748e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.40979881986931, -0.40876423083075e-1, -0.38169482469447,
                -0.10931956843993, -0.32073223327990e-1, 0.16877016216975e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP3,

        "Tmin": 238.0, "Tmax": 573.0, "Pmax": 30000.0, "rhomax": 9.410819, 
        "Pmin": 3.624503, "rhomin": 9.3861, 

        "nr1": [-0.117648900900e1, 0.163499095773e1, -0.366669005817,
                0.724947274043, -0.221919300269e1, 0.188671490348e1,
                -0.195774652096e1, 0.308440851184, 0.437424419722,
                -0.625853472351, 0.382868807091, -0.119467393955, 0.218631441082,
                0.485668874195e-1, -0.132198161379, 0.213549844850e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.117648900900e1, -0.163499095773e1, 0.366669005817,
                -0.363660829618e-2, 0.633672105685, -0.705792643982],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.968832]*6}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Starling (1973)",
        "__doc__":  u"""Starling, K.E., "Fluid Thermodynamic Properties for Light Petroleum Systems," Gulf Publishing Company, 1973.""",
        "R": 8.3159524,
        "cp": CP4,

        "Tmin": 177.0, "Tmax": 589.0, "Pmax": 55000.0, "rhomax": 10.2534, 
        "Pmin": 0.011064, "rhomin": 10.253, 

        "nr1": [0.175873733594e1, 0.485604047435, -0.111896446456e1,
                -0.685918143315, 0.368714111378e-1, -0.167498784887e-2,
                0.327765295239, -0.352742092747, -0.999487301826e-1,
                0.781999120830e-2, 0.221577806386e-2],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.175873733594e1, -0.411653507564],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.46812392]*2}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for pentane of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [2.20261753, 1.07797592, -3.82130221, 1.06627357e-1,
                3.07513215e-4, -2.84309667e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],
  
        "nr2": [-7.28441220e-2, -4.60943732e-1, 8.39360011e-2 , -1.50650444e-2,
                -2.03771872e-1, -7.90244277e-3, -5.68993564e-2, -2.99387974e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, GERG, helmholtz3, helmholtz4, helmholtz5

    _surface = {"sigma": [0.0562267, -0.0037496, -0.0029861],
                "exp": [1.25, 2.25, 3.25]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.10924],  "expt0": [-1.], "expd0": [1.],
                   "a1": [25.39, 0.025], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [78.39, 54.15, -12480, -4800.0],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.76322e-4,
                "Tmin": Tt, "Tmax": 2000.0,
                "a1": [-8647500000, 8647500001], "exp1": [0, 1.649],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73918e1, 0.31102e1, -0.22415e1, -0.31585e1, -0.90451],
        "exp": [1., 1.5, 1.74, 3.75, 8.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.10178e1, 0.42703, 0.11334e1, 0.41518, -0.47950e-1],
        "exp": [0.27, 0.44, 0.6, 4.0, 5.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.29389e1, -0.62784e1, -0.19941e2, -0.16709e2, -0.36543e2, -0.12799e3],
        "exp": [0.4, 1.18, 3.2, 6.6, 7.0, 15.0]}

    visco0 = {"eq": 2, "omega": 3,
              "__name__": "NIST14",
              "__doc__": """Coefficients are taken from NIST14, Version 9.08""",
              "ek": 341.10, "sigma": 0.5784,
              "n_chapman": 0.226720214/M**0.5,
              "F": [0, 0, 0, 100],
              "E": [-13.47938293, 1176.6275165, 14.2278439927, -21951.0293411,
                    0.03766867689, 70.1529173825, 21435.7720323],
              "rhoc": 3.215}

    visco1 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quiñones-Cisneros (2006)",
              "Tref": 469.7, "muref": 1.0,
              "ek": 341.1, "sigma": 0.5784, "n_chapman": 0,
              "n_ideal": [17.6805, -55.6942, 48.7177],
              "t_ideal": [0, 0.25, 0.5],

              "a": [1.08193e-5, -4.71699e-5, 0.0],
              "b": [1.21502e-4, -9.84766e-5, 0.0],
              "c": [5.08307e-5, -1.07e-5, 0.0],
              "A": [-2.10025e-10, -1.56583e-9, 0.0],
              "B": [1.98521e-8, 2.05972e-9, 0.0],
              "C": [-1.18487e-7, 1.69571e-7, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 1,
               "__name__": "NIST14",
               "__doc__": """Coefficients are taken from NIST14, Version 9.08""",

               "Tref": 341.1, "kref": 1e-3,
               "no": [1.35558587, -0.15569137, 1],
               "co": [0, -1, -96],

               "Trefb": 469.69, "rhorefb": 3.215, "krefb": 1e-3,
               "nb": [18.6089331038, -5.83657061299, 3.48987100529,
                      0.704467355508, -0.206501417728, -0.22307039402],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.9345e-9, "Tcref": 704.55}

    _thermal = thermo0,
