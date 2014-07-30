#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class H2S(MEoS):
    """Ecuación de estado de multiparametros para el sulfuro de hidrogeno

    >>> sulfuro=H2S(T=500, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (sulfuro.T, sulfuro.rho, sulfuro.u.kJkg, sulfuro.h.kJkg, sulfuro.s.kJkgK, sulfuro.cv.kJkgK, sulfuro.cp.kJkgK, sulfuro.w)
    500.0 0.82094 723.06 844.868 3.45347 0.8487 1.0944 396.06
    """
    name="hydrogen sulfide"
    CASNumber="7783-06-4"
    formula="H2S"
    synonym=""
    rhoc=unidades.Density(347.28)
    Tc=unidades.Temperature(373.1)
    Pc=unidades.Pressure(9000.0, "kPa")
    M=34.08088      #g/mol
    Tt=unidades.Temperature(187.7)
    Tb=unidades.Temperature(212.85)
    f_acent=0.1005
    momentoDipolar=unidades.DipoleMoment(0.97, "Debye")
    id=50

    CP1={  "ao": 4,
                "an": [0.14327e-5], "pow": [1.5],
                "ao_exp": [1.1364, 1.9721],
                "exp": [1823, 3965],
                "ao_hyp": [], "hyp": []}

    CP2={  "ao": 4,
                "an": [], "pow": [],
                "ao_exp": [0.9767422, 2.151898],
                "exp": [4.506266*Tc, -10.15526*Tc],
                "ao_hyp": [], "hyp": []}

    CP3={  "ao": 4.1012105,
                "an": [-0.16720073e-2, 0.75303152e-5, -0.62421053e-8, 0.18098453e-11],
                "pow": [1, 2, 3, 4],
                "ao_exp": [], "exp": [],
                "ao_hyp": [], "hyp": []}

    CP4={  "ao": 7.9468/8.3159524*4.184,
                "an": [], "pow": [],
                "ao_exp": [], "exp": [],
                "ao_hyp": [-1.5769761e4/8.3159524*4.184, 2.0329947e6/8.3159524*4.184, 1.3861204e7/8.3159524*4.184, -3.5044957e6/8.3159524*4.184],
                "hyp": [4.33801e2, 8.43792e2, 1.48143e3, 1.10223e3]}

    CP5={  "ao": 4.0,
                "an": [], "pow": [],
                "ao_exp": [], "exp": [],
                "ao_hyp": [3.11942, 1.00243, 0, 0],
                "hyp": [4.914580541*Tc, 2.270653980*Tc, 0, 0]}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for hydrogen sulfide of Lemmon and Span (2006).",
        "__doc__":  u"""Lemmon, E.W., Span, R. Short fundamental equations of state for 20 industrial fluids. J. Chem. Eng. Data 51 (2006), 785 – 850.""",
        "R": 8.314472,
        "cp": CP1,

        "nr1":  [0.87641, -2.0367, 0.21634, -0.050199, 0.066994, 0.00019076],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.20227, -0.0045348, -0.22230, -0.034714, -0.014885, 0.0074154],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz2={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of Sakoda and Uematsu (2004)",
        "__doc__":  u"""Sakoda, N., Uematsu, M. "A Thermodynamic Property Model for Fluid Phase Hydrogen Sulfide," Int. J. Thermophys., 25(3):709-737, 2004.""",
        "R": 8.314472,
        "cp": CP2,

        "nr1":  [0.1545780, -0.1717693e1, -0.1595211e1, 0.2046589e1, -0.1690358e1, 0.9483623, -0.6800772e-1, 0.4372273e-2, 0.3788552e-4, -0.3680980e-4, 0.8710726e-5],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4, 8, 9, 10],
        "t1": [0.241, 0.705, 1, 0.626, 1.15, 1.63, 0.21, 3.08, 0.827, 3.05, 3.05],

        "nr2": [0.6886876, 0.2751922e1, -0.1492558e1, 0.9202832, -0.2103469, 0.1084359e-2, 0.3754723e-1, -0.5885793e-1, -0.2329265e-1, -0.1272600e-3, -0.1336824e-1, 0.1053057e-1],
        "d2": [1, 1, 1, 2, 5, 1, 4, 4, 3, 8, 2, 3],
        "t2": [0.11, 1.07, 1.95, 0.142, 2.13, 4.92, 1.75, 3.97, 11.8, 10, 9.83, 14.2],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4],
        "gamma2": [1]*12}

    helmholtz3={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP3,

        "nr1": [0.135782366339e1, -0.153224981014e1, 0.329107661253, 0.195802782279e1, -0.301125182071e1, -0.126614059078e1, 0.129960331548e1, -0.185645977138, -0.160919744092e1, 0.234395817019e1, -0.378573094883, 0.758423219040, -0.973372615169, -0.120786235447, 0.209004959689, -0.919656385346e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.135782366339e1, 0.153224981014e1, -0.329107661253, 0.891427552242, -0.204776100441e1, 0.101366381241e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.9873538]*6}

    helmholtz4={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for hydrogen sulfide of Starling (1973)",
        "__doc__":  u"""Starling, K.E., "Fluid Thermodynamic Properties for Light Petroleum Systems," Gulf Publishing Company, 1973.""",
        "R": 8.3159524,
        "cp": CP4,

        "nr1":  [0.110928333109e1, 0.188834546108, -0.930906931583, -0.411249591635, 0.140676923412e-1, -0.169077883177e-4, 0.510265859853, -0.572402742986, -0.828859606622e-3, 0.971664064871e-2, 0.140700425434e-4],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.110928333109e1, -0.26913741657],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.48524558]*2}

    GERG={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for propane of Kunz and Wagner (2008).",
        "__doc__":  u"""Kunz, O.; Wagner, W. -- The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures- An Expansion of GERG-2004. J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091""",
        "R": 8.314472,
        "cp": CP5,

        "nr1":  [0.87641, -0.20367e1, 0.21634, -0.50199e-1, 0.66994e-1, 0.19076e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.20227, -0.45348e-2, -0.22230, -0.34714e-1, -0.14885e-1, 0.74154e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq=helmholtz1, helmholtz2, helmholtz3, helmholtz4, GERG

    _surface={"sigma": [0.082], "exp": [1.26]}
    _vapor_Pressure={ "eq": 5, "ao": [-6.5884, 2.1582, -1.6054, -2.3870 ], "exp": [1, 1.5, 2, 4.8]}
    _liquid_Density={ "eq": 2, "ao": [11.833, -17.019, 7.8047], "exp": [1.63, 1.95, 2.3]}
    _vapor_Density={ "eq": 3, "ao": [-3.9156, -7.7093, -19.543, -49.418], "exp": [0.49, 1.69, 4, 8]}

    visco0={"eq": 2, "omega": 1,
                "__name__": "NIST",
                "__doc__": """Coefficients are taken from NIST14, Version 9.08""",
                "ek": 301.1, "sigma": 0.36237,
                "n_chapman": 0.1558117/M**0.5,
                "F": [0, 0, 0, 100.],
                "E": [-12.328630418994, 782.29421491, 11.840322553, -10401.582791, -0.0482407464, 69.709031672, 256.31792390],
                "rhoc": 10.2}

    visco1={"eq": 4, "omega": 1,
                "__doc__": """K. A. G. Schmidt, J.J.Carroll, S.E.Quinones-Cisneros and B. Kvamme, "Hydrogen Sulphide Viscosity Model", proceedings of the 86th Annual GPA Convention, March 11-14,2007, San Antonio,TX.""",
                "__name__": "Schmidt (2007)",
                "Tref": 373.1, "muref": 1.0,
                "ek": 301.1, "sigma": 0.36237, "n_chapman": 0,
                "n_ideal": [4.36694e1, -12.1530e1, 9.35279e1],
                "t_ideal": [0, 0.25, 0.5],

                "a": [5.46919e-5, -7.32295e-6, -7.35622e-6],
                "b": [4.56159e-5, -1.82572e-5, -6.59654e-6 ],
                "c": [-4.33882e-6, 6.13716e-6, 0.0],
                "A": [6.67324e-9, -2.16365e-9, 0.0],
                "B": [-1.53973e-9, 2.17652e-9, 0.0],
                "C": [3.54228e-7, -4.76258e-8, 0.0],
                "D": [0.0, 0.0, 0.0 ]}

    _viscosity=visco0, visco1

    thermo0={"eq": 1,
                "__name__": "NIST14",
                "__doc__": """Coefficients are taken from NIST14, Version 9.08""",

                "Tref": 301.1, "kref": 1e-3,
                "no": [1.35558587, -0.1402163, 1],
                "co": [0, -1, -96],

                "Trefb": 373.4, "rhorefb": 10.2, "krefb": 1e-3,
                "nb": [21.7827447865, 10.8880930411, -7.45794247629, 3.65609005216, 1.89362258187, -1.10975687736],
                "tb": [0, 0, 0, -1, 0, -1],
                "db": [1, 3, 4, 4, 5, 5],
                "cb": [0]*6,

                "critical": 3,
                "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
                "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.3211e-9, "Tcref": 559.65}

    _thermal=thermo0,
    
