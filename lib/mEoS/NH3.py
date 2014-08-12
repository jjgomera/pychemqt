#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class NH3(MEoS):
    """Multiparameter equation of state for ammonia

    >>> amoniaco=NH3(T=500, P=50)
    >>> print "%0.1f %0.2f %0.1f %0.1f %0.4f %0.4f %0.4f %0.2f" % (amoniaco.T, amoniaco.rho, amoniaco.u.kJkg, amoniaco.h.kJkg, amoniaco.s.kJkgK, amoniaco.cv.kJkgK, amoniaco.cp.kJkgK, amoniaco.w)
    500.0 347.40 1330.2 1474.1 4.2087 2.6614 5.4776 712.88
    """
    name = "ammonia"
    CASNumber = "7664-41-7"
    formula = "NH3"
    synonym = "R-717"
    rhoc = unidades.Density(225.)
    Tc = unidades.Temperature(405.40)
    Pc = unidades.Pressure(11333.0, "kPa")
    M = 17.03026  # g/mol
    Tt = unidades.Temperature(195.495)
    Tb = unidades.Temperature(239.823)
    f_acent = 0.25601
    momentoDipolar = unidades.DipoleMoment(1.470, "Debye")
    id = 63

    CP1 = {"ao": 0,
           "an": [2.54985265683/405.40**-0.333, 4.86079045595/405.40**1.5,
                  -2.74637680305/405.40**1.75],
           "pow": [-1./3, 1.5, 1.75],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 5.111814,
           "an": [-0.42966650e2, -0.10243792e-1, 0.38750775e-4, -0.46406097e-7,
                  0.20268561e-10],
           "pow": [-1.001, 1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Tillner-Roth et al. (1993)",
        "__doc__":  u"""Tillner-Roth, R., Harms-Watzenberg, F., Baehr, H.D. Eine neue Fundamentalgleichung für Ammoniak. DKV-Tagungsbericht 20 (1993), 167 – 181.""",
        "R": 8.314471,
        "cp": CP1,

        "nr1": [-0.1858814e01, 0.4554431e-1, 0.7238548, 0.1229470e-1,
                0.2141882e-10],
        "d1": [1, 2, 1, 4, 15],
        "t1": [1.5, -0.5, 0.5, 1., 3.],

        "nr2": [-0.1430020e-1, 0.3441324, -0.2873571, 0.2352589e-4,
                -0.3497111e-1, 0.1831117e-2, 0.2397852e-1, -0.4085375e-1,
                0.2379275, -0.3548972e-1, -0.1823729, 0.2281556e-1,
                -0.6663444e-2, -0.8847486e-2, 0.2272635e-2, -0.5588655e-3],
        "d2": [3, 3, 1, 8, 2, 8, 1, 1, 2, 3, 2, 4, 3, 1, 2, 4],
        "t2": [0, 3, 4, 4, 5, 5, 3, 6, 8, 8, 10, 10, 5, 7.5, 15, 30],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3],
        "gamma2": [1]*16}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ammonia of Ahrendts and Baehr (1979)",
        "__doc__":  u"""Ahrendts, J. and Baehr, H.D. "The Thermodynamic Properties of Ammonia," VDI-Forsch., Number 596, 1979.""",
        "R": 8.31434,
        "cp": CP2,

        "nr1": [0.911447599671, -0.382129415537e1, 0.147730246416e1,
                0.580205129871e-1, -0.574413226616e-3, 0.153018094697,
                -0.256626062036, 0.445448838055, -0.1533210545,
                0.527996725202e-1, -0.484726581121e-1, 0.246579503330e-2,
                -0.107999941003e-3, -0.215298673010e-4, -0.306938893790e-4,
                0.839163613582e-5, 0.814833533876e-6, -0.314753664228e-7],
        "d1": [1, 1, 1, 1, 1, 2, 2, 3, 4, 5, 5, 7, 9, 9, 10, 11, 12, 14],
        "t1": [1, 2, 3, 6, 9, 0, 4, 2, 1, 1, 2, 3, 3, 5, 1, 1, 5, 5],

        "nr2": [0.642978802435, -0.139510669941e1, 0.956135683432,
                -0.272787386366, -0.189305337334e1, 0.479043603913e1,
                -0.245945016980e1, -0.121107723958e1, 0.500552271170e1,
                -0.615476024667e1, 0.210772481535e1, 0.298003513465,
                -0.152506723279, 0.115565883925e-2, -0.911244657201e-3,
                0.100587210000e-1, -0.120983155888e-1, 0.382694351151e-2],
        "d2": [1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 0, 0, 0, 0, 0],
        "t2": [2, 5, 6, 7, 5, 6, 7, 3, 4, 5, 6, 6, 7, 1, 2, 0, 1, 2],
        "c2": [2]*18,
        "gamma2": [0.86065403]*13+[506.2670781840292]*2+[50626.70781840292]*3}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for ammonia of Span and Wagner (2003)",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. III. Results for Polar Fluids," Int. J. Thermophys., 24(1):111-162, 2003.""",
        "R": 8.31451,
        "cp": CP1,

        "nr1": [0.7302272, -0.11879116e1, -0.68319136, 0.40028683e-1, 0.90801215e-4],
        "d1": [1, 1, 1, 3, 7],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.56216175e-1, 0.44935601, 0.29897121e-1, -0.18181684,
                -0.9841666e-1, -0.55083744e-1, -0.88983219e-2],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1, helmholtz2, helmholtz3

    _melting = {"eq": 1, "Tref": Tt, "Pref": 1000,
                "a1": [], "exp1": [], "a2": [], "exp2": [],
                "a3": [0.2533125e4], "exp3": [1]}
    _surface = {"sigma": [0.12190, 0.04015], "exp": [1.26, 2.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.70993e1, -0.24330e1, 0.87591e1, -0.64091e1, -0.21185e1],
        "exp": [1., 1.5, 1.7, 1.95, 4.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.34488e2, -0.12849e3, 0.17382e3, -0.10699e3, 0.30339e2],
        "exp": [0.58, 0.75, 0.9, 1.1, 1.3]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.38435, -0.40846e1, -0.66634e1, -0.31881e2, 0.21306e3, -0.24648e3],
        "exp": [0.218, 0.55, 1.5, 3.7, 5.5, 5.8]}

    visco0 = {"eq": 1, "omega": 1,
             "collision": [4.99318220, -0.61122364, 0.0, 0.18535124, -0.11160946],
             "__name__": "Fenghour (1995)",
             "__doc__": """Fenghour, A., Wakeham, W.A., Vesovic, V., Watson, J.T.R., Millat, J., and Vogel, E., "The viscosity of ammonia," J. Phys. Chem. Ref. Data, 24:1649-1667, 1995.""",
             "ek": 386., "sigma": 0.2957,
             "Tref": 1., "rhoref": 1.*M,
             "n_chapman": 8.8135503/M**0.5,

             "n_virial": [-0.17999496e1, 0.46692621e2, -0.53460794e3,
                          0.33604074e4, -0.13019164e5, 0.33414230e5,
                          -0.58711743e5, 0.71426686e5, -0.59834012e5,
                          0.33652741e5, -0.1202735e5, 0.24348205e4, -0.20807957e3],
             "t_virial": [0, -0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5, -5,
                          -5.5, -6],
             "Tref_virial": 386., "etaref_virial": 0.015570557,

             "Tref_res": 386., "rhoref_res": 1.*M, "etaref_res": 1,
             "n_poly": [2.19664285e-1, -0.83651107e-1, 0.17366936e-2,
                        -0.64250359e-2, 1.67668649e-4, -1.49710093e-4, 0.77012274e-4],
             "t_poly": [-2, -4, -0, -1, -2, -3, -4],
             "d_poly": [2, 2, 3, 3, 4, 4, 4],
             "g_poly": [0, 0, 0, 0, 0, 0, 0],
             "c_poly": [0, 0, 0, 0, 0, 0, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1, "critical": "NH3",
               "__name__": "Tufeu (1984)",
               "__doc__": """Tufeu, R., Ivanov, D.Y., Garrabos, Y., and Le Neindre, B., "Thermal conductivity of ammonia in a large temperature and pressure range including the critical region," Ber. Bunsenges. Phys. Chem., 88:422-427, 1984""",

               "Tref": 1., "kref": 1.,
               "no": [0.3589e-1, -0.1750e-3, 0.4551e-6, 0.1685e-9, -0.4828e-12],
               "co": [0, 1, 2, 3, 4],

               "Trefb": 1., "rhorefb": 0.05871901, "krefb": 1.,
               "nb": [0.16207e-3, 0.12038e-5, -0.23139e-8, 0.32749e-11],
               "tb": [0, 0, 0, 0],
               "db": [1, 2, 3, 4],
               "cb": [0, 0, 0, 0]}

    _thermal = thermo0,
