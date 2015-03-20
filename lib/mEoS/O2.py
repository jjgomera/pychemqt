#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class O2(MEoS):
    """Multiparameter equation of state for oxygen

    >>> oxigeno=O2(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (oxigeno.T, oxigeno.rho, oxigeno.u.kJkg, oxigeno.h.kJkg, oxigeno.s.kJkgK, oxigeno.cv.kJkgK, oxigeno.cp.kJkgK, oxigeno.w)
    300.0 1.28367 -77.28 0.617 -0.00121 0.6587 0.9199 329.72
    """
    name = "oxygen"
    CASNumber = "7782-44-7"
    formula = "O2"
    synonym = "R-732"
    rhoc = unidades.Density(436.14)
    Tc = unidades.Temperature(154.581)
    Pc = unidades.Pressure(5043.0, "kPa")
    M = 31.9988  # g/mol
    Tt = unidades.Temperature(54.361)
    Tb = unidades.Temperature(90.1878)
    f_acent = 0.0222
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 47
    _Tr = unidades.Temperature(150.875090)
    _rhor = unidades.Density(439.519141)
    _w = 0.023479051

    CP1 = {"ao": 3.51808732,
           "an": [], "pow": [],
           "ao_exp": [0.102323928e1, 0.784357918, 0.337183363e-2,
                      -.170864084e-1, 0.463751562e-1],
           "exp": [2246.32440, 11259.9763, 1201.26209, 69.0089445, 5328.05445],
           "ao_hyp": [], "hyp": []}

    Fi2 = {"ao_log": [1, 2.50146],
           "pow": [0, 1],
           "ao_pow": [10.001843586, -14.996095135],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [1.07558, 1.01334, 0, 0],
           "hyp": [14.461722565, 7.223325463, 0, 0]}

    CP2 = {"ao": 3.521876773671,
           "an": [-0.4981998537119e4, 0.2302477799952e3, -0.3455653235107e1,
                  0.3521876773671e1, -0.4354202160244e-4, 0.1346353450132e-7,
                  0.1620598259591e-10],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [1.031468515726], "exp": [2239.18105],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for oxygen of Schmidt and Wagner (1985).",
        "__doc__":  u"""Schmidt, R., Wagner, W. A new form of the equation of state for pure substances and its application to oxygen. Fluid Phase Equuilibria. 19 (1985), 175 – 200.""",
        "R": 8.31434,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 82000.0, "rhomax": 43.348, 
        "Pmin": 0.14628, "rhomin": 40.816, 

        "nr1": [0.39837687490, -0.1846157454e1, 0.4183473197, 0.2370620711e-1,
                0.9771730573e-1, 0.3017891294e-1, 0.2273353212e-1,
                0.1357254086e-1, -0.40526989430e-1, 0.54546285150e-3,
                0.51131822770e-3, 0.29534668830e-6, -0.86876450720e-4],
        "d1": [1, 1, 1, 2, 2, 2, 3, 3, 3, 6, 7, 7, 8],
        "t1": [0, 1.5, 2.5, -0.5, 1.5, 2, 0, 1, 2.5, 0, 2, 5, 2],

        "nr2": [-0.2127082589, 0.8735941958e-1, 0.127550919, -0.9067701064e-1,
                -0.3540084206e-1, -0.3623278059e-1, 0.132769929e-1,
                -0.3254111865e-3, -0.8313582932e-2, 0.2124570559e-2,
                -0.8325206232e-3, -0.2626173276e-4, 0.2599581482e-2,
                0.9984649663e-2, 0.2199923153e-2, -0.2591350486e-1,
                -0.12596308480, 0.14783556370, -0.10112510780e-1],
        "d2": [1, 1, 2, 2, 3, 3, 5, 6, 7, 8, 10, 2, 3, 3, 4, 4, 5, 5, 5],
        "t2": [5, 6, 3.5, 5.5, 3, 7, 6, 8.5, 4, 6.5, 5.5, 22, 11, 18, 11, 23,
               17, 18, 23],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*20}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for oxygen of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 1000.0, "Pmax": 82000.0, "rhomax": 43.348, 
#        "Pmin": 73.5, "rhomin": 29.2, 

        "nr1": [0.88878286369701, -0.24879433312148e1, 0.59750190775886,
                0.96501817061881e-2, 0.71970428712770e-1, 0.22337443000195e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.18558686391474, -0.38129368035760e-1, -0.15352245383006,
                -0.26726814910919e-1, -0.25675298677127e-1, 0.95714302123668e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.623, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*20,

        "nr3": [],
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for oxygen of Younglove (1982).",
        "__doc__": u"""Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
        "R": 8.31411,
        "cp": CP2,

        "Tmin": Tt, "Tmax": 400.0, "Pmax": 121000.0, "rhomax": 40.820, 
        "Pmin": 0.148, "rhomin": 40.820, 

        "b": [None, -0.4365859650e-3, 0.2005820677, -0.4197909916e1,
              0.1878215317e3, -0.1287473398e5, 0.1556745888e-4, 0.1343639359e-2,
              -0.2228415518e1, 0.4767792275e4, 0.4790846641e-6, 0.2462611107e-2,
              -0.1921891680, -0.6978320847e-5, -0.6214145909e-3, -0.1860852567,
              0.2609791417e-4, -0.2447611408e-6, 0.1457743352e-3, -0.1726492873e-5,
              -0.238489252e4, -0.2301807796e6, -0.2790303526e2, 0.9400577575e5,
              -0.4169449637e-1, 0.2008497853e1, -0.125607652e-3, -0.6406362964,
              -0.2475580168e-7, 0.1346309703e-4, -0.116150247e-9,
              -0.1034699798e-7, 0.2365936964e-6]}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for oxygen of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.", 
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"}, 
        "__test__": """
            >>> st=O2(T=700, rho=200, eq=3)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            1.0308 41.051 1.0969
            >>> st2=O2(T=750, rho=100, eq=3)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            50.89 0.26457
            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 43.348, 
        "Pmin": 0.14603, "rhomin": 40.885, 

        "nr1": [0.88878286, -0.24879433e1, 0.59750191, 0.96501817e-2,
                0.71970429e-1, 0.22337443e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],
 
        "nr2": [0.18558686, -0.38129368e-1, -0.15352245, -0.26726815e-1,
                -0.25675299e-1, 0.95714302e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    eq = helmholtz1, MBWR, GERG, helmholtz3
    _PR = -0.003157

    _surface = {"sigma": [0.038612652], "exp": [1.228]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [3.9578, 0.0065], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [0.575, 1.028, -8.96, -5.15],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.5, 2.5]}
    _melting = {"eq": 2, "Tref": Tt, "Pref": 0.14633,
                "Tmin": Tt, "Tmax": 300.0,
                "a1": [], "exp1": [],
                "a2": [-0.32463539e2, 0.14278011e3, -0.14702341e3, 0.520012e2],
                "exp2": [0.0625, 0.125, 0.1875, 0.25],
                "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 0.14633,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-20.714], "exp2": [1.06],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.60595e1, 0.13050e1, -0.54178, -0.19410e1, 0.35514],
        "exp": [1., 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.16622e1, 0.76846, -0.10041, 0.20480, 0.11551e-1],
        "exp": [0.345, 0.74, 1.2, 2.6, 7.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.22695e1, -0.46578e1, -0.99480e1, -0.22845e2, -0.45190e2,
               -0.25101e2],
        "exp": [0.3785, 1.07, 2.7, 5.5, 10., 20.]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Lemmon (2004)",
               "__doi__": {"autor": "Lemmon, E.W. and Jacobsen, R.T.",
                            "title": "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air", 
                            "ref": "Int. J. Thermophys., 25:21-69, 2004.",
                            "doi": "10.1023/B:IJOT.0000022327.04529.f3"}, 
               "__test__": """
                    >>> st=O2(T=100, rhom=0)
                    >>> print "%0.5f" % st.mu.muPas
                    7.70243
                    >>> st=O2(T=300, rhom=0)
                    >>> print "%0.4f" % st.mu.muPas
                    20.6307
                    >>> st=O2(T=100, rhom=35)
                    >>> print "%0.3f" % st.mu.muPas
                    172.136
                    >>> st=O2(T=200, rhom=10)
                    >>> print "%0.4f" % st.mu.muPas
                    22.4445
                    >>> st=O2(T=300, rhom=5)
                    >>> print "%0.4f" % st.mu.muPas
                    23.7577
                    >>> st=O2(T=154.6, rhom=13.6)
                    >>> print "%0.4f" % st.mu.muPas
                    24.7898
                    """, # Table V, Pag 28
                    
              "Tref": 1., "etaref": 1,
              "ek": 118.5, "sigma": 0.3428,
              "n_chapman": 0.151011418/M**0.5,

              "Tref_res": 154.581, "rhoref_res": 13.63*M, "etaref_res": 1,
              "n_poly": [17.67, 0.4042, 0.0001077, 0.3510, -13.67],
              "t_poly": [0.05, 0, 2.1, 0, 0.5],
              "d_poly": [1, 5, 12, 8, 1],
              "g_poly": [0, 0, 0, 1, 1],
              "c_poly": [0, 0, 0, 1, 2]}

    visco1 = {"eq": 2, "omega": 2,
              "collision": [-67.2093902106092, 277.148660965491, -399.192753863192,
                            166.828729537446, 143.163477478684, -191.767060368781,
                            98.4332230147836, -22.9410694301649, 2.12402264924749],
              "__name__": "Younglove (1982)",
              "__doc__": """Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
              "ek": 113., "sigma": 0.3437,
              "n_chapman": 0.15099557923496,
              "t_chapman": 0.0,
              "F": [1.39279625307e-2, -6.51536010579e-3, 1.4, 100],
              "E": [-14.45497211, 243.40689667, 12.9006761056004,
                    -1949.07966423848, -5.62078436742e-2,
                    21.3075467849104, 48.9965711691056],
              "rhoc": 13.5942597847419}

    visco2 = {"eq": 1, "omega": 1,
              "collision": [0.46649, -0.57015, 0.19164, -0.03708, 0.00241],
              "__name__": "Laesecke (1990)",
              "__doc__": """Laesecke, A., Krauss, R., Stephan, K., and Wagner, W., "Transport Properties of Fluid Oxygen," J. Phys. Chem. Ref. Data, 19(5):1089-1122, 1990.""",
              "ek": 116.2, "sigma": 0.34318867,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.151006/M**0.5,

              "Tref_res": 1, "rhoref_res": 13.63*M, "etaref_res": 18.8928,
              "n_poly": [-1.7993647, -0.397230772, 0.312536267, -0.0615559341],
              "t_poly": [0, 0, 0, 0],
              "d_poly": [0, 1, 2, 3],
              "g_poly": [0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0],
              "n_num": [-5.60288207],
              "t_num": [0],
              "d_num": [0],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1.0, -3.1138112],
              "t_den": [0, 0],
              "d_den": [1, 0],
              "g_den": [0, 0],
              "c_den": [0, 0]}

    _viscosity = visco0, visco1, visco2

    thermo0 = {"eq": 1,
               "__name__": "Lemmon (2004)",
               "__doi__": {"autor": "Lemmon, E.W. and Jacobsen, R.T.",
                            "title": "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air", 
                            "ref": "Int. J. Thermophys., 25:21-69, 2004.",
                            "doi": "10.1023/B:IJOT.0000022327.04529.f3"}, 
               "__test__": """
                    >>> st=O2(T=100, rhom=0)
                    >>> print "%0.5f" % st.k.mWmK
                    8.94334
                    >>> st=O2(T=300, rhom=0)
                    >>> print "%0.4f" % st.k.mWmK
                    26.4403
                    >>> st=O2(T=100, rhom=35)
                    >>> print "%0.3f" % st.k.mWmK
                    146.044
                    >>> st=O2(T=200, rhom=10)
                    >>> print "%0.4f" % st.k.mWmK
                    34.6124
                    >>> st=O2(T=300, rhom=5)
                    >>> print "%0.4f" % st.k.mWmK
                    32.5491
                    >>> st=O2(T=154.6, rhom=13.6)
                    >>> print "%0.4f" % st.k.mWmK
                    377.476
                    """, # Table V, Pag 28

               "Tref": 154.581, "kref": 1e-3,
               "no": [1.036, 6.283, -4.262],
               "co": [-97, 0.9, 0.6],

               "Trefb": 154.581, "rhorefb": 13.63, "krefb": 1e-3,
               "nb": [15.31, 8.898, -0.7336, 6.728, -4.374, -0.4747],
               "tb": [0, 0, 0.3, 4.3, 0.5, 1.8],
               "db": [1, 3, 4, 5, 7, 10],
               "cb": [0, 0, 0, 2, 2, 2],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.24e-9, "gam0": 0.055, "qd": 0.51e-9, "Tcref": 309.162}

    thermo1 = {"eq": 3,
               "__name__": "Younglove (1982)",
               "__doc__": """Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",

               "ek": 113, "sigma": 0.3437,
               "Nchapman": 0.15099557923496,
               "tchapman": 0,
               "b": [-1.41202117453516, 8.06267523869911, -19.44147946395,
                     25.78193316324, -20.5167203343277, 10.0087040966906,
                     -2.90450673487991, 0.459605807669332, -3.01906029521e-2],
               "F": [0.00097916328, 0.00089116658, 1.12, 100],
               "E": [-21.520741137, 473.50508788, 11.9072051301147,
                     -2122.44247203833, 0, 0, 0],
               "rhoc": 31.251171918947,
               "ff": 2.21064,
               "rm": 0.000000038896}

    thermo2 = {"eq": 1,
               "__name__": "Laesecke (1990)",
               "__doc__": """Laesecke, A., Krauss, R., Stephan, K., and Wagner, W., "Transport Properties of Fluid Oxygen," J. Phys. Chem. Ref. Data, 19(5):1089-1122, 1990.""",

               "Tref": 1, "kref": 1e-3,
               "no": [0.5825413, 0.0321266],
               "co": [-97, -98],

               "Trefb": 1, "rhorefb": 13.63, "krefb": 4.909e-3,
               "nb": [2.32825085, 4.23024231, -3.60798307, 2.01675631, -0.289731736],
               "tb": [0, 0, 0, 0, 0],
               "db": [1, 2, 3, 4, 5],
               "cb": [0, 0, 0, 0, 0],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 1.6e-10, "gam0": 0.08391, "qd": 0.4167e-9, "Tcref": 309.162}

    _thermal = thermo0, thermo1, thermo2
