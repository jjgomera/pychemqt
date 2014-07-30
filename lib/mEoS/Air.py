#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class Air(MEoS):
    """Ecuación de estado de multiparametros para el Aire

    >>> aire=Air(T=300, P=0.1)
    >>> print "%0.1f %0.4f %0.3f %0.3f %0.5f %0.4f %0.2f" % (aire.T, aire.rho, aire.h.kJkg, aire.s.kJkgK, aire.cv.kJkgK, aire.cp.kJkgK, aire.w)
    300.0 1.1613 -0.200 -0.001 0.71811 1.0064 347.35
    """
    name="air"
    CASNumber="1"
    formula="N2+Ar+O2"
    synonym="R-729"
    rhoc=unidades.Density(342.60456)
    Tc=unidades.Temperature(132.5306)
    Pc=unidades.Pressure(3786.0, "kPa")
    M=28.96546       #g/mol
    Tt=unidades.Temperature(59.75)
    Tb=unidades.Temperature(78.903)
    f_acent=0.0335
    momentoDipolar=unidades.DipoleMoment(0.0, "Debye")
    id=475

    CP1={  "ao": 0.34908880e1,
                "an": [0.23955256e-5, 0.71721112e-8, -0.31154131e-12, 0.22380669],
                "pow": [1, 2, 3, -1.5],
                "ao_exp": [0.79130951, 0.21223677],
                "exp": [3364.011, 2242.45],
                "ao_hyp": [],"hyp": []}

    CP2={  "ao": 0.34941563e1,
                "an": [-0.65392681e3, 0.29618973e2, 0.22380669, -0.47007760, -0.68351536e-5, 0.15136141e-7, -0.20027652e-11],
#                "pow": [-3, -2, -1.5, -1, 1, 2, 3],
                "pow": [-3, -2, -1.5, -1.01, 1, 2, 3],  #los coeficientes 0 y -1 producen overflow, el 0 se coloca en a0 pero desconozco que hacer con el -1
                "ao_exp": [0.78724442, 0.21223677],
                "exp": [3353.4061, 2242.45],
                "ao_hyp": [],"hyp": []}


    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for air of Lemmon et al. (2000)",
        "__doc__":  u"""Lemmon, E.W., Jacobsen, R.T, Penoncello, S.G., and Friend, D.G., "Thermodynamic Properties of Air and Mixtures of Nitrogen, Argon, and Oxygen from 60 to 2000 K at Pressures to 2000 MPa," J. Phys. Chem. Ref. Data, 29(3):331-385, 2000.""",
        "R": 8.31451, "Tref": 132.6312, "rhoref": 10.4477*M,
        "cp": CP1,

        "nr1": [0.118160747229, 0.713116392079, -0.161824192067e1, 0.714140178971e-1, -0.865421396646e-1, 0.134211176704, 0.112626704218e-1, -0.420533228842e-1,
                0.349008431982e-1, 0.164957183186e-3],
        "d1": [1, 1, 1, 2, 3, 3, 4, 4, 4, 6],
        "t1": [0, 0.33, 1.01, 0, 0, 0.15, 0, 0.2, 0.35, 1.35],

        "nr2": [-0.101365037912, -0.173813690970, -0.472103183731e-1, -0.122523554253e-1, -0.146629609713, -0.316055879821e-1, 0.233594806142e-3,
                0.148287891978e-1, -0.938782884667e-2],
        "d2": [1, 3, 5, 6, 1, 3, 11, 1, 3],
        "t2": [1.6, 0.8, 0.95, 1.25, 3.6, 6, 3.25, 3.5, 15],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3, 3],
        "gamma2": [1]*9}

    helmholtz2={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for air of Jacobsen et al. (1992)",
        "__doc__":  u"""Jacobsen, R.T, Penoncello, S.G., Beyerlein, S.W., Clarke, W.P., and Lemmon, E.W.,"A Thermodynamic Property Formulation for Air," Fluid Phase Equilibria, 79:113-124, 1992.""",
        "R": 8.31451,  "Tref": 132.6312, "rhoref": 10.4477*M,
        "cp": CP2,

        "nr1": [0.206604930965, 0.367099749382, -0.943192015369, 0.382519513142e-2, -0.865385542309e-1, 0.323019987452, 0.608695449299e-2, 0.128352106296e-3, -0.400058181940e-5],
        "d1": [1, 1, 1, 1, 2, 2, 4, 6, 7],
        "t1": [0, 0.25, 1, 3.5, 0, 0.25, 0.5, 2, 3],

        "nr2": [-0.544697915817, -0.526471065792, -0.608529300347, -0.124174277875, -0.595578533411e-2, -0.157523548353, -0.346463251040e-2, 0.837023084176e-2, -0.316701981142e-1, -0.721856676857e-2, 0.276838040645e-3, 0.160877459321e-4, 0.409235806738e-1, 0.652776125216e-3, -0.952903961290e-2, -0.100337820004e-1, 0.701111041628e-2, -0.472754336912e-2, 0.399257638569e-2, 0.968453675994e-2, -0.106826283630e-1, -0.489679885832e-2],
        "d2": [1, 2, 3, 5, 6, 1, 1, 2, 2, 3, 11, 11, 1, 1, 2, 3, 7, 8, 2, 4, 5, 2],
        "t2": [1.5, 1, 1, 1, 2, 3, 8, 0.5, 5.5, 9, 3, 6, 3, 9, 2, 13, 11, 11, 8, 22, 23, 11],
        "c2": [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5],
        "gamma2": [1]*22}

    eq=helmholtz1, helmholtz2

    _surface={"sigma": [0.03046], "exp": [1.28]}
    _melting={"eq": 1, "Tref": Tb, "Pref": 5.265, "a1": [1, 0.354935e5, -0.354935e5], "exp1": [0, 0.178963e1, 0], "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Density={ "eq": 3, "ao": [-0.20466e1, -0.4752e1, -0.13259e2, -0.47652e2], "exp": [0.41, 1, 2.8, 6.5]}
    _vapor_Pressure={ "eq": 5, "ao": [-0.1567266, -0.5539635e1, 0.7567212, -0.3514322e1], "exp": [0.5, 1, 2.5, 4]}

    visco0={"eq": 1, "omega": 1,
                "__name__": "Lemmon (2004)",
                "__doc__": """Lemmon, E.W. and Jacobsen, R.T, "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air," Int. J. Thermophys., 25:21-69, 2004.""",
                "ek": 103.3, "sigma": 0.36,
                "Tref": 1, "rhoref": 1.*M,

                "Tref_res": 132.6312, "rhoref_res": 10.4477*M,
                "n_poly": [10.72, 1.122, 0.002019, -8.876, -0.02916],
                "t_poly": [-.2, -.05, -2.4, -.6, -3.6],
                "d_poly": [1, 4, 9, 1, 8],
                "g_poly": [0, 0, 0, 1, 1],
                "c_poly": [0, 0, 0, 1, 1]}

    visco1={"eq": 1, "omega": 1,
                "collision": [0.5136, -0.5218, 0.8852e-1, 0.3445e-2, -0.2289e-2],
                "__name__": "Lemmon (2000)",
                "__doc__": """Lemmon, E.W. and Jacobsen, R.T, preliminary equation, 2000.""",
                "ek": 78.6, "sigma": 0.3711,
                "Tref": 1., "rhoref": 1.*M,

                "Tref_res": 132.6312, "rhoref_res": 10.4477*M,
                "n_poly": [0.3699e1, 0.2304e1, 0.2376e1, 0.6244e-3, 0.1616e-1],
                "t_poly": [0, 0, -0.1, -1.7, 0],
                "d_poly": [1, 2, 3, 10, 9],
                "g_poly": [0, 0, 0, 0, 1],
                "c_poly": [0, 0, 0, 0, 1]}

    _viscosity=visco0, visco1

    thermo0={"eq": 1,
                "__name__": "Lemmon (2004)",
                "__doc__": """Lemmon, E.W. and Jacobsen, R.T, "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air," Int. J. Thermophys., 25:21-69, 2004.""",

                "Tref": 132.6312, "kref": 1e-3,
                "no": [1.308, 1.405, -1.036],
                "co": [-97, 1.1, 0.3],

                "Trefb": 132.6312, "rhorefb": 10.4477, "krefb": 1e-3,
                "nb": [8.743, 14.76, -16.62, 3.793, -6.142, -0.3778],
                "tb": [-0.1, 0, -0.5, -2.7, -0.3, -1.3],
                "db": [1, 2, 3, 7, 7, 11],
                "cb": [0, 0, 2, 2, 2, 2],

                "critical": 3,
                "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
                "Xio": 0.11e-9, "gam0": 0.55e-1, "qd": 0.31e-9, "Tcref": 265.262}

    thermo1={"eq": 1,
                "__name__": "Lemmon (2000)",
                "__doc__": """Lemmon, E.W. and Jacobsen, R.T, preliminary equation, 2000.""",

                "Tref": 132.6312, "kref": 1e-3,
                "no": [1.1, 0],
                "co": [0, -96],

                "Trefb": 132.6312, "rhorefb": 10.4477, "krefb": 1e-3,
                "nb": [0.9759e1, 0.2259e2, -0.7995e1, -0.5714e2, 0.1324e2, 0.1456e2, 0.2577e1],
                "tb": [0, 0, -4, -0.15, -10.5, -0.5, -3],
                "db": [1, 2, 1, 3, 2, 4, 6],
                "cb": [0, 0, 1, 1, 2, 2, 2],

                "critical": 3,
                "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
                "Xio": 0.165e-9, "gam0": 0.55e-1, "qd": 0.386e-9, "Tcref": 265.2624}

    _thermal=thermo0, thermo1


