#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class nC8(MEoS):
    """Multiparameter equation of state for n-octane

    >>> octano=nC8(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (octano.T, octano.rho, octano.u.kJkg, octano.h.kJkg, octano.s.kJkgK, octano.cv.kJkgK, octano.cp.kJkgK, octano.w)
    300.0 696.78 11.33 11.47 0.99413 1.8108 2.2332 1139.3
    """
    name = "octane"
    CASNumber = "111-65-9"
    formula = "CH3-(CH2)6-CH3"
    synonym = ""
    rhoc = unidades.Density(234.9)
    Tc = unidades.Temperature(569.32)
    Pc = unidades.Pressure(2497.0, "kPa")
    M = 114.2285  # g/mol
    Tt = unidades.Temperature(216.37)
    Tb = unidades.Temperature(398.77)
    f_acent = 0.393
    momentoDipolar = unidades.DipoleMoment(0.07, "Debye")
    id = 12
    _Tr = unidades.Temperature(565.427917)
    _rhor = unidades.Density(234.605116)
    _w = 0.402698435
    
    Fi2 = {"ao_log": [1, 3.0],
           "pow": [0, 1],
           "ao_pow": [15.864687161, -97.370667555],
           "ao_exp": [], "titao": [], 
           "ao_hyp": [15.6865, 33.8029, 48.1731, 0],
           "hyp": [0.27914354, 1.431644769, 2.973845992, 0]}

    CP1 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.3961814e6, 0.2245626e8, 0.1380875e9, 0],
           "hyp": [0.1589220e3, 0.8150640e3, 0.1693070e4, 0]}

    CP2 = {"ao": 4,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [15.6865, 33.8029, 48.1731, 0],
           "hyp": [0.27914354*Tc, 1.431644769*Tc, 2.973845992*Tc, 0]}

    CP3 = {"ao": 3.018753,
           "an": [0.7297005e-1, -0.14171168e-4, -0.1225317e-7,  0.12912645e-11],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": 34.0847/8.3159524*4.184,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [2.603664e8/8.3159524*4.184, 4.1241363e7/8.3159524*4.184, 0, 0],
           "hyp": [1.6115500e3, 7.6884700e2, 0, 0]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for octane of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of state for technical applications. II. Results for nonpolar fluids.", 
                    "ref": "Int. J. Thermophys. 24 (2003), 41 – 109.",
                    "doi": "10.1023/A:1022310214958"}, 
        "__test__": """
            >>> st=nC8(T=700, rho=200)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            3.1537 6.363 3.8007
            >>> st2=nC8(T=750, rho=100)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            211.79 0.31183
            """, # Table III, Pag 46

        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 6.69, 
        "Pmin": 0.001989, "rhomin": 6.6864, 

        "nr1": [0.10722545e1, -0.24632951e1, 0.65386674, -0.36324974,
                0.1271327, 0.30713573e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.52656857, 0.19362863e-1, -0.58939427, -0.14069964,
                -0.78966331e-2, 0.33036598e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for octane of Kunz and Wagner (2004).",
        "__doi__": {"autor": "Kunz, O., Wagner, W.",
                    "title": "The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures: An Expansion of GERG-2004", 
                    "ref": "J. Chem. Eng. Data, 2012, 57 (11), pp 3032–3091",
                    "doi":  "10.1021/je300655b"}, 
        "R": 8.314472,
        "cp": Fi2,
        "ref": "OTO", 

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 6.69, 
#        "Pmin": 73.476, "rhomin": 29.249, 

        "nr1": [0.10722544875633e1, -0.24632951172003e1, 0.65386674054928,
                -0.36324974085628, 0.12713269626764, 0.30713572777930e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.5265685698754, 0.19362862857653e-1, -0.58939426849155,
                -0.14069963991934, -0.78966330500036e-2, 0.33036597968109e-2],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for octane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP3,

        "Tmin": 258.0, "Tmax": 500.0, "Pmax": 200000.0, "rhomax": 6.6355607, 
        "Pmin": 0.15134, "rhomin": 6.3907, 

        "nr1": [0.266117347782e1, -0.343810366899e1, 0.700476763325,
                0.573101545749e1, -0.411975339382e1, -0.771251551395e1,
                0.526137115388e1, -0.716144047789, -0.584632875151e1,
                0.736422551908e1, -0.100540027381e1, 0.158387242200e1,
                -0.153643650819e1, -0.142010818863, 0.333126039209e-1,
                0.271948869925e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [-0.266117347782e1, 0.343810366899e1, -0.700476763325,
                0.443217980268e1, -0.123858312597e2, 0.803373487925e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.9995725]*6}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for octane of Starling (1973)",
        "__doc__":  u"""Starling, K.E., "Fluid Thermodynamic Properties for Light Petroleum Systems," Gulf Publishing Company, 1973.""",
        "R": 8.3159524,
        "cp": CP4,

        "Tmin": 255.372, "Tmax": 644.0, "Pmax": 55000.0, "rhomax": 6.36203, 
        "Pmin": 0.099571, "rhomin": 6.3620, 

        "nr1": [0.253526486527e1, 0.616872653050, -0.941731168114,
                -0.109609729872e1, 0.849362892312e-1, -0.363538456997e-3,
                0.849748115039e-1, -0.961236603829e-1, -0.132591135067,
                0.269748328453e-2, 0.372085674947e-2],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.253526486527e1, -0.447291258549],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.35285564]*2}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [1.57750154, 1.15745614, -3.54867092, 1.18030671e-1,
                3.02753897e-4, -2.63074957e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [2.55299486e-2, -1.26632996e-1, 4.48343319e-1, -9.46702997e-3,
                -4.43927529e-1, -1.68224827e-2, -1.15864640e-1, -1.32417591e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, GERG, helmholtz3, helmholtz4, helmholtz5

    _surface = {"sigma": [0.0565399, -0.0104928, 0.0084723],
                "exp": [1.25, 2.25, 3.25]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.10924],  "expt0": [-1.], "expd0": [1.],
                   "a1": [39.74, 0.04], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [348.01, 494.18, -76838, -65772],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.79713e1, 0.17915e1, -0.34540e1, -0.82509e1, 0.53357e1],
        "exp": [1.0, 1.5, 2.64, 5.5, 6.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.56814e1, 0.38908e2, -0.75923e2, 0.59548e2, -0.19651e2],
        "exp": [0.1, 0.75, 0.9, 1.1, 1.25]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.16556, -0.59337e1, -0.18915e2, -0.36484e3, 0.72686e3, -0.50392e3],
        "exp": [0.09, 0.59, 2.4, 7.0, 8.0, 9.0]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.335103, -0.467898],
              "__name__": "Huber (2004)",
              "__doc__": """Huber, M.L., Laesecke, A. and Xiang, H.W. , "Viscosity correlations for minor constituent fluids in natural gas: n-octane, n-nonane and n-decane" Fluid Phase Equilibria 224(2004)263-270.""",
              "ek": 452.09, "sigma": 0.63617,
              "Tref": 1, "rhoref": 1.*M,
              "n_chapman": 0.228258776/M**0.5,

              "n_virial": [-0.19572881e2, 0.21973999e3, -0.10153226e4,
                           0.24710125e4, -0.33751717e4, 0.24916597e4,
                           -0.78726086e3, 0.14085455e2, -0.34664158],
              "t_virial": [0.0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 452.09, "etaref_virial": 0.1550494,

              "Tref_res": 569.32, "rhoref_res": 2.0564*M, "etaref_res": 1000,
              "n_packed": [0.20651e1, 0.307843e1, -0.879088],
              "t_packed": [0, 0.5, 1],
              "n_poly": [-0.103924, 0.113327e-1, 0.992302e-1, -0.322455e-1,
                         -0.606122],
              "t_poly": [-1, -1, -2, -2, 0],
              "d_poly": [2, 3, 2, 3, 1],
              "g_poly": [0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0],
              "n_num": [0.606122],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    visco1 = {"eq": 4, "omega": 1,
              "__doc__": """S.E.Quiñones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
              "__name__": "Quiñones-Cisneros (2006)",
              "Tref": 569.32, "etaref": 1.0,
              "ek": 452.09, "sigma": 0.636170, "n_chapman": 0,
              "n_ideal": [16.7562, -53.1705, 46.9105],
              "t_ideal": [0, 0.25, 0.5],

              "a": [8.68736376035937e-5, 0.0, -2.69591205491896e-5],
              "b": [1.46266597799792e-4, 0.0, -5.44584119633888e-5],
              "c": [1.28673387100000e-4, -1.76442029000000e-5, 0.0],
              "A": [-2.40884095261648e-9, 5.20715310859732e-11, 0.0],
              "B": [0.0, 6.62141302562572e-9, 1.60012396822086e-9],
              "C": [-9.50545390021906e-7, 1.03767490732769e-6, 0.0],
              "D": [0.0, 0.0, 0.0]}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 1,
               "__name__": "Huber (2005)",
               "__doc__": """Huber, M.L. and Perkins, R.A., "Thermal conductivity correlations for minor constituent fluids in natural gas: n-octane, n-nonane and n-decane" Fluid Phase Equilibria 227 (2005) 47-55""",

               "Tref": 569.32, "kref": 1,
               "no": [0.772930e-2, -0.371138e-1, 0.977580e-1, -0.288707e-1],
               "co": [0, 1, 2, 3],

               "Trefb": 569.32, "rhorefb": 2.0564, "krefb": 1,
               "nb": [0.285553e-1, -0.926155e-2, -0.171398e-1, 0.0,
                      0.659971e-2, 0.153496e-2, 0.0, 0.0, 0.0, 0.0],
               "tb": [0, 1]*5,
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.68628e-9, "Tcref": 853.98}

    _thermal = thermo0,
