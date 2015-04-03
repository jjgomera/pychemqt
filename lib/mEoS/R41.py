#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R41(MEoS):
    """Multiparameter equation of state for R41"""
    name = "fluoromethane"
    CASNumber = "593-53-3"
    formula = "CH3F"
    synonym = "R41"
    rhoc = unidades.Density(316.506)
    Tc = unidades.Temperature(317.28)
    Pc = unidades.Pressure(5897.0, "kPa")
    M = 34.03292  # g/mol
    Tt = unidades.Temperature(129.82)
    Tb = unidades.Temperature(194.84)
    f_acent = 0.2004
    momentoDipolar = unidades.DipoleMoment(1.851, "Debye")
    id = 225

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1, -1],
           "ao_pow": [-4.8676441160, 4.2527951258, -0.0268688568],
           "ao_exp": [5.6936, 2.9351],
           "titao": [1841/Tc, 4232/Tc]}

    CP1 = {"ao": 4.,
           "an": [0.00016937], "pow": [1.],
           "ao_exp": [5.6936, 2.9351], "exp": [1841, 4232],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 38.133739/8.314471,
           "an": [-7.88701e-2/8.314471, 3.29302e-4/8.314471, -2.37475e-7/8.314471],
           "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-41 of Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=R41(T=319, rho=9*34.03292)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            319 9 6129.100 13670.133 55.886 55.438 2796.224 189.549
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 
        
        "Tmin": Tt, "Tmax": 425.0, "Pmax": 70000.0, "rhomax": 29.66, 
        "Pmin": 0.345, "rhomin": 29.65, 

        "nr1": [1.6264, -2.8337, 0.0010932, 0.037136, 0.00018724],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.52, 1.12, 4, 0.03, 0.63],

        "nr2": [-.22189, .55021, .0461, -.056405, -.17005, -.032409, -.012276],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [3.4, 2.2, 1.5, 0.1, 4.8, 3.5, 15],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-41 of Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=R41(T=319, rho=9*34.03292, eq=1)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            319 9 6130.986 13699.382 55.982 53.388 2724.221 194.624
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 29.6, 
        "Pmin": 0.343, "rhomin": 29.6, 

        "nr1": [0.85316, -2.6366, 0.69129, 0.054681, 0.00012796],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.37093, 0.33920, -0.0017413, -0.095417, -0.078852, -0.030729,
                -0.011497],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-41 of Outcalt (1996).",
        "__doc__":  u"""Haynes, W.M., "Thermophysical properties of HCFC alternatives," National Institute of Standards and Technology, Boulder, Colorado, Final Report for ARTI MCLR Project Number 660-50800, 1996.""",
        "R": 8.314471,
        "cp": CP2,

        "Tmin": 175.0, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 27.1006, 
        "Pmin": 0.32, "rhomin": 29.66, 

        "b": [None, -0.326441485138e-1, 0.338620074694e1, -0.831696847103e2,
              0.139589938388e5, -0.156113972752e7, -0.165160386413e-2,
              0.118821153813e1, -0.137311604695e3, 0.176999573025e6,
              0.164945271187e-4, 0.595329911829e-1, -0.341969857376e2,
              -0.168552064750e-2, -0.758216269071e-2, -0.134800586220e2,
              0.311348265418e-2, -0.651499088798e-4, 0.184033192190e-1,
              -0.281459127843e-3, -0.186344956951e6, 0.110422095705e8,
              -0.147526754027e4, 0.261603025982e8, -0.744431617418e1,
              0.782355157170e3, -0.562784094508e-2, -0.843317187588e3,
              -0.600934897964e-4, 0.145050417148e-1, 0.222324172533e-7,
              -0.204419971811e-4, 0.245556593457e-3]}

    eq = helmholtz1, helmholtz2, MBWR

    _surface = {"sigma": [0.05049], "exp": [1.242]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.70970e1, 0.17409e1, -0.11668e1, -0.31830e1, 0.93827],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18181e2, -0.62193e2, 0.85171e2, -0.66958e2, 0.28790e2],
        "exp": [0.58, 0.8, 1.0, 1.3, 1.5]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.26966e2, .54303e2, -.36361e2, -.17816e2, -.48535e2, -.86727e2],
        "exp": [0.59, 0.72, 0.86, 3.2, 7.0, 15.0]}
