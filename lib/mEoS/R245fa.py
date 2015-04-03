#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R245fa(MEoS):
    """Multiparameter equation of state for R245fa"""
    name = "1,1,1,3,3-pentafluoropropane"
    CASNumber = "460-73-1"
    formula = "CF3CH2CHF2"
    synonym = "R245fa"
    rhoc = unidades.Density(516.0846)
    Tc = unidades.Temperature(427.16)
    Pc = unidades.Pressure(3651.0, "kPa")
    M = 134.04794  # g/mol
    Tt = unidades.Temperature(171.05)
    Tb = unidades.Temperature(288.29)
    f_acent = 0.3776
    momentoDipolar = unidades.DipoleMoment(1.549, "Debye")
    id = 671
    # id = 1817

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-13.4283638514, 9.87236538],
           "ao_exp": [5.5728, 10.385, 12.554],
           "titao": [222/Tc, 1010/Tc, 2450/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-245fa of Lemmon and Span (2006).",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=R245fa(T=429, rho=3*134.04794)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            429 3 3737.844 63909.823 235.875 172.283 1891.957 78.673
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 
        
        "Tmin": Tt, "Tmax": 440.0, "Pmax": 200000.0, "rhomax": 12.3, 
        "Pmin": 0.0125, "rhomin": 12.29, 

        "nr1": [1.2904, -3.2154, 0.50693, 0.093148, 0.00027638],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [.71458, .87252, -.015077, -.40645, -.11701, -.13062, -.022952],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1,

    _surface = {"sigma": [0.073586, 0.0103, -0.02663],
                "exp": [1.0983, 0.60033, 0.72765]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.79191e1, 0.24600e1, -0.33362e1, -0.42837e1, 0.41941],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.11301e2, -0.30007e2, 0.38458e2, -0.22241e2, 0.55179e1],
        "exp": [0.48, 0.68, 0.9, 1.2, 1.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.29132e1, -.67535e1, -.204e2, -.41608e2, -.11917e2, -.10261e3],
        "exp": [0.384, 1.17, 3.0, 6.0, 7.0, 12.0]}

    thermo0 = {"eq": 1,
               "__name__": "Marsh (2002)",
               "__doc__": """Marsh, K., Perkins, R., and Ramires, M.L.V., "Measurement and Correlation of the Thermal Conductivity of Propane from 86 to 600 K at Pressures to 70 MPa," J. Chem. Eng. Data (2002)47(4):932-940""",

               "Tref": 427.16, "kref": 1,
               "no": [0.300728e-1, -0.102742, 0.145703, -0.483106e-1],
               "co": [0, 1, 2, 3],

               "Trefb": 427.16, "rhorefb": 3.85, "krefb": 1,
               "nb": [-0.7391e-2, 0.887221e-2, -0.195198, 0.173498, 0.289485,
                      -.23028, -.126956, .892151e-1, .172567e-1, -.93749e-2],
               "tb": [0, 1]*5,
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.5e-9, "Tcref": 640.80}

    _thermal = thermo0,
