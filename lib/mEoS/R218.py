#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R218(MEoS):
    """Multiparameter equation of R218

#    >>> r218=R218(T=300, P=0.1)
#    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r218.T, r218.rho, r218.u.kJkg, r218.h.kJkg, r218.s.kJkgK, r218.cv.kJkgK, r218.cp.kJkgK, r218.w)
#    300.0 7.6962 302.73 315.72 1.4786 0.75013 0.79926 116.42
    """
    name = "octafluoropropane"
    CASNumber = "76-19-7"
    formula = "CF3CF2CF3"
    synonym = "R218"
    rhoc = unidades.Density(627.985)
    Tc = unidades.Temperature(345.02)
    Pc = unidades.Pressure(2640.0, "kPa")
    M = 188.01933  # g/mol
    Tt = unidades.Temperature(125.45)
    Tb = unidades.Temperature(236.36)
    f_acent = 0.3172
    momentoDipolar = unidades.DipoleMoment(0.14, "Debye")
    id = 671

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-15.6587335175, 11.4531412796],
           "ao_exp": [7.2198, 7.2692, 11.599],
           "titao": [326/Tc, 595/Tc, 1489/Tc]}

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [7.2198, 7.2692, 11.599], "exp": [326, 595, 1489],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-218 of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=R218(T=347, rho=3*188.01933)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            347 3 2742.100 58080.724 251.735 181.131 2375.958 57.554
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 
        
        "Tmin": Tt, "Tmax": 440.0, "Pmax": 20000.0, "rhomax": 10.69, 
        "Pmin": 0.00202, "rhomin": 10.69, 

        "nr1": [1.3270, -3.8433, 0.922, 0.1136, 0.00036195],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [1.1001, 1.1896, -.025147, -.65923, -.027969, -.1833, -.02163],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1,

    _surface = {"sigma": [0.04307], "exp": [1.22]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.78419e1, 0.28989e1, -0.33458e1, -0.33196e1, 0.25363],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-0.61027, 0.57453e1, -0.56835e1, 0.32137e1, 0.55194],
        "exp": [0.223, 0.39, 0.56, 0.75, 5.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.42658e1, -.69496e1, -.18099e2, -.4921e2, -.55945e2, -.74492e2],
        "exp": [0.481, 1.53, 3.2, 6.3, 12.0, 15.0]}
