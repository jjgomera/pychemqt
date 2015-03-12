#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Acetone(MEoS):
    """Multiparameter equation  of state for Acetone"""
    name = "acetone"
    CASNumber = "67-64-1"
    formula = "CH3COCH3"
    synonym = ""
    rhoc = unidades.Density(272.971958)
    Tc = unidades.Temperature(508.1)
    Pc = unidades.Pressure(4700.0, "kPa")
    M = 58.07914  # g/mol
    Tt = unidades.Temperature(178.5)
    Tb = unidades.Temperature(329.22)
    f_acent = 0.3071
    momentoDipolar = unidades.DipoleMoment(2.88, "Debye")
    id = 140

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-9.4883659997, 7.1422719708],
           "ao_exp": [3.7072, 7.0675, 11.012],
           "titao": [310/Tc, 3480/Tc, 1576/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for acetone of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=Acetone(T=510, rho=4*58.07914)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            510 4 4807.955 51782.004 157.331 138.449 3766.619 125.351
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 550.0, "Pmax": 700000.0, "rhomax": 15.73, 
        "Pmin": 0.0023, "rhomin": 15.72, 

        "nr1":  [0.90041, -2.1267, -0.083409, 0.065683, 0.00016527],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [-0.039663, 0.72085, 0.0092318, -0.17217, -0.14961, -0.076124,
                -0.018166],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1,

    _surface = {"sigma": [0.07], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.76214e1, 0.17441e1, -0.20514e1, -0.26644e1, -0.69437],
        "exp": [1, 1.5, 2.57, 4.43, 15.]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.11118e2, -0.29507e2, 0.35255e2, -0.14712e2, 0.95560],
        "exp": [0.456, 0.626, 0.8, 1.0, 2.47]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.25200e1, -0.66065e1, -0.25751e2, 0.78120e1, -0.53778e2, -0.11684e3],
        "exp": [0.36, 1.05, 3.2, 4.0, 6.5, 14.0]}
