#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R141b(MEoS):
    """Multiparameter equation of state for R141b"""
    name = "1,1-dichloro-1-fluoroethane"
    CASNumber = "1717-00-6"
    formula = "CCl2FCH3"
    synonym = "R141b"
    rhoc = unidades.Density(458.56)
    Tc = unidades.Temperature(477.5)
    Pc = unidades.Pressure(4212.0, "kPa")
    M = 116.94962  # g/mol
    Tt = unidades.Temperature(169.68)
    Tb = unidades.Temperature(305.20)
    f_acent = 0.2195
    momentoDipolar = unidades.DipoleMoment(2.014, "Debye")
    id = 236
    # id = 1633

    Fi1 = {"ao_log": [1, 3.],
           "pow": [0, 1],
           "ao_pow": [-15.5074814985, 9.1871858933],
           "ao_exp": [6.8978, 7.8157, 3.2039],
           "titao": [502/Tc, 1571/Tc, 4603/Tc]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-141b of Lemmon and Span (2006)",
        "__doi__": {"autor": "Lemmon, E.W., Span, R.",
                    "title": "Short Fundamental Equations of State for 20 Industrial Fluids", 
                    "ref": "J. Chem. Eng. Data, 2006, 51 (3), pp 785–850",
                    "doi":  "10.1021/je050186n"}, 
        "__test__": """
            >>> st=R141b(T=479, rho=3*116.94962)
            >>> print "%0.0f %0.0f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f" % (st.T, st.rhoM, st.P.kPa, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
            479 3 4267.596 61757.274 214.894 126.963 1485.482 93.674
            """, # Table 10, Pag 842
            
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 400000.0, "rhomax": 12.56, 
        "Pmin": 0.0065, "rhomin": 12.56, 

        "nr1": [1.1469, -3.6799, 1.3469, 0.083329, 0.00025137],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.32720, 0.46946, -0.029829, -0.31621, -0.026219, -0.078043,
                -0.020498],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1,

    _surface = {"sigma": [7.3958e-5, 0.059941], "exp": [0.066331, 1.2214]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73784e1, 0.52955e1, -0.46639e1, -0.31122e1, -0.18972e1],
        "exp": [1.0, 1.5, 1.7, 4.2, 9.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.10443e2, -0.24726e2, 0.27718e2, -0.11220e2, 0.75848],
        "exp": [0.49, 0.68, 0.88, 1.1, 2.9]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31177e1, -0.68872e1, -0.18566e2, -0.40311e2, -0.95472e1,
               -0.12482e3],
        "exp": [0.398, 1.33, 3.3, 6.7, 7.0, 14.0]}
