#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class Benzene(MEoS):
    """Ecuación de estado de multiparametros para el hexafluoruro de benceno

    >>> benceno=Benzene(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (benceno.T, benceno.rho, benceno.u.kJkg, benceno.h.kJkg, benceno.s.kJkgK, benceno.cv.kJkgK, benceno.cp.kJkgK, benceno.w)
    300.0 860.51 -155.13 -155.01 -0.45421 1.2679 1.7070 1295.0
    """
    name="benzene"
    CASNumber="71-43-2"
    formula="C6H6"
    synonym=""
    rhoc=unidades.Density(304.79239968)
    Tc=unidades.Temperature(562.02)
    Pc=unidades.Pressure(4906.3 , "kPa")
    M=78.11184      #g/mol
    Tt=unidades.Temperature(278.674)
    Tb=unidades.Temperature(353.22)
    f_acent=0.211
    momentoDipolar=unidades.DipoleMoment(0.0, "Debye")
    id=40

    CP1={  "ao": 3.94645,
                "an": [], "pow": [],
                "ao_exp": [7.36374, 18.6490, 4.01834], 
                "exp": [4116.0, 1511.0, 630.0],
                "ao_hyp": [],"hyp": []}
 
    CP2={  "ao": -0.478176/8.3143*78.108,
                "an": [0.618649e-2/8.3143*78.108, -0.380363e-5/8.3143*78.108, 0.699648e-9/8.3143*78.108, 0.42661e-13/8.3143*78.108],
                "pow": [1, 2, 3, 4],
                "ao_exp": [], "exp": [],
                "ao_hyp": [],"hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Thol et al. (2010).",
        "__doc__":  u"""Equation of state for benzene for temperatures from the melting line up to 750 K and pressures up to 500 MPa, High Temperatures -- High Pressures;2012, Vol. 41 Issue 2, p81""",
        "R": 8.314472,
        "cp": CP1,

        "nr1": [0.3513062e-1, 0.2229707e1, -0.3100459e1, -0.5763224, 0.2504179],
        "d1": [4, 1, 1, 2, 3],
        "t1": [1, 0.3, 0.744, 1.174, 0.68],

        "nr2": [-0.7049091, -0.1393433, 0.8319673, -0.3310741, -0.2793578e-1],
        "d2": [1, 3, 2, 2, 7],
        "t2": [2.5, 3.67, 1.26, 2.6, 0.95],
        "c2": [2, 2, 1, 2, 1],
        "gamma2": [1.]*5, 
        
        "nr3": [0.7087408, -0.3723906, -0.6267414e-1, -0.8629500],
        "d3": [1, 1, 3, 3],
        "t3": [1, 2.47, 3.35, 0.75],
        "alfa3": [1.032, 1.423, 1.071, 14.35],
        "beta3": [1.867, 1.766, 1.824, 297.5],
        "gamma3": [1.1180, 0.6392, 0.6536, 1.1640],
        "epsilon3": [0.7289, 0.9074, 0.7655, 0.8711],
        "nr4": []}
        
        
    helmholtz2={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Polt et al. (1992).",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP2,

        "nr1": [-0.918572178424, 0.155357491575e1, -0.356149241161, 0.817273664265, -0.331303917534e1, 0.335336626528e1, -0.256976312022e1, 0.427304812515, 0.406483484297, -0.329744378187, 0.208907540720, 0.777471199254e-1, -0.202621443063, -0.148580350700e-1, 0.503167715817e-1, 0.293012717053e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.918572178424, -0.155357491575e1, 0.356149241161, -0.447029533153e-1, 0.957712367542, -0.114688433057e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.95481]*6}

    eq=helmholtz1, helmholtz2

    _surface={"sigma": [0.0766304, -0.0157455, 0.0150819], "exp": [1.25, 2.25, 3.25]}
    _vapor_Pressure={ "eq": 5, "ao": [-0.71661e1, 0.21551e1, -0.20297e1, -0.40668e1, 0.38092], "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density={ "eq": 1, "ao": [0.18160e2, -0.56879e2, 0.87478e2, -0.64365e2, 0.18500e2], "exp": [0.534, 0.686, 0.84, 1.0, 1.2]}
    _vapor_Density={ "eq": 3, "ao": [-0.31147e1, -0.46689e1, -0.16161e2, -0.14650e3, 0.51887e3, -0.82772e3], "exp": [0.419, 1.12, 2.8, 7.3, 10., 12.]}

