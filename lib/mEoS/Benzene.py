#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Benzene(MEoS):
    """Multiparameter equation of state for benzene

    >>> benceno=Benzene(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (benceno.T, benceno.rho, benceno.u.kJkg, benceno.h.kJkg, benceno.s.kJkgK, benceno.cv.kJkgK, benceno.cp.kJkgK, benceno.w)
    300.0 860.51 -155.13 -155.01 -0.45421 1.2679 1.7070 1295.0
    """
    name = "benzene"
    CASNumber = "71-43-2"
    formula = "C6H6"
    synonym = ""
    rhoc = unidades.Density(304.79239968)
    Tc = unidades.Temperature(562.02)
    Pc = unidades.Pressure(4906.3, "kPa")
    M = 78.11184  # g/mol
    Tt = unidades.Temperature(278.674)
    Tb = unidades.Temperature(353.22)
    f_acent = 0.211
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 40

    Fi1 = {"ao_log": [1, 2.94645],
           "pow": [0, 1],
           "ao_pow": [-0.6740687105, 2.5560188958],
           "ao_exp": [7.36374, 18.6490, 4.01834],
           "titao": [4116/Tc, 1511/Tc, 630/Tc]}

    CP2 = {"ao": -0.478176/8.3143*78.108,
           "an": [0.618649e-2/8.3143*78.108, -0.380363e-5/8.3143*78.108,
                  0.699648e-9/8.3143*78.108, 0.42661e-13/8.3143*78.108],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Thol et al. (2010).",
        "__doi__": {"autor": "Thol M., Lemmon E.W., Span R.",
                    "title": "Equation of state for benzene for temperatures from the melting line up to 750 K and pressures up to 500 MPa", 
                    "ref": "High Temperatures-High Pressures 01/2012; 41:81.",
                    "doi": ""}, 
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 750., "Pmax": 500000.0, "rhomax": 11.45, 
        "Pmin": 4.78, "rhomin": 11.45, 

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

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Polt et al. (1992).",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe", 
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""}, 
        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP", 

        "Tmin": 278.7, "Tmax": 635.0, "Pmax": 78000.0, "rhomax": 11.45, 
        "Pmin": 6.0329, "rhomin": 11.385, 

        "nr1": [-0.918572178424, 0.155357491575e1, -0.356149241161,
                0.817273664265, -0.331303917534e1, 0.335336626528e1,
                -0.256976312022e1, 0.427304812515, 0.406483484297,
                -0.329744378187, 0.208907540720, 0.777471199254e-1,
                -0.202621443063, -0.148580350700e-1, 0.503167715817e-1,
                0.293012717053e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.918572178424, -0.155357491575e1, 0.356149241161,
                -0.447029533153e-1, 0.957712367542, -0.114688433057e1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.95481]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for benzene of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.314472,
        "cp": Fi1,
        "ref": "NBP", 

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [1.76284970, 1.02610647, -3.74263321, 9.57682041e-2,
                2.59179321e-4, -1.03082188e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [1.07359246e-1, -1.12562310e-1, 3.18737987e-1, -3.07549016e-2,
                -3.25082386e-1, 2.28099159e-2, -7.07431076e-2, -1.96809158e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, helmholtz2, helmholtz3

    _surface = {"sigma": [0.07298, -0.0007802, -0.0001756],
                "exp": [1.232, 0.8635, 0.3065]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.71661e1, 0.21551e1, -0.20297e1, -0.40668e1, 0.38092],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.18160e2, -0.56879e2, 0.87478e2, -0.64365e2, 0.18500e2],
        "exp": [0.534, 0.686, 0.84, 1.0, 1.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.31147e1, -0.46689e1, -0.16161e2, -0.14650e3, 0.51887e3, -0.82772e3],
        "exp": [0.419, 1.12, 2.8, 7.3, 10., 12.]}

if __name__ == "__main__":
    for eq in (0, 2):
        st=Benzene(T=353.22, P=101325., eq=eq)
        print "%0.6g %0.5g %0.1f %0.3f %0.3f %0.3f %0.3f %0.2f" % (\
            st.T, st.rhoM, st.uM.kJkmol, st.hM.kJkmol, st.sM.kJkmolK, st.cvM.kJkmolK, st.cpM.kJkmolK, st.w)
