#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class Cyclohexane(MEoS):
    """Multiparameter equation of state for cyclohexane

    >>> ciclohexano=Cyclohexane(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (ciclohexano.T, ciclohexano.rho, ciclohexano.u.kJkg, ciclohexano.h.kJkg, ciclohexano.s.kJkgK, ciclohexano.cv.kJkgK, ciclohexano.cp.kJkgK, ciclohexano.w)
    500.0 3.56 377.04 405.10 0.89052 2.4600 2.5333 166.4
    """
    name = "cyclohexane"
    CASNumber = "110-82-7"
    formula = "cyclo(CH2)6"
    synonym = ""
    rhoc = unidades.Density(273.)
    Tc = unidades.Temperature(553.64)
    Pc = unidades.Pressure(4075.0, "kPa")
    M = 84.1608  # g/mol
    Tt = unidades.Temperature(279.47)
    Tb = unidades.Temperature(353.886)
    f_acent = 0.20926
    momentoDipolar = unidades.DipoleMoment(0.3, "Debye")
    id = 38
    _Tr = unidades.Temperature(526.231121)
    _rhor = unidades.Density(274.647526)
    _w = 0.221837522

    CP1 = {"ao": 9.3683272,
           "an": [-0.56214088e8, 0.15261554e-1, -0.36352468e-5],
           "pow": [-3, 1, 2],
           "ao_exp": [.23766589e2],
           "exp": [2000],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclohexane of Penoncello et al. (1995)",
        "__doc__":  u"""Penoncello, S.G., Goodwin, A.R.H., and Jacobsen, R.T, "A Thermodynamic Property Formulation for Cyclohexane," Int. J. Thermophys., 16(2):519-531, 1995.""",
        "R": 8.31434,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 700.0, "Pmax": 80000.0, "rhomax": 9.77, 
        "Pmin": 5.2538, "rhomin": 9.4045, 

        "nr1": [0.8425412659, -0.3138388327e1, 0.1679072631e1, -0.153819249,
                0.1984911143, -0.144532594, 0.3746346428e-3, 0.1861479616e-3,
                0.1745721652e-3],
        "d1": [1, 1, 1, 2, 3, 3, 7, 6, 6],
        "t1": [0, 1.5, 2.5, 1.5, 1, 2.5, 2, 0.5, 3],

        "nr2": [-0.6427428062, 0.2280757615, -0.1868116802e1, -0.1028243711e1,
                0.5821457418, -0.255891152, 0.1276844113e-1, -0.5158613166e-2,
                0.6334794755e-1, -0.6014686589e-1, 0.4439056828, -0.6264920642,
                0.2132589969e1, -0.3620300991e-2, 0.2534453992,
                0.1669144715e-1, 0.3985052291e-2],
        "d2": [1, 1, 2, 3, 3, 5, 8, 10, 3, 4, 1, 1, 2, 2, 4, 4, 8],
        "t2": [5, 6, 5.5, 3, 7, 6, 6.5, 5.5, 11, 11, 0.5, 1, 4, 4, 1.5, 2, 0.5],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 3, 3, 2, 6, 2, 4, 2],
        "gamma2": [1]*17}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for cyclohexane of Span and Wagner (2003)",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. II. Results for Nonpolar Fluids," Int. J. Thermophys., 24(1):41-109, 2003.""",
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 9.77, 
        "Pmin": 5.2428, "rhomin":9.3999, 

        "nr1": [0.10232354e1, -0.29204964e1, 0.10736630e1, -0.19573985,
                0.12228111, 0.28943321e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.27231767, -0.4483332e-1, -0.38253334, -0.89835333e-1,
                -0.24874965e-1, 0.10836132e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for cyclohexane of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.31434,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [1.27436292, 1.15372124, -3.86726473, 8.84627298e-2,
                2.76478090e-4, 7.26682313e-2],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [7.10849914e-2, 4.46376742e-1, 7.64476190e-1, -4.23520282e-2,
                -3.96468623e-1, -1.41250071e-2, -1.08371284e-1, -2.50082884e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, helmholtz2, helmholtz3

    _surface = {"sigma": [0.0653], "exp": [1.26]}
    _melting = {"eq": 1, "Tref": 1, "Pref": 700,
                "Tmin": Tt, "Tmax": 370.0, 
                "a1": [0.1329969885, -374.255624], "exp1": [1.41, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.70408e1, 0.19110e1, -0.18308e1, -0.16218e2, 0.20237e2],
        "exp": [1., 1.5, 2.13, 5.9, 7.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.17177, 0.25797e3, -0.73723e3, 0.11669e4, -0.13041e4, 0.62129e3],
        "exp": [0.093, 1.14, 1.34, 1.67, 2.0, 2.2]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.12016, -0.57983e1, -0.95451e2, 0.22937e3, -0.22100e3, 0.34023e2],
        "exp": [0.102, 0.63, 3.0, 3.6, 4.2, 4.5]}


if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=Cyclohexane(T=300., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
