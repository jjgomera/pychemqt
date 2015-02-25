#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R161(MEoS):
    """Multiparameter equation of state for R161

    >>> r227ea=R161(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r227ea.T, r227ea.rho, r227ea.u.kJkg, r227ea.h.kJkg, r227ea.s.kJkgK, r227ea.cv.kJkgK, r227ea.cp.kJkgK, r227ea.w)
    300.0 7.0020 71.96 86.24 0.4638 0.75948 0.81568 122.13
    """
    name = "fluoroethane"
    CASNumber = "48.0595"
    formula = "C2H5F"
    synonym = "R227ea"
    rhoc = unidades.Density(301.81366)
    Tc = unidades.Temperature(375.3)
    Pc = unidades.Pressure(5091.0, "kPa")
    M = 48.0595  # g/mol
    Tt = unidades.Temperature(130.0)
    Tb = unidades.Temperature(235.6)
    f_acent = 0.217
    momentoDipolar = unidades.DipoleMoment(1.9397, "Debye")
    id = 247

    CP1 = {"ao": 3.985,
           "an": [], "pow": [],
           "ao_exp": [2.077, 9.265, 6.054], "exp": [420, 1548, 3882],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-161 of Lemmon (2005).",
        "__doc__":  u"""""",
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 400.0, "Pmax": 50000.0, "rhomax": 20.0, 
        "Pmin": 0.006, "rhomin": 19.95, 

        "nr1": [0.75688, -1.4110, -0.63922, 0.055685, 0.00028395],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [.73357, .67596, .011369, -.56406, -.094362, -.1678, .00034215],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2.0, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 2],
        "gamma2": [1]*7}

    eq = helmholtz1,

    _surface = {"sigma": [0.0589], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.75224e1, 0.29140e1, -0.30129e1, -0.44497e1, 0.24207e1],
        "exp": [1.0, 1.5, 2.3, 6.0, 7.0]}
    _liquid_Density = {
        "eq": 1,
        "ao": [-.22587e2, .13424e3, -.2671e3, .3389e3, -.31059e3, .13009e3],
        "exp": [0.56, 0.7, 0.9, 1.2, 1.5, 1.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.62548e1, 0.10499e2, -0.20353e2, -0.36709e2, -0.86781e2],
        "exp": [0.56, 1.3, 1.7, 5.0, 11.0]}

if __name__ == "__main__":
#    import doctest
#    doctest.testmod()

    cyc5=R161(T=300., P=0.1)
    print "%0.1f %0.2f %0.4f %0.6f %0.6f %0.6f %0.3f %0.5f %0.6f %0.9f" % (cyc5.T, cyc5.P.MPa, cyc5.rho, cyc5.cv.kJkgK, cyc5.cp.kJkgK, cyc5.cp0.kJkgK, cyc5.w, cyc5.joule.KMPa, cyc5.virialB, cyc5.virialC)
