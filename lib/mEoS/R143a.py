#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R143a(MEoS):
    """Multiparameter equation of state for R143a

    >>> r143a=R143a(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r143a.T, r143a.rho, r143a.u.kJkg, r143a.h.kJkg, r143a.s.kJkgK, r143a.cv.kJkgK, r143a.cp.kJkgK, r143a.w)
    300.0 3.4255 -6.63 22.57 0.1040 0.84195 0.94945 179.93
    """
    name = "1,1,1-trifluoroethane "
    CASNumber = "420-46-2"
    formula = "CF3CH3"
    synonym = "R143a"
    rhoc = unidades.Density(431.)
    Tc = unidades.Temperature(345.857)
    Pc = unidades.Pressure(3761.0, "kPa")
    M = 84.041  # g/mol
    Tt = unidades.Temperature(161.34)
    Tb = unidades.Temperature(225.909)
    f_acent = 0.2615
    momentoDipolar = unidades.DipoleMoment(2.34, "Debye")
    id = 243

    CP1 = {"ao": 0,
           "an": [1.0578], "pow": [0.33],
           "ao_exp": [4.4402, 3.7515], "exp": [1791, 823],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 0.10002060,
           "an": [-0.96337511e-3, 0.31822397e3, 0.46917620e-1],
           "pow": [1.5, -1.25, 1],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP3 = {"ao": 1.838736,
           "an": [3.01994e-2, -1.78455e-5, 4.42442e-9], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-143a of Lemmon and Jacobsen (2000).",
        "__doc__":  u"""Lemmon, E.W. and Jacobsen, R.T, "An International Standard Formulation for the Thermodynamic Properties of 1,1,1-Trifluoroethane (HFC-143a) for Temperatures from 161 to 450 K and Pressures to 50 MPa," J. Phys. Chem. Ref. Data, 29(4):521-552, 2000.""",
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 650.0, "Pmax": 100000.0, "rhomax": 15.85, 
        "Pmin": 1.0749, "rhomin": 15.832, 

        "nr1": [.77736443e1, -.870185e1, -.27779799, .1460922, .89581616e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.67, 0.833, 1.7, 1.82, 0.35],

        "nr2": [-0.20552116, 0.10653258, 0.23270816e-1, -0.13247542e-1,
                -0.42793870e-1, 0.36221685, -0.25671899, -0.92326113e-1,
                0.83774837e-1, 0.17128445e-1, -0.17256110e-1, 0.49080492e-2],
        "d2": [1, 3, 5, 7, 1, 2, 2, 3, 4, 2, 3, 5],
        "t2": [3.9, 0.95, 0, 1.19, 7.2, 5.9, 7.65, 7.5, 7.45, 15.5, 22, 19],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*12}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-143a of Outcalt and McLinden (1996)",
        "__doc__":  u"""Outcalt, S.L. and McLinden, M.O., "An equation of state for the thermodynamic properties of R143a (1,1,1-trifluoroethane)," Int. J. Thermophys., 18(6):1445-1463, 1997.""",
        "R": 8.314471,
        "cp": CP3,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 15.84, 
        "Pmin": 1.069, "rhomin": 15.8328, 

        "b": [None, -0.240561786316e-1, 0.262345913719e1, -0.650858041394e2,
              0.995952053681e4, -0.147536464961e7, 0.135498153308e-2,
              -0.281726617426e1, 0.134371062574e4, 0.850286316514e6,
              -0.180516636446e-3, 0.618889066246, -0.223083798271e3,
              -0.119095922349e-1, -0.173933336877e1, -0.420847601180e3,
              0.213502079796, -0.565708555185e-2, 0.185442296800e1,
              -0.520377059921e-1, -0.846735696108e6, -0.207964483848e8,
              -0.349977290513e5, 0.576427827667e9, -0.389131863941e3,
              0.103074054089e5, -0.447627052215e1, -0.106673161101e6,
              -0.219511369081e-1, 0.642186519493e1, -0.938317030843e-4,
              -0.478594713528e-1, -0.206555883874e1]}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-143a of Li et al. (1997).",
        "__doc__":  u"""Li, J., Tillner-Roth, R., Sato, H., and Watanabe, K., "An Equation of State for 1,1,1-Trifluoroethane (R-143a)," Int. J. Thermophys., 20(6):1639-1651, 1999.""",
        "R": 8.31451,
        "cp": CP2,
        
        "Tmin": Tt, "Tmax": 650.0, "Pmax": 50000.0, "rhomax": 15.84, 
        "Pmin": 1.0808, "rhomin": 15.819, 

        "nr1": [.1606645e-1, .4163515e1, -.5031058e1, -.1920208e-1, .1470093e-2],
        "d1": [5, 1, 1, 2, 4],
        "t1": [0, 0.5, 0.75, 2.5, 2.5],

        "nr2": [0.1775429, -0.7316069e-2, -0.9555916e-1, -0.5822518,
                -0.4211022e-3, -0.2059847e-1, 0.3711325e-1, 0.1799723e-3,
                -0.4145922e-1, 0.7682566e-4, -0.2089695e-2, 0.1958633e-2,
                -0.3198325e-5, -0.5376561e-2],
        "d2": [3, 8, 3, 1, 10, 1, 4, 8, 2, 12, 8, 2, 5, 3],
        "t2": [0.25, 0.25, 2, 3, 3, 8, 8, 8, 10, 8, 17, 20, 35, 27],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4],
        "gamma2": [1]*14}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-143a of Span and Wagner (2003)",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. III. Results for Polar Fluids," Int. J. Thermophys., 24(1):111-162, 2003.""",
        "R": 8.31451,
        "cp": CP2,
        
        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 15.82, 
        "Pmin": 1.072, "rhomin": 15.816, 

        "nr1": [.10306886e1, -.29497307e1, .6943523, .71552102e-1, .19155982e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.79764936e-1, 0.56859424, -0.90946566e-2, -0.24199452,
                -0.70610813e-1, -0.75041709e-1, -0.16411241e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq = helmholtz1, MBWR, helmholtz2, helmholtz3

    _surface = {"sigma": [0.0549], "exp": [1.26]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.73938e1, 0.19948e1, -0.18487e1, -0.41927e1, 0.14862e1],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.21135e1, 0.10200e2, -0.30836e2, 0.39909e2, -0.18557e2],
        "exp": [0.348, 1.6, 2.0, 2.4, 2.7]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.28673e1, -.63818e1, -.16314e2, -.45947e2, -.13956e1, -.24671e3],
        "exp": [0.384, 1.17, 3.0, 6.2, 7.0, 15.0]}
