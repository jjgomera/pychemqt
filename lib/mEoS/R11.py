#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class R11(MEoS):
    """Ecuación de estado de multiparametros para el R11

    >>> r11=R11(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r11.T, r11.rho, r11.u.kJkg, r11.h.kJkg, r11.s.kJkgK, r11.cv.kJkgK, r11.cp.kJkgK, r11.w)
    300.0 5.7039 13.47 31.00 0.2016 0.53664 0.61009 138.57
    """
    name="trichlorofluoromethane"
    CASNumber="75-69-4"
    formula="CCl3F"
    synonym="R11"
    rhoc=unidades.Density(554.)
    Tc=unidades.Temperature(471.11)
    Pc=unidades.Pressure(4407.638, "kPa")
    M=137.368      #g/mol
    Tt=unidades.Temperature(162.68)
    Tb=unidades.Temperature(296.858)
    f_acent=0.18875
    momentoDipolar=unidades.DipoleMoment(0.450, "Debye")
    id=217

    CP1={  "ao": 4.00564923,
                "an": [2.228875e-4], "pow": [1],
                "ao_exp": [1, 2, 1, 2, 1, 2],
                "exp": [1561.076, 1218.647, 770.035, 572.634, 502.854, 346.746],
                "ao_hyp": [],"hyp": []}

    CP2={  "ao": 4.0000024,
                "an": [], "pow": [],
                "ao_exp": [0.32960961e1, 0.28401126e1, 0.40350474, 0.30739271e1],
                "exp": [381.63168, 1368.22648, 3435.66931, 689.55053],
                "ao_hyp": [],"hyp": []}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-11 of Jacobsen et al. (1992)",
        "__doc__":  u"""Jacobsen, R.T, Penoncello, S.G., and Lemmon, E.W., "A fundamental equation for trichlorofluoromethane (R-11)," Fluid Phase Equilibria, 80:45-56, 1992.""",
        "R": 8.314510,
        "cp": CP1,

        "nr1": [0.125993633881e1, -0.260818574641e1, 0.982122542463e-2, -0.106085385839e1, 0.122820363510e1, 0.118000776439, -0.698956926463e-3, -0.355428373358e-1, 0.197169579643e-2, -0.848363012252e-2, 0.417997567653e-2, -0.242772533848e-3, 0.313371368974e-2, 0.396182646586e-5],
        "d1": [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 6, 8],
        "t1": [0.5, 1.5, 5, 1, 1.5, 0, 5, 2, 3, 1, 2, 4, 1, 4],

        "nr2": [0.339736319502, -0.203010634531, -0.106017859900, 0.451564882590, -0.339265767612, 0.114338523359, 0.319537833995e-1, 0.367908259780e-1, -0.961768948364e-5, 0.246717966418e-2, -0.167030256045e-2, 0.240710110806e-2, 0.156214678738e-2, -0.323352596704e-2],
        "d2": [1, 1, 2, 2, 2, 3, 4, 6, 10, 3, 5, 8, 9, 9],
        "t2": [5, 6, 3.5, 5.5, 7.5, 3, 2.5, 5, 1.5, 11, 9, 13, 5, 9],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 6, 6, 6, 6],
        "gamma2": [1]*14}

    helmholtz2={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-11 of Marx et al. (1992)",
        "__doc__":  u"""Marx, V., Pruss, A., and Wagner, W., "New Equation of State for R 12, R 22, R 11 and R 113," Fortschr.-Ber. VDI, Dusseldorf: VDI-Verlag, 19(57), 1992.""",
        "R": 8.31451,
        "cp": CP2,

        "nr1": [-0.219644325e1, 0.8562148696, 0.185864982e-1, 0.2807246052, -0.8526398864e-1, 0.1090334698e-1],
        "d1": [1, 1, 1, 2, 3, 5],
        "t1": [1.5, 2, 3, 0, 1.5, 1],

        "nr2": [0.4138515982, -0.3125498519, 0.1545749737, 0.1752299625, 0.2295443969e-1, -0.2094422944e-2, -0.1267942875e-8, 0.797272861e-2, -0.1520330549, 0.6448637628e-1, 0.2046144277e-3, -0.4100829613e-4, -0.123188575e-1, 0.6681486552e-2, -0.6742271171e-7],
        "d2": [1, 1, 2, 3, 5, 7, 14, 1, 2, 3, 11, 11, 4, 4, 10],
        "t2": [-0.5, 3.5, -0.5, 1, -0.5, 2, 4, 8, 8, 8, 4, 6, 18, 21, 33],
        "c2": [1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 4],
        "gamma2": [1]*15}

    helmholtz3={
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-11 of Span and Wagner (2003)",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. III. Results for Polar Fluids," Int. J. Thermophys., 24(1):111-162, 2003.""",
        "R": 8.31451,
        "cp": CP2,

        "nr1": [0.10656383e1, -0.32495206e1, 0.87823894, 0.87611569e-1, 0.29950049e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.42896949, 0.70828452, -0.17391823e-1, -0.37626521, 0.11605284e-1, -0.89550567e-1, -0.30063991e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*7}

    eq=helmholtz1, helmholtz2, helmholtz3

    _surface={"sigma": [0.062], "exp": [1.25]}
    _vapor_Pressure={ "eq": 5, "ao": [-0.70742e1, 0.38118e1, -0.32850e1, -0.76340e1, 0.50598e1], "exp": [1.0, 1.5, 1.73, 5.2, 6.0]}
    _liquid_Density={ "eq": 1, "ao": [0.20368e1, 0.12850e2, -0.22521e2, 0.11340e2, -0.94375], "exp": [0.357, 1.5, 1.7, 2.0, 3.0]}
    _vapor_Density={ "eq": 3, "ao": [-0.30296e1, -0.60723e1, -0.15890e2, -0.63024e2, 0.87167e2, -0.15715e3], "exp": [0.417, 1.25, 3.1, 6.8, 10.0, 12.0]}
