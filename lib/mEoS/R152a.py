#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class R152a(MEoS):
    """Multiparameter equation of state for R152a

    >>> r152a=R152a(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.2f %0.3f %0.5f %0.4f %0.4f %0.2f" % (r152a.T, r152a.rho, r152a.u.kJkg, r152a.h.kJkg, r152a.s.kJkgK, r152a.cv.kJkgK, r152a.cp.kJkgK, r152a.w)
    300.0 3.4255 -6.63 22.57 0.1040 0.84195 0.94945 179.93
    """
    name = "1,1-difluoroethane"
    CASNumber = "75-37-6"
    formula = "CHF2CH3"
    synonym = "R152a"
    rhoc = unidades.Density(368.)
    Tc = unidades.Temperature(386.411)
    Pc = unidades.Pressure(4516.75, "kPa")
    M = 66.051  # g/mol
    Tt = unidades.Temperature(154.56)
    Tb = unidades.Temperature(249.127)
    f_acent = 0.27521
    momentoDipolar = unidades.DipoleMoment(2.262, "Debye")
    id = 245

    CP1 = {"ao": 27.89465,
           "an": [9.134686e-2, 2.079961e-4, -2.317613e-7], "pow": [1, 2, 3],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 0,
           "an": [1.4652739, 0.2627677e-4, -0.29988241e-10],
           "pow": [0.25, 2, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for R-152a of Outcalt and McLinden (1996)",
        "__doc__":  u"""Outcalt, S.L. and McLinden, M.O., "A modified Benedict-Webb-Rubin equation of state for the thermodynamic properties of R152a (1,1-difluoroethane)," J. Phys. Chem. Ref. Data, 25(2):605-636, 1996.""",
        "R": 8.314471,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 18.07, 
        "Pmin": 0.0641, "rhomin": 18.061, 

        "b": [None, -0.101623317192e-1, 0.215677129618e1, -0.648581254334e2,
              0.122535596303e5, -0.206805988259e7, -0.379836507323e-3,
              -0.441333232984, 0.158248874708e3, 0.564062216256e6,
              -0.124115350431e-3, 0.494972178825, -0.208058039834e3,
              -0.131403187106e-1, 0.212083848812, -0.151263785082e3,
              0.311108025395e-1, -0.115280979645e-2, 0.437040025765,
              -0.965596535032e-2, -0.242705525346e6, -0.518042519989e8,
              -0.119070545681e5, 0.459333195257e9, -0.719317286511e2,
              -0.840102861460e4, -0.102910957390e1, -0.325913880841e5,
              -0.412362182230e-2, 0.175102808144e1, -0.198636624640e-4,
              -0.421363036104e-2, -0.198696760653e1]}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz transform of MBWR EOS for R-152a of Outcalt and McLinden (1996).",
        "__doc__":  u"""Outcalt, S.L. and McLinden, M.O., "A modified Benedict-Webb-Rubin equation of state for the thermodynamic properties of R152a (1,1-difluoroethane)," J. Phys. Chem. Ref. Data, 25(2):605-636, 1996.""",
        "R": 8.314471,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 500.0, "Pmax": 60000.0, "rhomax": 18.07, 
        "Pmin": 0.0641, "rhomin": 18.061, 

        "nr1": [-0.354657949982e1, -0.364631280620, 0.333233335558e-1,
                -0.6809684351170, 0.735212646801e1, -0.112473063838e2,
                0.549916715657e1, -0.240186327322e1, -0.709036447042e-1,
                -0.213200886814, 0.197839736368, 0.182494769909e1,
                -0.860546479693e-1, 0.888137366540, -0.966127346370,
                -0.985223479324e-1, 0.183419368472e-1, -0.338550204252e-1,
                0.124921101016e-1, -0.221056706423e-2, 0.216879133161e-2,
                -0.233597690478e-3],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 6, 7, 7, 8],
        "t1": [3, 4, 5, 0, 0.5, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 1, 2, 3, 2, 2, 3, 3],

        "nr2": [0.354657949982e1, 0.364631280620, -0.333233335558e-1,
                0.276133830254e1, -0.691185711880e-1, -0.333233335558e-1,
                0.782761327717, -0.345592855940e-1, 0.137813531906,
                0.186173126153, -0.341119393297e-1, 0.459378439687e-1,
                0.216470012607e-1, -0.852798483242e-2, 0.620394038634e-2,
                0.185210290813e-2, 0.101674662734e-2, 0.124078807727e-2],
        "d2": [0, 0, 0, 2, 2, 2, 4, 4, 4, 6, 6, 6, 8, 8, 8, 10, 10, 10],
        "t2": [3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5, 3, 4, 5],
        "c2": [2]*18,
        "gamma2": [1]*18}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-152a of Tillner-Roth (1995).",
        "__doc__":  u"""Tillner-Roth, R., "A Fundamental Equation of State for 1,1-Difluoroethane (HFC-152a)," Int. J. Thermophys., 16(1):91-100, 1995.""",
        "R": 8.314471,
        "cp": CP2,
        
        "Tmin": Tt, "Tmax": 435.0, "Pmax": 30000.0, "rhomax": 18.03, 
        "Pmin": 0.065395176, "rhomin": 18.020671, 

        "nr1": [0.3552260, -0.1425660e1, -0.4631621e-1, 0.6903546e-1,
                0.1975710e-1, 0.7486977e-3, 0.4642204e-3],
        "d1": [1, 1, 1, 1.5, 3, 6, 6],
        "t1": [0, 1.5, 3, -0.5, -0.5, -0.5, 1.5],

        "nr2": [-0.2603396, -0.7624212e-1, 0.2233522, 0.1992515e-1, 0.3449040,
                -0.4963849, 0.1290719, 0.9760790e-3, 0.5066545e-2,
                -0.1402020e-1, 0.5169918e-2, 0.2679087e-3],
        "d2": [1, 1, 3, 4, 1, 1, 1, 8, 2, 3, 5, 6],
        "t2": [3, 4, 3, 2, 4, 5, 6, 5, 12.5, 25, 20, 25],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
        "gamma2": [1]*12}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for R-152a of Span and Wagner (2003)",
        "__doi__": {"autor": "Span, R., Wagner, W.",
                    "title": "Equations of State for Technical Applications. III. Results for Polar Fluids", 
                    "ref": "Int. J. Thermophys., 24(1):111-162, 2003.",
                    "doi": "10.1023/A:1022362231796"}, 
        "__test__": """
            >>> st=R152a(T=700, rho=200, eq=3)
            >>> print "%0.4f %0.3f %0.4f" % (st.cp0.kJkgK, st.P.MPa, st.cp.kJkgK)
            1.4632 21.594 2.1580
            >>> st2=R152a(T=750, rho=100, eq=3)
            >>> print "%0.2f %0.5f" % (st2.h.kJkg-st.h.kJkg, st2.s.kJkgK-st.s.kJkgK)
            270.60 0.60934
            """, # Table III, Pag 117
            
        "R": 8.31451,
        "cp": CP2,
        
        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 18.1, 
        "Pmin": 0.064093, "rhomin": 18.031, 

        "nr1": [.95702326, -.23707196e1, .18748463, .63800843e-1, .16625977e-3],
        "d1": [1, 1, 1, 3, 7],
        "t1": [0.25, 1.25, 1.5, 0.25, 0.875],

        "nr2": [0.82208165e-1, 0.57243518, 0.39476701e-2, -0.23848654,
                -0.80711618e-1, -0.73103558e-1, -0.15538724e-1],
        "d2": [1, 2, 5, 1, 1, 4, 2],
        "t2": [2.375, 2, 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*12}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-152a of Astina and Sato (2004)",
        "__doc__":  u"""Astina, I.M. and Sato, H., "A Rigorous Thermodynamic Property Model for Fluid - Phase 1,1-Difluoroethane (R-152a)," Int. J. Thermophys., 25(6):1713-1733, 2004.""",
        "R": 8.314472,
        "cp": CP1,
        
        "Tmin": Tt, "Tmax": 450.0, "Pmax": 60000.0, "rhomax": 18.04 , 
        "Pmin": 0.064, "rhomin": 18.04 , 

        "nr1": [1.753847317, -4.049760759, -2.277389257e-1, 7.087751950e-1,
                -5.528619502e-1, -3.025046686e-2, 1.396289974e-1, 1.121238954e-4],
        "d1": [1, 1, 1, 2, 2, 3, 3, 4],
        "t1": [0.5, 0.125, 2.875, 0.875, 1.875, 0.5, 1.875, 4],

        "nr2": [1.181005890, 1.535785579, 7.468363045e-1, -1.252266405e-1,
                -3.898223986e-2, -7.260588801e-2, -2.659302250e-3,
                4.210849329e-3, 2.015953966e-4],
        "d2": [1, 2, 3, 1, 2, 3, 3, 4, 5],
        "t2": [1.25, 2, 2.75, 6, 9, 6, 22, 20, 32],
        "c2": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*9}

    eq = MBWR, helmholtz1, helmholtz2, helmholtz3, helmholtz4

    _surface = {"sigma": [0.05906], "exp": [1.221]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.74821e1, 0.21105e1, -0.20761e1, -0.35539e1, 0.58004],
        "exp": [1.0, 1.5, 2.2, 4.8, 6.2]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.19914e2, -0.68624e2, 0.99821e2, -0.77984e2, 0.29913e2],
        "exp": [0.56, 0.76, 0.95, 1.2, 1.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-.33621e1, -.85985e1, -.2683e1, -.2414e2, -.43159e2, -.28045e2],
        "exp": [0.406, 1.42, 3.6, 3.9, 8.0, 9.0]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.4425728, -0.5138403, 0.1547566, -0.2821844e-1,
                            0.1578286e-2],
              "__name__": "Krauss (1996)",
              "__doc__": """Krauss, R., Weiss, V.C., Edison, T.A., Sengers, J.V., and Stephan, K., "Transport properties of 1,1-Difluoroethane (R152a)," Int. J. Thermophysics 17:731-757, 1996.""",
              "ek": 354.84, "sigma": 0.46115,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.2169614/M**0.5,

              "Tref_res": 1., "rhoref_res": 5.571537*M, "etaref_res": 51.12,
              "n_poly": [-.139987, -.737927e-1, .517924, -.308875, .108049],
              "t_poly": [0, 0, 0, 0, 0],
              "d_poly": [0, 1, 2, 3, 4],
              "g_poly": [0, 0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0, 0],
              "n_num": [-0.408387],
              "t_num": [0],
              "d_num": [0],
              "g_num": [0],
              "c_num": [0],
              "n_den": [-2.91733, 1.0],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [0, 0],
              "c_den": [0, 0]}

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Krauss (1996)",
               "__doc__": """Krauss, R., Weiss, V.C., Edison, T.A., Sengers, J.V., and Stephan, K., "Transport properties of 1,1-Difluoroethane (R152a)," Int. J. Thermophysics 17:731-757, 1996""",

               "Tref": 1., "kref": 1e-3,
               "no": [-1.49420e1, 9.73283-2],
               "co": [0, 1],

               "Trefb": 1., "rhorefb": 5.57145, "krefb": 1.115e-3,
               "nb": [9.1809, 1.18577e1, -5.44730, 1.71379],
               "tb": [0]*4,
               "db": [1, 2, 3, 4],
               "cb": [0]*4,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 1.894e-10, "gam0": 0.0487, "qd": 4.37e-10, "Tcref": 579.617}

    _thermal = thermo0,
