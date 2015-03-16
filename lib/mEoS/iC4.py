#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class iC4(MEoS):
    """Multiparameter equation of state for isobutane

    >>> butano=iC4(T=500, P=0.1)
    >>> print "%0.1f %0.4f %0.2f %0.2f %0.4f %0.4f %0.4f %0.2f" % (butano.T, butano.rho, butano.u.kJkg, butano.h.kJkg, butano.s.kJkgK, butano.cv.kJkgK, butano.cp.kJkgK, butano.w)
    500.0 1.4048 357.41 428.59 1.0807 2.4261 2.5726 274.09

    >>> butano=iC4(u=-300.4, s=-1.1807)
    >>> print "%0.1f %0.4f %0.2f %0.2f %0.4f %0.4f %0.4f %0.2f" % (butano.T, butano.rho, butano.u.kJkg, butano.h.kJkg, butano.s.kJkgK, butano.cv.kJkgK, butano.cp.kJkgK, butano.w)
    330.0 573.4842 -300.40 -239.17 -1.1807 1.8387 2.3936 1031.02
    """
    name = "isobutane"
    CASNumber = "75-28-5"
    formula = "CH(CH3)3"
    synonym = "R-600a"
    rhoc = unidades.Density(225.5)
    Tc = unidades.Temperature(407.81)
    Pc = unidades.Pressure(3629.0 , "kPa")
    M = 58.1222  # g/mol
    Tt = unidades.Temperature(113.73)
    Tb = unidades.Temperature(261.401)
    f_acent = 0.184
    momentoDipolar = unidades.DipoleMoment(0.132, "Debye")
    id = 5
    _Tr = unidades.Temperature(390.355535)
    _rhor = unidades.Density(228.302484)
    _w = 0.178714317

    CP1 = {"ao": 4.05956619,
           "an": [], "pow": [],
           "ao_exp": [4.94641014, 4.09475197, 15.6632824, 9.73918122],
           "exp": [387.94064, 973.80782, 1772.71103, 4228.52424],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 4.06714,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [8.97575, 5.25156, 25.1423, 16.1388],
           "hyp": [1.074673199*Tc, 0.485556021*Tc, 4.671261865*Tc, 2.19158348*Tc]}

    CP3 = {"ao": 4.059347,
           "an": [], "pow": [],
           "ao_exp": [4.940314, 4.090139, 15.68832, 9.739581],
           "exp": [387.75987, 972.01102, 1772.81924, 4235.81166],
           "ao_hyp": [], "hyp": []}

    CP4 = {"ao": -1.7231723278e1,
           "an": [1.7027919006e7, -4.7269724737e5, 4.7301406581e3,
                  5.8491344291e-2, 8.9440351886e-6, -1.8274599197e-8],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-1.9283021962e1], "exp": [3000],
           "ao_hyp": [], "hyp": []}

    CP5 = {"ao": 4.06714,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.1724067e7, 0.2059196e6, 0.9124395e8, 0.1289193e8],
           "hyp": [0.4382700e3, 0.1980180e3, 0.1905020e4, 0.8937650e3]}

    CP6 = {"ao": 0.397893/8.3143*58.124,
           "an": [0.412501e-2/8.3143*58.124, -0.196195e-6/8.3143*58.124,
                  0.380185e-8/8.3143*58.124, -0.523950e-11/8.3143*58.124],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Buecker and Wagner (2006)",
        "__doc__":  u"""Bücker, D., Wagner, W. Reference equations of state for the thermodynamic properties of fluid phase n-butane and isobutane. J. Phys. Chem. Ref. Data 35 (2006), 929 – 1020.""",
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 35000.0, "rhomax": 12.9, 
        "Pmin": 0.0000219, "rhomin": 12.74, 

        "nr1":  [0.20686820727966e1, -0.36400098615204e1, 0.51968754427244,
                 0.17745845870123, -0.12361807851599, 0.45145314010528e-1,
                 0.30476479965980e-1],
        "d1": [1, 1, 1, 2, 3, 4, 4],
        "t1": [0.50, 1.00, 1.50, 0.00, 0.50, 0.50, 0.75],
        "nr2": [0.75508387706302, -0.85885381015629, 0.36324009830684e-1,
                -0.19548799450550e-1, -0.44452392904960e-2, 0.46410763666460e-2,
                -0.71444097992825e-1, -0.80765060030713e-1, 0.15560460945053,
                0.20318752160332e-2, -0.10624883571689, 0.39807690546305e-1,
                0.16371431292386e-1, 0.53212200682628e-3, -0.78681561156387e-2,
                -0.30981191888963e-2],
        "d2": [1, 1, 2, 7, 8, 8, 1, 2, 3, 3, 4, 5, 5, 10, 2, 6],
        "t2": [2.00, 2.50, 2.50, 1.50, 1.00, 1.50, 4.00, 7.00, 3.00, 7.00,
               3.00, 1.00, 6.00, 0.00, 6.00, 13.00],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3],
        "gamma2": [1]*16,

        "nr3": [-0.42276036810382e-1, -0.53001044558079e-2],
        "d3": [1, 2],
        "t3": [2., 0.],
        "alfa3": [10, 10],
        "beta3": [150, 200],
        "gamma3": [1.16, 1.13],
        "epsilon3": [0.85, 1.]}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for isobutane of Younglove and Ely (1987)",
        "__doc__":  u"""Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
        "R": 8.31434,
        "cp": CP4,

        "Tmin": 113.55, "Tmax": 600.0, "Pmax": 35000.0, "rhomax": 12.89, 
        "Pmin": 1.948e-5, "rhomin": 12.755, 

        "b": [None, 0.1307325972e-1, 0.3927802742, -0.3185427394e2,
              0.7608825192e4, -0.1753919859e7, -0.2090019755e-2, 0.8959557971e1,
              -0.6816710130e4, -0.1111271045e7, 0.3248737572e-3, -0.1046526456e1,
              0.6536598969e3, 0.3726503734e-1, 0.8553649395e1, 0.2109987236e4,
              -0.1401267363e1, 0.5213089327e-1, -0.1925026382e2, 0.7640067895,
              0.3425854273e7, -0.3373475924e9, 0.1180683444e6, 0.1529683738e10,
              0.3323837416e4, 0.6423169487e5, 0.3891706042e2, -0.1494755736e7,
              -0.1720240173e-1, 0.2894195375e3, 0.2005086329e-2, -0.4448393005,
              0.8028488415e2]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Kunz and Wagner (2004).",
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP2,

        "Tmin": Tt, "Tmax": 575.0, "Pmax": 35000.0, "rhomax": 12.9, 
#        "Pmin": 7.36, "rhomin": 38.2, 

        "nr1":  [0.10429331589100e1, -0.28184272548892e1, 0.86176232397850,
                 -0.10613619452487, 0.98615749302134e-1, 0.23948208682322e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.30330004856950, -0.41598156135099e-1, -0.29991937470058,
                -0.80369342764109e-1, -0.29761373251151e-1, 0.13059630303140e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Miyamoto and Watanabe (2001)",
        "__doc__":  u"""Miyamoto, H. and Watanabe, K. "A Thermodynamic Property Model for Fluid-Phase Isobutane," Int. J. Thermophys., 23(2):477-499, 2002.""",
        "R": 8.314472,
        "cp": CP3,

        "Tmin": 113.56, "Tmax": 573.0, "Pmax": 35000.0, "rhomax": 12.9, 
        "Pmin": 0.000021, "rhomin": 12.738, 

        "nr1":  [2.892737e-1, -1.342570, -7.976713e-3, 2.025793e-1,
                 -4.241612e-2, 2.617971e-3, 5.068955e-5, -1.144596e-6],
        "d1": [1, 1, 2, 2, 3, 5, 8, 8],
        "t1": [-0.25, 1.5, -0.75, 0, 1.25, 1.5, 0.5, 2.5],

        "nr2": [-1.930153, 1.982609, 2.076533e-3, -4.958752e-3, 1.377372e-3,
                -1.582662e-1, -4.961892e-2, 9.451030e-4, -3.037276e-2,
                -1.382675e-2, 8.876254e-5],
        "d2": [3, 3, 8, 5, 6, 1, 5, 7, 2, 3, 15],
        "t2": [1.5, 1.75, -0.25, 3, 3, 4, 2, -1, 2, 19, 5],
        "c2": [21, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3],
        "gamma2": [1]*11}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for isobutane of Span and Wagner (2003)",
        "__doc__":  u"""Span, R., Wagner, W. Equations of state for technical applications. II. Results for nonpolar fluids. Int. J. Thermophys. 24 (2003), 41 – 109.""",
        "R": 8.31451,
        "cp": CP5,

        "Tmin": 113.55, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 12.89, 
        "Pmin": 0.000020860, "rhomin": 12.784, 

        "nr1":  [0.10429332e1, -0.28184273e1, 0.86176232, -0.10613619,
                 0.986157490e-1, 0.23948209e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.30330005, -0.41598156e-1, -0.29991937, -0.80369343e-1,
                -0.29761373e-1, 0.1305963e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz5 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP6,

        "Tmin": 120.0, "Tmax": 498.0, "Pmax": 35000.0, "rhomax": 12.89, 
        "Pmin": 0.46491e-4, "rhomin": 12.649, 

        "nr1":  [-0.958589873652, 0.818846326211, -0.115814967179,
                 0.345513148715, -0.168751721524e1, 0.936693300209,
                 -0.106644545724e1, 0.980958295776e-1, 0.495941129005,
                 -0.261313404262, 0.485109471188, -0.177275820736,
                 -0.209415485311e-1, 0.788178884079e-1, -0.102751671767,
                 0.178645875838e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.958589873652, -0.818846326211, 0.115814967179,
                0.537585249054, -0.71942446879, 0.245830118086],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1.0071072]*6}

    helmholtz6 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for isobutane of Sun and Ely (2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering application: Algorithm and  application to non-polar and polar fluids", 
                    "ref": "Fluid Phase Equilib., 222-223:107-118, 2004.",
                    "doi": "10.1016/j.fluid.2004.06.028"}, 
        "R": 8.314472,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40., 
        "Pmin": 0.1, "rhomin": 40., 

        "nr1": [1.18083775, 9.46903331e-1, -2.90618044, 8.51346220e-2,
                2.79868503e-4, -1.68266335e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],
  
        "nr2": [-2.01202825e-1, -3.32570120e-2, 2.42967225e-1, -4.20931100e-3,
                -2.24528572e-1, -1.41307663e-2, -5.93401702e-2, -2.27862942e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    eq = helmholtz1, MBWR, GERG, helmholtz3, helmholtz4, helmholtz5, helmholtz6

    _surface = {"sigma": [0.05756, -0.009554], "exp": [1.29, 2.29]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [0.388417],  "expt0": [-1.], "expd0": [1.],
                   "a1": [20.534, 0.02], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [126.25, 52.91, -7501.4, -2672.9],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.9, 2.9]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.000022891,
                "Tmin": Tt, "Tmax": 575.0,
                "a1": [-1953637129., 1953637130.], "exp1": [0, 6.12],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.85093103, 1.36543198, -1.32542691, -2.56190994],
        "exp": [1, 1.5, 2.5, 4.5]}
    _liquid_Density = {
        "eq": 2,
        "ao": [2.04025104, 0.850874089, -0.479052281, 0.348201252],
        "exp": [1.065, 3, 4, 7]}
    _vapor_Density = {
        "eq": 6,
        "ao": [-2.12933323, -2.93790085, -0.89441086, -3.46343707],
        "exp": [1.065, 2.5, 9.5, 13]}

    visco0 = {"eq": 1, "omega": 1,
              "collision": [0.53583008, -0.45629630, 0.049911282],
              "__name__": "Vogel (2000)",
              "__doc__": """Vogel, E., Kuechenmeister, C., and Bich, E., "Viscosity Correlation for Isobutane over Wide Ranges of the Fluid Region," Int. J. Thermophys, 21(2):343-356, 2000.""",
              "ek": 307.55, "sigma": 0.46445,
              "Tref": 1., "rhoref": 1.*M,
              "n_chapman": 0.1628213/M**0.5,

              "n_virial": [-19.572881, 219.73999, -1015.3226, 2471.01251,
                           -3375.1717, 2491.6597, -787.26086, 14.085455,
                           -0.34664158],
              "t_virial": [0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -2.5, -5.5],
              "Tref_virial": 307.55, "etaref_virial": 0.0603345,

              "Tref_res": 407.817, "rhoref_res": 3.86*M, "etaref_res": 1,
              "n_packed": [0.233859774637e1, 0.235255150838e1],
              "t_packed": [0, 0.5],
              "n_poly": [0.103511763411e3, -0.312670896234e3, 0.145253750239e3,
                         -0.210649894193e3, 0.386269696509e3, -0.214963015527e3,
                         0.112580360920e3, -0.223242033154e3, 0.119114788598e3,
                         -0.181909745900e2, 0.360438957232e2, -0.213960184050e2,
                         -0.194037606990e4],
              "t_poly": [0, -1, -2, -0, -1, -2, 0, -1, -2, 0-1, -2, 0],
              "d_poly": [2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 1],
              "g_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1],
              "c_poly": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              "n_num": [0.194037606990e4],
              "t_num": [0],
              "d_num": [1],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1, -1],
              "t_den": [0, 0],
              "d_den": [0, 1],
              "g_den": [1, 0],
              "c_den": [0, 0]}

    visco1 = {"eq": 2, "omega": 2,
              "__name__": "Younglove (1987)",
              "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
              "ek": 418.0, "sigma": 0.509217,
              "n_chapman": 0.203525266/M**0.5,
              "F": [1.687838652, 0.0, 1.40, 407.85],
              "E": [-0.2055498053e2, 0.1357076181e4, 0.1893774336e2,
                    -0.1822277344e5, -0.4599387773e-2, 0.6305247065e2,
                    0.1282253921e5],
              "rhoc": 3.86}

    _viscosity = visco0, visco1

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2002)",
               "__doc__": """Perkins, R.A., "Measurement and Correlation of the Thermal Conductivity of Isobutane from 114 K to 600 K at Pressures to 70 MPa," J. Chem. Eng. Data, 47(5):1272-1279, 2002.""",

               "Tref": 407.85, "kref": 1,
               "no": [-2.37901e-3, 1.06601e-2, 2.15811e-2],
               "co": [0, 1, 2],

               "Trefb": 407.85, "rhorefb": 3.86, "krefb": 1,
               "nb": [-4.11789e-2, 4.76346e-2, 1.46805e-1, -1.28445e-1,
                      -1.19190e-1, 1.07565e-1, 4.10226e-2, -3.85968e-2,
                      -4.88704e-3, 5.20901e-3],
               "tb": [0, 1]*5,
               "db": [1, 1, 2, 2, 3, 3, 4, 4, 5, 5],
               "cb": [0]*10,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.657661e-9, "Tcref": 611.73}

    thermo1 = {"eq": 2, "omega": 2,
               "__name__": "Younglove (1987)",
               "__doc__": """Younglove, B.A. and Ely, J.F.,"Thermophysical properties of fluids. II. Methane, ethane, propane, isobutane and normal butane," J. Phys. Chem. Ref. Data, 16:577-798, 1987.""",
               "visco": visco1,
               "n_chapman": 2.0352526600e-1,
               "G": [0.1449797353e1, -0.1685643887],
               "E": [0.4307008989e-2, -0.1509010974e1, 0.4693712392e3,
                     -0.3554280979e-3, 0.1841552874, -0.3892338766e2,
                     -0.9354624917e-1, 0.7114330590e1],

               "critical": 2,
               "X": [0.0034718, 10.1207, 0.466392, 1.00344],
               "Z": 9.10218e-10}

    _thermal = thermo0, thermo1
