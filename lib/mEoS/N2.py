#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class N2(MEoS):
    """Multiparamente equation of state for nitrogen

    >>> nitrogeno=N2(T=300, P=.1)
    >>> print "%0.1f %0.4f %0.3f %0.2f %0.4f %0.5f %0.4f %0.2f" % (nitrogeno.T, nitrogeno.rho, nitrogeno.u.kJkg, nitrogeno.h.kJkg, nitrogeno.s.kJkgK, nitrogeno.cv.kJkgK, nitrogeno.cp.kJkgK, nitrogeno.w)
    300.0 1.1233 222.171 311.20 6.8457 0.74316 1.0413 353.16
    """
    name = "nitrogen"
    CASNumber = "7727-37-9"
    formula = "N2"
    synonym = "R-728"
    rhoc = unidades.Density(313.3)
    Tc = unidades.Temperature(126.192)
    Pc = unidades.Pressure(3395.8, "kPa")
    M = 28.01348  # g/mol
    Tt = unidades.Temperature(63.151)
    Tb = unidades.Temperature(77.355)
    f_acent = 0.0372
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 46
    _Tr = unidades.Temperature(122.520245)
    _rhor = unidades.Density(316.134310)
    _w = 0.043553140

    CP1 = {"ao": 3.5,
           "an": [3.066469e-6, 4.70124e-9, -3.987984e-13], "pow": [1, 2, 3],
           "ao_exp": [1.012941], "exp": [3364.011],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 3.50418363823,
           "an": [-0.837079888737e3, 0.379147114487e2, -0.601737844275,
                  -0.874955653028e-5, 0.148958507239e-7, -0.256370354277e-11],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [1.00773735767], "exp": [3353.4061],
           "ao_hyp": [], "hyp": []}

    CP3 = {"ao": 3.50031,
           "an": [], "pow": [],
           "ao_exp": [], "exp": [],
           "ao_hyp": [0.13732, -0.1466, 0.90066, 0],
           "hyp": [5.251822620*Tc, -5.393067706*Tc, 13.788988208*Tc, 0],
           "R": 8.31451}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Span et al. (2000).",
        "__doc__":  u"""Span, R., Lemmon, E.W., Jacobsen, R.T, Wagner, W., Yokozeki, A. A reference equation of state for the thermodynamic properties of nitrogen for temperatures from 63.151 to 1000 K and pressures to 2200 MPa. J. Phys. Chem. Ref. Data 29 (2000), 1361 – 1433.""",
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 2200000.0, "rhomax": 53.15, 
        "Pmin": 12.5198, "rhomin": 30.957, 

        "nr1": [0.924803575275, -0.492448489428, 0.661883336938,
                -0.192902649201e1, -0.622469309629e-1, 0.349943957581],
        "d1": [1, 1, 2, 2, 3, 3],
        "t1": [0.25, 0.875, 0.5, 0.875, 0.375, 0.75],

        "nr2": [0.564857472498, -0.161720005987e1, -0.481395031883,
                0.421150636384, -0.161962230825e-1, 0.172100994165,
                0.735448924933e-2, 0.168077305479e-1, -0.107626664179e-2,
                -0.137318088513e-1, 0.635466899859e-3, 0.304432279419e-2,
                -0.435762336045e-1, -0.723174889316e-1, 0.389644315272e-1,
                -0.212201363910e-1, 0.408822981509e-2, -0.551990017984e-4,
                -0.462016716479e-1, -0.300311716011e-2, 0.368825891208e-1,
                -0.255856846220e-2, 0.896915264558e-2, -0.441513370350e-2,
                0.133722924858e-2, 0.264832491957e-3],
        "d2": [1, 1, 1, 3, 3, 4, 6, 6, 7, 7, 8, 8, 1, 2, 3, 4, 5, 8, 4, 5, 5,
               8, 3, 5, 6, 9],
        "t2": [0.5, 0.75, 2., 1.25, 3.5, 1., 0.5, 3., 0., 2.75, 0.75, 2.5, 4.,
               6., 6., 3., 3., 6., 16., 11., 15., 12., 12., 7., 4., 16.],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3,
               3, 4, 4, 4, 4],
        "gamma2": [1]*26,

        "nr3": [0.196688194015e2, -0.209115600730e2, 0.167788306989e-1,
                0.262767566274e4],
        "d3": [1, 1, 3, 2],
        "t3": [0., 1., 2., 3.],
        "alfa3": [20, 20, 15, 25],
        "beta3": [325, 325, 300, 275],
        "gamma3": [1.16, 1.16, 1.13, 1.25],
        "epsilon3": [1]*4,
        "nr4": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for deuterium of McCarty (1989)",
        "__doc__": u"""McCarty, R.D., "Correlations for the Thermophysical Properties of Deuterium," National Institute of Standards and Technology, Boulder, CO, 1989.""",
        "R": 8.31434,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 1900.0, "Pmax": 1013000.0, "rhomax": 30.977, 
        "Pmin": 12.463, "rhomin": 30.977, 

        "b": [None, 0.1380297474657e-2, 0.1084506501349, -0.2471324064362e1,
              0.3455257980807e2, -0.4279707690666e4, 0.1064911566998e-3,
              -0.1140867079735e-1, 0.1444902497287e-3, 0.1871457567553e5,
              0.8218876886831e-7, 0.2360990493348e-2, -0.5144803081201,
              0.4914545013668e-4, -0.1151627162399e-2, -0.7168037246650,
              0.7616667619500e-4, -0.1130930066213e-5, 0.3736831166831e-3,
              -0.2039851507581e-5, -0.1719662008990e5, -0.1213055199748e6,
              -0.9881399141428e2, 0.5619886893511e5, -0.1823043964118,
              -0.2599826498477e1, -0.4191893423157e-3, -0.2596406670530,
              -0.1258683201921e-6, 0.1049286599400e-4, -0.5458369305152e-9,
              -0.7674511670597e-8, 0.5931232870994e-7]}

    GERG = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Kunz and Wagner (2004).",
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP3,

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 2200000.0, "rhomax": 53.15, 
#        "Pmin": 73.476, "rhomin": 29.249, 

        "nr1": [0.59889711801201, -0.16941557480731e1, 0.24579736191718,
                -0.23722456755175, 0.17954918715141e-1, 0.14592875720215e-1],
        "d1": [1, 1, 2, 2, 4, 4],
        "t1": [0.125, 1.125, 0.375, 1.125, 0.625, 1.5],

        "nr2": [0.10008065936206, 0.73157115385532, -0.88372272336366,
                0.31887660246708, 0.20766491728799, -0.19379315454158e-1,
                -0.16936641554983, 0.13546846041701, -0.33066712095307e-1,
                -0.60690817018557e-1, 0.12797548292871e-1, 0.58743664107299e-2,
                -0.18451951971969e-1, 0.47226622042472e-2, -0.52024079680599e-2,
                0.43563505956635e-1, -0.36251690750939e-1, -0.28974026866543e-2],
        "d2": [1, 1, 1, 2, 3, 6, 2, 3, 3, 4, 4, 2, 3, 4, 5, 6, 6, 7],
        "t2": [0.625, 2.625, 2.75, 2.125, 2, 1.75, 4.5, 4.75, 5, 4, 4.5, 7.5,
               14, 11.5, 26, 28, 30, 16],
        "c2": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 6, 6, 6, 6],
        "gamma2": [1]*18,

        "nr3": [],
        "nr4": []}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for nitrogen of Jacobsen et al. (1986).",
        "__doc__":  u"""Jacobsen, R.T, Stewart, R.B., and Jahangiri, M., "Thermodynamic properties of nitrogen from the freezing line to 2000 K at pressures to 1000 MPa," J. Phys. Chem. Ref. Data, 15(2):735-909, 1986.""",
        "R": 8.31434,
        "cp": CP2,

        "Tmin": Tt, "Tmax": 2000.0, "Pmax": 1000000.0, "rhomax": 30.96, 
        "Pmin": 12.52, "rhomin": 31.046, 

        "nr1": [0.9499541827, 0.2481718513, -0.2046287122, -0.1748429008,
                0.6387017148, -0.5272986168, -0.2049741504e1, 0.5551383553e-1,
                -0.8191106396e-3, -0.5032519699e-1, 0.2650110798, 0.7311459372e-1,
                -0.2813080718e-1, 0.1659823569e-2, 0.6012817812e-1,
                -0.3785445194, 0.1895290433, -0.7001895093e-2],
        "d1": [1, 2, 3, 2, 3, 3, 1, 4, 6, 2, 1, 2, 4, 6, 2, 1, 2, 4],
        "t1": [0.25, 0.25, 0.25, 0.5, 0.5, 0.75, 1, 1, 1, 1, 1.5, 2, 2, 2, 2,
               3, 3, 3],

        "nr2": [-0.4927710927e-1, 0.6512013679e-1, 0.113812194200,
                -0.955140963197e-1, 0.2118354140e-1, -0.1100721771e-1,
                0.1284432210e-1, -0.1054474910e-1, -0.1484600538e-3,
                -0.5806483467e-2],
        "d2": [1, 4, 1, 2, 4, 2, 4, 4, 2, 3],
        "t2": [3, 4, 4, 5, 6, 8, 14, 18, 20, 22],
        "c2": [3, 2, 3, 2, 2, 4, 4, 4, 4, 3],
        "gamma2": [1]*10,

        "nr3": [],
        "nr4": []}

    helmholtz4 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for nitrogen of Span and Wagner (2003).",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. II. Results for Nonpolar Fluids," Int. J. Thermophys., 24(1):41-109, 2003.""",
        "R": 8.31451,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 600.0, "Pmax": 100000.0, "rhomax": 53.15, 
        "Pmin": 12.566, "rhomin": 30.935, 

        "nr1": [0.92296567, -0.25575012e1, 0.64482463, 0.1083102e-1,
                0.73924167e-1, 0.23532962e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.18024854, -0.45660299e-1, -0.1552106, -0.3811149e-1,
                -0.31962422e-1, 0.15513532e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6,

        "nr3": [],
        "nr4": []}

    eq = helmholtz1, MBWR, GERG, helmholtz3, helmholtz4
    _PR = -0.004032

    _surface = {"sigma": [0.029324108], "exp": [1.259]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [4.3872, 0.00226], "expt1": [0, 1], "expd1": [1, 1],
                   "a2": [2.206, 1.135, -169., -35.83],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3.1, 3.1]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 12.523,
                "Tmin": Tt, "Tmax": 2000.0,
                "a1": [1, 12798.61, -12798.61], "exp1": [0, 1.78963, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _sublimation = {"eq": 3, "Tref": Tt, "Pref": 12.523,
                    "Tmin": Tt, "Tmax": Tt,
                    "a1": [], "exp1": [],
                    "a2": [-13.088692], "exp2": [1],
                    "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 6,
        "ao": [-0.612445284e1, 0.126327220e1, -0.765910082, -0.177570564e1],
        "exp": [2, 3, 5, 10]}
    _liquid_Density = {
        "eq": 4,
        "ao": [0.148654237e1, -0.280476066, 0.894143085e-1, -0.119879866],
        "exp": [0.9882, 2, 8, 17.5]}
    _vapor_Density = {
        "eq": 6,
        "ao": [-0.170127164e1, -0.370402649e1, 0.129859383e1, -0.561424977,
               -0.268505381e1],
        "exp": [1.02, 2.5, 3.5, 6.5, 14]}

    visco0 = {"eq": 1, "omega": 1,
              "__name__": "Lemmon (2004)",
              "__doc__": """Lemmon, E.W. and Jacobsen, R.T, "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air," Int. J. Thermophys., 25:21-69, 2004.""",
              "Tref": 1., "etaref": 1,
              "ek": 98.94, "sigma": 0.3656,
              "n_chapman": 0.141294895/M**0.5,

              "Tref_res": 126.192, "rhoref_res": 11.1839*M, "etaref_res": 1,
              "n_poly": [10.72, 0.03989, 0.001208, -7.402, 4.62],
              "t_poly": [-.1, -.25, -3.2, -.9, -0.3],
              "d_poly": [2, 10, 12, 2, 1],
              "g_poly": [0, 0, 0, 1, 1],
              "c_poly": [0, 1, 1, 2, 3]}

    visco1 = {"eq": 2, "omega": 2,
              "collision": [-136.985150760851, 734.241371453542, -1655.39131952744,
                            2062.67809686969, -1579.52439123889, 777.942880032361,
                            -232.996787901831, 40.0691427576552, -2.99482706239363],
              "__name__": "Younglove (1982)",
              "__doc__": """Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
              "ek": 118., "sigma": 0.354,
              "n_chapman": 0.141286429751707,
              "t_chapman": 0.0,
              "F": [-3.14276193277e-3, 9.22071479907e-4, 1.4, 118],
              "E": [-12.128154129, 68.46443564, 11.2569594404402, -565.76279020055,
                    9.56677570672e-2, -.355533724265011, 618.536783201947],
              "rhoc": 11.2435750999429}

    visco2 = {"eq": 1, "omega": 1,
              "collision": [0.46649, -0.57015,  0.19164, -0.03708,  0.00241],
              "__name__": "Stephan (2004)",
              "__doc__": """Stephan, K., Krauss, R., and Laesecke, A., "Viscosity and Thermal Conductivity of Nitrogen for a Wide Range of Fluid States," J. Phys. Chem. Ref. Data, 16(4):993-1023, 1987""",
              "Tref": 1., "etaref": 1,
              "ek": 100.01654, "sigma": 0.36502496,
              "n_chapman": 0.141290/M**0.5,

              "Tref_res": 1, "rhoref_res": 11.2088889*M, "etaref_res": 14.,
              "n_poly": [-5.8470232, -1.4470051, -0.27766561e-1, -0.21662362],
              "t_poly": [0, 0, 0, 0],
              "d_poly": [0, 1, 2, 3],
              "g_poly": [0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0],
              "n_num": [-20.09997],
              "t_num": [0],
              "d_num": [0],
              "g_num": [0],
              "c_num": [0],
              "n_den": [1.0, -3.4376416],
              "t_den": [0, 0],
              "d_den": [1, 0],
              "g_den": [0, 0],
              "c_den": [0, 0]}

    visco3 = {"eq": 1, "omega": 1,
              "collision": [0.5136, -0.5218, 0.8852e-1, 0.3445e-2, -0.2289e-2],
              "__name__": "Lemmon (2001)",
              "__doc__": """Lemmon, E.W. and Jacobsen, R.T, unpublished equation, 2001.""",
              "Tref": 1., "etaref": 1,
              "ek": 71.4, "sigma": 0.3798,
              "n_chapman": 0.141296/M**0.5,

              "Tref_res": 126.192, "rhoref_res": 11.1839*M, "etaref_res": 1,
              "n_poly": [0.4310e1, 0.3135e1, 0.2695e1, 0.1273e-2, 0.5795e-1],
              "t_poly": [0, 0, -0.2, -2.3, -0.1],
              "d_poly": [1, 2, 3, 10, 9],
              "g_poly": [0, 0, 0, 0, 0],
              "c_poly": [0, 0, 0, 0, 1]}

    _viscosity = visco0, visco1, visco2, visco3

    thermo0 = {"eq": 1,
               "__name__": "Lemmon (2004)",
               "__doc__": """Lemmon, E.W. and Jacobsen, R.T, "Viscosity and Thermal Conductivity Equations for Nitrogen, Oxygen, Argon, and Air," Int. J. Thermophys., 25:21-69, 2004.""",

               "Tref": 126.192, "kref": 1e-3,
               "no": [1.511, 2.117, -3.332],
               "co": [-97, 1, 0.7],

               "Trefb": 126.192, "rhorefb": 11.1839, "krefb": 1e-3,
               "nb": [8.862, 31.11, -73.13, 20.03, -0.7096, 0.2672],
               "tb": [0, -0.03, -0.2, -0.8, -0.6, -1.9],
               "db": [1, 2, 3, 4, 8, 10],
               "cb": [0, 0, 1, 2, 2, 2],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.17e-9, "gam0": 0.055, "qd": 0.40e-9, "Tcref": 252.384}

    thermo1 = {"eq": 3,
               "__name__": "Younglove (1982)",
               "__doc__": """Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",

               "ek": 118, "sigma": 0.354,
               "Nchapman": 0.141286429751707,
               "tchapman": 0,
               "b": [-.15055520615565, 0.183477124982509, 1.45008451566007,
                     -4.88031780663869, 6.68390592664363, -4.90242883649539,
                     2.02630917877999, -.439826733340102, 3.91906706514e-2],
               "F": [1.50938067650e-3, 1.70975795748e-4, 1.2, 118],
               "E": [-38.613291627, -31.826109485, 26.0197970589236,
                     -27.2869897441495, 0, 0, 0],
               "rhoc": 35.6938892061679,
               "ff": 1.67108,
               "rm": 0.00000003933}

    thermo2 = {"eq": 1, "critical": 0,
               "__name__": "Stephan (2004)",
               "__doc__": """Stephan, K., Krauss, R., and Laesecke, A., "Viscosity and Thermal Conductivity of Nitrogen for a Wide Range of Fluid States," J. Phys. Chem. Ref. Data, 16(4):993-1023, 1987""",

               "Tref": 1, "kref": 1e-3,
               "no": [0.6950401, 0.03643102],
               "co": [-97, -98],

               "Trefb": 1, "rhorefb": 11.2088889, "krefb": 4.17e-3,
               "nb": [3.3373542, 0.37098251, 0.89913456, 0.16972505],
               "tb": [0, 0, 0, 0],
               "db": [1, 2, 3, 4],
               "cb": [0, 0, 0, 0]}

    thermo3 = {"eq": 1,
               "__name__": "Lemmon (2001)",
               "__doc__": """Lemmon, E.W. and Jacobsen, R.T, unpublished equation, 2001.""",
               "Tref": 126.192, "kref": 1e-3,
               "no": [1.113, 0],
               "co": [0, -96],

               "Trefb": 126.192, "rhorefb": 11.1839, "krefb": 1e-3,
               "nb": [0.9153e1, 0.6094, 0.3466e2, -0.1138e2, -0.9117e2,
                      0.4004e1, 0.1886e2, 0.6146e1, 0.1587e2],
               "tb": [0, -2.5, -.14, -3, -.4, -.75, -.75, -3.85, -9],
               "db": [1, 1, 2, 1, 3, 2, 4, 6, 2],
               "cb": [0, 0, 0, 1, 1, 2, 2, 2, 3],

               "critical": 3,
               "gnu": 0.63, "gamma": 1.2415, "R0": 1.01,
               "Xio": 0.165e-9, "gam0": 0.055, "qd": 0.386e-9, "Tcref": 252.384}

    _thermal = thermo0, thermo1, thermo2, thermo3
