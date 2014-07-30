#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades

class nC7(MEoS):
    """Ecuación de estado de multiparametros para el n-heptano

    >>> heptano=nC7(T=300, P=0.1)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.5f %0.4f %0.4f %0.1f" % (heptano.T, heptano.rho, heptano.u.kJkg, heptano.h.kJkg, heptano.s.kJkgK, heptano.cv.kJkgK, heptano.cp.kJkgK, heptano.w)
    300.0 678.03 13.82 13.97 0.86061 1.7816 2.2475 1120.8
    """
    name="heptane"
    CASNumber="142-82-5"
    formula="CH3-(CH2)5-CH3"
    synonym=""
    rhoc=unidades.Density(231.9977)
    Tc=unidades.Temperature(540.13)
    Pc=unidades.Pressure(2736.0, "kPa")
    M=100.202      #g/mol
    Tt=unidades.Temperature(182.55)
    Tb=unidades.Temperature(371.53)
    f_acent=0.349
    momentoDipolar=unidades.DipoleMoment(0.07, "Debye")
    id=11
    _Tr=unidades.Temperature(525.389862)
    _rhor=unidades.Density(235.977855)
    _w=0.350780196

    CP1={  "ao": 4,
                "an": [], "pow": [],
                "ao_exp": [], "exp": [],
                "ao_hyp": [0.3957146e6, 0.2130579e8, 0.1349899e9, 0],
                "hyp": [0.1697890e3, 0.8361950e3, 0.1760460e4, 0]}

    CP2={  "ao": 4,
                "an": [], "pow": [],
                "ao_exp": [], "exp": [],
                "ao_hyp": [13.7266, 30.4707, 43.5561, 0],
                "hyp": [0.314348398*Tc, 1.54813656*Tc, 3.259326458*Tc, 0]}

    CP3={  "ao": 1.157528,
                "an": [0.70489617e-1, -0.23419686e-4, -0.14768221e-8, -0.20117611e-11],
                "pow": [1, 2, 3, 4],
                "ao_exp": [], "exp": [],
                "ao_hyp": [], "hyp": []}


    CP4={  "ao": 30.4029/8.3159524*4.184,
                "an": [], "pow": [],
                "ao_exp": [], "exp": [],
                "ao_hyp": [2.5273083e8/8.3159524*4.184, 3.9046536e7/8.3159524*4.184, 0, 0],
                "hyp": [1.6693200e3, 7.86001e2, 0, 0]}

    helmholtz1={
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for heptane of Span and Wagner (2003)",
        "__doc__":  u"""Span, R., Wagner, W. Equations of state for technical applications. II. Results for nonpolar fluids. Int. J. Thermophys. 24 (2003), 41 – 109.""",
        "R": 8.31451,
        "cp": CP1,

        "nr1":  [0.10543747645262e1, -0.26500681506144e1, 0.81730047827543, -0.30451391253428, 0.12253868710800, 0.27266472743928e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.49865825681670, -0.71432815084176e-3, -0.54236895525450, -0.13801821610756, -0.61595287380011e-2, 0.48602510393022e-3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    GERG={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heptane of Kunz and Wagner (2004).",
        "__doc__":  u"""Kunz, O., Klimeck, R., Wagner, W., Jaeschke, M. "The GERG-2004 Wide-Range Reference Equation of State for Natural Gases and Other Mixtures," to be published as a GERG Technical Monograph, Fortschr.-Ber. VDI, VDI-Verlag, Düsseldorf, 2006.""",
        "R": 8.314472,
        "cp": CP2,

        "nr1":  [0.10543747645262e1, -0.26500681506144e1, 0.81730047827543, -0.30451391253428, 0.122538687108, 0.27266472743928e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.49865825681670, -0.71432815084176e-3, -0.54236895525450, -0.13801821610756, -0.61595287380011e-2, 0.48602510393022e-3],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12.],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    helmholtz3={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heptane of Polt et al. (1992)",
        "__doc__":  u"""Polt, A., Platzer, B., and Maurer, G., "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe," Chem. Tech. (Leipzig), 44(6):216-224, 1992.""",
        "R": 8.3143,
        "cp": CP3,

        "nr1":  [-0.52030538102, 0.338196304523, -0.491117643215e-2, 0.200594802481, -0.260824422526e-1, -0.191516844204e1, 0.364407895089, -0.142523250539, -0.160069782510, 0.578283584822, 0.476898816887, 0.937511885529e-1, -0.442185898133, 0.553661375084e-1, -0.303420126133e-1, 0.138649129298e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.520305381020, -0.338196304523, 0.491117643215e-2, 0.256518106995e1, -0.528051955217e1, 0.266827442122e1, ],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [1]*6}

    helmholtz4={
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for heptane of Starling (1973)",
        "__doc__":  u"""Starling, K.E., "Fluid Thermodynamic Properties for Light Petroleum Systems," Gulf Publishing Company, 1973.""",
        "R": 8.3159524,
        "cp": CP4,

        "nr1":  [0.153471579811e1, 0.521386289098, -0.107860953728e1, -0.902616154206, 0.117182735038, -0.986768914864e-4, 0.287014205217, -0.359887681359, -0.860848441514e-2, 0.952855119365e-2, 0.227922178775e-3],
        "d1": [0, 1, 1, 1, 1, 1, 2, 2, 2, 5, 5],
        "t1": [3, 0, 1, 3, 4, 5, 0, 1, 2, 1, 2],

        "nr2": [-0.153471579811e1, -0.397448776976],
        "d2": [0, 2],
        "t2": [3, 3],
        "c2": [2]*2,
        "gamma2": [0.51794447]*2}

    eq=helmholtz1, GERG, helmholtz3, helmholtz4

    _surface={"sigma": [0.0541778, -0.00075856, 0.0039897], "exp": [1.25, 2.25, 3.25]}
    _dielectric={"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                            "a0": [0.10924],  "expt0": [-1.], "expd0": [1.],
                            "a1": [34.96, 0.035], "expt1": [0, 1], "expd1": [1, 1],
                            "a2": [162.24, 308.9, -37446, -39684], "expt2": [0, 1, 0, 1], "expd2": [2, 2, 3, 3]}
    _vapor_Pressure={ "eq": 5, "ao": [-0.76995e1, 0.14238e1, -0.20583e2, -0.50767e1, 0.19986e2], "exp": [1.0, 1.5, 3.4, 4.7, 3.6]}
    _liquid_Density={ "eq": 1, "ao": [-0.26395e1, 0.21806e2, -0.28896e2, 0.12609e2, 0.40749], "exp": [0.322, 0.504, 0.651, 0.816, 6.4]}
    _vapor_Density={ "eq": 3, "ao": [-0.24571, -0.63004e1, -0.19144e2, -0.96970e2, 0.21643e3, -0.27953e3], "exp": [0.097, 0.646, 2.56, 6.6, 9.3, 10.7]}

    visco0={"eq": 2, "omega": 3,
                "__name__": "NIST",
                "__doc__": """Coefficients are taken from NIST14, Version 9.08""",
                "ek": 400., "sigma": 0.64947,
                "n_chapman": 0.26718615/M**0.5,
                "F": [0, 0, 0, 100.],
                "E": [-17.168627495994, 3387.5906558, 16.943704644, -54960.940794, -0.24749641622, 163.37738185, 46932.568528],
                "rhoc": 2.315}

    visco1={"eq": 4, "omega": 1,
                "__doc__": """S.E.Quinones-Cisneros and U.K. Deiters, "Generalization of the Friction Theory for Viscosity Modeling," J. Phys. Chem. B 2006, 110,12820-12834.""",
                "__name__": "Quinones-Cisneros (2006)",
                "Tref": 540.13, "muref": 1.0,
                "ek": 400., "sigma": 0.64947, "n_chapman": 0,
                "n_ideal": [19.6036, -59.7839, 50.7528],
                "t_ideal": [0, 0.25, 0.5],

                "a": [3.76297120152080e-5, 0.0, -4.40242197269552e-5],
                "b": [1.38067766234763e-4, 0.0, -9.11095867363485e-5],
                "c": [9.93870811000000e-5, -6.36532780000000e-6, 0.0],
                "A": [-3.76786095828018e-9, 1.92499718242368e-9, 0.0],
                "B": [0.0, 9.75462662440927e-9, 2.71873666825660e-9],
                "C": [-1.24466129111157e-6, 8.83260990875321e-7, 0.0],
                "D": [0.0, 0.0, 0.0 ]}

    _viscosity=visco0, visco1

    thermo0={"eq": 1,
                "__name__": "NIST14",
                "__doc__": """Coefficients are taken from NIST14, Version 9.08""",

                "Tref": 400.0, "kref": 1e-3,
                "no": [1.35558587, -0.152682526035, 1],
                "co": [0, -1, -96],

                "Trefb": 540.15, "rhorefb": 2.315, "krefb": 1e-3,
                "nb": [15.900635275, 3.96318667803, -1.72336149946, 0.437228619593, 0.490514843565, -0.163256898944],
                "tb": [0, 0, 0, -1, 0, -1],
                "db": [1, 3, 4, 4, 5, 5],
                "cb": [0]*6,

                "critical": 3,
                "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
                "Xio": 0.194e-9, "gam0": 0.0496, "qd": 1.1246e-9, "Tcref": 810.195}

    _thermal=thermo0,



if __name__ == "__main__":
    import doctest
    doctest.testmod()
