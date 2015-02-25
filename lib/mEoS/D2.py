#!/usr/bin/python
# -*- coding: utf-8 -*-

from lib.meos import MEoS
from lib import unidades


class D2(MEoS):
    """Multiparameter equation of state for deuterium

    >>> deuterio=D2(T=300, P=0.1)
    >>> print "%0.1f %0.5f %0.3f %0.3f %0.5f %0.4f %0.4f %0.2f" % (deuterio.T, deuterio.rho, deuterio.u.kJkg, deuterio.h.kJkg, deuterio.s.kJkgK, deuterio.cv.kJkgK, deuterio.cp.kJkgK, deuterio.w)
    300.0 0.16039 -613.532 9.947 0.05958 3.1162 5.1932 1019.58
    """
    name = "deuterium"
    CASNumber = "7782-39-0"
    formula = "D2"
    synonym = ""
    rhoc = unidades.Density(69.7966)
    Tc = unidades.Temperature(38.34)
    Pc = unidades.Pressure(1665.3, "kPa")
    M = 4.0282  # g/mol
    Tt = unidades.Temperature(18.71)
    Tb = unidades.Temperature(23.3097)
    f_acent = -0.175
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 1

    CP1 = {"ao": 2.4512991,
           "an": [0.43563077e-02, -0.53169470e-03, 0.17067184e-04,
                  -0.53819932e-08, 0.89310438e-12],
           "pow": [1, 1.5, 2, 3, 4],
           "ao_exp": [0.18403263e2, -0.21257617e2, 0.41091635e1],
           "exp": [319, 361, 518],
           "ao_hyp": [], "hyp": []}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for deuterium of McCarty (1989)",
        "__doc__": u"""McCarty, R.D., "Correlations for the Thermophysical Properties of Deuterium," National Institute of Standards and Technology, Boulder, CO, 1989.""",
        "R": 8.31434,
        "cp": CP1,

        "Tmin": Tt, "Tmax": 423.0, "Pmax": 320000.0, "rhomax": 43.38, 
        "Pmin": 19.462, "rhomin": 43.365, 

        "b": [None, 0.4894244053982e-4, 0.5600164604601e-1, -0.6301493491211,
              0.2538329946038e1, 0.1723475985309e3, 0.2956238369436e-4,
              -0.3926317169317e-2, 0.1195764193293e-1, 0.1136916678824e5,
              -0.1916378195727e-6, 0.3153535946452e-3, 0.2122937335070e-1,
              -0.1057999371607e-5, -0.6722062598854e-4, -0.3030166828627,
              0.1980817195099e-5, -0.1453922641871e-7, 0.1780919116891e-3,
              -0.1823145348424e-5, -0.1135358616578e5, -0.1943542941899e4,
              -0.3632847669580e2, 0.1087745118380e3, -0.4078276062687e-1,
              0.6460021864005e-2, -0.4480242189217e-4, -0.2475011206216e-3,
              -0.8834384656760e-8, -0.1081622159862e-8, -0.1478159334303e-10,
              0.7926922356112e-11, 0.5721547329378e-11]}

    eq = MBWR,

    _surface = {"sigma": [0.00795], "exp": [1.]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.31333e1, -0.53100e1, 0.26450e2, -0.51890e2, 0.31468e2],
        "exp": [1., 1.5, 2.7, 3.5, 4.1]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.88968e2, -0.32365e3, 0.55793e3, -0.46395e3, 0.38200e3, -0.24647e3],
        "exp": [0.98, 1.25, 1.6, 3.0, 3.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.10716e1, -0.19044e3, 0.53964e3, -0.93931e3, 0.36302e3, 0.21020e3],
        "exp": [0.1, 2.25, 2.75, 3.4, 3.6, 3.8]}
