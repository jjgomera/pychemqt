#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp

from lib.meos import MEoS
from lib import unidades


class Ethylene(MEoS):
    """Multiparameter equation of state for ehylene

    >>> etileno=Ethylene(T=315, P=200.)
    >>> print "%0.1f %0.2f %0.2f %0.2f %0.4f %0.4f %0.4f %0.1f" % (etileno.T, etileno.rho, etileno.u.kJkg, etileno.h.kJkg, etileno.s.kJkgK, etileno.cv.kJkgK, etileno.cp.kJkgK, etileno.w)
    315.0 579.15 -453.64 -108.31 -2.7619 1.6569 2.2448 1679.8
    """
    name = "ethylene"
    CASNumber = "74-85-1"
    formula = "CH2=CH2"
    synonym = "R-1150"
    rhoc = unidades.Density(214.24656512)
    Tc = unidades.Temperature(282.35)
    Pc = unidades.Pressure(5041.8, "kPa")
    M = 28.05376  # g/mol
    Tt = unidades.Temperature(103.986)
    Tb = unidades.Temperature(169.379)
    f_acent = 0.0866
    momentoDipolar = unidades.DipoleMoment(0.0, "Debye")
    id = 22
    _Tr = unidades.Temperature(273.316763)
    _rhor = unidades.Density(216.108926)
    _w = 0.085703183

    CP1 = {"ao": 4.,
           "an": [], "pow": [],
           "ao_exp": [1]*12,
           "exp": [4353.90715, 2335.22515, 1930.91322, 1471.92565, 4464.69725,
                   1778.39697, 1365.45204, 1356.81905, 4469.01375, 1188.47564,
                   4300.67034, 2077.67413],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 3.554495281,
           "an": [0.5603615762e6, -0.2141069802e5, 0.2532008897e3,
                  -0.9951927478e-2, 0.5108931070e-4, -0.1928667482e-7],
           "pow": [-3, -2, -1.001, 1, 2, 3],
           "ao_exp": [-0.2061703241e2],
           "exp": [3000],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Smukala et al. (2000)",
        "__doc__":  u"""Smukala, J., Span, R., Wagner, W. New equation of state for ethylene covering the fluid region from the melting line to 450 K at pressures up to 300 MPa. J. Phys. Chem. Ref. Data 29 (2000), 1053 – 1121.""",
        "R": 8.31451,
        "cp": CP1,

        "nr1": [0.18617429100670e1, -0.30913708460844e1, -0.17384817095516,
                0.80370985692840e-1, 0.23682707317354, 0.21922786610247e-1],
        "d1": [1, 1, 1, 2, 2, 4],
        "t1": [0.5, 1., 2.5, 0.0, 2.0, 0.5],

        "nr2": [0.11827885813193, -0.21736384396776e-1, 0.44007990661139e-1,
                0.12554058863881, -0.13167945577241, -0.52116984575897e-2,
                0.15236081265419e-3, -0.24505335342756e-4, 0.28970524924022,
                -0.18075836674288, 0.15057272878461, -0.14093151754458,
                0.22755109070253e-1, 0.14026070529061e-1, 0.61697454296214e-2,
                -0.41286083451333e-3, 0.12885388714785e-1, -0.69128692157093e-1,
                0.10936225568483, -0.81818875271794e-2, -0.56418472117170e-1,
                0.16517867750633e-2, 0.95904006517001e-2, -0.26236572984886e-2],
        "d2": [1, 1, 3, 4, 5, 7, 10, 11, 1, 1, 2, 2, 4, 4, 6, 7, 4, 5, 6, 6, 7,
               8, 9, 10],
        "t2": [1., 4., 1.25, 2.75, 2.25, 1., 0.75, 0.5, 2.5, 3.5, 4., 6., 1.5,
               5., 4.5, 15., 20., 23., 22., 29., 19., 15., 13., 10.],
        "c2": [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4,
               4, 4, 4],
        "gamma2": [1]*24,

        "nr3": [-0.50242414011355e2, 0.74846420119299e4, -0.68734299232625e4,
                -0.93577982814338e3, 0.94133024786113e3],
        "d3": [2, 2, 2, 3, 3],
        "t3": [1., 0., 1., 2., 3.],
        "alfa3": [25.]*5,
        "beta3": [325, 300, 300, 300, 300],
        "gamma3": [1.16, 1.19, 1.19, 1.19, 1.19],
        "epsilon3": [1.]*5}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for ethylene of Jahangiri et al. (1986)",
        "__doc__":  u"""Jahangiri, M., Jacobsen, R.T, Stewart, R.B., and McCarty, R.D., "Thermodynamic properties of ethylene from the freezing line to 450 K at pressures to 260 MPa," J. Phys. Chem. Ref. Data, 15(2):593-734, 1986.""",
        "R": 8.31434,
        "cp": CP1,

        "nr1": [0.324893703388e1, -0.101727886161e2, 0.738660405252e1,
                -0.156891635862e1, -0.888451428662e-1, 0.602106814262e-1,
                0.107832458846, -0.200402521069e-1, 0.195049141244e-2,
                0.671800640346e-1, -0.420045146918e-1, -0.162050762577e-2,
                0.555515679497e-3, 0.758367114630e-3, -0.287854402074e-3],
        "d1": [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 6, 6, 6],
        "t1": [0.5, 1, 1.25, 1.75, 4, 2, 4, 5, 6, 0.25, 3, 0.25, 0.5, 2.5, 3],

        "nr2": [0.6258987063e-1, -0.641843116e-1, -0.1368693752, 0.517920766,
                -0.3026331319, 0.7757213872, -0.2639890864e1, 0.2927563554e1,
                -0.1066267599e1, -0.538047154e-1, 0.127792108, -0.745015231e-1,
                -0.1624304356e-1, 0.1476032429, -0.2003910489, 0.2926905618,
                -0.1389040901, 0.5913513541e1, -0.380037013e2, 0.969194057e2,
                -0.1226256839e3, 0.7702379476e2, -0.1922684672e2,
                -0.3800045701e-2, 0.1118003813e-1, 0.2945841426e-2],
        "d2": [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4,
               4, 4, 8, 8, 8],
        "t2": [0.5, 1, 0.5, 2, 4, 3, 4, 5, 6, 2, 3, 4, 1.5, 0.5, 1.5, 4, 5, 1,
               2, 3, 4, 5, 6, 0.5, 1, 5],
        "c2": [3, 3, 2, 2, 2, 4, 4, 4, 4, 6, 6, 6, 3, 2, 2, 2, 2, 4, 4, 4, 4,
               4, 4, 2, 2, 2],
        "gamma2": [1]*26}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "short Helmholtz equation of state for ethylene of Span and Wagner (2003)",
        "__doc__":  u"""Span, R. and Wagner, W. "Equations of State for Technical Applications. II. Results for Nonpolar Fluids," Int. J. Thermophys., 24(1):41-109, 2003.""",
        "R": 8.31451,
        "cp": CP1,

        "nr1": [0.9096223, -0.24641015e1, 0.56175311, -0.19688013e-1,
                0.78831145e-1, 0.21478776e-3],
        "d1": [1, 1, 1, 2, 3, 7],
        "t1": [0.25, 1.125, 1.5, 1.375, 0.25, 0.875],

        "nr2": [0.23151337, -0.37804454e-1, -0.20122739, -0.44960157e-1,
                -0.2834296e-1, 0.12652824e-1],
        "d2": [2, 5, 1, 4, 3, 4],
        "t2": [0.625, 1.75, 3.625, 3.625, 14.5, 12],
        "c2": [1, 1, 2, 2, 3, 3],
        "gamma2": [1]*6}

    MBWR = {
        "__type__": "MBWR",
        "__name__": "MBWR equation of state for ethylene of McCarty and Jacobsen (1981)",
        "__doc__":  u"""Younglove, B.A., "Thermophysical Properties of Fluids.  I. Argon, Ethylene, Parahydrogen, Nitrogen, Nitrogen Trifluoride, and Oxygen," J. Phys. Chem. Ref. Data, Vol. 11, Suppl. 1, pp. 1-11, 1982.""",
        "R": 8.31434,
        "cp": CP2,

        "b": [None, -0.2146684366683e-1, 0.1791433722534e1, -0.3675315603930e2,
              0.3707178934669e4, -0.3198282566709e6, 0.5809379774732e-3,
              -0.7895570824899, 0.1148620375835e3, 0.2713774629193e6,
              -0.8647124319107e-4, 0.1617727266385, -0.2731527496271e2,
              -0.2672283641459e-2, -0.4752381331990e-1, -0.6255637346217e2,
              0.4576234964434e-2, -0.7534839269320e-4, 0.1638171982209,
              -0.3563090740740e-2, -0.1833000783170e6, -0.1805074209985e8,
              -0.4794587918874e4, 0.3531948274957e8, -0.2562571039155e2,
              0.1044308253292e4, -0.1695303363659, -0.1710334224958e4,
              -0.2054114462372e-3, 0.6727558766661e-1, -0.1557168403328e-5,
              -0.1229814736077e-3, 0.4234325938573e-3]}

    eq = helmholtz1, MBWR, helmholtz2, helmholtz3

    _surface = {"sigma": [0.050195], "exp": [1.26]}
    _dielectric = {"eq": 3, "Tref": 273.16, "rhoref": 1000.,
                   "a0": [],  "expt0": [], "expd0": [],
                   "a1": [10.725], "expt1": [0], "expd1": [1],
                   "a2": [55.19, 49.5, -2045, -1154.],
                   "expt2": [0, 1, 0, 1], "expd2": [2, 2, 2.9, 2.9]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 1000,
                "a1": [0.1225e-3, 0.357924e3, -0.357924e3],
                "exp1": [0, 0.20645e1, 0],
                "a2": [], "exp2": [], "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-6.3905741, 1.4060338, -1.6589923, 1.0278028, -2.5071716],
        "exp": [1.0, 1.5, 2.5, 3.0, 4.5]}
    _liquid_Density = {
        "eq": 4,
        "ao": [1.8673079, -0.61533892, -0.058973772, 0.10744720],
        "exp": [1.029, 1.5, 4.0, 6.0]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-1.9034556, -0.75837929, -3.7717969, -8.7478586, -23.885296,
               -54.197979],
        "exp": [1.047, 2.0, 3.0, 7.0, 14.5, 28.0]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "Holland (1983)"}

    visco1 = {"eq": 2, "omega": 2,
              "__name__": "NIST",
              "__doc__": """Coefficients are taken from NIST14, Version 9.08""",
              "ek": 224.7, "sigma": 0.4163,
              "n_chapman": 0.141374566253583/M**0.5,
              "F": [0, 0, 0, 100.],
              "E": [-8.03553028329404, -439.8962514, 8.69536237617, 5773.08496161,
                    0.267589139152, -34.39391627, 66.4795135739],
              "rhoc": 7.63299886259}

    _viscosity = visco0, visco1

    def _visco0(self):
        """Holland, P.M., Eaton, B.E., and Hanley, H.J.M., "A Correlation of the Viscosity and Thermal Conductivity Data of Gaseous and Liquid Ethylene," J. Phys. Chem. Ref. Data, 12(4):917-932, 1983."""
        GV = [-3.5098225018e6, 2.5008406184e6, -5.8365540744e5, 4.5549146583e3,
              2.2881683403e4, -4.7318682077e3, 4.5022249258e2, -2.1490688088e1,
              4.1649263233e-1]
        muo = 0
        for i in range(-3, 6):
            muo += GV[i+3]*self.T**(i/3.)

        mu1 = 0
        tita = (self.rho-self.rhoc)/self.rhoc
        j = [0, -4.8544486732, 1.3033585236e1, 2.7808928908e4, -1.8241971308e3,
             1.5913024509, -2.0513573927e2, -3.9478454708e4]
        deltamu = exp(j[1]+j[4]/self.T)*(exp(self.rho.gcc**0.1*(j[2]+j[3]/self.T**1.5)+tita*self.rho.gcc**0.5*(j[5]+j[6]/self.T+j[7]/self.T**2))-1.)

        return unidades.Viscosity((muo+mu1+deltamu)*1e-7, "Pas")

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "Holland (1983)"}

    thermo1 = {"eq": 1, "critical": 0,
               "__name__": "NIST14",
               "__doc__": """Coefficients are taken from NIST14, Version 9.08""",

               "Tref": 224.7, "kref": 1e-3,
               "no": [1.35558587, -0.14207565869509, 1],
               "co": [0, -1, -96],

               "Trefb": 282.350007277, "rhorefb": 7.63299886259, "krefb": 1e-3,
               "nb": [15.3064493136, 25.0280721432, -15.4526955192,
                      0.8590418672, 3.32700049633, -0.333048907849],
               "tb": [0, 0, 0, -1, 0, -1],
               "db": [1, 3, 4, 4, 5, 5],
               "cb": [0]*6}

    _thermal = thermo0, thermo1

    def _thermo0(self):
        """Holland, P.M., Eaton, B.E., and Hanley, H.J.M., "A Correlation of the Viscosity and Thermal Conductivity Data of Gaseous and Liquid Ethylene," J. Phys. Chem. Ref. Data, 12(4):917-932, 1983."""
        GT = [-2.903423528e5, 4.680624952e5, -1.8954783215e5, -4.8262235392e3,
              2.243409372e4, -6.6206354818e3, 8.9937717078e2, -6.0559143718e1,
              1.6370306422]
        lo = 0
        for i in range(-3, 6):
            lo += GT[i+3]*self.T**(i/3.)

        tita = (self.rho.gcc-0.221)/0.221
        j = [0, -1.304503323e1, 1.8214616599e1, 9.903022496e3, 7.420521631e2,
             -3.0083271933e-1, 9.6456068829e1, 1.350256962e4]
        l1 = exp(j[1]+j[4]/self.T)*(exp(self.rho.gcc**0.1*(j[2]+j[3]/self.T**1.5)+tita*self.rho.gcc**0.5*(j[5]+j[6]/self.T+j[7]/self.T**2))-1.)

        lc = 0
        # FIXME: no sale
#        deltarho=(self.rho/self.M-0.221)/0.221
#        deltaT=(self.T-282.34)/282.34
#        xkt=(1.0/self.rho/self.M/self.derivative("P", "rho", "T")*1e3)**0.5
#        b=abs(deltarho)/abs(deltaT)**1.19
#        xts=(self.rho/self.M)**2*xkt*5.039/.221**2
#        g=xts*abs(deltaT)**1.19
#        xi=0.69/(b**2*5.039/g/Boltzmann/282.34)**0.5
#        f=exp(-18.66*deltaT**2-4.25*deltarho**4)
#        c=(self.M/self.rho.gcc/Avogadro/Boltzmann/self.T)**0.5
#        lc=c*Boltzmann*self.T**2/6.0/pi/self.mu.muPas/xi*self.dpdT**2*self.kappa**0.5*f
#        print lo, l1
        return unidades.ThermalConductivity(lo+l1+lc, "mWmK")

