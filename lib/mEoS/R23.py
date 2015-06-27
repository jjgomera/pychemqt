#!/usr/bin/python
# -*- coding: utf-8 -*-

from scipy import exp

from lib.meos import MEoS
from lib import unidades


class R23(MEoS):
    """Multiparameter equation of state for R23"""
    name = "trifluoromethane"
    CASNumber = "75-46-7"
    formula = "CHF3"
    synonym = "R23"
    rhoc = unidades.Density(526.504)
    Tc = unidades.Temperature(299.293)
    Pc = unidades.Pressure(4832.0, "kPa")
    M = 70.01385  # g/mol
    Tt = unidades.Temperature(118.02)
    Tb = unidades.Temperature(191.132)
    f_acent = 0.263
    momentoDipolar = unidades.DipoleMoment(1.649, "Debye")
    id = 643

    CP1 = {"ao": 3.999,
           "an": [], "pow": [],
           "ao_exp": [2.371, 3.237, 2.61, 0.8274],
           "exp": [744, 1459, 2135, 4911],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 3.999509244,
           "an": [], "pow": [],
           "ao_exp": [1.070326018, 1.566866769, 0.848051597, 1.847243699,
                      1.649657530, 2.043965290],
           "exp": [4368.102594, 1607.104940, 1007.138279, 1973.991027,
                   1657.461854, 729.455868],
           "ao_hyp": [], "hyp": []}

    CP3 = {"ao": 4.0101431,
           "an": [-.55274742e-2, .74008258e-4, -.12590943e-6, .69472178e-10],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    helmholtz1 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-23 of Penoncello et al. (2003)",
        "__doi__": {"autor": "Penoncello, S.G., Lemmon, E.W., Jacobsen, R.T, Shan, Z.",
                    "title": "A Fundamental Equation for Trifluoromethane (R-23)", 
                    "ref": "J. Phys. Chem. Ref. Data 32, 1473 (2003).",
                    "doi":  "10.1063/1.1559671"}, 
        # TODO: Find paper to search test
        
        "R": 8.314472,
        "cp": CP1,
        "ref": "IIR", 
        
        "Tmin": Tt, "Tmax": 475.0, "Pmax": 120000.0, "rhomax": 24.31, 
        "Pmin": 0.058, "rhomin": 24.31, 

        "nr1": [.7041529e1, -.8259512e1, .805304e-2, -.8617615e-1, .633341e-2],
        "d1": [1, 1, 1, 2, 5],
        "t1": [0.744, 0.94, 4.3, 1.46, 0.68],

        "nr2": [-0.1863285, 0.3280510, 0.5191023, 0.6916144e-1, -0.5045875e-2,
                -0.1744221e-1, -0.5003972e-1, 0.4729813e-1, -0.6164031e-1,
                0.1583585e-1, -0.1795790e-2, -0.1099007e-2],
        "d2": [1, 2, 3, 5, 1, 2, 2, 4, 4, 4, 2, 2],
        "t2": [4.8, 1.5, 2.07, 0.09, 9.6, 0.19, 11.2, 0.27, 1.6, 10.3, 14., 15.],
        "c2": [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4],
        "gamma2": [1]*12}

    helmholtz2 = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for R-23 of Penoncello et al. (2000)",
        "__doi__": {"autor": "Penoncello, S.G., Shan, Z., and Jacobsen, R.T.",
                    "title": "A fundamental equation for the calculation of the thermodynamic properties of trifluoromethane (R23)", 
                    "ref": "ASHRAE Trans. 106(Part 1), 2000.",
                    "doi":  ""}, 

        "R": 8.31451,
        "cp": CP2,
        "ref": "IIR", 
        
        "Tmin": Tt, "Tmax": 473.15, "Pmax": 120000.0, "rhomax": 23.0, 
        "Pmin": 0.05888, "rhomin": 22.851535, 

        "nr1": [0.350093635099, -0.131185838025e1, -0.254118065769,
                .104275296122, -.205326997924, .256040993750, .118078220087e-1,
                0.532850915621e-3,  0.956700157221e-3, -0.118990410423e-5],
        "d1": [1, 1, 1, 2, 2, 2, 3, 4, 6, 8],
        "t1": [-0.14, 1.49, 2.41, 0.05, 1.59, 2.04, -0.27, 2.76, -0.06, 3.25],

        "nr2": [-0.180609172794, 0.138077199166, 0.507828500811e-1,
                0.439772083175e-1, -0.723557234469e-1, 0.256500006055e-2,
                0.263213487134e-1, 0.139266509424e-1, -0.105325247813e-1,
                0.136475671500e-2, -0.592653649931e-2, -0.644925101471e-1,
                -0.227635186710e-1, 0.122367812706, 0.318153208563e-1,
                0.146725272055e-1, -0.923639585566e-1],
        "d2": [1, 2, 2, 3, 3, 6, 6, 7, 7, 10, 2, 3, 4, 4, 5, 5, 5],
        "t2": [5.36, 5.28, 4.23, 3.35, 6.93, 8.48, 6.01, 3.34, 7.1, 5.46,
               16.06, 19.37, 10.81, 22.79, 34.95, 9.94, 29.16],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4],
        "gamma2": [1]*17}

    helmholtz3 = {
        "__type__": "Helmholtz",
        "__name__": "Bender equation of state for R-23 of Platzer et al. (1990).",
        "__doi__": {"autor": "Platzer, B., Polt, A., and Maurer, G.",
                    "title": "Thermophysical properties of refrigerants", 
                    "ref": "Berlin:  Springer-Verlag, 1990.",
                    "doi": ""}, 
        "R": 8.31451,
        "cp": CP3,
        "ref": "NBP", 
        
        "Tmin": 90.0, "Tmax": 475.0, "Pmax": 60000.0, "rhomax": 16.65, 
        "Pmin": 2.5664104, "rhomin": 22.851535, 

        "nr1": [-0.133234251368e1, 0.210373595421e1, -0.376198728030,
                0.881622087335, -0.272053790906e1, 0.247468024356e1,
                -0.234010064393e1, 0.303959507238, 0.317372750273e-1,
                .329392142221e-1, .20583853186, .133550139894, -.181698216766,
                -0.245123269882e-1, 0.247477874180e-1, 0.589916583383e-2],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.133234251368e1, -0.210373595421e1, 0.376198728030,
                0.574267667948, -0.762218931280, 0.472710395636e-1],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.70304082]*6}

    eq = helmholtz1, helmholtz2, helmholtz3

    _surface = {"sigma": [-0.32359, 0.37702], "exp": [1.6055, 1.5232]}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-7.2631,  1.3140, -0.78507, -3.1991],
        "exp": [1, 1.5, 2.4, 3.9]}
    _liquid_Density = {
        "eq": 1,
        "ao": [2.2636, 0.47007, 0.28660],
        "exp": [0.37, 0.94, 3.1]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-3.5136, -7.7491, -24.871, -65.637],
        "exp": [0.43, 1.4, 3.7, 8]}

    visco0 = {"eq": 0, "omega": 1,
              "method": "_visco0",
              "__name__": "Shan (2000)",
              "__doi__": {"autor": "Shan, Z., Penoncello, S.G., and Jacobsen, R.T.",
                          "title": "A Generalized Model for Viscosity and Thermal Conductivity of Trifluoromethane (R-23)", 
                          "ref": "ASHRAE Transactions, Volume 106:1-11, 2000.",
                          "doi": ""}, 

              "ek": 243.91, "sigma": 0.4278,
              "n_chapman": 0.2233755/M**0.5,
              "collision": [0.4425728, -0.5138403, 0.1547566, -0.02821844,
                            0.001578286]}

    _viscosity = visco0,

    def _visco0(self):
        rhol = 32.174
        C1 = 1.3163
        C2 = 0.1832
        deltaG = 771.23
        nmax = 3.967

        Drho = rhol-self.rho/self.M
        delta = self.rho/self.M-7.5114
        tau = self.T-299.28

        no = self._Visco0()
        ng = no*(Drho/rhol)**C1
        nr = (self.rho/self.M/rhol)**C1*C2*rhol**2/Drho*self.T**0.5*exp(self.rho/self.M/Drho*deltaG/self.R.kJkgK/self.M/self.T)
        nc = 4*nmax/(exp(delta)+exp(-delta))/(exp(tau)+exp(-tau))
        return unidades.Viscosity(ng+nr+nc, "muPas")

    thermo0 = {"eq": 0,
               "method": "_thermo0",
               "__name__": "Shan (2000)", 
               "__doi__": {"autor": "Shan, Z., Penoncello, S.G., and Jacobsen, R.T.",
                           "title": "A Generalized Model for Viscosity and Thermal Conductivity of Trifluoromethane (R-23)", 
                           "ref": "ASHRAE Transactions, Volume 106:1-11, 2000.",
                           "doi": ""}}

    _thermal = thermo0,

    def _thermo0(self):
        rhol = 68.345
        B1 = -2.5370
        B2 = 0.05366
        C1 = 0.94215
        C2 = 0.14914
        deltaG = 2508.58
        lmax = 25.

        Drho = rhol-self.rho/self.M
        delta = self.rho/self.M-7.5114
        tau = self.T-299.28

        lg = (B1+B2*self.T)*(Drho/rhol)**C1
        lr = (self.rho/self.M/rhol)**C1*C2*rhol**2/Drho*self.T**0.5*exp(self.rho/self.M/Drho*deltaG/self.R.kJkgK/self.M/self.T)
        lc = 4*lmax/(exp(delta)+exp(-delta))/(exp(tau)+exp(-tau))
        return unidades.ThermalConductivity(lg+lr+lc, "mWmK")
