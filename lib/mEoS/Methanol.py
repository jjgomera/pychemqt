#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from scipy import exp
from scipy import __version__
from scipy.constants import pi, Avogadro
from scipy.constants import Boltzmann

from lib.meos import MEoS
from lib import unidades


class Methanol(MEoS):
    """Multiparameter equation of state for methanol"""
    name = "methanol"
    CASNumber = "67-56-1"
    formula = "CH3OH"
    synonym = ""
    rhoc = unidades.Density(275.5626)
    Tc = unidades.Temperature(512.6)
    Pc = unidades.Pressure(8103.5, "kPa")
    M = 32.04216  # g/mol
    Tt = unidades.Temperature(175.61)
    Tb = unidades.Temperature(337.632)
    f_acent = 0.5625
    momentoDipolar = unidades.DipoleMoment(1.7, "Debye")
    id = 117

    CP1 = {"ao": 3.9007912,
           "an": [], "pow": [],
           "ao_exp": [0.10992677e2, 0.18336830e2, -0.16366004e2, -0.62332348e1,
                      0.28035363e1, 0.10778099e1, 0.96965697],
           "exp": [2115.01542, 1676.18569, 1935.16717, 1504.97016, 4222.83691,
                   5296.17127, 273.36934],
           "ao_hyp": [], "hyp": []}

    CP2 = {"ao": 0.964220/8.3143*32.,
           "an": [0.532325e-4/8.3143*32., 0.672819e-5/8.3143*32.,
                  -0.768411e-8/8.3143*32., 0.275220e-11/8.3143*32.],
           "pow": [1, 2, 3, 4],
           "ao_exp": [], "exp": [],
           "ao_hyp": [], "hyp": []}

    reuck = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of de Reuck and Craven (1993)",
        "__doi__": {"autor": "de Reuck, K.M. and Craven, R.J.B.",
                    "title": "Methanol, International Thermodynamic Tables of the Fluid State - 12",
                    "ref": "IUPAC, Blackwell Scientific Publications, London, 1993.",
                    "doi": ""},
        "R": 8.31448,
        "cp": CP1,
        "ref": "NBP",

        "Tref": 513.38, "rhoref": 8.78517*M,

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 35.57,
        "Pmin": 0.18629e-3, "rhomin": 28.23,

        "nr1": [-0.280062505988e1, 0.125636372418e2, -0.130310563173e2,
                0.326593134060e1, -0.411425343805e1, 0.346397741254e1,
                -0.836443967590e-1, -0.369240098923, 0.313180842152e-2,
                0.603201474111, -0.231158593638, 0.106114844945,
                -0.792228164995e-1, -0.422419150975e-4, 0.758196739214e-2,
                -0.244617434701e-4, 0.115080328802e-5],
        "d1": [1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 6, 7],
        "t1": [0, 1, 2, 3, 1, 2, 3, 4, 6, 0, 3, 4, 0, 7, 1, 6, 7],

        "nr2": [-0.125099747447e2, 0.270392835391e2, -0.212070717086e2,
                0.632799472270e1, 0.143687921636e2, -0.287450766617e2,
                0.185397216068e2, -0.388720372879e1, -0.416602487963e1,
                0.529665875982e1, 0.509360272812, -0.330257604839e1,
                -0.311045210826, 0.273460830583, 0.518916583979,
                -0.227570803104e-2, 0.211658196182e-1, -0.114335123221e-1,
                0.249860798459e-2],
        "d2": [1, 1, 1, 1, 2, 2, 2, 2, 3, 4, 5, 5, 5, 5, 6, 9, 6, 6, 4],
        "t2": [1, 2, 3, 4, 1, 2, 3, 5, 1, 2, 1, 2, 4, 5, 2, 5, 9, 14, 19],
        "c2": [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 6],
        "gamma2": [1.01733510223052]*16+[1.03497071023039]*2+[1.05291203329783],

        "nr3": [-0.819291988442e1, 0.478601004557, -0.444161392885,
                0.179621810410, -0.687602278259, 0.240459848295e1,
                -0.688463987466e1, 0.113992982501e1],
        "d3": [1, 1, 1, 1, 1, 3, 3, 3],
        "t3": [0]*8,
        "alfa3": [4.06934040892209, 8.20892015621185, 9.15601592007471,
                  83.8326275286616, 16.2773616356884, 27.705105527215,
                  16.2773616356884, 264.95250181898],
        "beta3": [-3.8940745646517, -3.8940745646517, -3.8940745646517,
                  -3.8940745646517, -3.8940745646517, -23.0649031906293,
                  -23.0649031906293, -23.0649031906293],
        "gamma3": [1.54080254509371, 1.54080254509371, 1.54080254509371,
                   1.54080254509371, 1.54080254509371, 1.08389789427588,
                   1.08389789427588, 1.08389789427588],
        "epsilon3": [0]*8,
        "exp1": [2, 3, 2, 4, 2, 3, 2, 4],
        "exp2": [1]*8}

    sun = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of Sun and Ely "
                    "(2004)",
        "__doi__": {"autor": "Sun, L. and Ely, J.F.",
                    "title": "Universal equation of state for engineering "
                             "application: Algorithm and  application to "
                             "non-polar and polar fluids",
                    "ref": "Fluid Phase Equilib., 222-223 (2004) 107-118",
                    "doi": "10.1016/j.fluid.2004.06.028"},

        "R": 8.3143,
        "cp": CP1,
        "ref": "NBP",

        "Tmin": Tt, "Tmax": 620.0, "Pmax": 800000.0, "rhomax": 40.,
        "Pmin": 0.1, "rhomin": 40.,

        "nr1": [-2.4578394, 1.39060027, 8.56114069e-1, -4.20843418e-2,
                3.63682442e-5, 7.05598662e-1],
        "d1": [1, 1, 1, 3, 7, 2],
        "t1": [1.5, 0.25, 1.25, 0.25, 0.875, 1.375],

        "nr2": [3.70573369e-1, 2.46303468, 1.50253790, 7.47553687e-2,
                -3.06417876e-1, -7.48402758e-1, -1.01432849e-1, 8.06830693e-2],
        "d2": [1, 1, 2, 5, 1, 1, 4, 2],
        "t2": [0, 2.375, 2., 2.125, 3.5, 6.5, 4.75, 12.5],
        "c2": [1, 1, 1, 1, 2, 2, 2, 3],
        "gamma2": [1]*8}

    polt = {
        "__type__": "Helmholtz",
        "__name__": "Helmholtz equation of state for methanol of Polt et al. (1992)",
        "__doi__": {"autor": "Polt, A., Platzer, B., and Maurer, G.",
                    "title": "Parameter der thermischen Zustandsgleichung von Bender fuer 14 mehratomige reine Stoffe",
                    "ref": "Chem. Technik 22(1992)6 , 216/224",
                    "doi": ""},
        "R": 8.3143,
        "cp": CP2,
        "ref": "NBP",

        "Tmin": 298., "Tmax": 703.0, "Pmax": 63000.0, "rhomax": 26.0625,
        "Pmin": 16.803, "rhomin": 24.576,

        "nr1": [-0.412043979985e1, 0.541210456547e1, -0.974639417666,
                -0.909437999343, -0.143467597275, 0.557052459597e1,
                -0.697445416557e1, 0.860535902136, 0.244117735035e1,
                -0.449073510921e1, 0.223855290012e1, -0.71733653794, 0.876135006507,
                0.151777405466, -0.233178058896, 0.140022534721e-1],
        "d1": [0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5],
        "t1": [3, 4, 5, 0, 1, 2, 3, 4, 0, 1, 2, 0, 1, 0, 1, 1],

        "nr2": [0.412043979985e1, -0.541210456547e1, 0.974639417666,
                -0.4642672133, 0.944015617353, -0.449348200461],
        "d2": [0, 0, 0, 2, 2, 2],
        "t2": [3, 4, 5, 3, 4, 5],
        "c2": [2]*6,
        "gamma2": [0.591872]*6}

    eq = reuck, sun, polt

    _surface = {"sigma": [0.22421, -0.21408, 0.083233],
                "exp": [1.3355, 1.677, 4.4402]}
    _melting = {"eq": 1, "Tref": Tt, "Pref": 0.187e-3,
                "Tmin": Tt, "Tmax": 620.,
                "a1": [1], "exp1": [0],
                "a2": [5.320770e9, 4.524780e9, 3.888861e10], "exp2": [1, 1.5, 4],
                "a3": [], "exp3": []}
    _vapor_Pressure = {
        "eq": 5,
        "ao": [-0.87414e1, 0.15035e1, -0.28720e1, -0.51345],
        "exp": [1., 1.5, 2.5, 5.]}
    _liquid_Density = {
        "eq": 1,
        "ao": [0.60230e-1, 0.18855e2, -0.27626e2, 0.11213e2, 0.69039],
        "exp": [0.1, 0.65, 0.79, 0.95, 4.4]}
    _vapor_Density = {
        "eq": 3,
        "ao": [-0.81104, -0.55661e1, -0.79326e3, 0.19234e4, -0.29219e4,
               0.29660e4, -0.13210e4],
        "exp": [0.25, 0.6, 3.5, 4.0, 5.0, 6.0, 7.0]}

    visco0 = {"eq": 0,
              "method": "_visco0",
              "__name__": "Xiang (2006)",
              "__doi__": {"autor": "Xiang, H.W., Huber, M.L. and Laesecke, A.",
                          "title": "A New Reference Correlation for the Viscosity of Methanol",
                          "ref": "J. Phys. Chem. Ref. Data 35, 1597 (2006)",
                          "doi": "10.1063/1.2360605"},
              "__test__":
                  # Table 5, Pag 15
                  """
                  >>> st=Methanol(T=175.63, x=0.5)
                  >>> print "%0.2f %0.3g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  175.63 1.86e-7 4.09e-6 0.005822 904.56 12.80
                  >>> st=Methanol(T=200, x=0.5)
                  >>> print "%0.0f %0.4g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  200 6.098e-6 1.1754e-4 0.006563 880.28 4.506
                  >>> st=Methanol(T=250, x=0.5)
                  >>> print "%0.0f %0.4g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  250 0.0008103 0.012577 0.008112 831.52 1.236
                  >>> st=Methanol(T=300, x=0.5)
                  >>> print "%0.0f %0.5g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  300 0.018682 0.24623 0.009678 784.51 0.5291
                  >>> st=Methanol(T=350, x=0.5)
                  >>> print "%0.0f %0.5g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  350 0.16172 1.9053 0.01118 735.84 0.2838
                  >>> st=Methanol(T=400, x=0.5)
                  >>> print "%0.0f %0.5g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  400 0.77374 8.7343 0.01251 678.59 0.1714
                  >>> st=Methanol(T=450, x=0.5)
                  >>> print "%0.0f %0.5g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  450 2.5433 30.831 0.01388 600.49 0.1058
                  >>> st=Methanol(T=500, x=0.5)
                  >>> print "%0.0f %0.5g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  500 6.5250 109.88 0.01891 451.53 0.05748
                  >>> st=Methanol(T=512, x=0.5)
                  >>> print "%0.0f %0.5g %0.3g %0.4g %0.5g %0.4g" % (\
                      st.T, st.P.MPa, st.Liquido.rho, st.Liquido.mu.muPas, st.Gas.rho, st.Gas.mu.muPas)
                  512 8.0195 202.99 0.02838 341.17 0.04174
                  """
                  # Table 6, Pag 16
                  """
                  >>> st=Methanol(T=180, P=1e4)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  0.01 900.27 10.44
                  >>> st=Methanol(T=200, P=1e5)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  0.10 880.34 4.510
                  >>> st=Methanol(T=220, P=4e5)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  0.40 860.76 2.460
                  >>> st=Methanol(T=280, P=1e6)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  1.00 804.12 0.7228
                  >>> st=Methanol(T=300, P=1e4)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  0.01 0.12955 0.009696
                  >>> st=Methanol(T=400, P=1e8)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  100 784.82 0.2846
                  >>> st=Methanol(T=500, P=8e8)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  800 954.52 0.3908
                  >>> st=Methanol(T=600, P=3e6)
                  >>> print "%0.2f %0.5g %0.4g" % (st.P.MPa, st.rho, st.mu.muPas)
                  3.00 20.717 0.01971
                  """}

    def _visco0(self, rho, T, fase):
        # FIXME: No sale
        rhoc = 273.
        ek = 577.87
        sigma0 = 0.3408e-9
        delta = 0.4575
        sigmac = 0.7193422e-9
        a = [1.16145, -0.14874, 0.52487, -0.77320, 2.16178, -2.43787,
             0.95976e-3, 0.10225, -0.97346, 0.10657, -0.34528, -0.44557, -2.58055]
        b = [-19.572881, 219.73999, -1015.3226, 2471.0125, -3375.1717,
             2491.6597, -787.26086, 14.085455, -0.34664158]
        d = [-1.181909, 0.5031030, -0.6268461, 0.5169312, -0.2351349,
             5.3980235e-2, -4.9069617e-3]
        e = [0, 4.018368, -4.239180, 2.245110, -0.5750698, 2.3021026e-2,
             2.5696775e-2, -6.8372749e-3, 7.2707189e-4, -2.9255711e-5]

        T_ = self.T/ek
        OmegaLJ = a[0]*T_**a[1]+a[2]*exp(a[3]*T_)+a[4]*exp(a[5]*T_)
        OmegaD = a[7]*T_**a[8]+a[9]*exp(a[10]*T_)+a[11]*exp(a[12]*T_)
        OmegaSM = OmegaLJ*(1+delta**2*OmegaD/(1+a[6]*delta**6))
        no = 5.*(self.M/Avogadro*Boltzmann*self.T/pi)**0.5/(16*sigma0**2*OmegaSM)
        B = (sum([b[i]/T_**(0.25*i) for i in range(7)])+b[7]/T_**2.5+b[8]/T_**5.5)*Avogadro*sigma0**3
        C = 1.86222085e-3*T_**3*exp(9.990338/T_**0.5)*(Avogadro*sigma0**3)**2
        ng = 1+B*rho/self.M*1000+C*(rho/self.M*1000)**2

        Tr = self.T/self.Tc
        rhor = rho/rhoc
        sigmaHS = sigmac*(sum([d[i]/Tr**i for i in range(7)])+sum([e[i]*rhor**(i) for i in range(1, 10)]))
        b = 2*pi*Avogadro*sigmaHS**3/3
        Xi = b*rho/self.M*1000/4
        g = (1-0.5*Xi)/(1-Xi)**3
        ne = 1./g+0.8*b*rho/self.M*1000+0.761*g*sigmaHS*b**2*(rho/self.M*1000)**2

        f = 1/(1+exp(5*(rhor-1)))
        n = no * (f*ng + (1-f)*ne)
        return unidades.Viscosity(n)

    _viscosity = visco0,

    thermo0 = {"eq": 1,
               "__name__": "Perkins (2002)",
               "__doi__": {"autor": "",
                           "title": "",
                           "ref": "unpublished preliminary correlation, NIST, MLH, Aug. 2006",
                           "doi": ""},

               "Tref": 1., "kref": 1,
               "no": [5.7992e-7],
               "co": [1.7862],

               "Trefb": 513.38, "rhorefb": 8.78517, "krefb": 1.,
               "nb": [0.405435, -0.293791, -0.289002, 0.226890, 0.579019e-1,
                      -0.399576e-1],
               "tb": [0, 1]*3,
               "db": [1, 1, 2, 2, 3, 3],
               "cb": [0]*6,

               "critical": 3,
               "gnu": 0.63, "gamma": 1.239, "R0": 1.03,
               "Xio": 0.194e-9, "gam0": 0.0496, "qd": 0.342e-9, "Tcref": 768.9}

    _thermal = thermo0,
