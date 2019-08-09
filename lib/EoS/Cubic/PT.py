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


from scipy import roots
from scipy.constants import R

from lib.EoS.cubic import Cubic, CubicHelmholtz


# Table I in [1]_ and Table III, IV in [3]_
dat = {
       98 : (0.328, 0.450751),
       46 : (0.329, 0.516798),
       47 : (0.327, 0.487035),
       2 : (0.324, 0.455336),
       3 : (0.317, 0.561567),
       22 : (0.313, 0.554369),
       4 : (0.317, 0.648049),
       23 : (0.324, 0.661305),
       65 : (0.310, 0.664179),
       6 : (0.309, 0.678389),
       5 : (0.315, 0.683133),
       24 : (0.315, 0.696423),
       8 : (0.308, 0.746470),
       7 : (0.314, 0.741095),
       10 : (0.305, 0.801605),
       11 : (0.305, 0.868856),
       12 : (0.301, 0.918544),
       13 : (0.301, 0.982750),
       14 : (0.297, 1.021919),
       15 : (0.297, 1.080416),
       16 : (0.294, 1.115585),
       17 : (0.295, 1.179982),
       18 : (0.291, 1.188785),
       21 : (0.283, 1.297054),
       90 : (0.276, 1.276058),
       92 : (0.277, 1.409671),
       49 : (0.309, 0.707727),
       48 : (0.328, 0.535060),
       51 : (0.310, 0.797391),
       50 : (0.320, 0.583165),
       62 : (0.269, 0.689803),
       38 : (0.303, 0.665434),
        # : (0.310, 0.859036),  Quinoline
       346 : (0.300, 1.000087),
       406 : (0.305, 1.082667),
       191 : (0.297, 0.827417),
       63 : (0.283, 0.642740),
       40 : (0.311, 0.698911),
       41 : (0.306, 0.753893),
       42 : (0.305, 0.812845),
       43 : (0.301, 0.816962),
       44 : (0.300, 0.807023),
       117 : (0.274, 0.965347),
       134 : (0.292, 1.171714),
       146 : (0.302, 1.211304),
       160 : (0.305, 1.221182),
       313 : (0.308, 1.240459),
       335 : (0.330, 1.433586),
       357 : (0.301, 1.215380),
       360 : (0.308, 1.270267),
       130 : (0.258, 0.762043),
       143 : (0.295, 1.146553),
       154 : (0.329, 1.395151),
       306 : (0.292, 1.174746),
       510 : (0.291, 1.272986),
       540 : (0.292, 1.393678),
       545 : (0.290, 1.496554),
       140 : (0.283, 0.701112),
       162 : (0.308, 0.787322),
       100 : (0.314, 0.694866),
       155 : (0.296, 0.842965),
       166 : (0.295, 0.882502),
       165 : (0.294, 0.826046)}


# Table I in [2]_
PT2 = {
    135: [-1.11765, -1.81779, 0.47892, 3],
    129: [0.82082, -2.80514, 0, 0],
    62: [0.60462, -2.56713, 0, 0],
    22: [1.94572, -3.59956, -0.37410, 2],
    49: [0.63199, -2.69935, 0, 0],
    50: [0.66433, -2.39792, -0.00669, 10],
    2: [0.32274, -1.47606, -0.02025, 6],
    140: [0.19454, -1.45357, 0.32485, -0.5],
    46: [0.09339, -1.26573, 0, 0],
    1: [-0.72258, 1.08363, -1.4928e-6, -8]}

# TODO: Add compound specific parameters from
# Georgeton, G.K., Smith, R.L. Jr.. and Teja, A.S: pp. 434-451, 1985. in
# Chao, K.C. and Robinson, R.L. (editors)
# Equations of State. Theories and Applications
# ACS Svmposium series Nº 300


class PT(Cubic):
    r"""Patel-Teja cubic equation of state implementation
    
    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+c\left(V-b\right)}\\
        a = \Omega_a\frac{R^2T_c^2}{P_c}\alpha\\
        b = \Omega_b\frac{RT_c}{P_c}\\
        c = \Omega_c\frac{RT_c}{P_c}\\
        \Omega_c = 1 - 3\zeta_c\\
        \Omega_a = 3\zeta_c^2 + 3\left(1-2\zeta_c\right)\Omega_b + \Omega_b^2
        + 1 - 3\zeta_c\\
        \alpha^{0.5} = 1 + F\left(1-Tr^{0.5}\right)\\
        \end{array}

    :math:`\Omega_b` is the smallest positive root or the equation::

    .. math::
        \Omega_b^3 + \left(2-3\zeta_c\right)\Omega_b^2 + 3\zeta_c^2\Omega_b -
        \zeta_c^3 = 0

    The paper give generation correlation for F and ζc, valid only for nonpolar
    compounds.

    .. math::
        \begin{array}[t]{l}
        F = 0.452413 + 1.30982\omega - 0.295937\omega^2\\
        \zeta_c = 0.329032 - 0.076799\omega + 0.0211947\omega^2\\
        \end{array}
    
    This library implement the generalized correlation using this correlation.
    The paper give values for these parameters for several compounds.

    Examples
    --------
    Example 4.3 from 1_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = PT(300, 9.9742e5, mix)
    >>> '%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol)
    '2050 91.4'
    >>> eq = PT(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '88.5'
    """

    __title__="Patel-Teja (1982)"
    __status__="PT"
    __doi__ = (
            {
        "autor": "Patel, N.C., Teja, A.S.",
        "title": "A New Cubic Equation of State for Fluids and Fluid Mixtures",
        "ref": "Chem. Eng. Sci. 37(3) (1982) 463-473",
        "doi": "10.1016/0009-2509(82)80099-7"},
            {
        "autor": "Patel, N.C.",
        "title": "Improvements of the Patel-Teja Equation of State",
        "ref": "Int. J. Thermophysics 17(3) (1996) 673-682",
        "doi": "10.1007/bf01441513"},
            {
        "autor": "Georgeton, G.K., Smith, R.L.Jr., Teja, A.S",
        "title": "Application of Cubic Equations of State to Polar Fluids "
                 "and Fluid Mixtures",
        "ref": "in Chao, K.C., Robinson, R.L. Equations of State. Theories "
               "and Applications, 1985, ACS Svmposium 300, pp. 434-451",
        "doi": ""})

    def __init__(self, T, P, mezcla):
        """Initialization procedure
        
        Parameters
        ----------
        
        The library implement too the enhancement alfa function for several
        compounds given in [2]_

        """

        self.T = T
        self.P = P
        self.mezcla = mezcla

        ai = []
        bi = []
        ci = []
        for cmp in mezcla.componente:
            a, b, c = self._lib(cmp, T)
            ai.append(a)
            bi.append(b)
            ci.append(c)

        am, bm, cm = self._mixture(None, mezcla.ids, [ai, bi, ci])

        self.ai = ai
        self.bi = bi
        self.b = bm
        self.tita = am
        self.delta = bm+cm
        self.epsilon = -bm*cm

        super(PT, self).__init__(T, P, mezcla)

    def _lib(self, cmp, T):
        if cmp.id in dat:
            # Use the compound specific parameters values
            xic, f = dat[cmp.id]
        else:
            # Use the generalization correlations, Eq 20-21
            f = 0.452413 + 1.30982*cmp.f_acent - 0.295937*cmp.f_acent**2
            xic = 0.329032 - 0.076799*cmp.f_acent + 0.0211947*cmp.f_acent**2

        # Eq 8
        c = (1-3*xic)*R*cmp.Tc/cmp.Pc

        # Eq 10
        b = roots([1, 2-3*xic, 3*xic**2, -xic**3])
        Bpositivos=[]
        for i in b:
            if i>0:
                Bpositivos.append(i)
        Omegab = min(Bpositivos)
        b = Omegab*R*cmp.Tc/cmp.Pc

        # Eq 9
        Omegaa = 3*xic**2 + 3*(1-2*xic)*Omegab + Omegab**2 + 1 - 3*xic

        if cmp.id in PT2:
            # Using improved alpha correlation from [2]_
            c1, c2, c3, n = PT2[cmp.id]
            alfa = 1 + c1*(T/cmp.Tc-1) + c2*((T/cmp.Tc)**0.5-1) + \
                c3*((T/cmp.Tc)**n-1)
        else:
            alfa = (1+f*(1-(T/cmp.Tc)**0.5))**2
        a = Omegaa*alfa*R**2*cmp.Tc**2/cmp.Pc

        return a, b, c


class PTV(PT):
    r"""Patel-Teja-Valderrama cubic equation of state implementation

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+c\left(V-b\right)}\\
        a = \Omega_a\frac{R^2T_c^2}{P_c}\alpha\\
        b = \Omega_b\frac{RT_c}{P_c}\\
        c = \Omega_c\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + F\left(1-Tr^{0.5}\right)\\
        \Omega_a = 0.66121 - 0.76105 Z_c\\
        \Omega_b = 0.02207 + 0.20868 Z_c\\
        \Omega_c = 0.57765 - 1.87080 Z_c\\
        F = 0.46283 + 3.5823\left(\omega·Z_c\right)
        + 8.19417\left(\omega·Z_c\right)^2\\
        \end{array}
    """

    __title__="Patel-Teja-Valderrama (1990)"
    __status__="PTV"
    __doi__ = {
        "autor": "Valderrama, J.O.",
        "title": "A Generalized Patel-Teja Equation of Stte for Polar and "
                 "Nonpolar Fluids and their Mixtures",
        "ref": "J. Chem. Eng. Jap. 23(1) (1990) 87-91",
        "doi": "10.1252/jcej.23.87"},

    def _lib(self, cmp, T):
        if cmp.Tc != 0 and cmp.Pc != 0 and cmp.Vc != 0:
            Zc = cmp.Pc.kPa*cmp.Vc*cmp.M/R/cmp.Tc
        else:
            Zc = 0.329032 - 0.076799*cmp.f_acent + 0.0211947*cmp.f_acent**2

        # Eq 5
        F = 0.46283 + 3.5823*cmp.f_acent*Zc + 8.19417*cmp.f_acent**2*Zc**2

        # Eq 4
        Omegaa = 0.66121 - 0.76105*Zc
        Omegab = 0.02207 + 0.20868*Zc
        Omegac = 0.57765 - 1.8708*Zc

        # Eq 3
        alfa = (1+F*(1-(T/cmp.Tc)**0.5))**2

        # Eq 2
        a = Omegaa*R**2*cmp.Tc**2/cmp.Pc
        b = Omegab*R*cmp.Tc/cmp.Pc
        c = Omegac*R*cmp.Tc/cmp.Pc

        return a*alfa, b, c


class PTVC(PT):
    r"""Patel-Teja-Valderrama-Cisternas cubic equation of state implementation

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+c\left(V-b\right)}\\
        a = \Omega_a\frac{R^2T_c^2}{P_c}\alpha\\
        b = \Omega_b\frac{RT_c}{P_c}\\
        c = \Omega_c\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + F\left(1-Tr^{0.5}\right)\\
        F = -6.608 + 70.43Z_c - 159.0Z_c^2\\
        \Omega_a = 0.69368018 - 1.0634424 Z_c + 0.68289995 Z_c^2 -
            0.21044403 Z_c^3 + 0.003752658 Z_c^4\\
        \Omega_b = 0.025987178 + 0.180754784 Z_c + 0.061258949 Z_c^2\\
        \Omega_c = 0.577500514 - 1.898414283 Z_c\\
        \end{array}
    """

    __title__="Patel-Teja-Valderrama-Cisternas (1986)"
    __status__="PTVC"
    __doi__ = {
        "autor": "Valderrama, J.O., Cisternas, L.A.",
        "title": "A Cubic Equation of State for Polar and Other Complex "
                 "Mixtures",
        "ref": "Fluid Phase Equilibria 29 (1986) 431-438",
        "doi": "10.1016/0378-3812(86)85041-5"},

    def _lib(self, cmp, T):
        if cmp.Tc != 0 and cmp.Pc != 0 and cmp.Vc != 0:
            Zc = cmp.Pc.kPa*cmp.Vc*cmp.M/R/cmp.Tc
        else:
            Zc = 0.329032 - 0.076799*cmp.f_acent + 0.0211947*cmp.f_acent**2

        F = -6.608 + 70.43*Zc - 159.0*Zc**2                             # Eq 6
        alfa = (1+F*(1-(T/cmp.Tc)**0.5))**2                             # Eq 5

        # Eq 4
        Omegaa = 0.69368018 - 1.0634424*Zc + 0.68289995*Zc**2 - \
            0.21044403*Zc**3 + 0.003752658*Zc**4
        Omegab = 0.025987178 + 0.180754784*Zc + 0.061258949*Zc**2
        Omegac = 0.577500514 - 1.898414283*Zc

        # Eq 3
        a = Omegaa*R**2*cmp.Tc**2/cmp.Pc
        b = Omegab*R*cmp.Tc/cmp.Pc
        c = Omegac*R*cmp.Tc/cmp.Pc

        return a*alfa, b, c
