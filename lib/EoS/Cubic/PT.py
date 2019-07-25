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


from csv import reader
import os

from scipy import roots
from scipy.constants import R

from lib.EoS.cubic import Cubic, CubicHelmholtz


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

    The paper give generation correlation for F and ζc

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
    __doi__ = {
        "autor": "Patel, N.C., Teja, A.S.",
        "title": "A New Cubic Equation of State for Fluids and Fluid Mixtures",
        "ref": "Chem. Eng. Sci. 37(3) (1982) 463-473",
        "doi": "10.1016/0009-2509(82)80099-7"},

    def __init__(self, T, P, mezcla):
        """Initialization procedure
        
        Parameters
        ----------

        """

        self.T = T
        self.P = P
        self.mezcla = mezcla

        bi = []
        ci = []
        ai = []
        for cmp in mezcla.componente:
            # Generalization Eq 20-21
            f = 0.452413 + 1.30982*cmp.f_acent - 0.295937*cmp.f_acent**2
            xic = 0.329032 - 0.076799*cmp.f_acent + 0.0211947*cmp.f_acent**2

            # Eq 8
            ci.append((1-3*xic)*R*cmp.Tc/cmp.Pc)

            # Eq 10
            b = roots([1, 2-3*xic, 3*xic**2, -xic**3])
            Bpositivos=[]
            for i in b:
                if i>0:
                    Bpositivos.append(i)
            Omegab = min(Bpositivos)
            bi.append(Omegab*R*cmp.Tc/cmp.Pc)

            # Eq 9
            Omegaa = 3*xic**2 + 3*(1-2*xic)*Omegab + Omegab**2 + 1 - 3*xic

            alfa = (1+f*(1-(T/cmp.Tc)**0.5))**2
            ai.append(Omegaa*alfa*R**2*cmp.Tc**2/cmp.Pc)

        am, bm, cm = self._mixture(None, mezcla.ids, [ai, bi, ci])

        self.ai = ai
        self.bi = bi
        self.b = bm
        self.tita = am
        self.delta = bm+cm
        self.epsilon = -bm*cm

        # tdadt=0
        # for i in range(len(mezcla.componente)):
            # for j in range(len(mezcla.componente)):
                # tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])
        # self.dTitadT=tdadt
        # self.u=2
        # self.w=-1

        super(PT, self).__init__(T, P, mezcla)
