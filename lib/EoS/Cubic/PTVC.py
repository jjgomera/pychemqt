#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

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


from scipy.constants import R

from lib.EoS.Cubic.PT import PT


class PTVC(PT):
    r"""Patel-Teja-Valderrama-Cisternas cubic equation of state implementation
    as explain in [1]_

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

    __title__ = "Patel-Teja-Valderrama-Cisternas (1986)"
    __status__ = "PTVC"
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
