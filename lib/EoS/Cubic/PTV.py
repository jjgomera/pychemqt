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


class PTV(PT):
    r"""Patel-Teja-Valderrama cubic equation of state implementation as explain
    in [1]_

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

    __title__ = "Patel-Teja-Valderrama (1990)"
    __status__ = "PTV"
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


# Aznar, M., Silva-Telles, A., Valderrama, J.O.
# Parameters for the Attractive Coefficient of the Patel-Teja-Valderrama
# Equation of State
# Chem. Eng. Comm. 190(10) (2003) 1411-1426
# 10.1080/00986440302156
# This paper give selected parameters for the three-parameter form of the
# temperature dependence of the attractive term applied to PTV cubic EoS. The
# paper give only values for 50 selected compounds, to get the complete list
# is necessary ask to author
