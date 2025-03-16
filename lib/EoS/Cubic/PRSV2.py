#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from lib.EoS.Cubic.PRSV import PRSV, dat


class PRSV2(PRSV):
    r"""Peng-Robinson cubic equation of state with a modified dependence of
    temperature by Stryjek-Vera v.2 [1]_

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + k\left(1-Tr^{0.5}\right)\\
        k = k_0+\left[k_1+k_2\left(k_3-T_r\right)\left(1-\sqrt{T_r}\right)
        \right]\left(1+\sqrt{T_r}\right)\left(0.7-T_r\right)\\
        k_0 = 0.378893+1.4897153\omega-0.17131848\omega^2+0.0196554\omega^3\\
        \end{array}

    :math:`k_1`, :math:`k_2` and :math:`k_3` are parameters characteristic
    compound specific
    """

    __title__ = "PRSV2 (1986)"
    __status__ = "PRSV2"
    __doi__ = {
        "autor": "Stryjek, R., Vera, J.H.",
        "title": "PRSV2: A Cubic Equation of State for Accurate Vapor—Liquid "
                 "Equilibria calculations",
        "ref": "Can. J. Chem. Eng. 64 (1986) 820–826",
        "doi": "10.1002/cjce.5450640516"},

    def _k(self, cmp, Tr):
        # Eq 11
        ko = 0.378893 + 1.4897153*cmp.f_acent - \
            0.17131848*cmp.f_acent**2 + 0.0196554*cmp.f_acent**3

        if cmp.id in dat and Tr < 1:
            k1, k2, k3 = dat[cmp.id]
        else:
            k1, k2, k3 = 0, 0, 0
        k = ko + (k1+k2*(k3-Tr)*(1+Tr**0.5))*(1+Tr**0.5)*(0.7-Tr)
        return k


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRSV2(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRSV2(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
