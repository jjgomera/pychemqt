#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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


from scipy.constants import R

from lib.EoS.cubic import Cubic


class PRYuLu(Cubic):
    r"""Peng-Robinson equation of state modified by Yu-Lu

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+c\right)+b\left(3V+c\right)}\\
        a = \Omega_{ac}\frac{R^2T_c^2}{P_c}\alpha\\
        b = \Omega_{bc}\frac{RT_c}{P_c}\\
        c = \Omega_{cc}\frac{RT_c}{P_c}\\
        \Omega_{ac} = 0.468630-0.0379304\omega+0.00751969\omega^2\\
        \Omega_{bc} = 0.0892828-0.0340903\omega-0.00518289\omega^2\\
        \Omega_{cc} = w \Omega_{bc}
        c = wb\\
        w+3 = u = 1.70083+0.648463\omega+0.895926\omega^2\\
        \alpha = 10^{M\left(A_0+A_1T_r+A_2T_r"2\right)\left(1-T_r\right)}\\
        \end{array}

    The alpha parameter depend of acentric factor. For ω ≤ 0.49:

    .. math::
        \begin{array}[t]{l}
        M = 0.406846+1.87907\omega-0.792636\omega^2+0.737519\omega^3\\
        A_0 = 0.536843\\
        A_1 = -0.39244\\
        A_2 = 0.26507\\
        \end{array}

    For 0.49 < ω ≤ 1:

    .. math::
        \begin{array}[t]{l}
        M = 0.581981-0.171416\omega-1.84441\omega^2+1.19047\omega^3\\
        A_0 = 0.79355\\
        A_1 = -0.53409\\
        A_2 = 0.37273\\
        \end{array}
    """

    __title__ = "PR-Yu-Lu (1987)"
    __status__ = "PRYuLu"
    __doi__ = {
        "autor": "Yu, J.-M., Lu, B.C.-Y.",
        "title": "A Three-Parameter Cubic Equation of State for Asymmetric "
                 "Mixture Density Calculations",
        "ref": "Fluid Phase Equilibria 34 (1987) 1-19",
        "doi": "10.1016/0378-3812(87)85047-1"},

    def _cubicDefinition(self, T):
        """Definition of individual components coefficients"""
        ai = []
        bi = []
        ci = []
        for cmp in self.componente:
            a, b, c = self._lib(cmp, T)
            ai.append(a)
            bi.append(b)
            ci.append(c)

        self.ai = ai
        self.bi = bi
        self.ci = ci

    def _GEOS(self, xi):
        am, bm, cm = self._mixture(xi, [self.ai, self.bi, self.ci])

        delta = 3*bm + cm
        epsilon = bm*cm
        return am, bm, delta, epsilon

    def _lib(self, cmp, T):
        Tr = T/cmp.Tc

        # Eq 14
        Omegaa = 0.468630 - 0.0378304*cmp.f_acent + 0.00751969*cmp.f_acent**2

        # Eq 15
        Omegab = 0.0892828 - 0.0340903*cmp.f_acent - 0.00518289*cmp.f_acent**2

        # Eq 16
        u = 1.70083 + 0.648463*cmp.f_acent + 0.895926*cmp.f_acent**2
        w = u-3

        # Eq 13
        Omegac = w*Omegab

        if cmp.f_acent <= 0.49:
            # Eq 18
            M = 0.406846 + 1.87907*cmp.f_acent - 0.792636*cmp.f_acent**2 + \
                0.737519*cmp.f_acent**3
            A0 = 0.536843
            A1 = -0.39244
            A2 = 0.26507
        else:
            # Eq 19
            M = 0.581981 - 0.171416*cmp.f_acent + 1.84441*cmp.f_acent**2 - \
                1.19047*cmp.f_acent**3
            A0 = 0.79355
            A1 = -0.53409
            A2 = 0.37273

        if Tr > 1:
            alfa = 10**(M*(A0 + A1 + A2)*(1-Tr))
        else:
            # Eq 17
            alfa = 10**(M*(A0 + A1*Tr + A2*Tr**2)*(1-Tr))

        a = Omegaa*alfa*R**2*cmp.Tc**2/cmp.Pc                           # Eq 10
        b = Omegab*R*cmp.Tc/cmp.Pc                                      # Eq 11
        c = Omegac*R*cmp.Tc/cmp.Pc                                      # Eq 12

        return a, b, c


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PRYuLu(300, 9.9742e5, mix)
    print('%0.1f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PRYuLu(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
