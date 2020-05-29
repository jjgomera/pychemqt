#!/usr/bin/python3
# -*- coding: utf-8 -*-

r"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


from scipy.constants import R

from lib.EoS.cubic import Cubic


class ALS1983(Cubic):
    r"""Adachi modification to SRK cubic equation of state [1]_

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b_1}-\frac{a}{\left(V-b_2\right)\left(V+b_3\right)}\\
        a = A\frac{R^2T_c^2}{P_c}\alpha\\
        b_1 = B_1\frac{RT_c}{P_c}\\
        b_2 = B_2\frac{RT_c}{P_c}\\
        b_3 = B_3\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + m\left(1-\sqrt{T_r}\right)\\
        m = 0.4070 + 1.3787\omega - 0.2933\omega^2\\
        \end{array}

    The paper give generation correlation for A, B₁, B₂ and B₃

    .. math::
        \begin{array}[t]{l}
        B_1 = 0.08974 - 0.03452\omega + 0.0033\omega^2\\
        B_2 = 0.03686 + 0.00405\omega - 0.01073\omega^2 + 0.00157\omega^3\\
        B_3 = 0.154 + 0.14122\omega - 0.00272\omega^2 - 0.00484\omega^3\\
        A = 0.44869 + 0.04024\omega + 0.01111\omega^2 - 0.00576\omega^3\\
        \end{array}
    """

    __title__ = "Adachi (1983)"
    __status__ = "Adachi83"
    __doi__ = {
        "autor": "Adachi, Y., Lu, B.C.-Y., Sugie, H.",
        "title": "A Four-Parameter Equation of State",
        "ref": "Fluid Phase Equilibria 11 (1983) 29-48",
        "doi": "10.1016/0378-3812(83)85004-3"},

    def _cubicDefinition(self, T):
        """Definition of individual components coefficients"""

        ai = []
        b1i = []
        b2i = []
        b3i = []
        for cmp in self.componente:

            B1 = 0.08974 - 0.03452*cmp.f_acent + 0.0033*cmp.f_acent**2
            B2 = 0.03686 + 0.00405*cmp.f_acent - 0.01073*cmp.f_acent**2 + \
                0.00157*cmp.f_acent**3
            B3 = 0.154 + 0.14122*cmp.f_acent - 0.00272*cmp.f_acent**2 - \
                0.00484*cmp.f_acent**3
            A = 0.44869 + 0.04024*cmp.f_acent + 0.01111*cmp.f_acent**2 - \
                0.00576*cmp.f_acent**3

            m = 0.4070 + 1.3787*cmp.f_acent - 0.2933*cmp.f_acent**2
            alfa = 1+m*(1-(T/cmp.Tc)**0.5)**2

            ai.append(A*alfa*R**2*cmp.Tc**2/cmp.Pc)
            b1i.append(B1*R*cmp.Tc/cmp.Pc)
            b2i.append(B2*R*cmp.Tc/cmp.Pc)
            b3i.append(B3*R*cmp.Tc/cmp.Pc)

        self.ai = ai
        self.bi = b1i
        self.b2i = b2i
        self.b3i = b3i

    def _GEOS(self, xi):
        coef = [self.ai, self.bi, self.b2i, self.b3i]
        am, b1m, b2m, b3m = self._mixture(None, xi, coef)

        delta = b3m-b2m
        epsilon = -b2m*b3m

        return am, b1m, delta, epsilon


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = ALS1983(300, 9.9742e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
    eq = ALS1983(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vg.ccmol))
