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


class Nasrifar(Cubic):
    r"""Nasrifar-Moshfeguian cubic equation of state, [1]_

    This is a two parameter cubic equation of state with both parameters
    temperature dependent

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.497926\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.094451\frac{RT_c}{P_c}\\
        \\
        \theta = \frac{T-T_{pt}}{T_c-T_{pt}}\\
        a = a_c\left(1+m_a\left(1-\sqrt{\theta}\right)\right)^2\\
        m_a = \sqrt{\frac{a_{pt}}{a_c}}-1\\
        b = b_c\left(1+m_b\left(1-\theta\right)\right)\\
        m_b = \frac{b_{pt}}{b_c}-1\\
        \\
        \frac{a_{pt}}{b_{pt}RT_{pt}} = 29.7056\\
        \frac{T_{pt}}{T_c} = 0.2498 + 0.3359\omega - 0.1037\omega^2\\
        \frac{b_{pt}}{b_c} = 1 - 0.1519\omega - 3.9462\omega^2 +
        7.0538\omega^3\\
        \end{array}
    """

    __title__ = "Nasrifar-Moshfeghian (2001)"
    __status__ = "Nasrifar"
    __doi__ = {
        "autor": "Nasrifar, Kh., Moshfeghian, M.",
        "title": "A New Cubic Equation of State for Simple Fluids: Pure and "
                 "Mixture",
        "ref": "Fluid Phase Equilibria 190 (2001) 73-88",
        "doi": "10.1016/s0378-3812(01)00592-1"},

    def _cubicDefinition(self, T):
        """Definition of individual components coefficients"""

        # Schmidt-Wenzel factorization of terms
        self.u = 1+3**0.5
        self.w = 1-3**0.5

        ai = []
        bi = []
        for cmp in self.componente:
            a, b = self._lib(cmp, T)
            ai.append(a)
            bi.append(b)

        self.ai = ai
        self.bi = bi

    def _GEOS(self, xi):
        am, bm = self._mixture(None, xi, [self.ai, self.bi])

        delta = 2*bm
        epsilon = -2*bm**2
        return am, bm, delta, epsilon

    def _lib(self, cmp, T):

        ac = 0.497926*R**2*cmp.Tc**2/cmp.Pc                             # Eq 3
        bc = 0.094451*R*cmp.Tc/cmp.Pc                                   # Eq 4

        # Eq 20
        Tpt = cmp.Tc*(0.2498 + 0.3359*cmp.f_acent - 0.1037*cmp.f_acent**2)

        Tita = (T-Tpt)/(cmp.Tc-Tpt)                                     # Eq 6

        # Eq 21
        bptr = 1 - 0.1519*cmp.f_acent - 3.9462*cmp.f_acent**2 + \
            7.0539*cmp.f_acent**3
        bpt = bptr*bc

        mb = bptr-1                                                     # Eq 17
        apt = 29.7056*bpt*R*Tpt                                         # Eq 19

        b = bc*(1+mb*(1-Tita))                                          # Eq 16

        ma = (apt/ac)**0.5-1                                            # Eq 11
        a = ac*(1+ma*(1-Tita**0.5))**2                                  # Eq 10

        return a, b


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = Nasrifar(300, 9.9742e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
    eq = Nasrifar(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vg.ccmol))

    mix = Mezcla(5, ids=[46, 2], caudalMolar=1, fraccionMolar=[0.04752, 0.95248])
    eq = Nasrifar(105, 5e7, mix)
    print(eq.rhoL.kmolm3)

