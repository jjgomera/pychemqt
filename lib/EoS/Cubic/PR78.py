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


from lib.EoS.Cubic.PR import PR


class PR78(PR):
    r"""Peng-Robinson cubic equation of state

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + m\left(1-Tr^{0.5}\right)\\
        m = 0.37464 + 1.54226\omega-0.26992\omega^2\\
        \end{array}
    """

    __title__ = "Peng-Robinson (1978)"
    __status__ = "PR78"
    __doi__ = (
      {
        "autor": "Peng, D.-Y., Robinson, D.B.",
        "title": "A New Two-Constant Equation of State",
        "ref": "Ind. Eng. Chem. Fund. 15(1) (1976) 59-64",
        "doi": "10.1021/i160057a011"},
      {
        "autor": "Peng, D.-Y., Robinson, D.B.",
        "title": "Two- and Three-Phase Equilibrium Calculations for Coal "
                 "Gasification and Related Processes",
        "ref": "in Newman, Barner, Klein, Sandler (Eds.). Thermodynamic of "
               "Aqueous Systems with Industrial Applications. ACS (1980), pag."
               " 393-414",
        "doi": "10.1021/bk-1980-0133.ch020"})

    def _alfa(self, cmp, T):
        """Custom expression for α"""
        if cmp.f_acent <= 0.491:
            m = 0.3746 + 1.54226*cmp.f_acent - 0.26992*cmp.f_acent**2
        else:
            m = 0.379642 + 1.48503*cmp.f_acent - 0.164423*cmp.f_acent**2 + \
                0.016666*cmp.f_acent**3
        alfa = (1+m*(1-(T/cmp.Tc)**0.5))**2                         # Eq 17
        return m, alfa


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PR78(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PR78(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
