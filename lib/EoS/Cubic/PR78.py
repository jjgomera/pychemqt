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
        m = 0.37464 + 1.54226\omega-0.26992\omega^2 if \omega < 0.491\\
        m = 0.379642 + 1.48503\omega - 0.164423*\omega^2 + 0.016666*\omega^3\\
        \end{array}

    Examples
    --------
    Helmholtz energy formulation example for supplementary documentatión from
    [4]_, the critical parameter are override for the valued used in paper to
    get the values of test with high precision

    >>> from lib.mezcla import Mezcla
    >>> from lib import unidades
    >>> from lib.compuestos import Componente
    >>> ch4 = Componente(2)
    >>> ch4.Tc, ch4.Pc, ch4.f_acent = 190.564, 4599200, 0.011
    >>> o2 = Componente(47)
    >>> o2.Tc, o2.Pc, o2.f_acent = 154.581, 5042800, 0.022
    >>> ar = Componente(98)
    >>> ar.Tc, ar.Pc, ar.f_acent = 150.687, 4863000, -0.002
    >>> mix = Mezcla(5, customCmp=[ch4, o2, ar], caudalMolar=1,
    ...              fraccionMolar=[0.5, 0.3, 0.2])
    >>> eq = PR78(800, 36451227.52066596, mix, R=8.3144598)
    >>> fir = eq._phir(800, 5000, eq.yi)
    >>> delta = 5000
    >>> tau = 1/800
    >>> print("fir: %0.15f" % (fir["fir"]))
    fir: 0.084339749584296
    >>> print("fird: %0.15f" % (fir["fird"]*delta))
    fird: 0.096019116018396
    >>> print("firt: %0.14f" % (fir["firt"]*tau))
    firt: -0.10134978074971
    >>> print("firdd: %0.15f" % (fir["firdd"]*delta**2))
    firdd: 0.023611667278971
    >>> print("firdt: %0.15f" % (fir["firdt"]*delta*tau))
    firdt: -0.092099683110520
    >>> print("firtt: %0.15f" % (fir["firtt"]*tau**2))
    firtt: -0.078186052271240
    >>> print("firddd: %0.16f" % (fir["firddd"]*delta**3))
    firddd: 0.0017433108161805
    >>> print("firddt: %0.15f" % (fir["firddt"]*delta**2*tau))
    firddt: 0.015574974734224
    >>> print("firdtt: %0.15f" % (fir["firdtt"]*delta*tau**2))
    firdtt: -0.071050085995025
    >>> print("firttt: %0.14f" % (fir["firttt"]*tau**3))
    firttt: 0.11727907840686
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
            m = 0.37464 + 1.54226*cmp.f_acent - 0.26992*cmp.f_acent**2
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
