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


from math import exp

from scipy.constants import R

from lib.EoS.cubic import Cubic


class SRK(Cubic):
    r"""Equation of state of Soave-Redlich-Kwong (1972)

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)}\\
        a = 0.42747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.08664\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + m\left(1-Tr^{0.5}\right)\\
        m = 0.48 + 1.574\omega - 0.176\omega^2\\
        \end{array}

    In supercritical states, the α temperature dependence can use different
    extrapolation correlation:

    * Boston-Mathias expression, [3]_

    .. math::
        \begin{array}[t]{l}
        \alpha = \exp\left(c\left(1-T_r^d\right)\right)\\
        d = 1+\frac{m}{2}\\
        c = \frac{m}{d}\\
        \end{array}

    * Nasrifar-Bolland expression, [4]_

    .. math::
        \begin{array}[t]{l}
        \alpha = \frac{b_1}{T_r} + \frac{b_2}{T_r^2} + \frac{b_3}{T_r^3}\\
        b_1 = 0.25\left(12 - 11m + m^2\right)\\
        b_2 = 0.5\left(-6 + 9m - m^2\right)\\
        b_3 = 0.25\left(4 - 7m + m^2\right)\\
        \end{array}


    Parameters
    ----------
    alpha : int
        Correlation index for alpha expresion at supercritical temperatures:
            * 0 - Boston
            * 1 - Nasrifar

    Examples
    --------
    Example 4.3 from [2]_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = SRK(300, 9.9742e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '98.4'
    >>> eq = SRK(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vg.ccmol)
    '95.1'
    """

    __title__ = "Soave-Redlich-Kwong (1972)"
    __status__ = "SRK72"
    __doi__ = (
      {
        "autor": "Soave, G.",
        "title": "Equilibrium Constants from a modified Redlich-Kwong "
                 "Equation of State",
        "ref": "Chem. Eng. Sci. 27 (1972) 1197-1203",
        "doi": "10.1016/0009-2509(72)80096-4"},
      {
         "autor": "Poling, B.E, Prausnitz, J.M, O'Connell, J.P",
         "title": "The Properties of Gases and Liquids 5th Edition",
         "ref": "McGraw-Hill, New York, 2001",
         "doi": ""},
      {
         "autor": "Boston, J.F., Mathias, P.M.",
         "title": "Phase Equilibria in a Third-Generation Process Simulator",
         "ref": "Presented at: 'Phase Equilibria and Fluid Properties in the "
                "Chemical Industries', Berlin, March 17-21, 1980.",
         "doi": ""},
      {
         "autor": "Nasrifar, Kh., Bolland, O.",
         "title": "Square-Well Potential and a New α Function for the Soave-"
                  "Redlich-Kwong Equation of State",
         "ref": "Ind. Eng. Chem. Res. 43(21) (2004) 6901-6909",
         "doi": "10.1021/ie049545i"},
      )

    def _cubicDefinition(self):
        """Definition of coefficients for generic cubic equation of state"""
        # Schmidt-Wenzel factorization of terms
        self.u = 1
        self.w = 0

        ai = []
        bi = []
        for cmp in self.componente:
            a, b = self._lib(cmp, self.T)
            ai.append(a)
            bi.append(b)

        self.ai = ai
        self.bi = bi
        self.Bi = [bi*self.P/R/self.T for bi in self.bi]
        self.Ai = [ai*self.P/(R*self.T)**2 for ai in self.ai]

    def _GEOS(self, xi):
        am, bm = self._mixture("SRK", xi, [self.ai, self.bi])

        delta = bm
        epsilon = 0

        return am, bm, delta, epsilon

    def _lib(self, cmp, T):
        ao = 0.42747*R**2*cmp.Tc**2/cmp.Pc                              # Eq 5
        b = 0.08664*R*cmp.Tc/cmp.Pc                                     # Eq 6
        alfa = self._alfa(cmp, T)
        return ao*alfa, b

    def _alfa(self, cmp, T):
        """α parameter calculation procedure, separate of general procedure
        to let define derived equation where only change this term.

        This method use the original alpha formulation for temperatures below
        the critical temperature and can choose by configuration between:

         * Boston-Mathias formulation
         * Nasrifar-Bolland formulation

        Parameters
        ----------
        cmp : componente.Componente
            Componente instance
        T : float
            Temperature, [K]

        Returns
        -------
        alpha : float
            alpha parameter of equation, [-]
        """
        Tr = T/cmp.Tc
        m = 0.48 + 1.574*cmp.f_acent - 0.175*cmp.f_acent**2             # Eq 15

        if Tr > 1:
            alfa = self.kwargs.get("alpha", 0)
            if alfa == 0:
                # Use the Boston-Mathias supercritical extrapolation, ref [3]_
                d = 1+m/2                                               # Eq 10
                c = m/d                                                 # Eq 11
                alfa = exp(c*(1-Tr**d))                                 # Eq 9
            elif alfa == 1:
                # Use the Nasrifar-Bolland supercritical extrapolation, ref [4]
                b1 = 0.25*(12-11*m+m**2)                                # Eq 17
                b2 = 0.5*(-6+9*m-m**2)                                  # Eq 18
                b3 = 0.25*(4-7*m+m**2)                                  # Eq 19
                alfa = b1/Tr + b2/Tr**2 + b3/Tr**3                      # Eq 16

        else:
            alfa = (1+m*(1-Tr**0.5))**2                                 # Eq 13
        return alfa


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    # eq = SRK(300, 9.9742e5, mix, alpha=1)
    # print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = SRK(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))
