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
            * 0 - Original
            * 1 - Boston
            * 2 - Nasrifar

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
    >>> eq = SRK(800, 36451227.52066596, mix, R=8.3144598)
    >>> fir = eq._phir(800, 5000, eq.yi)
    >>> delta = 5000
    >>> tau = 1/800
    >>> print("fir: %0.14f" % (fir["fir"]))
    fir: 0.11586323513845
    >>> print("fird: %0.14f" % (fir["fird"]*delta))
    fird: 0.12741566551477
    >>> print("firt: %0.15f" % (fir["firt"]*tau))
    firt: -0.082603152680518
    >>> print("firdd: %0.15f" % (fir["firdd"]*delta**2))
    firdd: 0.024895937945147
    >>> print("firdt: %0.15f" % (fir["firdt"]*delta*tau))
    firdt: -0.077752734990782
    >>> print("firtt: %0.14f" % (fir["firtt"]*tau**2))
    firtt: -0.10404751064185
    >>> print("firddd: %0.16f" % (fir["firddd"]*delta**3))
    firddd: 0.0060986538256190
    >>> print("firddt: %0.16f" % (fir["firddt"]*delta**2*tau))
    firddt: 0.0089488831000362
    >>> print("firdtt: %0.15f" % (fir["firdtt"]*delta*tau**2))
    firdtt: -0.097937890490398
    >>> print("firttt: %0.14f" % (fir["firttt"]*tau**3))
    firttt: 0.15607126596277
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

    def _cubicDefinition(self, T):
        """Definition of coefficients for generic cubic equation of state"""
        # Schmidt-Wenzel factorization of terms
        self.u = 1
        self.w = 0

        ao = []
        ai = []
        bi = []
        mi = []
        for cmp in self.componente:
            a0, b = self._lib(cmp)
            alfa = self._alfa(cmp, T)
            m = self._m(cmp)
            ao.append(a0)
            ai.append(a0*alfa)
            bi.append(b)
            mi.append(m)

        self.ao = ao
        self.ai = ai
        self.bi = bi
        self.mi = mi

    def _lib(self, cmp):
        ao = 0.42747*self.R**2*cmp.Tc**2/cmp.Pc                          # Eq 5
        b = 0.08664*self.R*cmp.Tc/cmp.Pc                                 # Eq 6
        return ao, b

    def _GEOS(self, xi):
        am, bm = self._mixture("SRK", xi, [self.ai, self.bi])

        delta = bm
        epsilon = 0

        return am, bm, delta, epsilon

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
        m = self._m(cmp)

        if Tr > 1:
            alfa = self.kwargs.get("alpha", 0)
            if alfa == 0:
                alfa = (1+m*(1-Tr**0.5))**2                             # Eq 13
            elif alfa == 1:
                # Use the Boston-Mathias supercritical extrapolation, ref [3]_
                d = 1+m/2                                               # Eq 10
                c = m/d                                                 # Eq 11
                alfa = exp(c*(1-Tr**d))                                 # Eq 9
            elif alfa == 2:
                # Use the Nasrifar-Bolland supercritical extrapolation, ref [4]
                b1 = 0.25*(12-11*m+m**2)                                # Eq 17
                b2 = 0.5*(-6+9*m-m**2)                                  # Eq 18
                b3 = 0.25*(4-7*m+m**2)                                  # Eq 19
                alfa = b1/Tr + b2/Tr**2 + b3/Tr**3                      # Eq 16

        else:
            alfa = (1+m*(1-Tr**0.5))**2                                 # Eq 13
        return alfa

    def _m(self, cmp):
        """Calculate the intermediate parameter for alpha expression"""
        # Eq 15
        return 0.48 + 1.574*cmp.f_acent - 0.176*cmp.f_acent**2

    def _da(self, tau, x):
        """Calculate the derivatives of α, this procedure is used for Helmholtz
        energy formulation of EoS for calculation of properties, alternate alfa
        formulation must define this procedure for any change of formulation
        """
        Tr, rhor = self._Tr()

        # Eq 64-67
        Di = []
        Dt = []
        Dtt = []
        Dttt = []
        for cmp in self.componente:
            Di.append(1-(Tr/cmp.Tc)**0.5/tau**0.5)
            Dt.append((Tr/cmp.Tc)**0.5/2/tau**1.5)
            Dtt.append(-3*(Tr/cmp.Tc)**0.5/4/tau**2.5)
            Dttt.append(15*(Tr/cmp.Tc)**0.5/8/tau**3.5)

        # Eq 63
        Bi = []
        for c1, d in zip(self.mi, Di):
            Bi.append(1+c1*d)

        # Eq 69-71
        Bt = []
        Btt = []
        Bttt = []
        for c1, d, dt, dtt, dttt in zip(self.mi, Di, Dt, Dtt, Dttt):
            Bt.append(c1*dt)
            Btt.append(c1*d*dtt*d**-1)
            Bttt.append(c1*d**2*dttt*d**-2)

        # Eq 73-75
        dait = []
        daitt = []
        daittt = []
        for a, B, bt, btt, bttt in zip(self.ao, Bi, Bt, Btt, Bttt):
            dait.append(2*a*B*bt)
            daitt.append(2*a*(B*btt+bt**2))
            daittt.append(2*a*(B*bttt+3*bt*btt))

        # Eq 52
        uij = []
        for aii in self.ai:
            uiji = []
            for ajj in self.ai:
                uiji.append(aii*ajj)
            uij.append(uiji)

        # Eq 59-61
        duijt = []
        duijtt = []
        duijttt = []
        for aii, diit, diitt, diittt in zip(self.ai, dait, daitt, daittt):
            duijit = []
            duijitt = []
            duijittt = []
            for ajj, djjt, djjtt, djjttt in zip(self.ai, dait, daitt, daittt):
                duijit.append(aii*djjt + ajj*diit)
                duijitt.append(aii*djjtt + 2*diit*djjt + ajj*diitt)
                duijittt.append(
                    aii*djjttt + 3*diit*djjtt + 3*diitt*djjt + ajj*diittt)
            duijt.append(duijit)
            duijtt.append(duijitt)
            duijttt.append(duijittt)

        # Eq 54-56
        daijt = []
        daijtt = []
        daijttt = []
        for uiji, duijit, duijitt, duijittt, kiji in zip(
                uij, duijt, duijtt, duijttt, self.kij):
            daijit = []
            daijitt = []
            daijittt = []
            for u, ut, utt, uttt, k in zip(
                    uiji, duijit, duijitt, duijittt, kiji):
                daijit.append((1-k)/2/u**0.5*ut)
                daijitt.append((1-k)/4/u**1.5*(2*u*utt-ut**2))
                daijittt.append(
                    (1-k)/8/u**2.5*(4*u**2*uttt - 6*u*ut*utt + 3*ut**3))
            daijt.append(daijit)
            daijtt.append(daijitt)
            daijttt.append(daijittt)

        # Eq 51
        damt = 0
        damtt = 0
        damttt = 0
        for xi, daijit, daijitt, daijittt in zip(x, daijt, daijtt, daijttt):
            for xj, dat, datt, dattt in zip(x, daijit, daijitt, daijittt):
                damt += xi*xj*dat
                damtt += xi*xj*datt
                damttt += xi*xj*dattt

        # Eq 126
        aij = []
        for a_i in self.ai:
            aiji = []
            for a_j in self.ai:
                aiji.append((a_i*a_j)**0.5)
            aij.append(aiji)

        daxi = []
        for i, (xi, aiji) in enumerate(zip(x, aij)):
            daxij = 0
            for xj, a in zip(x, aiji):
                daxij += 2*xj*a
            daxi.append(daxij)

        kw = {}
        kw["dat"] = damt
        kw["datt"] = damtt
        kw["dattt"] = damttt
        kw["daxi"] = daxi

        return kw

if __name__ == "__main__":
    from lib.mezcla import Mezcla
    from lib import unidades

    # # mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    # # eq = SRK(300, 9.9742e5, mix, alpha=1)
    # # print('%0.1f' % (eq.Vl.ccmol))
    # # eq = SRK(300, 42.477e5, mix)
    # # print('%0.1f' % (eq.Vg.ccmol))

    # mix = Mezcla(5, ids=[46, 2], caudalMolar=1, fraccionMolar=[0.2152, 0.7848])
    # eq = SRK(144.26, 2.0684e6, mix)
    # print(eq.rhoL.kmolm3)

    # # Ejemplo 6.6, Wallas, pag 340




    # mezcla = Mezcla(2, ids=[1, 2, 40, 41], caudalUnitarioMolar=[0.3177, 0.5894, 0.0715, 0.0214])
    # P = unidades.Pressure(500, "psi")
    # T = unidades.Temperature(120, "F")
    # eq = SRK(T, P, mezcla)
    # print(T)
    # print(eq.x)
    # print([x*(1-eq.x)*5597 for x in eq.xi])
    # print([y*eq.x*5597 for y in eq.yi])
    # print(eq.Ki)

    # Example 6.6 wallas
    P = unidades.Pressure(20, "atm")
    mix = Mezcla(5, ids=[23, 5], caudalMolar=1, fraccionMolar=[0.607, 0.393])
    eq1 = SRK(300, P, mix)
    eq2 = SRK(400, P, mix)
    print(eq1._Dew_T(P))
    print(eq2._Dew_T(P))

    # eq = SRK(500, P, mezcla)
    # print(eq._Dew_T(P))

