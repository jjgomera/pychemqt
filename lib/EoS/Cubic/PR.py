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


class PR(Cubic):
    r"""Peng-Robinson cubic equation of state, [1]_

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + m\left(1-Tr^{0.5}\right)\\
        m = 0.37464 + 1.54226\omega-0.26992\omega^2\\
        \end{array}

    This EoS include too a special alpha temperature dependence for water as
    described in [2]_

    In supercritical states, the α temperature dependence use the
    Boston-Mathias expression, [4]_.

    .. math::
        \begin{array}[t]{l}
        \alpha = \exp\left(c\left(1-T_r^d\right)\right)\\
        d = 1+\frac{m}{2}\\
        c = \frac{m}{d}\\
        \end{array}

    Examples
    --------
    Example 4.3 from [3]_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = PR(300, 9.9742e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '86.7'
    >>> eq = PR(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vg.ccmol)
    '84.1'
    """

    __title__ = "Peng-Robinson (1976)"
    __status__ = "PR76"
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
               "Aqueous Systems with Industrial Applications. ACS (1980), "
               "pag. 393-414",
        "doi": "10.1021/bk-1980-0133.ch020"},
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
         "doi": ""})

    OmegaA = 0.45724
    OmegaB = 0.0778

    def _cubicDefinition(self):
        """Definition of individual components coefficients"""

        # Schmidt-Wenzel factorization of terms
        self.u = 1+2**0.5
        self.w = 1-2**0.5

        ao = []
        ai = []
        bi = []
        mi = []
        for cmp in self.componente:
            a0 = self.OmegaA*R**2*cmp.Tc**2/cmp.Pc                      # Eq 9
            b = self.OmegaB*R*cmp.Tc/cmp.Pc                             # Eq 10

            m, alfa = self._alfa(cmp, self.T)

            ao.append(a0)
            ai.append(a0*alfa)
            bi.append(b)
            mi.append(m)

        self.ao = ao
        self.mi = mi
        self.ai = ai
        self.bi = bi
        self.Bi = [bi*self.P/R/self.T for bi in self.bi]
        self.Ai = [ai*self.P/(R*self.T)**2 for ai in self.ai]

    def _alfa(self, cmp, T):
        """α parameter calculation procedure, separate of general procedure
        to let define derived equation where only change this term like the
        1978 version.

        This method use the original alpha formulation for temperatures below
        the critical temperature and the Boston-Mathias formulation for
        supercritical states.

        The procedure return too the m parameters because it's used in
        Helmholtz formulation

        Parameters
        ----------
        cmp : componente.Componente
            Componente instance
        T : float
            Temperature, [K]

        Returns
        -------
        m : float
            m parameter of equation, [-]
        alpha : float
            alpha parameter of equation, [-]
        """
        Tr = T/cmp.Tc
        m = 0.37464 + 1.54226*cmp.f_acent - 0.2699*cmp.f_acent**2   # Eq 18

        if cmp.id == 62 and Tr < 0.85:
            # Special case for water from [2]_
            alfa = (1.0085677 + 0.82154*(1-(T/cmp.Tc)**0.5))**2         # Eq 6
            m = 0

        elif Tr > 1:
            # Use the Boston-Mathias supercritical extrapolation, ref [4]_
            d = 1+m/2                                                   # Eq 10
            c = m/d                                                     # Eq 11
            alfa = exp(c*(1-Tr**d))                                     # Eq 9

        else:
            alfa = (1+m*(1-(T/cmp.Tc)**0.5))**2                         # Eq 17
        return m, alfa

    def _GEOS(self, xi):
        am, bm = self._mixture("PR", xi, [self.ai, self.bi])

        delta = 2*bm
        epsilon = -bm**2
        return am, bm, delta, epsilon

    def _da(self, tau, zi):
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
        for cmp in self.mezcla.componente:
            tc = cmp.Tc
            Di.append(1-(Tr/tc)**0.5/tau**0.5)
            Dt.append((Tr/tc)**0.5/2/tau**1.5)
            Dtt.append(-3*(Tr/tc)**0.5/4/tau**2.5)
            Dttt.append(15*(Tr/tc)**0.5/8/tau**3.5)

        # Eq 63
        Bi = []
        for c1, d in zip(self.mi, Di):
            Bi.append(1+c1*d)

        # Eq 69-71
        Bt = []
        Btt = []
        Bttt = []
        for c1, d, dt, dtt, dttt in zip(self.mi, Di, Dt, Dtt, Dttt):
            Bt.append(c1**dt)
            Btt.append(c1*(d*dtt)/d)
            Bttt.append(c1*(d**2*dttt)/d**2)

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
                duijittt.append(aii*djjttt + 3*diit*djjtt + 3*diitt*djjt +
                                ajj*diittt)
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
                daijittt.append((1-k)/8/u**2.5 *
                                (4*u**2*uttt - 6*u*ut*utt + 3*ut**3))
            daijt.append(daijit)
            daijtt.append(daijitt)
            daijttt.append(daijittt)

        # Eq 51
        damt = 0
        damtt = 0
        damttt = 0
        for xi, daijit, daijitt, daijittt in zip(zi, daijt, daijtt, daijttt):
            for xj, dat, datt, dattt in zip(zi, daijit, daijitt, daijittt):
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
        for i, (xi, aiji) in enumerate(zip(zi, aij)):
            daxij = 0
            for xj, a in zip(zi, aiji):
                daxij += 2*xj*a
            daxi.append(daxij)

        kw = {}
        kw["dat"] = damt
        kw["datt"] = damtt
        kw["dattt"] = damttt
        kw["daxi"] = daxi
        return kw

    # def _fug(self, Zl, xi, Zv, yi, A, B):
        # tital = self._fugl(Zl, xi, A, B)
        # titav = self._fugl(Zv, yi, A, B)

        # tital2 = self._fugl2(self.Zl, xi)
        # titav2 = self._fugl2(self.Zg, yi)
        # from pprint import pprint
        # print("tital:")
        # pprint(tital)
        # pprint(tital2)
        # print("titav:")
        # pprint(titav)
        # pprint(titav2)
        # return tital, titav

    # def _fug2(self, Z, xi):
    def _fugl2(self, Z, zi, a, b):

        Tr, rhor = self._Tr()
        print(Tr, rhor)

        V = Z*R*self.T/self.P*1e6  # l/mol
        rho = 1/V
        tau = Tr/self.T
        delta = rho/rhor

        fir = self._phir(tau, delta, zi, a, b)
        # Derivation in GERG-2008 paper, Table B4
        firxi = fir["firxi"]
        fird = fir["fird"]
        fir = fir["fir"]

        # Both reducing parameters are fixed, critical values for pure
        # component and unity for mixtures, so their derivatives with molar
        # fraction are zero

        dfxm = sum(x*fx for x, fx in zip(zi, firxi))
        nfirni = []
        for fx in firxi:
            nfirni.append(delta*fird + fx - dfxm)
        nfirni = [fir + ndfni for ndfni in nfirni]
        fi = [x*rho*R*self.T*exp(nfni) for x, nfni in zip(zi, nfirni)]
        # print("fug_helm: ", fi)
        return fi


if __name__ == "__main__":
    from lib.corriente import Mezcla
    from lib import unidades
    # from lib.compuestos import Componente
    # from lib.mEoS import C3

    # c3 = Componente(4)
    # c3.M = C3.M
    # c3.Pc = C3.Pc
    # c3.Tc = C3.Tc
    # c3.f_acent = C3.f_acent

    # mix = Mezcla(5, customCmp=[c3], caudalMolar=1, fraccionMolar=[1])
    # eq = PR(300, 42.477e5, mix)
    # print(eq.V)

    # mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    # eq = PR(300, 42.477e5, mix)
    # print(eq.V)

    # st = C3(T=300, P=42.477e5)
    # print(st.v*st.M)

    # mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    # eq = PR(300, 9.9742e5, mix)
    # print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    # eq = PR(300, 42.477e5, mix)
    # print('%0.1f' % (eq.Vl.ccmol))
    # print('V = %0.1f' % (1/eq.Vl.m3kmol*eq.mezcla.M))

    # mezcla = Mezcla(1, ids=[4], caudalUnitarioMasico=[1])
    # eq = PR(300, 101325, mezcla)
    # print('V = %0.4f' % (1/eq.Vg.m3kmol*eq.mezcla.M))

    # mix = Mezcla(5, ids=[2, 47, 98], caudalMolar=1, fraccionMolar=[0.5, 0.3, 0.2])
    # eq = PR(800, 34.933e6, mix)
    # mix = Mezcla(2, ids=[10, 38, 22, 61], caudalUnitarioMolar=[0.3, 0.5, 0.05, 0.15])
    # eq = PR(340, 101325, mix)

    # mix = Mezcla(2, ids=[2, 4], caudalUnitarioMolar=[0.234, 0.766])
    # eq = PR(260, 1.7e6, mix)
    # mix = Mezcla(2, ids=[2, 4], caudalUnitarioMolar=[0.234, 0.766])
    # eq = PR(420, 1e7, mix)

    # Example 2.3
    mezcla = Mezcla(2, ids=[1, 2, 40, 41], caudalUnitarioMolar=[0.31767, 0.58942, 0.07147, 0.02144])
    P = unidades.Pressure(485, "psi")
    T = unidades.Temperature(120, "F")
    eq = PR(T, P, mezcla)

    # mix = Mezcla(2, ids=[2, 3, 4, 6, 5, 8, 46, 49, 50, 22], caudalUnitarioMolar=[1]*10)
    # eq = PR(293.15, 5e6, mix)
