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

from lib.bip import Kij, Mixing_Rule
from lib.EoS.cubic import Cubic, CubicHelmholtz


class PR(Cubic):
    """Peng-Robinson cubic equation of state
    
    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^{2.5}}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + k\left(1-Tr^{0.5}\right)\\
        m = 0.37464 + 1.54226\omega-0.26992\omega^2\\
        \end{array}

    Examples
    --------
    Example 4.3 from 1_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = PR(300, 9.9742e5, mix)
    >>> '%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol)
    '2039 86.7'

    # Tiny desviation
    # '2038 86.8'
    >>> eq = PR(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '84.1'
    """

    __title__ = "Peng-Robinson (1976)"
    __status__ = "PR76"
    __doi__ = {
        "autor": "Peng, D.-Y., Robinson, D.B.",
        "title": "A New Two-Constant Equation of State",
        "ref": "Ind. Eng. Chem. Fund. 15(1) (1976) 59-64",
        "doi": "10.1021/i160057a011"}

    def __init__(self, T, P, mezcla):

        x = mezcla.fraccion

        ao = []
        ai = []
        bi = []
        C1 = []
        C2 = []
        C3 = []
        for componente in mezcla.componente:
            a0, a, b, m = self.__lib(componente, T)
            ao.append(a0)
            ai.append(a)
            bi.append(b)
            C1.append(m)
            C2.append(0)
            C3.append(0)

        self.kij = Kij(mezcla.ids, PR)
        am, bm = Mixing_Rule(x, [ai, bi], self.kij)

        self.ao = ao
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3
        self.ai = ai
        self.bi = bi
        self.b = bm
        self.tita = am
        self.delta = 2*bm
        self.epsilon = -bm**2
        # print("tita", am)
        # print("kij", self.kij)
        # print("ai", ai)

        # self.u = 2
        # self.w = -1

        super(PR, self).__init__(T, P, mezcla)

        # Tr = mezcla.Tc
        # tau = Tr/T
        # rhor = 1/mezcla.Vc/1000  # m3/mol
        # # rhor = 1
        # print("V", self.V)
        # rho = 1/self.V[0]
        # delta = rho/rhor
        # print(delta, rho, rhor)
        # print((self.Z[0]-1))

        # self._phir(mezcla, tau, delta, T, rho, ao, ai, C1, C2, C3)

    def _phir(self, mezcla, tau, delta, T, rho, ao, ai, C1, C2, C3):

        Tr = mezcla.Tc

        # Eq 64-67
        Di = []
        Dt = []
        Dtt = []
        Dttt = []
        for cmp in mezcla.componente:
            tc = cmp.Tc
            Di.append(1-(Tr/tc)**0.5/tau**0.5)
            Dt.append((Tr/tc)**0.5/2/tau**1.5)
            Dtt.append(-3*(Tr/tc)**0.5/4/tau**2.5)
            Dttt.append(15*(Tr/tc)**0.5/8/tau**3.5)

        # Eq 63
        Bi = []
        for c1, c2, c3, d in zip(C1, C2, C3, Di):
            Bi.append(1+c1*d+c2*d**2+c3*d**3)

        # Eq 69-71
        Bt = []
        Btt = []
        Bttt = []
        for c1, c2, c3, d, dt, dtt, dttt in zip(C1, C2, C3, Di, Dt, Dtt, Dttt):
            cs = (c1, c2, c3)
            bt = 0
            btt = 0
            bttt = 0
            for n, c in enumerate(cs):
                n += 1
                bt += n*c*d**(n-1)*dt
                btt += n*c*((n-1)*dt**2+d*dtt)*d**(n-2)
                bttt += n*c*(3*(n-1)*d*dt*dtt+(n**2-3*n+2)*dt**3+d**2*dttt)*d**(n-3)
            Bt.append(bt)
            Btt.append(btt)
            Bttt.append(bttt)

        # Eq 73-75
        dait = []
        daitt = []
        daittt = []
        for a, B, bt, btt, bttt in zip(ao, Bi, Bt, Btt, Bttt):
            dait.append(2*a*B*bt)
            daitt.append(2*a*(B*btt+bt**2))
            daittt.append(2*a*(B*bttt+3*bt*btt))

        # Eq 52
        uij = []
        for aii in ai:
            uiji = []
            for ajj in ai:
                uiji.append(aii*ajj)
            uij.append(uiji)

        # Eq 59-61
        duijt = []
        duijtt = []
        duijttt = []
        for aii, diit, diitt, diittt in zip(ai, dait, daitt, daittt):
            duijit = []
            duijitt = []
            duijittt = []
            for ajj, djjt, djjtt, djjttt in zip(ai, dait, daitt, daittt):
                duijit.append(aii*djjt + ajj*diit)
                duijitt.append(aii*djjtt + 2*diit*djjt + ajj*diitt)
                duijittt.append(aii*djjttt + 3*diit*djjtt + 3*diitt*djjt + ajj*diittt)
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
            for u, ut, utt, uttt, k in zip(uiji, duijit, duijitt, duijittt, kiji):
                daijit.append((1-k)/2/u**0.5*ut)
                daijitt.append((1-k)/4/u**1.5*(2*u*utt-ut**2))
                daijittt.append((1-k)/8/u**2.5*(4*u**2*uttt - 6*u*ut*utt + 3*ut**3))
            daijt.append(daijit)
            daijtt.append(daijitt)
            daijttt.append(daijittt)

        # Eq 51
        damt = 0
        damtt = 0
        damttt = 0
        for xi, daijit, daijitt, daijittt in zip(mezcla.fraccion, daijt, daijtt, daijttt):
            for xj, dat, datt, dattt in zip(mezcla.fraccion, daijit, daijitt, daijittt):
                damt += xi*xj*dat
                damtt += xi*xj*datt
                damttt += xi*xj*dattt

        kw = {}
        kw["rhoc"] = 1/mezcla.Vc
        kw["Tc"] = mezcla.Tc
        kw["Delta1"] = 1+2**0.5
        kw["Delta2"] = 1-2**0.5
        kw["b"] = self.b
        kw["a"] = self.tita
        kw["dat"] = damt
        kw["datt"] = damtt
        kw["dattt"] = damttt

        print(tau, delta, kw)
        fir = CubicHelmholtz(tau, delta, **kw)
        # print(fir)
        print(delta, fir["fird"], R, T, rho)
        print("P", (1+delta*fir["fird"])*R*T*rho*1000)

    def __lib(self, cmp, T):
        """Calculate compound specific properties"""
        ao = 0.45724*R**2*cmp.Tc**2/cmp.Pc                              # Eq 9
        b = 0.0778*R*cmp.Tc/cmp.Pc                                      # Eq 10
        m = 0.37464 + 1.54226*cmp.f_acent - 0.2699*cmp.f_acent**2       # Eq 18
        alfa = (1+m*(1-(T/cmp.Tc)**0.5))**2                             # Eq 17
        return ao, ao*alfa, b, m


if __name__ == "__main__":
    # from lib.corriente import Mezcla
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
    
    from lib.mezcla import Mezcla
    mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    eq = PR(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = PR(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))

