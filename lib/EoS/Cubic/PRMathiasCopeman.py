#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from csv import reader
import os

from lib.EoS.cubic import Cubic


# Load parameters from database file
dat = {}
fp = os.path.join(os.environ["pychemqt"], "dat", "MathiasCopeman.csv")
with open(fp) as file:
    my_data = reader(file)
    for string in my_data:
        data = string[0].split(":")

        if data[0] == "CAS-nr.":
            continue

        cas = data[0]
        c1 = float(data[1])
        c2 = float(data[2])
        c3 = float(data[3])
        dat[cas] = (c1, c2, c3)


class PRMathiasCopeman(Cubic):
    r"""Mathias-Copeman alpha temperature dependency modification of
    Peng-Robinson cubic equation of state, [1]_

    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^2}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + c_1\left(1-Tr^{0.5}\right) +
        c_2\left(1-Tr^{0.5}\right)^2 + c_3\left(1-Tr^{0.5}\right)^3\\
        \end{array}

    In case T > Tc it use a special alpha temperature dependence give in [3]_

    .. math::
        \alpha^{0.5} = 1 + c_1\left(1-Tr^{0.5}\right)

    The C1, C2 and C3 parameters are those given in [2]_

    Parameters
    ----------
    alpha : int
        Correlation index for alpha expresion at supercritical temperatures:
            * 0 - Original
            * 1 - Coquelet formulation

    Examples
    --------
    Helmholtz energy formulation example for supplementary documentatión from
    [4]_, the critical parameter are override for the valued used in paper
    to get the values of test

    >>> from lib.mezcla import Mezcla
    >>> from lib.compuestos import Componente
    >>> ch4 = Componente(2)
    >>> ch4.Tc, ch4.Pc, ch4.f_acent = 190.564, 4599200, 0.011
    >>> o2 = Componente(47)
    >>> o2.Tc, o2.Pc, o2.f_acent = 154.581, 5042800, 0.022
    >>> ar = Componente(98)
    >>> ar.Tc, ar.Pc, ar.f_acent = 150.687, 4863000, -0.002
    >>> zi = [0.5, 0.3, 0.2]
    >>> cmpList = [ch4, o2, ar]
    >>> mix = Mezcla(5, customCmp=cmpList, caudalMolar=1, fraccionMolar=zi)
    >>> eq = PRMathiasCopeman(800, 34933409.8798343, mix, R=8.3144598)
    >>> fir = eq._phir(800, 5000, eq.yi)
    >>> delta = 5000
    >>> tau = 1/800
    >>> print("fir: %0.15f" % (fir["fir"]))
    fir: 0.034118184296355
    >>> print("fird: %0.15f" % (fir["fird"]*delta))
    fird: 0.050381225002564
    >>> print("firt: %0.14f" % (fir["firt"]*tau))
    firt: 0.10841024634867
    >>> print("firdd: %0.15f" % (fir["firdd"]*delta**2))
    firdd: 0.031329489702333
    >>> print("firdt: %0.15f" % (fir["firdt"]*delta*tau))
    firdt: 0.098515746761245
    >>> print("firtt: %0.14f" % (fir["firtt"]*tau**2))
    firtt: -0.55088266208097
    >>> print("firddd: %0.16f" % (fir["firddd"]*delta**3))
    firddd: -0.0018875965519497
    >>> print("firddt: %0.15f" % (fir["firddt"]*delta**2*tau))
    firddt: -0.016659995071735
    >>> print("firdtt: %0.14f" % (fir["firdtt"]*delta*tau**2))
    firdtt: -0.50060412793624
    >>> print("firttt: %0.13f" % (fir["firttt"]*tau**3))
    firttt: 2.5911592464473
    """
    __title__ = "PR-Mathias-Copeman (1983)"
    __status__ = "PR-MC"
    __doi__ = (
        {"autor": "Mathias, P.M., Copeman, T.W.",
         "title": "Extension of the Peng-Robinson Equation of State to Complex"
                  " Mixtures: Evaluation of the Various Forms of the Local "
                  "Composition Concept",
         "ref": "Fluid Phase Equilibria 13 (1983) 91-108",
         "doi": "10.1016/0378-3812(83)80084-3"},
        {"autor": "Horstmann, S., Jabloniec, A., Krafczyk, J., Fischer, K., "
                  "Gmehling, J.",
         "title": "PSRK group contribution equation of state: comprehensive "
                  "revision and extension IV, including critical constants "
                  "and α-function parameters for 1000 components",
         "ref": "Fluid Phase Equilibria 227 (2005) 157-164",
         "doi": "10.1016/j.fluid.2004.11.002"},
        {"autor": "Coquelet, C., Chapoy, A., Richon, D.",
         "title": "Development of a New Alpha Function for the Peng-Robinson "
                  "Equation of State: Comparative Study of Alpha Function "
                  "Models for Pure Gases (Natural Gas Components) and "
                  "Water-Gas Systems",
         "ref": "Int. J. Thermophys. 25(1) (2004) 133-158",
         "doi": "10.1023/b_ijot.0000022331.46865.2f"},
        {"autor": "Bell, I.H., Jäger, A.",
         "title": "Helmholtz Energy Transformations of Common Cubic Equations "
                  "of State for Use with Pure Fluids and Mixtures",
         "ref": "J. Res. of NIST 121 (2016) 236-263",
         "doi": "10.6028/jres.121.011"})

    def _cubicDefinition(self, T):
        """Definition of coefficients for generic cubic equation of state"""

        # Schmidt-Wenzel factorization of terms
        self.u = 1+2**0.5
        self.w = 1-2**0.5

        ao = []
        bi = []
        ai = []
        C1 = []
        C2 = []
        C3 = []
        for cmp in self.componente:
            if cmp.CASNumber in dat:
                c1, c2, c3 = dat[cmp.CASNumber]
            else:
                # Use the generalized correlation from [2]_
                # Eq 17-19
                c1 = 0.1316*cmp.f_acent**2 + 1.4031*cmp.f_acent + 0.3906
                c2 = -1.3127*cmp.f_acent**2 + 0.3015*cmp.f_acent - 0.1213
                c3 = 0.7661*cmp.f_acent + 0.3041

            ao.append(0.45724*self.R**2*cmp.Tc**2/cmp.Pc)
            bi.append(0.0778*self.R*cmp.Tc/cmp.Pc)

            term = 1-(T/cmp.Tc)**0.5

            alfa = self.kwargs.get("alpha", 0)
            if T/cmp.Tc > 1 and alfa == 1:
                alfa = (1 + c1*term)**2
            else:
                alfa = (1 + c1*term + c2*term**2 + c3*term**3)**2
            ai.append(ao[-1]*alfa)

            C1.append(c1)
            C2.append(c2)
            C3.append(c3)

        self.ai = ai
        self.ao = ao
        self.bi = bi
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3

    def _GEOS(self, xi):
        am, bm = self._mixture(None, xi, [self.ai, self.bi])

        delta = 2*bm
        epsilon = -bm**2

        return am, bm, delta, epsilon

    def _da(self, tau, x):
        """Calculate the derivatives of α, this procedure is used for Helmholtz
        energy formulation of EoS for calculation of properties, alternate alfa
        formulation must define this procedure for any change of formulation
        """
        C1 = self.C1
        C2 = self.C2
        C3 = self.C3

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
                bttt += n*c*(3*(n-1)*d*dt*dtt+(n**2-3*n+2)*dt**3+d**2*dttt) * \
                    d**(n-3)
            Bt.append(bt)
            Btt.append(btt)
            Bttt.append(bttt)

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
    from lib.compuestos import Componente

    ch4 = Componente(2)
    ch4.Tc, ch4.Pc, ch4.f_acent = 190.564, 4599200, 0.011

    o2 = Componente(47)
    o2.Tc, o2.Pc, o2.f_acent = 154.581, 5042800, 0.022

    ar = Componente(98)
    ar.Tc, ar.Pc, ar.f_acent = 150.687, 4863000, -0.002

    mix = Mezcla(5, customCmp=[ch4, o2, ar], caudalMolar=1,
                 fraccionMolar=[0.5, 0.3, 0.2])
    eq = PRMathiasCopeman(800, 34933409.8798343, mix, R=8.3144598)
    phir = eq._phir(800, 5000, eq.yi)
    print(phir["fir"])
    print(0.034118184296355)

    # from lib.mezcla import Mezcla
    # from lib.compuestos import Componente
    # ch4 = Componente(2)
    # ch4.Tc, ch4.Pc, ch4.f_acent = 190.564, 4599200, 0.011
    # o2 = Componente(47)
    # o2.Tc, o2.Pc, o2.f_acent = 154.581, 5042800, 0.022
    # ar = Componente(98)
    # ar.Tc, ar.Pc, ar.f_acent = 150.687, 4863000, -0.002
    # mix = Mezcla(5, customCmp=[ch4, o2, ar], caudalMolar=1,
                 # fraccionMolar=[0.5, 0.3, 0.2])
    # eq = PRMathiasCopeman(800, 34933409.8798343, mix)

    # eq = PRMathiasCopeman(800, 34933409.8798343, mix)
    # print(eq._phir(800, 5000, eq.yi))
