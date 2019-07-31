#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from csv import reader
import os

from scipy.constants import R

from lib.EoS.cubic import Cubic, CubicHelmholtz


# TODO: Upgrade CAS-number for oxygen in database and in MathiasCopeman.csv
# Correct value in mEoS.O2, mising last term in other place

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
    """Mathias-Copeman alpha temperature dependency modification of
    Peng-Robinson cubic equation of state
    
    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V\left(V+b\right)+b\left(V-b\right)}\\
        a = 0.45747\frac{R^2T_c^{2.5}}{P_c}\alpha\\
        b = 0.0778\frac{RT_c}{P_c}\\
        \alpha^{0.5} = 1 + c_1\left(1-Tr^{0.5}\right) +
        c_2\left(1-Tr^{0.5}\right)^2 + c_3\left(1-Tr^{0.5}\right)^3\\
        \end{array}

    The C1, C2 and C3 parameters are those given in [2]_
    """
    __title__="PR-Mathias-Copeman (1983)"
    __status__="PR-MC"
    __doi__ = (
            {
        "autor": "Mathias, P.M., Copeman, T.W.",
        "title": "Extension of the Peng-Robinson Equation of State to Complex "
                 "Mixtures: Evaluation of the Various Forms of the Local "
                 "Composition Concept",
        "ref": "Fluid Phase Equilibria 13 (1983) 91-108",
        "doi": "10.1016/0378-3812(83)80084-3"},
            {
        "autor": "Horstmann, S., Jabloniec, A., Krafczyk, J., Fischer, K., "
                 "Gmehling, J.",
        "title": "PSRK group contribution equation of state: comprehensive "
                 "revision and extension IV, including critical constants "
                 "and α-function parameters for 1000 components",
        "ref": "Fluid Phase Equilibria 227 (2005) 157-164",
        "doi": "10.1016/j.fluid.2004.11.002"},
            {
        "autor": "Coquelet, C., Chapoy, A., Richon, D.",
        "title": "Development of a New Alpha Function for the Peng-Robinson "
                 "Equation of State: Comparative Study of Alpha Function "
                 "Models for Pure Gases (Natural Gas Components) and "
                 "Water-Gas Systems",
        "ref": "Int. J. Thermophys. 25(1) (2004) 133-158",
        "doi": "10.1023/b_ijot.0000022331.46865.2f"})

    def __init__(self, T, P, mezcla):
        """Initialization procedure
        
        Parameters
        ----------

        """

        self.T = T
        self.P = P
        self.mezcla = mezcla
        
        ao = []
        bi = []
        ai = []
        C1 = []
        C2 = []
        C3 = []
        for cmp in mezcla.componente:
            if cmp.CASNumber in dat:
                c1, c2, c3 = dat[cmp.CASNumber]
            else:
                # Use the generalized correlation from [2]_
                # Eq 17-19
                c1 = 0.1316*cmp.f_acent**2 + 1.4031*cmp.f_acent + 0.3906
                c2 = -1.3127*cmp.f_acent**2 + 0.3015*cmp.f_acent - 0.1213
                c3 = 0.7661*cmp.f_acent + 0.3041

            ao.append(0.45724*R**2*cmp.Tc**2/cmp.Pc)
            bi.append(0.0778*R*cmp.Tc/cmp.Pc)

            term = 1-(T/cmp.Tc)**0.5
            # if T > cmp.Tc:
                # alfa = (1 + c1*term)**2
            # else:
            alfa = (1 + c1*term + c2*term**2 + c3*term**3)**2
            ai.append(ao[-1]*alfa)

            C1.append(c1)
            C2.append(c2)
            C3.append(c3)

        am, bm = self._mixture(None, mezcla.ids, [ai, bi])

        self.ai = ai
        self.bi = bi
        self.b = bm
        self.tita = am
        self.delta = 2*bm
        self.epsilon = -bm**2

        # tdadt=0
        # for i in range(len(mezcla.componente)):
            # for j in range(len(mezcla.componente)):
                # tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])
        # self.dTitadT=tdadt
        # self.u=2
        # self.w=-1

        super(PRMathiasCopeman, self).__init__(T, P, mezcla)
        print(1/self.Vl, 1/self.Vg)
        print(self.rho)

        Tr = 1
        tau = Tr/T
        self.rho = self.rho[0]*1000  # m3/mol
        rhor = 1  # m3/mol
        delta = self.rho/rhor

        self._phir(tau, delta, ao, ai, C1, C2, C3)
        excess = self._excess(tau, delta)
        H_exc = excess["H"]*R*T/mezcla.M
        print(self.mezcla._Ho(T), H_exc)

    def _phir(self, tau, delta, ao, ai, C1, C2, C3):

        # Tr = self.mezcla.Tc
        Tr = 1

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
        for xi, daijit, daijitt, daijittt in zip(self.mezcla.fraccion, daijt, daijtt, daijttt):
            for xj, dat, datt, dattt in zip(self.mezcla.fraccion, daijit, daijitt, daijittt):
                damt += xi*xj*dat
                damtt += xi*xj*datt
                damttt += xi*xj*dattt

        kw = {}
        # kw["rhoc"] = 1/self.mezcla.Vc
        # kw["Tc"] = self.mezcla.Tc
        kw["Delta1"] = 1+2**0.5
        kw["Delta2"] = 1-2**0.5
        kw["b"] = self.b
        kw["a"] = self.tita
        kw["dat"] = damt
        kw["datt"] = damtt
        kw["dattt"] = damttt

        print(tau, delta, kw)
        p = CubicHelmholtz(tau, delta, **kw)
        p["tau"] = tau
        p["delta"] = delta
        print("fir: ", p["fir"])
        print("fird: ", p["fird"]*delta)
        print("firt: ", p["firt"]*tau)
        print("firdd: ", p["firdd"]*delta**2)
        print("firdt: ", p["firdt"]*delta*tau)
        print("firtt: ", p["firtt"]*tau**2)
        print("firddd: ", p["firddd"]*delta**3)
        print("firddt: ", p["firddt"]*delta**2*tau)
        print("firdtt: ", p["firdtt"]*delta*tau**2)
        print("firttt: ", p["firttt"]*tau**3)
        print("P", (1+delta*p["fird"])*R*self.T*self.rho-self.P)
        self.fir = p


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mix = Mezcla(tipo=5, caudalMolar=1, ids=[2, 47, 98], fraccionMolar=[0.5, 0.3, 0.2])
    eq = PRMathiasCopeman(800, 34937532, mix)
