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


from math import exp, log

from scipy import roots
from scipy.constants import R

from lib.EoS.cubic import Cubic, CubicHelmholtz


dat = {
    98 : (0.27693, 0),
    105 : (0.35482, 0),
    # : (0.31265, 0),  # phosphine
    208 : (0.35541, 0.12274),
    # 951 : (0.39401, 0.18031),
    104 : (0.38358, 0),
    1 : (0.03067, 0),
    62 : (0.46195, 0.23002),
    50 : (0.37505, 0),
    63 : (0.48939, 0.29235),
    46 : (0.44178, 0),
    110 : (0.42432, 0),
    107 : (0.24137, 0),
    47 : (0.31734, 0.10963),
    51 : (0.56481, 0),
    111: (0.68490, 0),
    215 : (0.49592, 0),
    216 : (0.49195, 0),
    101 : (0.54091, 0),
    217 : (0.48867, 0),
    100 : (0.47175, 0),
    218 : (0.48214, 0),
    48 : (0.41227, 0),
    49 : (0.44277, 0),
    220 : (0.49946, 0),
    642 : (0.42062, 0.24053),
    115 : (0.45529, 0),
    2 : (0.36438, 0.05927),
    117 : (0.66130, 0),
    # : (0.42504, 0),  # C1Mercaptan
    65 : (0.49620, 0),
    125 : (0.36837, 0),
    22 : (0.39732, 0.15798),
    130 : (0.45330, 0),
    131 : (0.48028, 0),
    132 : (0.45081, 0),
    3 : (0.43271, 0.13725),
    133 : (0.48354, 0),
    134 : (0.86691, 0),
    137 : (0.47451, 0),
    136 : (0.45143, 0),
    66 : (0.46765, 0),
    258 : (0.43029, 0),
    23 : (0.37929, 0),
    140: (0.49283, 0),
    141 : (0.51114, 0),
    142 : (0.54989, 0),
    4 : (0.47576, 0.16608),
    146 : (0.84083, 0),
    24 : (0.45998, 0),
    155 : (0.58129, 0),
    156 : (0.57988, 0),
    157 : (0.55656, 0),
    6 : (0.52441, 0),
    5 : (0.54526, 0),
    162 : (0.52362, 0),
    290 : (0.54596, 0),
    294 : (0.52911, 0),
    166 : (0.60904, 0),
    309 : (0.61855, 0),
    310 : (0.61344, 0),
    311 : (0.60160, 0),
    8 : (0.56773, 0),
    7 : (0.50036, 0),
    318 : (0.57911, 0),
    171 : (0.49460, 0),
    172 : (0.51065, 0),
    319 : (0.50023, 0),
    173 : (0.50353, 0),
    40 : (0.50078, 0.21953),
    38 : (0.51281, 0),
    10 : (0.60943, 0.24797),
    11 : (0.59452, 0),
    12 : (0.63794, 0)}


class TB(Cubic):
    r"""Trebble-Bishnoi cubic equation of state
    
    .. math::
        \begin{array}[t]{l}
        P = \frac{RT}{V-b}-\frac{a}{V^2+(b+c)V+bc+d^2}\\
        a = a_c\exp{q_1\left(1-T_r\right)\\
        b = b_c\left(1+q_2\left(1-T_r+\ln{T_r}\right)\right)\\
        c = \frac{RT_c}{P_c}\left(1-3\zetta\right)\\
        q_1 = f\left(\omega, Z_c\right)\\
        q_2 = f\left(\omega\right)\\
        \zetta = 1.075Z_c\\
        d(cc/mol) = 0.241V_c(cc/mol)-5\\
        \end{array}

    This four parameter cubic equation include a b parameter with temperature
    dependence. The model need two compound specific parameters but include too
    a generalized correlation for both parameters.

    Examples
    --------
    Example 4.3 from 1_, Propane saturated at 300K

    >>> from lib.mezcla import Mezcla
    >>> mix = Mezcla(5, ids=[4], caudalMolar=1, fraccionMolar=[1])
    >>> eq = TB(300, 9.9742e5, mix)
    >>> '%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol)
    '2019 90.6'
    >>> eq = TB(300, 42.477e5, mix)
    >>> '%0.1f' % (eq.Vl.ccmol)
    '88.9'

    # Tiny desviation
    # '2025 89.4'
    # '88.3'
    """

    __title__ = "Trebble-Bishnoi (1987)"
    __status__ = "TB"
    __doi__ = {
        "autor": "Trebble, M.A., Bishnoi, P.R.",
        "title": "Development of a New Four-Parameter Cubic Equation of State",
        "ref": "Fluid Phase Equilibria 35 (1987) 1-8",
        "doi": "10.1016/0378-3812(87)80001-8"},

    def __init__(self, T, P, mezcla):

        self.T = T
        self.P = P
        self.mezcla = mezcla

        ao = []
        ai = []
        bi = []
        ci = []
        di = []
        for cmp in mezcla.componente:
            Tr = T/cmp.Tc

            TcPci = cmp.Tc*cmp.Pc.MPa
            if cmp.f_acent < 0.225:
                # Eq 36
                TcPcH = 775.9 + 12003*cmp.f_acent - 57335*cmp.f_acent**2 + \
                    91393*cmp.f_acent**3
            elif cmp.f_acent <=1:
                TcPcH = 1876 - 1160*cmp.f_acent                         # Eq 37
            else:
                TcPcH = TcPci

            if cmp.f_acent >= -0.14:
                # Eq 34
                Zc = 0.29 - 0.0885*cmp.f_acent - 0.0005/(TcPci**0.5-TcPcH**0.5)
            else:
                Zc = 0.3024

            Xc = 1.075*Zc                                               # Eq 28
            d = 0.341*cmp.Vc.ccg*cmp.M - 0.005                          # Eq 29

            # Convert d from cc/mol to m3/mol
            d /= 1e6

            Dc = d*cmp.Pc/R/cmp.Tc                                      # Eq 12
            Cc = 1-3*Xc                                                 # Eq 6

            # Bc calculated ad the smallest positive root of Eq 8
            B = roots([1, 2-3*Xc, 3*Xc**2, -Dc**2-Xc**3])
            Bpositivos=[]
            for Bi in B:
                if Bi > 0:
                    Bpositivos.append(Bi)
            Bc = min(Bpositivos)

            Ac = 3*Xc**2 + 2*Bc*Cc + Bc + Cc + Bc**2 + Dc**2            # Eq 7
            
            if cmp.id in dat:
                q1, q2 = dat[cmp.id]
            else:
                if cmp.id == 212 and Tr <= 1:
                    q1 = -0.31913
                elif cmp.f_acent < -0.1:
                    # Eq 30
                    q1 = 0.66208 + 4.63961*cmp.f_acent + 7.45183*cmp.f_acent**2
                elif cmp.f_acent <= 0.4:
                    # Eq 32
                    q1 = 0.35 + 0.7924*cmp.f_acent + 0.1875*cmp.f_acent**2 - \
                        28.93*(0.3-Zc)**2
                elif cmp.f_acent > 0.4:
                    # Eq 33
                    q1 = 0.32 + 0.9424*cmp.f_acent - 28.93*(0.3-Zc)**2

                if cmp.f_acent < -0.0423:
                    q2 = 0
                elif cmp.f_acent <= 0.3:
                    # Eq 26
                    q2 = 0.05246 + 1.15058*cmp.f_acent - \
                        1.99348*cmp.f_acent**2 + 1.5949*cmp.f_acent**3 - \
                        1.39267*cmp.f_acent**4
                else:
                    q2 = 0.17959 + 0.23471*cmp.f_acent                  # Eq 27

            alfa = exp(q1*(1-Tr))
            ac = Ac*R**2*cmp.Tc**2/cmp.Pc
            a = ac*alfa                                                 # Eq 23

            if T <= cmp.Tc:
                beta = 1 + q2*(1 - Tr + log(Tr))                        # Eq 19
            else:
                beta = 1
            bc = Bc*R*cmp.Tc/cmp.Pc
            b = bc*beta                                                 # Eq 18

            c = Cc*R*cmp.Tc/cmp.Pc

            ao.append(ac)
            ai.append(a)
            bi.append(b)
            ci.append(c)
            di.append(d)

        am, bm, cm, dm = self._mixture(None, mezcla.ids, [ai, bi, ci, di])

        self.ao = ao
        self.ai = ai
        self.bi = bi
        self.b = bm
        self.tita = am
        self.delta = bm+cm
        self.epsilon = -bm*cm - dm**2

        # self.u = 2
        # self.w = -1

        super(TB, self).__init__(T, P, mezcla)

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


    def TB_Fugacidad(self, T, P):
        """Método de cálculo de la fugacidad mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        z=self.TB_Z(T, P)
        A=a*P/R_atml**2/T**2
        B=b*P/R_atml/T
        u=1+c/b
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        tita=abs(t)**0.5
        if t>=0:
            lamda=log((2*z+B*(u-tita))/(2*z+B*(u+tita)))
        else:
            lamda=2*arctan((2*z+u*B)/B/tita)-pi
        fi=z-1-log(z-B)+A/B/tita*lamda
        return unidades.Pressure(P*exp(fi), "atm")

    def TB_U_exceso(self, T, P):
        """Método de cálculo de la energía interna de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        v=self.TB_V(T, P)
        z=P*v/R_atml/T
        A=a*P/R_atml**2/T**2
        B=b*P/R_atml/T
        u=1+c/b
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        tita=abs(t)**0.5
        if t>=0:
            lamda=log((2*z+B*(u-tita))/(2*z+B*(u+tita)))
        else:
            lamda=2*arctan((2*z+u*B)/B/tita)-pi

        delta=v**2+(b+c)*v-b*c-d**2
        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        da=-q1*a/self.Tc
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
        else:
            db=0
        U=lamda/b/tita*(a-da*T)+db*T*(-R_atml*T/(v-b)+a/b**2/t*((v*(3*c+b)-b*c+c**2-2*d**2)/delta+(3*c+b)*lamda/b/tita)) #atm*l/mol
        return unidades.Enthalpy(U*101325/1000/self.peso_molecular, "Jkg")

    def TB_H_exceso(self, T, P):
        """Método de cálculo de la entalpía de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        p=unidades.Pressure(P, "atm")
        U=self.TB_U_exceso(T, P)
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        if t>=0:
            v=self.TB_V(T, P).m3g*self.peso_molecular
            return unidades.Enthalpy(p*v-R/T/self.peso_molecular+U.Jg, "Jg")
        else:
            Z=self.TB_Z(T, P)
            return unidades.Enthalpy(R/self.peso_molecular*T*(Z-1)+U.Jg, "Jg")

    def TB_Entalpia(self, T, P):
        """Método de cálculo de la entalpía mediante la ecuación de estado de Trebble-Bishnoi"""
        Ho=self._Ho(T)
        Delta=self.TB_H_exceso(T, P)
        return unidades.Enthalpy(Delta+Ho)

    def TB_S_exceso(self, T, P):
        """Método de cálculo de la entropía de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        H=self.TB_H_exceso(T, P)
        f=self.TB_Fugacidad(T, P)
        return unidades.SpecificHeat(H.Jg/T-R/self.peso_molecular*log(f.atm/P), "JgK")

    def TB_Entropia(self, T, P):
        """Método de cálculo de la entropía mediante la ecuación de estado de Trebble-Bishnoi"""
        So=self._so(T)
        Delta=self.TB_S_exceso(T, P)
        return unidades.SpecificHeat(Delta+So)

    def TB_Cv_exceso(self, T, P):
        """Método de cálculo de la capacidad calorífica a volumen constante de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        v=self.TB_V(T, P)
        z=P*v/R_atml/T
        t=1+6*c/b+c**2/b**2+4*d**2/b**2
        tita=abs(t)**0.5
        A=a*P/R_atml**2/T**2
        B=b*P/R_atml/T
        u=1+c/b
        delta=v**2+(b+c)*v-b*c-d**2
        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        da=-q1*a/self.Tc
        dda=q1**2*a/self.Tc**2
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
            ddb=-q2*b/beta/T**2
        else:
            db=0
            ddb=0

        dt=-db/b**2*(6*c+2*c**2/b+8*d**2/b)
        dtita=abs(dt)/20
        if t>=0:
            lamda=log((2*z+B*(u-tita))/(2*z+B*(u+tita)))
            dlamda=(db-db*tita-b*dtita)/(2*v+b+c-b*tita)-(db+db*tita+b*dtita)/((2*v+b+c+b*tita))
        else:
            lamda=2*arctan((2*z+u*B)/B/tita)-pi
            dlamda=2/(1+((2*v+b+c)/b/tita)**2)*(db/b/tita-(2*v+b+c)*(db/b**2/tita+dtita/b/tita**2))

        Cv=1/b/tita*(dlamda*(a-da*T)-lamda*dda*T-lamda*(a-da*T)*(db/b+dtita/tita))+(ddb*T+db)*(-R_atml*T/(v-b)+a/b**2/t*((v*(3*c+b)-b*c+c**2-2*d**2)/delta+(3*c+b)*lamda/b/tita))+db*T*(-R_atml/(v-b)-R_atml*T*db/(v-b)**2+1/b**2/t*(da-2*a*db/b-a*dt/t)*((v*(3*c+b)-b*c+c**2-2*d**2)/delta+(3*c+b)*lamda/b/tita)+a/b**2/t*(db*(v-c)*(v**2-2*c*v-c**2+d**2)/delta**2+db*lamda/b/tita+(3*c+b)/b/tita*(dlamda-lamda*(db/b+dtita/tita))))
        return unidades.SpecificHeat(Cv*101325/1000/self.peso_molecular, "JkgK")

    def TB_Cv(self, T, P):
        """Método de cálculo de la capacidad calorifica a volumen constante mediante la ecuación de estado de Trebble-Bishnoi"""
        Cvo=self.Cv_ideal(T)
        Delta=self.TB_Cv_exceso(T, P)
        return unidades.SpecificHeat(Delta+Cvo)

    def TB_Cp_exceso(self, T, P):
        """Método de cálculo de la capacidad calorífica a presión constante de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        a, b, c, d, q1, q2=self.TB_lib(T, P)
        v=self.TB_V(T, P)
        Cv=self.TB_Cv_exceso(T, P)

        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        delta=v**2+(b+c)*v-b*c-d**2
        da=-q1*a/self.Tc
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
        else:
            db=0

        dpdt=R_atml/(v-b)+R*T*db/(v-b)**2-da/delta+a*(v-c)*db/delta**2
        dpdv=-R_atml*T/(v-b)**2+a*(2*v+b+c)/delta**2
        Cp=-R-T*dpdt**2/dpdv*101325/1000+Cv.JgK*self.peso_molecular

        return unidades.SpecificHeat(Cp/self.peso_molecular, "JgK")

    def TB_Cp(self, T, P):
        """Método de cálculo de la capacidad calorifica a presión constante mediante la ecuación de estado de Trebble-Bishnoi"""
        Cpo=self.Cp_ideal(T)
        Delta=self.TB_Cp_exceso(T, P)
        return unidades.SpecificHeat(Delta+Cpo)

    def TB_Joule_Thomson(self, T, P):
        """Método de cálculo del coeficiente de Joule-Thomson mediante la ecuación de estado de Trebble-Bishnoi"""
        v=self.TB_V(T, P)
        Cp=self.TB_Cp(T, P)
        a, b, c, d, q1, q2=self.TB_lib(T, P)

        beta=1.+q2*(1-self.tr(T)+log(self.tr(T)))
        delta=v**2+(b+c)*v-b*c-d**2
        da=-q1*a/self.Tc
        if self.tr(T)<=1.0:
            db=b/beta*(1/T-1/self.Tc)
        else:
            db=0
        dpdt=R_atml/(v-b)+R*T*db/(v-b)**2-da/delta+a*(v-c)*db/delta**2
        dpdv=-R_atml*T/(v-b)**2+a*(2*v+b+c)/delta**2

        return -(T*dpdt+v*dpdv)/Cp/dpdv




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

    def _alfa(self, cmp, T):
        if cmp.id == 62 and T/cmp.Tc < 0.85:
            # Special case for water from [2]_
            alfa = (1.0085677 + 0.82154*(1-(T/cmp.Tc)**0.5))**2         # Eq 6
            m = 0
        else:
            m = 0.37464 + 1.54226*cmp.f_acent - 0.2699*cmp.f_acent**2   # Eq 18
            alfa = (1+m*(1-(T/cmp.Tc)**0.5))**2                         # Eq 17
        return m, alfa



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
    eq = TB(300, 9.9742e5, mix)
    print('%0.0f %0.1f' % (eq.Vg.ccmol, eq.Vl.ccmol))
    eq = TB(300, 42.477e5, mix)
    print('%0.1f' % (eq.Vl.ccmol))

