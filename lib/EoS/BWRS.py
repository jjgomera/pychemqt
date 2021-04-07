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



###############################################################################
# Library to add EoS common functionality
###############################################################################

from numpy import exp
from numpy.lib.scimath import log
from scipy import sqrt
from scipy import r_
from scipy.optimize import fsolve

from lib import unidades
from lib.eos import EoS
from lib.EoS.Cubic import SRK
from lib.physics import R_atml
from lib.bip import Kij


class BWRS(EoS):
    """Ecuación de estado de Benedict-Webb-Rubin-Starling"""
    __title__ = "Benedict-Webb-Rubin-Starling"
    __status__ = "BWRS"

    def __init__(self, T, P, mezcla):
        self.T = unidades.Temperature(T)
        self.P = unidades.Pressure(P)
        self.componente = mezcla.componente
        self.zi = mezcla.fraccion
        self.kij = Kij("BWRS")

        Aoi = []
        Boi = []
        Coi = []
        Doi = []
        Eoi = []
        ai = []
        bi = []
        ci = []
        di = []
        alfai = []
        gammai = []
        for compuesto in self.componente:
            Ao_, Bo_, Co_, Do_, Eo_, a_, b_, c_, d_, alfa_, gamma_ = self._lib(compuesto)
            Aoi.append(Ao_)
            Boi.append(Bo_)
            Coi.append(Co_)
            Doi.append(Do_)
            Eoi.append(Eo_)
            ai.append(a_)
            bi.append(b_)
            ci.append(c_)
            di.append(d_)
            alfai.append(alfa_)
            gammai.append(gamma_)

        Ao = Co = Do = Eo = Bo = a = b = c = d = alfa = gamma = 0
        for i in range(len(self.componente)):
            Bo += self.zi[i]*Boi[i]
            a += self.zi[i]*ai[i]**(1./3)
            b += self.zi[i]*bi[i]**(1./3)
            c += self.zi[i]*ci[i]**(1./3)
            d += self.zi[i]*di[i]**(1./3)
            alfa += self.zi[i]*alfai[i]**(1./3)
            gamma += self.zi[i]*gammai[i]**0.5
        a = a**3
        b = b**3
        c = c**3
        d = d**3
        alfa = alfa**3
        gamma = gamma**2

        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
                Ao += self.zi[i]*self.zi[j]*Aoi[i]**0.5*Aoi[j]**0.5*(1-self.kij[i][j])
                Co += self.zi[i]*self.zi[j]*Coi[i]**0.5*Coi[j]**0.5*(1-self.kij[i][j])**3
                Do += self.zi[i]*self.zi[j]*Doi[i]**0.5*Doi[j]**0.5*(1-self.kij[i][j])**4
                Eo += self.zi[i]*self.zi[j]*Eoi[i]**0.5*Eoi[j]**0.5*(1-self.kij[i][j])**5

        self.Aoi = Aoi
        self.Boi = Boi
        self.Coi = Coi
        self.Doi = Doi
        self.Eoi = Eoi
        self.ai = ai
        self.bi = bi
        self.ci = ci
        self.di = di
        self.alfai = alfai
        self.gammai = gammai
        self.Ao = Ao
        self.Co = Co
        self.Do = Do
        self.Eo = Eo
        self.Bo = Bo
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.alfa = alfa
        self.gamma = gamma

        Vm = lambda V: self.P.atm-R_atml*self.T/V - (Bo*R_atml*self.T-Ao-Co/self.T**2+Do/self.T**3-Eo/self.T**4)/V**2 - (b*R_atml*self.T-a-d/self.T)/V**3 - alfa*(a+d/self.T)/V**6 - c/self.T**2/V**3*(1+gamma/V**2)*exp(-gamma/V**2)

        # Usamos SRK para estimar los volumenes de ambas fases usados como valores iniciales en la iteeración
        srk=SRK(T, P, mezcla)
        Z_srk=srk.Z
        Vgo=Z_srk[0]*R_atml*T/P
        Vlo=Z_srk[1]*R_atml*T/P
        Vg=fsolve(Vm, Vgo)
        Vl=fsolve(Vm, Vlo)
        self.V=r_[Vg, Vl]       #mol/l
        self.Z=P*self.V/R_atml/T

        self.H_exc=(Bo*R_atml*self.T-2*Ao-4*Co/self.T**2+5*Do/self.T**3-6*Eo/self.T**4)/self.V+(2*b*R_atml*self.T-3*a-4*d/self.T)/2/self.V**2+alfa/5*(6*a+7*d/self.T)/self.V**5+c/gamma/self.T**2*(3-(3+gamma/self.V**2/2-gamma**2/self.V**4)*exp(-gamma/self.V**2))

        self.x, self.xi, self.yi, self.Ki=self._Flash()

    def _fug(self, Z, xi):
        rho=self.P.atm/Z/R_atml/self.T
        tita=[]
        for i in range(len(self.componente)):
            suma=0
            for j in range(len(self.componente)):
                suma+=xi[j]*(-self.Ao**0.5*self.Aoi[i]**0.5*(1-self.kij[i][j]) \
                                    -  self.Co**0.5*self.Coi[i]**0.5*(1-self.kij[i][j])**3/self.T**2 \
                                    + self.Do**0.5*self.Doi[i]**0.5*(1-self.kij[i][j])**4/self.T**3 \
                                    -  self.Eo**0.5*self.Eoi[i]**0.5*(1-self.kij[i][j])**5/self.T**4)
#                print suma
            lo=R_atml*self.T*log(rho*R_atml*self.T*xi[i]) \
                    + rho*(self.Bo+self.Boi[i])*R_atml*self.T \
                    + 2*rho*suma \
                    + rho**2/2*(3*(self.b**2*self.bi[i])**(1./3)*R_atml*self.T-3*(self.a**2*self.ai[i])**(1./3)-3*(self.d**2*self.di[i])**(1./3)/self.T) \
                    + self.alfa*rho**5/5*(3*(self.a**2*self.ai[i])**(1./3)+3*(self.d**2*self.di[i])**(1./3)/self.T) \
                    + 3*rho**5/5*(self.a+self.d/self.T)*(self.alfa**2*self.alfai[i])**(1./3) \
                    + 3*(self.c**2*self.ci[i])**(1./3)*rho**2/self.T**2*((1-exp(-self.gamma*rho**2))/self.gamma/rho**2-exp(-self.gamma*rho**2)/2) \
                    - (2*self.c*sqrt(self.gammai[i]/self.gamma)**0.5/self.gamma/self.T**2)*((1-exp(-self.gamma*rho**2))*(1+self.gamma*rho**2+self.gamma**2*rho**4/2))
            tita.append(exp(lo/R_atml/self.T))
        return tita

    def _lib(self, cmp):
        """Library for parameter calculation in Benedict-Webb-Rubin-Starling"""
        a1, b1 = 0.443690, 0.115449
        a2, b2 = 1.28438, -0.920731
        a3, b3 = 0.356306, 1.70871
        a4, b4 = 0.544979, -0.270896
        a5, b5 = 0.528629, 0.349261
        a6, b6 = 0.484011, 0.754130
        a7, b7 = 0.0705233, -0.044448
        a8, b8 = 0.504087, 1.32245
        a9, b9 = 0.0307452, 0.179433
        a10, b10 = 0.0732828, 0.463492
        a11, b11 = 0.006450, -0.022143

        Bo = (a1+b1*cmp.f_acent)*cmp.Vc
        Ao = (a2+b2*cmp.f_acent)*R_atml*cmp.Tc*cmp.Vc
        Co = (a3+b3*cmp.f_acent)*R_atml*cmp.Tc**3*cmp.Vc
        gamma = (a4+b4*cmp.f_acent)*(cmp.Vc)**2
        b = (a5+b5*cmp.f_acent)*(cmp.Vc)**2
        a = (a6+b6*cmp.f_acent)*R_atml*cmp.Tc*(cmp.Vc)**2
        alfa = (a7+b7*cmp.f_acent)*(cmp.Vc)**3
        c = (a8+b8*cmp.f_acent)*R_atml*cmp.Tc**3*(cmp.Vc)**2
        Do = (a9+b9*cmp.f_acent)*R_atml*cmp.Tc**4*cmp.Vc
        d = (a10+b10*cmp.f_acent)*R_atml*cmp.Tc**2*(cmp.Vc)**2
        Eo = (a11+b11*cmp.f_acent*exp(-3.8*cmp.f_acent))*R_atml*cmp.Tc**5*cmp.Vc

        return Ao, Bo, Co, Do, Eo, a, b, c, d, alfa, gamma


_all = [BWRS]


if __name__ == "__main__":
    from lib.mezcla import Mezcla
    mezcla = Mezcla(2, ids=[10, 38, 22, 61], caudalUnitarioMolar=[0.3, 0.5, 0.05, 0.15])
    eq = BWRS(340, 101325, mezcla)
    print(eq.x)
