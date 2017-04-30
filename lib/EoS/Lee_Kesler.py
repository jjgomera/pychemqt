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
# Virial equation of state implementation
###############################################################################

from scipy import zeros, log, exp
from scipy.constants import R

from lib.eos import EoS
from lib.physics import R_atml

class Lee_Kesler(EoS):
    """Ecuación de estado de Lee-Kesler"""
    __title__="Lee Kesler"
    __status__="LK"

    def __init__(self, T, P, mezcla):
        self.T=unidades.Temperature(T)
        self.P=unidades.Pressure(P, "atm")
        self.componente=mezcla.componente
        self.fraccion=mezcla.fraccion

        zci=[]
        Vci=[]
        for componente in self.componente:
            zci.append(0.2905-0.085*componente.f_acent)
            Vci.append(zci[-1]*R_atml*componente.Tc/componente.Pc.atm)

        sumaV1=sumaV2=sumaV3=sumaT1=sumaT2=sumaT3=0
        for i, componente in enumerate(self.componente):
            sumaV1+=self.fraccion[i]*Vci[i]
            sumaV2+=self.fraccion[i]*Vci[i]**(2./3)
            sumaV3+=self.fraccion[i]*Vci[i]**(1./3)
            sumaT1+=self.fraccion[i]*Vci[i]*componente.Tc
            sumaT2+=self.fraccion[i]*Vci[i]**(2./3)*componente.Tc**(1./2)
            sumaT3+=self.fraccion[i]*Vci[i]**(1./3)*componente.Tc**(1./2)

        Vmc=(sumaV1+3*sumaV2*sumaV3)/4.
        Tmc=(sumaT1+3*sumaT2*sumaT3)/4/Vmc
        Pmc=(0.2905-0.085*mezcla.f_acent)*R_atml*Tmc/Vmc

        Tr=T/Tmc
        Pr=P/Pmc

        b1=0.1181193, 0.2026579
        b2=0.265728, 0.331511
        b3=0.154790, 0.027655
        b4=0.030323, 0.203488
        c1=0.0236744, 0.0313385
        c2=0.0186984, 0.0503618
        c3=0.0, 0.016901
        c4=0.042724, 0.041577
        d1=0.155488e-4, 0.48736e-4
        d2=0.623689e-4, 0.0740336e-4
        beta=0.65392, 1.226
        gamma=0.060167, 0.03754

        Bo=b1[0]-b2[0]/Tr-b3[0]/Tr**2-b4[0]/Tr**3
        Co=c1[0]-c2[0]/Tr+c3[0]/Tr**3
        Do=d1[0]+d2[0]/Tr
        Vr=lambda V: 1+Bo/V+Co/V**2+Do/V**5+c4[0]/Tr**3/V**2*(beta[0]+gamma[0]/V**2)*exp(-gamma[0]/V**2)-Pr*V/Tr

        Bh=b1[1]-b2[1]/Tr-b3[1]/Tr**2-b4[1]/Tr**3
        Ch=c1[1]-c2[1]/Tr+c3[1]/Tr**3
        Dh=d1[1]+d2[1]/Tr
        Vrh=lambda V: 1+Bh/V+Ch/V**2+Dh/V**5+c4[1]/Tr**3/V**2*(beta[1]+gamma[1]/V**2)*exp(-gamma[1]/V**2)-Pr*V/Tr

        #Usamos SRK para estimar los volumenes de ambas fases usados como valores iniciales en la iteración
        srk=SRK(T, P, mezcla)
        Z_srk=srk.Z
        Vgo=Z_srk[0]*R_atml*T/P
        Vlo=Z_srk[1]*R_atml*T/P

        vr0v=fsolve(Vr, Vgo)
        vrhv=fsolve(Vrh, Vgo)
        vr0l=fsolve(Vr, Vlo)
        vrhl=fsolve(Vrh, Vlo)

        z0l=Pr*vr0l/Tr
        zhl=Pr*vrhl/Tr
        z0v=Pr*vr0v/Tr
        zhv=Pr*vrhv/Tr
        self.Z=r_[z0v+mezcla.f_acent/factor_acentrico_octano*(zhv-z0v), z0l+mezcla.f_acent/factor_acentrico_octano*(zhl-z0l)]
        self.V=self.Z*R_atml*self.T/self.P.atm  #mol/l

        E=c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0v**2)*exp(-gamma[0]/vr0v**2))
        H0=-Tr*(z0v-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0v-c2[0]/Tr/2/vr0v**2+d2[0]/5/Tr/vr0v**5+3*E)
        E=c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhv**2)*exp(-gamma[1]/vrhv**2))
        Hh=-Tr*(zhv-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhv-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhv**2+d2[1]/5/Tr/vrhv**5+3*E)
        Hv=H0+mezcla.f_acent/factor_acentrico_octano*(Hh-H0)

        E=c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0l**2)*exp(-gamma[0]/vr0l**2))
        H0=-Tr*(z0l-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0l-c2[0]/Tr/2/vr0l**2+d2[0]/5/Tr/vr0l**5+3*E)
        E=c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhl**2)*exp(-gamma[1]/vrhl**2))
        Hh=-Tr*(zhl-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhl-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhl**2+d2[1]/5/Tr/vrhl**5+3*E)
        Hl=H0+mezcla.f_acent/factor_acentrico_octano*(Hh-H0)
        self.H_exc=r_[Hv, Hl]
        self.x, self.xi, self.yi, self.Ki=srk._Flash()


_all = [Lee_Kesler]

if __name__ == "__main__":
    from lib.corriente import Mezcla
    mezcla = Mezcla(1, ids=[98], caudalUnitarioMasico=[1.])
    for T in [125, 135, 145, 165, 185, 205]:
        eq = Virial(T, 1, mezcla)
        print(eq.H_exc)
