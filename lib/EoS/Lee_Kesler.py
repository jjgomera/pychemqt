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


from numpy import exp
from numpy.lib.scimath import log
from scipy import zeros
from scipy.constants import R

from lib import unidades
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



    # def Cp_Lee_Kesler(self, T, P, fase=None):
        # """Método alternativo para el cálculo de la capacidad calorífica
        # Procedure API 7D3.6 Pag.711"""
        # Tr=self.tr(T)
        # if fase==None:
            # fase=self.Fase(T, P)
        # Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(T, P, fase)

        # B=0.1181193-0.265728/Tr-0.154790/Tr**2-0.030323/Tr**3
        # C=0.0236744-0.0186984/Tr
        # D=0.155488e-4+0.623689e-4/Tr
        # dpdt_0=1/vr0*(1+(0.1181193+0.154790/Tr**2+2*0.030323/Tr**3)/vr0+0.0236744/vr0**2+0.155488e-4/vr0**5-2*0.042724/Tr**3/vr0**2*((0.65392+0.060167/vr0**2)*exp(-0.060167/vr0**2)))
        # dpdv_0=-Tr/vr0**2*(1+2*B/vr0+3*C/vr0**2+6*D/vr0**5+0.042724/Tr**3/vr0**2*(3*0.65392+(5-2*(0.65392+0.060167/vr0**2))*0.060167/vr0**2)*exp(-0.060167/vr0**2))
        # Cp0=1+Tr*dpdt_0**2/dpdv_0+Cv0

        # B=0.2026579-0.331511/Tr-0.027655/Tr**2-0.203488/Tr**3
        # C=0.0313385-0.0503618/Tr+0.016901/Tr**3
        # D=0.48736e-4+0.0740336e-4/Tr
        # dpdt_h=1/vrh*(1+(0.2026579+0.027655/Tr**2+2*0.203488/Tr**3)/vrh+(0.0313385-2*0.016901/Tr**3)/vrh**2+0.48736e-4/vrh**5-2*0.041577/Tr**3/vrh**2*((1.226+0.03754/vrh**2)*exp(-0.03754/vrh**2)))
        # dpdv_h=-Tr/vrh**2*(1+2*B/vrh+3*C/vrh**2+6*D/vrh**5+0.041577/Tr**3/vrh**2*(3*1.226+(5-2*(1.226+0.03754/vrh**2))*0.03754/vrh**2)*exp(-0.03754/vrh**2))
        # Cph=1+Tr*dpdt_h**2/dpdv_h+Cvh

        # Cp_adimensional=Cp0+self.f_acent/factor_acentrico_octano*(Cph-Cp0)
        # return unidades.SpecificHeat(self._Cpo(T).JgK-R/self.M*Cp_adimensional, "JgK")

    # def Cv_Lee_Kesler(self, T, P, fase=None):
        # """Método de cálculo de la capacidad calorífica a volumen constante
        # Procedure API 7E1.6 Pag.726"""
        # #FIXME: No sale, un factor de 100 tengo que añadir no sé de donde
        # Pr=P/self.Pc
        # Tr=T/self.Tc
        # if fase==None:
            # fase=self.Fase(T, P)
        # Cpo=self._Cpo(T)
        # Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(Tr, Pr, fase)
        # Cv_adimensional=Cv0+self.f_acent/factor_acentrico_octano*(Cvh-Cv0)
        # return unidades.SpecificHeat(100*(Cpo.JgK-R/self.M*(1+Cv_adimensional)), "JgK")


    # def Cp_Cv_Lee_Kesler(self, T, P):
        # """Método de cálculo de la capacidad calorífica a volumen constante
        # Procedure API 7E1.6 Pag.726"""
        # Cv=self.Cv_Lee_Kesler(T, P.atm)
        # Cp=self.Cp_Lee_Kesler(T, P.atm)
# #        print Cp.BtulbF, Cv
        # return Cp/Cv


    # def Fugacidad_Lee_Kesler(self, T, P):
        # """Método de cálculo de la fugacidad
        # Procedure API 7G1.8 Pag.752"""
        # Tr=T/self.Tc
        # Pr=P/self.Pc
        # f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P))
        # return unidades.Pressure(P*exp(f), "atm")

    # def Entropia_Lee_Kesler(self, T, P):
        # """Método de cálculo de la entropia
        # Procedure API 7F1.7 Pag.739"""
        # Tr=T/self.Tc
        # Pr=P/self.Pc
        # S0=self._so(T)
        # H_adimensional=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        # f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        # S=H_adimensional+f+log(P/101325)

        # return unidades.SpecificHeat(S0.JgK-R*S/self.M, "JgK")

    # def Hv_Lee_Kesler(self, T):
        # """Método alternativo para el cálculo del calor de vaporización haciendo uso de las propiedades críticas
        # Procedure API 7C1.16 Pag.680
        # Valor en J/mol"""
        # Pv=self.Pv_DIPPR(T)
        # Tr=T/self.Tc
        # Pr=Pv/self.Pc
        # H_adimensional_vapor=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 1)
        # H_adimensional_liquido=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 0)
        # return unidades.Enthalpy(R*self.Tc/self.M*(H_adimensional_vapor-H_adimensional_liquido), "Jg")


    # def RhoG_Lee_Kesler(self, T, P):
        # a, b=eos.SRK_lib(self, T)
        # Z_srk=eos.Z_Cubic_EoS(T, P, b, a, b, 0, b)
        # Vvo=Z_srk[0]*R_atml*T/P

        # vr0v, vrhv, vr0l, vrhl=eos.Lee_Kesler_lib(T/self.Tc, P/self.Pc.atm, fase=1, Vvo=Vvo)
        # z0v=P/self.Pc.atm*vr0v/T*self.Tc
        # zhv=P/self.Pc.atm*vrhv/T*self.Tc
        # z=z0v+self.f_acent/factor_acentrico_octano*(zhv-z0v)
        # return P/z/R_atml/T


_all = [Lee_Kesler]

if __name__ == "__main__":
    from lib.corriente import Mezcla
    mezcla = Mezcla(1, ids=[98], caudalUnitarioMasico=[1.])
    for T in [125, 135, 145, 165, 185, 205]:
        eq = Virial(T, 1, mezcla)
        print(eq.H_exc)
