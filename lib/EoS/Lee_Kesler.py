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


###############################################################################
# Lee-Kesler equation of state implementation
###############################################################################


from numpy import exp, r_
from numpy.lib.scimath import log
from scipy.constants import R
from scipy.optimize import fsolve

from lib.compuestos import RhoL_Costald

from lib.eos import EoS
from lib.physics import R_atml, factor_acentrico_octano


class Lee_Kesler(EoS):
    r"""
    Corresponding state equation of state of Lee-Kesler

    .. math::
        \begin{array}[t]{l}
        Z = Z^{(0)} + \omega Z^{(1)}\\
        Z = \frac{P_rV_r}{T_r} = 1 + \frac{B}{V_r} + \frac{C}{V_r^2} +
        \frac{D}{V_r^5} + \frac{c_4}{T_r^3V_r^2}\left(\beta+\frac{\gamma}
        {V_r^2}\right)\exp{\left(-\frac{\gamma}{V_r^2}\right)}\\
        B = b_1 - \frac{b_2}{T_r} - \frac{b_3}{T_r^2}-\frac{b_4}{T_r^3}\\
        C = c_1 - \frac{c_2}{T_r} - \frac{c_3}{T_r^3}\\
        D = d_1 - \frac{d_2}{T_r}\\
        \end{array}

    """
    __title__ = "Lee Kesler"
    __status__ = "LK"

    __doi__ = (
        {"autor": "Lee, B.I., Kesler, M.G.",
         "title": "A Generalized Thermodynamic Correlation Based on "
                  "Three-Parameter Corresponding States",
         "ref": "AIChE Journal 21(3) (1975) 510-527",
         "doi": "10.1002/aic.690210313"},
        {"autor": "Joffe, J.",
         "title": "Vapor-Liquid Equilibria by the Pseudocritical Method",
         "ref": "Ind. Eng. Chem. Fundam. 15(4) (1976) 298-303",
         "doi": "10.1021/i160060a013"})

    def __init__(self, T, P, mezcla):
        EoS.__init__(self, T, P, mezcla)

        # Mixing rules, Eq 20-25
        Zci = [0.2905-0.085*cmp.f_acent for cmp in self.componente]
        Vci = [zc*R*cmp.Tc/cmp.Pc.kPa for zc, cmp in zip(Zci, self.componente)]
        Tci = [cmp.Tc for cmp in self.componente]
        sumV = 0
        sumT = 0
        for xj, Vcj, Tcj in zip(self.zi, Vci, Tci):
            for xk, Vck, Tck in zip(self.zi, Vci, Tci):
                sumV += xj*xk*(Vcj**(1/3)+Vck**(1/3))**3
                sumT += xj*xk*(Vcj**(1/3)+Vck**(1/3))**3*(Tcj*Tck)**0.5
        Vc = sumV/8
        Tc = sumT/8/Vc
        Pc = (0.2905-0.085*mezcla.f_acent)*R_atml*Tc/Vc

        Tr = T/Tc
        Pr = P/Pc

        # Table 1
        b1 = 0.1181193, 0.2026579
        b2 = 0.265728, 0.331511
        b3 = 0.154790, 0.027655
        b4 = 0.030323, 0.203488
        c1 = 0.0236744, 0.0313385
        c2 = 0.0186984, 0.0503618
        c3 = 0.0, 0.016901
        c4 = 0.042724, 0.041577
        d1 = 0.155488e-4, 0.48736e-4
        d2 = 0.623689e-4, 0.0740336e-4
        beta = 0.65392, 1.226
        gamma = 0.060167, 0.03754

        Bo = b1[0]-b2[0]/Tr-b3[0]/Tr**2-b4[0]/Tr**3
        Co = c1[0]-c2[0]/Tr+c3[0]/Tr**3
        Do = d1[0]+d2[0]/Tr

        def Vr(V):
            Vr = 1 + Bo/V + Co/V**2 + Do/V**5 + c4[0]/Tr**3/V**2 * \
                (beta[0]+gamma[0]/V**2) * exp(-gamma[0]/V**2)-Pr*V/Tr
            return Vr

        Bh = b1[1]-b2[1]/Tr-b3[1]/Tr**2-b4[1]/Tr**3
        Ch = c1[1]-c2[1]/Tr+c3[1]/Tr**3
        Dh = d1[1]+d2[1]/Tr

        def Vrh(V):
            Vrh = 1 + Bh/V + Ch/V**2 + Dh/V**5 + c4[1]/Tr**3/V**2 * \
                (beta[1]+gamma[1]/V**2) * exp(-gamma[1]/V**2)-Pr*V/Tr
            return Vrh

        # Used initial values for iteration
        Vlo = RhoL_Costald(T, Tc, mezcla.f_acent, Vc)
        Vgo = R_atml*T/P

        vr0v = fsolve(Vr, Vgo)
        vrhv = fsolve(Vrh, Vgo)
        vr0l = fsolve(Vr, Vlo)
        vrhl = fsolve(Vrh, Vlo)

        z0l = Pr*vr0l/Tr
        zhl = Pr*vrhl/Tr
        z0v = Pr*vr0v/Tr
        zhv = Pr*vrhv/Tr
        self.Z = r_[z0v+mezcla.f_acent/factor_acentrico_octano*(zhv-z0v), z0l+mezcla.f_acent/factor_acentrico_octano*(zhl-z0l)]
        self.V = self.Z*R_atml*self.T/self.P.atm  #mol/l

        E = c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0v**2)*exp(-gamma[0]/vr0v**2))
        H0 = -Tr*(z0v-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0v-c2[0]/Tr/2/vr0v**2+d2[0]/5/Tr/vr0v**5+3*E)
        E = c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhv**2)*exp(-gamma[1]/vrhv**2))
        Hh = -Tr*(zhv-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhv-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhv**2+d2[1]/5/Tr/vrhv**5+3*E)
        Hv = H0+mezcla.f_acent/factor_acentrico_octano*(Hh-H0)

        E = c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0l**2)*exp(-gamma[0]/vr0l**2))
        H0 = -Tr*(z0l-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0l-c2[0]/Tr/2/vr0l**2+d2[0]/5/Tr/vr0l**5+3*E)
        E = c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhl**2)*exp(-gamma[1]/vrhl**2))
        Hh = -Tr*(zhl-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhl-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhl**2+d2[1]/5/Tr/vrhl**5+3*E)
        Hl = H0+mezcla.f_acent/factor_acentrico_octano*(Hh-H0)
        self.H_exc = r_[Hv, Hl]
#         self.x, self.xi, self.yi, self.Ki = srk._Flash()


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

def Lee_Kesler_lib(Tr, Pr, fase=2, Vvo=0.0001, Vlo=5):
    """Librería para el cálculo de la EoS de Lee-Kesler
    Procedure API 6B1.8 pag 518
    Perry pag 2-358
    fase: fase para la que se realiza el cálculo
        0   -   Liquido
        1   -   Vapor
        2   -   Ambas"""

    Bo=b1[0]-b2[0]/Tr-b3[0]/Tr**2-b4[0]/Tr**3
    Co=c1[0]-c2[0]/Tr+c3[0]/Tr**3
    Do=d1[0]+d2[0]/Tr
    def Vr(V):
        return 1+Bo/V+Co/V**2+Do/V**5+c4[0]/Tr**3/V**2*(beta[0]+gamma[0]/V**2)*exp(-gamma[0]/V**2)-Pr*V/Tr

    Bh=b1[1]-b2[1]/Tr-b3[1]/Tr**2-b4[1]/Tr**3
    Ch=c1[1]-c2[1]/Tr+c3[1]/Tr**3
    Dh=d1[1]+d2[1]/Tr
    def Vrh(V):
        return 1+Bh/V+Ch/V**2+Dh/V**5+c4[1]/Tr**3/V**2*(beta[1]+gamma[1]/V**2)*exp(-gamma[1]/V**2)-Pr*V/Tr

    vr0v=vr0l=vrhv=vrhl=None
    if fase!=0:
        vr0v=fsolve(Vr, Vvo)
        vrhv=fsolve(Vrh, Vvo)
    elif fase!=1:
        vr0l=fsolve(Vr, Vlo)
        vrhl=fsolve(Vrh, Vlo)

    return vr0v, vrhv, vr0l, vrhl


def Z_Lee_Kesler(T, P, mezcla):
    """Factor de compresibilidad según la ecuación de estado de Lee-Kesler"""
    a, b, ai, bi, tdadt=SRK_lib(mezcla)
    Z_srk=Z_Cubic_EoS(b, a, b, 0, b)
#    self.titail, self.titaiv=self.Fugacidad_Cubic_EoS(Z_srk, b, a, ai, bi, 1, 0)

    Vvo=Z_srk[0]*R_atml*self.T/self.P.atm
    Vlo=Z_srk[1]*R_atml*self.T/self.P.atm

    Tmc, Pmc, f_acent, Vmc=Mix_Lee_Kesler(mezcla.fraccion, mezcla.componentes, mezcla.f_acent)
    vr0v, vrhv, vr0l, vrhl=Lee_Kesler_lib(T/Tmc, P/Pmc, Vvo, Vlo)
    z0l=self.P.atm/Pmc*vr0l/self.T*Tmc
    zhl=self.P.atm/Pmc*vrhl/self.T*Tmc
    z0v=self.P.atm/Pmc*vr0v/self.T*Tmc
    zhv=self.P.atm/Pmc*vrhv/self.T*Tmc

    return z0v+f_acent/factor_acentrico_octano*(zhv-z0v), z0l+f_acent/factor_acentrico_octano*(zhl-z0l)


def Lee_Kesler_lib_Cp(Tr, Pr, fase=1):
    """Librería para el cálculo de capacidades calorificas, usada a continuación en diferentes funciones
    Procedure API 7E1.6 Pag.726"""
    #FIXME: No concuerdan mucho los valores de cp y cv con los valores por DIPPR
    vr0v, vrhv, vr0l, vrhl=Lee_Kesler_lib(Tr, Pr, fase)
    if fase:
        E=c4[0]/2/Tr**3/gamma[0]*(beta[0]+1-(beta[0]+1+gamma[0]/vr0v**2)*exp(-gamma[0]/vr0v**2))
        Cv0=-2*(b3[0]+3*b4[0]/Tr)/Tr**2/vr0v+3*c3[0]/Tr**3/vr0v**2+6*E
        E=c4[1]/2/Tr**3/gamma[1]*(beta[1]+1-(beta[1]+1+gamma[1]/vrhv**2)*exp(-gamma[1]/vrhv**2))
        Cvh=-2*(b3[1]+3*b4[1]/Tr)/Tr**2/vrhv+3*c3[1]/Tr**3/vrhv**2+6*E
        vr0=vr0v
        vrh=vrhv
    else:
        E=c4[0]/2/Tr**3/gamma[0]*(beta[0]+1-(beta[0]+1+gamma[0]/vr0l**2)*exp(-gamma[0]/vr0l**2))
        Cv0=-2*(b3[0]+3*b4[0]/Tr)/Tr**2/vr0l+3*c3[0]/Tr**3/vr0l**2+6*E
        E=c4[1]/2/Tr**3/gamma[1]*(beta[1]+1-(beta[1]+1+gamma[1]/vrhl**2)*exp(-gamma[1]/vrhl**2))
        Cvh=-2*(b3[1]+3*b4[1]/Tr)/Tr**2/vrhl+3*c3[1]/Tr**3/vrhl**2+6*E
        vr0=vr0l
        vrh=vrhl

    return Cv0, Cvh, vr0, vrh


def Lee_Kesler_Entalpia_lib(Tr, Pr, w, fase=1):
    """Librería para el cálculo del factor adimensional de influencia de la presión sobre la temperatura.
    Usado en diversos métodos a continuación
    eq 7B3.7-1 pag 643"""

    vr0v, vrhv, vr0l, vrhl=Lee_Kesler_lib(Tr, Pr, fase)
    if fase:
        z0=Pr*vr0v/Tr
        zh=Pr*vrhv/Tr
        E=c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0v**2)*exp(-gamma[0]/vr0v**2))
        H0=-Tr*(z0-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0v-c2[0]/Tr/2/vr0v**2+d2[0]/5/Tr/vr0v**5+3*E)

        E=c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhv**2)*exp(-gamma[1]/vrhv**2))
        Hh=-Tr*(zh-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhv-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhv**2+d2[1]/5/Tr/vrhv**5+3*E)
        H=H0+w/factor_acentrico_octano*(Hh-H0)
    else:
        z0=Pr*vr0l/Tr
        zh=Pr*vrhl/Tr
        E=c4[0]/(2*Tr**3*gamma[0])*(beta[0]+1-(beta[0]+1+gamma[0]/vr0l**2)*exp(-gamma[0]/vr0l**2))
        H0=-Tr*(z0-1-(b2[0]+2*b3[0]/Tr+3*b4[0]/Tr**2)/Tr/vr0l-c2[0]/Tr/2/vr0l**2+d2[0]/5/Tr/vr0l**5+3*E)

        E=c4[1]/(2*Tr**3*gamma[1])*(beta[1]+1-(beta[1]+1+gamma[1]/vrhl**2)*exp(-gamma[1]/vrhl**2))
        Hh=-Tr*(zh-1-(b2[1]+2*b3[1]/Tr+3*b4[1]/Tr**2)/Tr/vrhl-(c2[1]-3*c3[1]/Tr**2)/Tr/2/vrhl**2+d2[1]/5/Tr/vrhl**5+3*E)
        H=H0+w/factor_acentrico_octano*(Hh-H0)

    return H

def Lee_Kesler_Fugacidad_lib(Tr, Pr, w, fase=1):
    """Librería para el cálculo de la fugacidad, entropia...
    Procedure API 7G1.8 Pag.752"""
    vr0, vrh=self.Lee_Kesler_lib(Tr, Pr, fase)
    z0=Pr*vr0/Tr
    B=b1[0]-b2[0]/Tr-b3[0]/Tr**2-b4[0]/Tr**3
    C=c1[0]-c2[0]/Tr+c3[0]/Tr**3
    D=d1[0]+d2[0]/Tr
    E=c4[0]/2/Tr**3/gamma[0]*(beta[0]+1-(beta[0]+1+gamma[0]/vr0**2)*exp(-gamma[0]/vr0**2))
    f0=z0-1-log(z0)+B/vr0+C/2/vr0**2+D/5/vr0**5+E

    zh=Pr*vrh/Tr
    B=b1[1]-b2[1]/Tr-b3[1]/Tr**2-b4[1]/Tr**3
    C=c1[1]-c2[1]/Tr+c3[1]/Tr**3
    D=d1[1]+d2[1]/Tr
    E=c4[1]/2/Tr**3/gamma[1]*(beta[1]+1-(beta[1]+1+gamma[1]/vrh**2)*exp(-gamma[1]/vrh**2))
    fh=zh-1-log(zh)+B/vrh+C/2/vrh**2+D/5/vrh**5+E

    return f0+w/factor_acentrico_octano*(fh-f0)

def Entalpia_Lee_Kesler(self):
    """Librería para el cálculo del factor adimensional de influencia de la presión sobre la temperatura.
    Usado en diversos métodos a continuación
    eq 7B3.7-1 pag 643"""
    Tmc, Pmc, f_acent, Vmc=self.mezcla.Mix_Lee_Kesler()

    a, b, ai, bi, tdadt=self.SRK_lib()
    Z_srk=self.Z_Cubic_EoS(b, a, b, 0, b)
    Vvo=Z_srk[0]*R_atml*self.T/self.P.atm
    Vlo=Z_srk[1]*R_atml*self.T/self.P.atm

    vr0, vrh, vr0l, vrhl=self.Lee_Kesler_lib(Tmc, Pmc, Vmc, Vvo, Vlo)
    Tr=self.T/Tmc
    Pr=self.P.atm/Pmc
    z0=Pr*vr0/self.T*Tmc
    E=0.041724/(2*Tr**3*0.060167)*(0.65392+1-(0.65392+1+0.060167/vr0**2)*exp(-0.060167/vr0**2))
    H0=-Tr*(z0-1-(0.2658728+2*0.154790/Tr+3*0.030323/Tr**2)/Tr/vr0-0.0186984/Tr/2/vr0**2+0.623689e-4/5/Tr/vr0**5+3*E)

    zh=Pr*vrh/Tr
    E=0.041577/(2*Tr**3*0.03754)*(1.226+1-(1.226+1+0.03754/vrh**2)*exp(-0.03754/vrh**2))
    Hh=-Tr*(zh-1-(0.331511+2*0.027655/Tr+3*0.203488/Tr**2)/Tr/vrh-(0.0503618-3*0.016901/Tr**2)/Tr/2/vrh**2+0.0740336e-4/5/Tr/vrh**5+3*E)
    return H0+f_acent/factor_acentrico_octano*(Hh-H0)


_all = [Lee_Kesler]

if __name__ == "__main__":
    from lib.corriente import Mezcla
    mezcla = Mezcla(1, ids=[98], caudalUnitarioMasico=[1.])
    for T in [125, 135, 145, 165, 185, 205]:
        eq = Lee_Kesler(T, 1, mezcla)
#         print(eq.H_exc)
        print(eq.Z)
