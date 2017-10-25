#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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


###############################################################################
# Virial equation of state implementation
###############################################################################

from scipy import roots, r_, log, exp, sqrt

from PyQt5.QtWidgets import QApplication

from lib import unidades, config
from lib.eos import EoS
from lib.physics import R_atml
from lib.bip import Kij


# # TODO: Añadir parametros S1,S2 a la base de datos, API databook, pag 823
# self.SRKGraboski = [0, 0]

# # TODO: Añadir parámetros, archivo /media/datos/Biblioteca/archivos/Melhem, Almeida - A data Bank of Parameters for the Attractive-Aznar Telles.pdf
# self.Melhem = [0, 0]          #Alcoholes en archivo de abajo
# self.Almeida = [0, 0]

# # TODO: Añadir parámetros, archivo /media/datos/Biblioteca/archivos/alfas.pdf
# self.Mathias = 0
# self.MathiasCopeman = [0, 0, 0]
# self.Adachi = [0, 0]
# self.Andoulakis = [0, 0, 0]
# self.Yu_Lu = [0, 0, 0]

alfa = (QApplication.translate("pychemqt", "Original"),
        "Boston-Mathias",
        "Twu",
        "Doridon")


class Cubic(EoS):
    """Clase que modela de manera generalizada las ecuaciones de estado cúbicas
    ref. Prausnick  Propiedades de gases y liquidos, pag 203"""
    def __init__(self, T, P, mezcla):
        self.T=unidades.Temperature(T)
        self.P=unidades.Pressure(P, "atm")
        self.mezcla=mezcla
        self.componente=mezcla.componente
        self.fraccion=mezcla.fraccion

        self.B=self.b*self.P.atm/R_atml/self.T
        self.Tita=self.tita*self.P.atm/(R_atml*self.T)**2

        delta=self.delta*self.P.atm/R_atml/self.T
        epsilon=self.epsilon*(self.P.atm/R_atml/self.T)**2
        eta=self.eta*self.P.atm/R_atml/self.T
        Z=roots([1, delta-self.B-1, self.Tita+epsilon-delta*(self.B+1), -epsilon*(self.B+1)-self.Tita*eta])
        self.Z=r_[Z[0].real, Z[2].real]

        self.V=self.Z*R_atml*self.T/self.P.atm  #mol/l
        self.x, self.xi, self.yi, self.Ki=self._Flash()
        self.H_exc=-(self.tita+self.dTitadT)/R_atml/self.T/(self.delta**2-4*self.epsilon)**0.5*log((2*self.V+self.delta-(self.delta**2-4*self.epsilon)**0.5)/(2*self.V+self.delta+(self.delta**2-4*self.epsilon)**0.5))+1-self.Z

    def _fug(self, Z, xi):
        Ai=[]
        for i in range(len(self.componente)):
            suma=0
            for j in range(len(self.componente)):
                suma+=self.fraccion[j]*self.ai[j]**0.5*(1-self.kij[i][j])
            Ai.append(1/self.tita*2*self.ai[i]**0.5*suma)
        tita=[]
        for i in range(len(self.componente)):
            tita.append(exp(self.bi[i]/self.b*(Z-1)-log(Z-self.B)-self.Tita/self.B/sqrt(self.u**2-4*self.w)*(Ai[i]-self.bi[i]/self.b)*log((Z+self.B/2*(self.u+sqrt(self.u**2-4*self.w)))/(Z+self.B/2*(self.u-sqrt(self.u**2-4*self.w))))).real)
        return tita


class _2ParameterCubic(Cubic):
    pass

class van_Waals(Cubic):
    """Ecuación de estado de van der Waals
        van der Waals, J.D. Over de continuiteit van den gas- en vloestof-toestand. Dissertation, Leiden University, Leiden, Niederlande, 1873."""
    __title__="van der Waals (1890)"
    __status__="vdW"
    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        for componente in mezcla.componente:
            a, b=self.__lib(componente)
            ai.append(a)
            bi.append(b)
        self.kij=Kij(None)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=0
        self.epsilon=0
        self.eta=b

        #TODO: Find relation between u,w y tita,delta, epsilon...
        self.u=0
        self.w=0

        self.dTitadT=0
        super(van_Waals, self).__init__(T, P, mezcla)


    def __lib(self, compuesto):
        a=0.421875*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.125*R_atml*compuesto.Tc/compuesto.Pc.atm
        return  a, b


class RK(Cubic):
    """Ecuación de estado de Redlich-Kwong
    Redlich, O.; Kwong, J.N.S., On The Thermodynamics of Solutions. Chem. Rev. 1949, 44, 233."""
    __title__="Redlich-Kwong (1949)"
    __status__="RK"
    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        for componente in mezcla.componente:
            a, b=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
        self.kij=Kij(None)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        #TODO: Find relation between u,w y tita,delta, epsilon...
        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(RK, self).__init__(T, P, mezcla)

    def __lib(self, compuesto, T):
        a=0.42747*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        alfa=compuesto.tr(T)**-0.5
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        return  a*alfa, b


class Wilson(Cubic):
    """Ecuación de estado de Wilson 1964
    Wilson, G. M.: Adv. Cryogenic Eng., 9: 168 (1964)."""
    __title__="Wilson (1964)"
    __status__="Wilson"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        for componente in mezcla.componente:
            a, b=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
        self.kij=Kij(None)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(Wilson, self).__init__(T, P, mezcla)

    def __lib(self, compuesto, T):
        """Librería de cálculo de la ecuación de estado de Wilson"""
        a=0.42747*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        alfa=1+(1.57+1.62*compuesto.f_acent/(compuesto.tr(T)-1))*compuesto.tr(T)
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        return  a*alfa, b


class Fuller(Cubic):
    """Ecuación de estado de Fuller 1976
    Fuller, G. G.: Ind. Eng. Chem. Fundam., 15: 254 (1976)."""
    __title__="Fuller (1976)"
    __status__="Fuller"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        ci=[]
        for componente in mezcla.componente:
            a, b, c=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            ci.append(ac)
        self.kij=Kij(SRK)
        a, b=mezcla.Mix_van_der_Waals([ai, bi, ci], self.kij)
        tdadt=0

        self.ai=ai
        self.bi=bi
        self.ci=ci
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=c
        self.w=0

        self.dTitadT=tdadt
        super(Fuller, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        #FIXME: La expresión para calcular beta no es correcta, necesita el valor de b que se culcula a partir de el, algún fallo habra en la bibliografia
        b=compuesto.vc*compuesto.peso_molecular
        beta=b/compuesto.vc*compuesto.peso_molecular
        c=1/beta*(sqrt(1/beta-0.75)-1.5)
        Wb=beta*((1-beta)*(2+c*beta)-(1+c*beta))/((2+c*beta)*(1-beta)**2)
        b=Wb*R_atml*compuesto.Tc/compuesto.Pc.atm
        Wa=(1+c*beta)**2*Wb/beta/(1-beta)**2/(2+c*beta)
        m=0.48+1.574*compuesto.f_acent-0.176*compuesto.f_acent**2
        q=(beta/0.26)**0.25*m
        alfa=(1+q*(1-compuesto.tr(T)**0.5))**2
        a=Wa*R_atml**2*compuesto.Tc*alfa/compuesto.Pc.atm
        return  a, b, c


class SRK(Cubic):
    """Ecuación de estado de Soave-Redlich-Kwong
    Soave, G. Equilibrium constants from a modified Redlich-Kwong equation of state. Chem. Eng. Sci. 1972, 27, 1197."""
    __title__="SRK (1972)"
    __status__="SRK"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(SRK)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(SRK, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        """Librería de cálculo de la ecuación de estado de Soave-Redlich-Kwong,"""
        Tr=T/compuesto.Tc
        ac=0.42748*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        m=0.48+1.574*compuesto.f_acent-0.176*compuesto.f_acent**2

        Config=config.getMainWindowConfig()
        Alpha_Mathias=Config.getint("Thermo","Alfa")
        if Alpha_Mathias==1 and Tr>1:
            d=1.+m/2.
            c=1.-1./d
            alfa=exp(c*(1-Tr**d))**2
        else:
            alfa=(1+m*(1-Tr**0.5))**2
        return ac*alfa, b, ac, m


class SRK_API(Cubic):
    """Ecuación de estado de Soave-Redlich-Kwong modificada publicada en el API Technical Databook
    Soave, G.: Inst. Chem. Eng. Symp. Ser., 56(1.2): 1 (1979)."""
    __title__="SRK-API (1979)"
    __status__="SRK-API"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(SRK)
        a, b=mezcla.Mix_van_der_Waals([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(SRK_API, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        """Librería de cálculo de la ecuación de estado de Soave-Redlich-Kwong,"""
        Tr=T/compuesto.Tc
        ac=0.42748*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        m=0.48505+1.55171*compuesto.f_acent-0.15613*compuesto.f_acent**2

        Config=config.getMainWindowConfig()
        Alpha_Mathias=Config.getint("Thermo","Alfa")
        if Alpha_Mathias==1 and Tr>1:
            d=1.+m/2.
            c=1.-1./d
            alfa=exp(c*(1-Tr**d))**2
        else:
            alfa=(1+m*(1-Tr**0.5))**2
        return ac*alfa, b, ac, m


class MSRK(Cubic):
    """Ecuación de estado de Soave-Redlich-Kwong modificada de dos parámetros
    Soave, G.: Chem. Eng. Sci., 39: 357 (1984)."""
    __title__="M-SRK (1984)"
    __status__="MSRK"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(SRK)
        a, b=mezcla.Mix_van_der_Waals([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(MSRK, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        ac=0.42748*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        if compuesto.MSRK[0]==0 and compuesto.MSRK[1]==0:
            m=0.48+1.574*compuesto.f_acent-0.176*compuesto.f_acent**2
            alf=alfa=(1+m*(1-compuesto.tr(T)**0.5))**2
        else:
            m=0.48+1.574*compuesto.f_acent-0.176*compuesto.f_acent**2
            alf=1.+(1-compuesto.tr(T))*(compuesto.MSRK[0]+compuesto.MSRK[1]/compuesto.tr(T))
        return ac*alf, b, ac, m


class SRK_Graboski(Cubic):
    """Ecuación de estado de Soave-Redlich-Kwong Graboski Daubert API 8D4.1, par 820
       Graboski, M. S., Daubert, T. E., “A Modified Soave Equation of State for Phase Equilibrium Calculations-II. Systems Containing CO,, H,S, N2, and C0,”Ind. Eng. Chem. ProcessDes. Develop. 17 (1978)."""
    __title__="SRK-Graboski-Daubert (1978)"
    __status__="SRK-GD"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        for componente in mezcla.componente:
            a, b=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
        self.kij=Kij(SRK)
        a, b=mezcla.Mix_van_der_Waals([ai, bi], self.kij)
        tdadt=0

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(SRK_Thorwart, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        a=0.42748*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        if not compuesto.SRKGraboski[1]:
            m=0.48505+1.55171*compuesto.f_acent-0.15613*compuesto.f_acent**2

            Config=config.getMainWindowConfig()
            Alpha_Mathias=Config.getint("Thermo","Alfa")
            if Alpha_Mathias==1 and Tr>1:
                d=1.+m/2.
                c=1.-1./d
                alfa=exp(c*(1-Tr**d))**2
            else:
                alfa=(1+m*(1-Tr**0.5))**2
        elif not compuesto.SRKGraboski[0]:
            S1=0.48508+1.55171*compuesto.f_acent-0.15613*compuesto.f_acent**2
            S2=compuesto.SRKGraboski[1]
            alfa=(1+S1*(1-Tr**0.5)+S2*(1-Tr**0.5)/Tr**0.5)**2
        return  a*alfa, b


class SRK_Mathias(Cubic):
    """Ecuación de estado de Soave-Redlich-Kwong modificada por Mathias
    Mathias, P.M.: A versatile phase equilibrium equation of state. Industrial and Engineering Chemistry PRocess Design and Development 22, 385-391 (1983)"""
    __title__="SRK-Mathias (1983)"
    __status__="SRK-Math"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(SRK)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(SRK, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        """Librería de cálculo de la ecuación de estado de Soave-Redlich-Kwong,"""
        Tr=T/compuesto.Tc
        ac=0.42748*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        m=0.48508+1.55191*compuesto.f_acent-0.15613*compuesto.f_acent**2

        Config=config.getMainWindowConfig()
        Alpha_Mathias=Config.getint("Thermo","Alfa")
        if Alpha_Mathias==1 and Tr>1:
            d=1.+m/2.+0.3*compuesto.Mathias
            c=1.-1./d
            alfa=exp(c*(1-Tr**d))**2
        else:
            alfa=(1+m*(1-Tr**0.5)*compuesto.Mathias*(1-Tr)*(0.7-Tr))**2
        return ac*alfa, b, ac, m


class SRK_Adachi(Cubic):
    """Ecuación de estado de Soave-Redlich-Kwong modificada por Adachi-Lu
   Adachi, Y., Lu, B.C.Y.: Simplest equation of state for vapor-liquid equilibrium calculation: a modification of the van der Walls equation. Journal of the American Institute of Chemical Engineers 30, 991-993 (1984)"""
    __title__="SRK-Adachi-Lu (1984)"
    __status__="SRK-Adachi"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(SRK)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(SRK, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        """Librería de cálculo de la ecuación de estado de Soave-Redlich-Kwong,"""
        Tr=T/compuesto.Tc
        ac=0.42748*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        m=0.48508+1.55191*compuesto.f_acent-0.15613*compuesto.f_acent**2
        alfa=compuesto.Adachi[0]*10**(compuesto.Adachi[1]*(1-Tr))
        return ac*alfa, b, ac, m


class SRK_Androulakis(Cubic):
    """Ecuación de estado de Soave-Redlich-Kwong modificada por Androulakis
    Andoulakis I.P., Kalospiros, N.S., Tassios, D.P.: Thermophysical properties of pure polar and nonpolar compounds with a modified vdW-711 equation of state. Fluid Phase Equilibria 45, 135-163 (1989)"""
    __title__="SRK-Androulakis (1984)"
    __status__="SRK-And"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(SRK)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=b
        self.epsilon=0
        self.eta=b

        self.u=1
        self.w=0

        self.dTitadT=tdadt
        super(SRK, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        """Librería de cálculo de la ecuación de estado de Soave-Redlich-Kwong,"""
        Tr=T/compuesto.Tc
        ac=0.42748*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
        m=0.48508+1.55191*compuesto.f_acent-0.15613*compuesto.f_acent**2

        Config=config.getMainWindowConfig()
        Alpha_Mathias=Config.getint("Thermo","Alfa")
        if Alpha_Mathias==1 and Tr>1:
            alfa=exp(compuesto.Androulakis[0]*(1-Tr**(2./3)))
        else:
            alfa=(1+compuesto.Androulakis[0]*(1-Tr**(2./3))+compuesto.Androulakis[1]*(1-Tr**(2./3))**2+compuesto.Androulakis[2]*(1-Tr**(2./3))**3)**2
        return ac*alfa, b, ac, m


class PR(Cubic):
    """Ecuación de estado de Peng Robinson
    Peng, D.-Y.; Robinson, D.B. A New Two-Constant Equation of State. I&EC Fundam. 1976, 15(1), 59."""
    __title__="Peng-Robinson (1976)"
    __status__="PR"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        a=0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm
        m=0.37464+1.54226*compuesto.f_acent-0.26992*compuesto.f_acent**2

        Config=config.getMainWindowConfig()
        Alpha_Mathias=Config.getint("Thermo","Alfa")
        if Alpha_Mathias==1 and Tr>1:
            d=1.+m/2.
            c=1.-1./d
            alfa=exp(c*(1-Tr**d))**2
        else:
            alfa=(1+m*(1-Tr**0.5))**2
        return a*alfa, b, a, m



class PRSV(Cubic):
    """Ecuación de estado de Peng Robinson modificada por Stryjek y Vera, v1"""
    __title__="PR-SV (1986)"
    __status__="PR-SV"
    __doi__ = {"autor": "Stryjek, R.; Vera, J.H.",
               "title": "PRSV: An improved peng—Robinson equation of state for pure compounds and mixtures",
               "ref": "Can. J. Chem. Eng. 1986, 64: 323–333",
               "doi":  "10.1002/cjce.5450640224"},

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, k=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(k)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR_SV, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        w = compuesto.f_acent
        if w>=0.49:
            ko=0.378893+1.4897153*w-0.17131848*w**2+0.0196554*w**3
        else:
            ko=0.37464+1.54226*w-0.26992*w**2

        if compuesto.PRSV_k1:
            # TODO: Add data from journal to database
            k1 = compuesto.PRSV_k1
        elif Tr>=0.7:
            k1=0
        elif 1<compuesto.C<=18:
            k1=[-0.00159, 0.02669, 0.03136, 0.03443, 0.03946, 0.05104, 0.04648,
                0.04464, 0.04104, 0.04510, 0.02919, 0.05426, 0.04157, 0.02686,
                0.01892, 0.02665, 0.04048, 0.08291][compuesto.C]
        else:
            k1=0
        k=ko+k1*(1.+sqrt(Tr))*(0.7-Tr)
        alfa=(1+k*(1-Tr**0.5))**2
        a=0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm

        return a*alfa, b, a, k


class PRSV2(Cubic):
    """Ecuación de estado de Peng Robinson modificada por Stryjek y Vera, v2"""
    __title__="PR-SV2 (1986)"
    __status__="PR-SV2"
    __doi__ = {"autor": "Stryjek, R.; Vera, J.H.",
               "title": "PRSV2: A cubic equation of state for accurate vapor—liquid equilibria calculations",
               "ref": "Can. J. Chem. Eng., 64: 820–826",
               "doi":  "10.1002/cjce.5450640516"},

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, k=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(k)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR_SV, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        w = compuesto.f_acent
        if w>=0.49:
            ko=0.378893+1.4897153*w-0.17131848*w**2+0.0196554*w**3
        else:
            ko=0.37464+1.54226*w-0.26992*w**2

        if compuesto.PRSV_k1 and compuesto.PRSV_k2:
            # TODO: Add data from journal to database
            k1 = compuesto.PRSV_k1
            k2 = compuesto.PRSV_k2
            k3 = compuesto.PRSV_k3
        else:
            # Use PRSV V1
            k2 = 0
            K3 = 0
            if Tr>=0.7:
                k1=0
            elif 1<compuesto.C<=18:
                k1=[-0.00159, 0.02669, 0.03136, 0.03443, 0.03946, 0.05104,
                    0.04648, 0.04464, 0.04104, 0.04510, 0.02919, 0.05426,
                    0.04157, 0.02686, 0.01892, 0.02665, 0.04048, 0.08291][compuesto.C]
            else:
                k1=0
        k = ko+(k1+k2*(k3-Tr)*(1-Tr**0.5))*(1+Tr**0.5)*(0.7-Tr)
        alfa = (1+k*(1-Tr**0.5))**2
        a = 0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b = 0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm

        return a*alfa, b, a, k


class PR_Gasem(Cubic):
    """Ecuación de estado de Peng Robinson modificada por Gasem (2001)
    Gasem, Gao, Pan & Robinson: Fluid Phase Equilibria, 181, 113-125 (2001)"""
    __title__="PR Gassem (2001)"
    __status__="PR-Gas"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR_Gasem, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        m=0.134+0.508*compuesto.f_acent-0.0467*compuesto.f_acent**2
        alfa=exp((2.+0.836*Tr)*(1-Tr**m))
        a=0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm
        return a*alfa, b, a, m


class PR_Melhem(Cubic):
    """Ecuación de estado de Peng Robinson modificada por Melhem
    Melhem, G.A.; Saini, R.; Goodwin, B.M. A Modified Peng-Robinson Equation of State. Fluid Phase Eq. 1989, 47, 189."""
    __title__="PR Melhem (1989)"
    __status__="PR-Mel"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        for componente in mezcla.componente:
            a, b, ac=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR_Melhem, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        alfa=exp(compuesto.Melhem[0]*(1-Tr)+compuesto.Melhem[1]*(1-Tr**0.5)**2)
        a=0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm
        return a*alfa, b, a


class PR_Almeida(Cubic):
    """Ecuación de estado de Peng Robinson modificada por Almeida
    Almeida, G.S.; Aznar, M. and Silva Telles, A., Uma Nova Forma de Dependência com a Temperatura do Termo Atrativo de Equaçöes de Estado Cúbicas, RBE, Cad. Eng. Quim., 8, 95-123, (1991)"""
    __title__="PR Almeida (1991)"
    __status__="PR-Alm"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        for componente in mezcla.componente:
            a, b, ac=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR_Almeida, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        alfa=exp(compuesto.Almeida[0]*(1-Tr)*abs(1-Tr)**(compuesto.Almeida[2]-1)+compuesto.Almeida[1]*(Tr**-1-1))
        a=0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm
        return a*alfa, b, a


class PR_Mathias_Copeman(Cubic):
    """Ecuación de estado de Peng-Robinson modificada por Mathias-Copeman
    Mathias, P.M., Copeman, T.W.: Extension of the Peng-Robinson equation of the various forms of the local composition concept. Fluid Phase Equilibria 13, 91-108."""
    __title__="PR-Mathias-Copeman (1983)"
    __status__="PR-MC"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR_Mathias_Copeman, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        a=0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm

        Config=config.getMainWindowConfig()
        Alpha_Mathias=Config.getint("Thermo","Alfa")
        if Alpha_Mathias==1 and Tr>1:
            alfa=(1+compuesto.MathiasCopeman[0]*(1-Tr**0.5))**2
        else:
            alfa=(1+compuesto.MathiasCopeman[0]*(1-Tr**0.5)+compuesto.MathiasCopeman[1]*(1-Tr**0.5)**2+compuesto.MathiasCopeman[2]*(1-Tr**0.5)**3)**2
        return a*alfa, b, a, m


class PR_Yu_Lu(Cubic):
    """Ecuación de estado de Peng-Robinson modificada por Yu Lu
    Yu, J.-M.; Lu, B.C.-Y. A three-parameter cubic equation of state for asymmetric mixture density calculations. Fluid Phase Eq. 1987, 34, 1."""
    __title__="PR-Yu Lu (1987)"
    __status__="PR-YL"

    def __init__(self, T, P, mezcla):
        ai=[]
        bi=[]
        aci=[]
        mi=[]
        for componente in mezcla.componente:
            a, b, ac, m=self.__lib(componente, T)
            ai.append(a)
            bi.append(b)
            aci.append(ac)
            mi.append(m)
        self.kij=Kij(PR)
        a, b=mezcla.Mixing_Rule([ai, bi], self.kij)
        tdadt=0
        for i in range(len(mezcla.componente)):
            for j in range(len(mezcla.componente)):
                tdadt-=mezcla.fraccion[i]*mezcla.fraccion[j]*mi[j]*(aci[i]*aci[j]*mezcla.componente[j].tr(T))**0.5*(1-self.kij[i][j])

        self.ai=ai
        self.bi=bi
        self.b=b
        self.tita=a
        self.delta=2*b
        self.epsilon=-b**2
        self.eta=b

        self.u=2
        self.w=-1

        self.dTitadT=tdadt
        super(PR_Yu_Lu, self).__init__(T, P, mezcla)


    def __lib(self, compuesto, T):
        Tr=T/compuesto.Tc
        ac=(0.46863-0.0378304*compuesto.f_acent-0.00751969*compuesto.f_acent**2)*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        if compuesto.Yu_Lu==[0, 0, 0]:
            if compuesto.f_acent<=0.49:
                m=0.406846+1.87907*compuesto.f_acent-0.792636*compuesto.f_acent**2+0.737519*compuesto.f_acent**3
                A0=0.535843
                A1=-0.39244
                A2=0.26507
            else:
                m=0.581981+0.17141*compuesto.f_acent-1.84441*compuesto.f_acent**2+1.19047*compuesto.f_acent**3
                A0=0.79355
                A1=-0.53409
                A2=0.37273
        else:
            A0, A1, A2=compuesto.Yu_Lu
        if Tr<1:
            alfa=10**(m*(A0+A1*Tr+A2*Tr**2)*(1-Tr))
        else:
            alfa=10**(m*(A0+A1+A2)*(1-Tr))
        b=(0.0892828-0.0640903*compuesto.f_acent-0.00518289*compuesto.f_acent**2)*R_atml*compuesto.Tc/compuesto.Pc.atm
        c=b*(-1.29917+0.648463*compuesto.f_acent+0.895926*compuesto.f_acent**2)
        return ac*alfa, b, c



_all=[van_Waals, RK, Wilson, Fuller, SRK, SRK_API, MSRK, SRK_Graboski, PR, PRSV, PR_Gasem, PR_Melhem, PR_Almeida]

if __name__ == "__main__":
    from lib.corriente import Mezcla
    mezcla = Mezcla(1, ids=[98], caudalUnitarioMasico=[1.])
    for T in [125, 135, 145, 165, 185, 205]:
        eq = Virial(T, 1, mezcla)
        print(eq.H_exc)
