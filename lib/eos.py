#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library to add EoS common functionality
###############################################################################

from scipy import exp, log, log10, tan, sinh, tanh, arctan, sqrt
from scipy import roots, r_
from scipy.constants import pi, Avogadro, R
from scipy.optimize import fsolve

import unidades
import config
from physics import R_atml, factor_acentrico_octano
from bip import SRK, PR, BWRS



def Alfa(Tr, m, f_acent, EOS="SRK", alfa=None):
    """Diferentes expresiónes de cálculo de la alfa de la ecuaciones de estado de SRK y PR
    alfa: Indica el método de cálculo de alfa
        0 - Original
        1 - Boston Mathias
        2 - Twu
        3 - Doridon
        Boston, J.F.; Mathias, P.M. Phase Equilibria in a Third-Generation Process Simulator. Proc. 2nd. Int. Conf. On Phase Equilibria and Fluid Properties in the Chemical Process Industries, Berlin, Germany 17.-21.3.1980, p. 823.
        """
    if not alfa:
        Config=config.getMainWindowConfig()
        alfa=Config.getint("Thermo","Alfa")

    if alfa==2: #Función alfa de Twu et Alt.
        if EOS=="PR":
            if Tr>1:
                L0=0.401219
                M0=4.963075
                N0=-0.2
                L1=0.024655
                M1=1.248088
                N1=-8.
            else:
                L0=0.125283
                M0=0.911807
                N0=1.948153
                L1=0.511614
                M1=0.784054
                N1=2.812522
        else:
            if Tr>1:
                L0=0.441411
                M0=6.500018
                N0=-0.2
                L1=0.032580
                M1=1.289098
                N1=-8.
            else:
                L0=0.141599
                M0=0.919422
                N0=2.496441
                L1=0.500315
                M1=0.799457
                N1=3.29179
        alfa0=Tr**(N0*(M0-1))*exp(L0*(1-Tr**(N0*M0)))
        alfa1=Tr**(N1*(M1-1))*exp(L1*(1-Tr**(N1*M1)))
        alfa=alfa0+f_acent*(alfa1-alfa0)
    elif alfa==3:
        if f_acent<0.4:
            m=0.418+1.58*f_acent-0.58*f_acent**2
        else:
            m=0.212+2.2*f_acent-0.831*f_acent**2
        alfa=exp(m*(1-Tr))
    elif alfa==1 and Tr>1:
        d=1.+m/2.
        c=1.-1./d
        alfa=exp(c*(1-Tr**d))**2
    else:
        alfa=(1+m*(1-Tr**0.5))**2
    return alfa

def volume_Correction_Peneloux(compuesto):
    """Peneloux, A.; Rauzy, E.; Freze, R. A consistent correction for Redlich-Kwong-Soave volumes. Fluid Phase Eq. 1982, 8, 7."""
    return 0.40768*(0.29441-compuesto.rackett)*R_atml*compuesto.Tc/compuesto.Pc.atm


class EoS(object):

    def _Flash(self):
        """Cálculo de los coeficientes de reparto entre fases, Ref Naji - Conventional and rapid flash claculations"""
        #Estimación inicial de K mediante correlación wilson Eq 19
        Ki=[]
        for i in self.componente:
            Ki.append(i.Pc/self.P*exp(5.37*(1.+i.f_acent)*(1.-i.Tc/self.T)))

        Rachford=lambda x: sum([zi*(ki-1.)/(1.-x+x*ki) for zi, ki in zip(self.fraccion, Ki)])

        if Rachford(0)>0 and Rachford(1)>0: #x>1, por tanto solo fase vapor
            xi=self.fraccion
            yi=self.fraccion
            x=1.
        elif Rachford(0)<0 and Rachford(1)<0: #x<0, por tanto solo fase líquida
            xi=self.fraccion
            yi=self.fraccion
            x=0.
        else:
            x=0.5
            while True:
                xo=x
                solucion=fsolve(Rachford, x, full_output=True)
                if solucion[2]!=1:
                    print solucion
                    break
                else:
                    x=solucion[0][0]

                    xi=[]
                    yi=[]
                    for zi, ki in zip(self.fraccion, Ki):
                        xi.append(float(zi/(1-x+x*ki)))
                        yi.append(float(zi*ki/(1-x+x*ki)))

                    tital=self._fug(self.Z[1], xi)
                    titav=self._fug(self.Z[0], yi)
                    fiv=[z*t*self.P for z, t in zip(yi, titav)]
                    fil=[z*t*self.P for z, t in zip(xi, tital)]
                    #criterio de convergencia Eq 21
                    if sum([abs(l/v-1) for l, v in zip(fil, fiv)])< 1e-14 and abs(x-xo) < 1e-10:
                        break
                    else:
                        Ki=[l/v for l, v in zip(tital, titav)]

        return x, xi, yi, Ki


    def _Bubble_T(self):
        def f(T):
            eq=self.__class__(T, self.P.atm, self.mezcla)
            return sum([k*x for k, x in zip(eq.Ki, self.fraccion)])-1.

        T=fsolve(f, self.T)
        return unidades.Temperature(T)

    def _Bubble_P(self):
        def f(P):
            eq=self.__class__(self.T, P, self.mezcla)
            return sum([k*x for k, x in zip(eq.Ki, self.fraccion)])-1.

        P=fsolve(f, self.P.atm)
        return unidades.Pressure(P, "atm")

    def _Dew_T(self):
        def f(T):
            eq=self.__class__(T, self.P.atm, self.mezcla)
            return 1./sum([x/k for k, x in zip(eq.Ki, self.fraccion)])-1.

        T=fsolve(f, self.T)
        return unidades.Temperature(T)

    def _Dew_P(self):
        def f(P):
            eq=self.__class__(self.T, P, self.mezcla)
            return sum([x/k for k, x in zip(eq.Ki, self.fraccion)])-1.

        P=fsolve(f, self.P.atm)
        return unidades.Pressure(P, "atm")


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
        self.kij=mezcla.Kij(None)
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
        self.kij=mezcla.Kij(None)
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
        self.kij=mezcla.Kij(None)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(SRK)
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
        self.kij=mezcla.Kij(PR)
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



class PR_SV(Cubic):
    """Ecuación de estado de Peng Robinson modificada por Stryjek y Vera
    Stryjek, R.; Vera, J.H. PRSV2: A Cubic Equation of State for Accurate Vapor-Liquid Equilibria Calculations. Can. J. Chem. Eng. 1986b, 64, 820."""
    __title__="PR-SV (1986)"
    __status__="PR-SV"

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
        self.kij=mezcla.Kij(PR)
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
        """Stryjek, R.; Vera, J.H. PRSV2: A Cubic Equation of State for Accurate Vapor-Liquid Equilibria Calculations. Can. J. Chem. Eng. 1986b, 64, 820."""
        Tr=T/compuesto.Tc
        if compuesto.f_acent>=0.49:
            mo=0.378893+1.4897153*compuesto.f_acent-0.17131848*compuesto.f_acent**2+0.0196554*compuesto.f_acent**3
        else:
            mo=0.37464+1.54226*compuesto.f_acent-0.26992*compuesto.f_acent**2
        if Tr>=0.7:
            m1=0
        else:
            if compuesto.indice==46:
                m1=0.01996
            elif compuesto.indice==47:
                m1=0.01512
            elif compuesto.indice==49:
                m1=0.04285
            elif compuesto.indice==63:
                m1=0.001
            elif compuesto.indice==62:
                m1=-0.06635
            elif compuesto.indice==104:
                m1=0.01989
            else:
                if 1<compuesto.C<=18:
                    m1=[-0.00159, 0.02669, 0.03136, 0.03443, 0.03946, 0.05104, 0.04648, 0.04464, 0.04104, 0.04510, 0.02919, 0.05426, 0.04157, 0.02686, 0.01892, 0.02665, 0.04048, 0.08291][compuesto.C]
                else:
                    m1=0
        m=mo+m1*(1.+sqrt(Tr))*(0.7-Tr)
        alfa=(1+m*(1-Tr**0.5))**2
        a=0.457235*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
        b=0.077796*R_atml*compuesto.Tc/compuesto.Pc.atm

        return a*alfa, b, a, m


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
        self.kij=mezcla.Kij(PR)
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
        self.kij=mezcla.Kij(PR)
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
        self.kij=mezcla.Kij(PR)
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
        self.kij=mezcla.Kij(PR)
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
        self.kij=mezcla.Kij(PR)
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



_cubic=[van_Waals, RK, Wilson, Fuller, SRK, SRK_API, MSRK, SRK_Graboski, PR, PR_SV, PR_Gasem, PR_Melhem, PR_Almeida]






def PT_lib(compuesto, T):
    """Librería de cálculo de la ecuación de estado de Patel-Teja"""
    if compuesto.Tc<>0 and compuesto.Pc<>0 and compuesto.vc<>0:
        Zc=compuesto.Pc.atm*compuesto.vc*compuesto.peso_molecular/R_atml/compuesto.Tc
    else:
        Zc=0.329032+0.076799*compuesto.f_acent-0.0211947*compuesto.f_acent**2

    c=(1-3*Zc)*R_atml*compuesto.Tc/compuesto.Pc.atm
    omega=roots([1, 2-3*Zc, 3*Zc**2, -Zc**3])
    Bpositivos=[]
    for i in omega:
        if i>0:
            Bpositivos.append(i)
    omegab=min(Bpositivos)
    b=omegab*R_atml*compuesto.Tc/compuesto.Pc.atm

    f=0.452413+1.30982*compuesto.f_acent-0.295937*compuesto.f_acent**2
    alfa=(1+f*(1-sqrt(T/compuesto.Tc)))**2
    omegaa=3*Zc**2+3*(1-2*Zc)*omegab+omegab**2+1-3*Zc
    a=omegaa*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm**2*alfa
    return a, b, c

def PTC_lib(compuesto, T):
    """Librería de cálculo de la ecuación de estado de Patel-Teja-Crause"""
    if compuesto.Tc<>0 and compuesto.Pc<>0 and compuesto.vc<>0:
        Zc=compuesto.Pc.atm*compuesto.vc*compuesto.peso_molecular/R_atml/compuesto.Tc
    else:
        Zc=0.253168556+0.09253329*exp(-0.0048018*compuesto.peso_molecular)

    c=(1-3*Zc)*R_atml*compuesto.Tc/compuesto.Pc.atm
    omega=roots([1, 2-3*Zc, 3*Zc**2, -Zc**3])
    Bpositivos=[]
    for i in omega:
        if i>0:
            Bpositivos.append(i)
    omegab=min(Bpositivos)
    b=omegab*R_atml*compuesto.Tc/compuesto.Pc.atm

    f=0.002519*compuesto.peso_molecular+0.70647
    alfa=(1+f*(1-sqrt(T/compuesto.Tc)))**2
    omegaa=3*Zc**2+3*(1-2*Zc)*omegab+omegab**2+1-3*Zc
    a=omegaa*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm**2*alfa
    return a, b, c

def PTV_lib(compuesto, T):
    """Librería de cálculo de la ecuación de estado de Patel-Teja-Valderrama (1990)"""
    if compuesto.Tc<>0 and compuesto.Pc<>0 and compuesto.vc<>0:
        Zc=compuesto.Pc.atm*compuesto.vc*compuesto.peso_molecular/R_atml/compuesto.Tc
    else:
        Zc=0.329032-0.076799*compuesto.f_acent+0.0211947*compuesto.f_acent**2
    Tr=T/compuesto.Tc
    m=0.4628+3.5823*compuesto.f_acent*Zc+8.1942*compuesto.f_acent**2*Zc**2
    alfa=(1+m*(1-Tr**0.5))**2
    omegaa=0.6612-0.7616*Zc
    a=omegaa*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
    omegab=0.0221+0.2087*Zc
    b=omegab*R_atml*compuesto.Tc/compuesto.Pc.atm
    omegac=0.5777-1.8718*Zc
    c=omegac*R_atml*compuesto.Tc/compuesto.Pc.atm
    return a*alfa, b, c

def PTVC_lib(compuesto, T):
    """Librería de cálculo de la ecuación de estado de Patel-Teja-Valderrama-Cisternas (1986)"""
    if compuesto.Tc<>0 and compuesto.Pc<>0 and compuesto.vc<>0:
        Zc=compuesto.Pc.atm*compuesto.vc*compuesto.peso_molecular/R_atml/compuesto.Tc
    else:
        Zc=0.329032-0.076799*compuesto.f_acent+0.0211947*compuesto.f_acent**2
    Tr=T/compuesto.Tc
    m=-6.608+70.43*Zc-159.0*Zc**2
    alfa=(1+m*(1-Tr**0.5))**2
    omegaa=0.694-1.063*Zc+0.683*Zc**2-0.21*Zc**3+0.00375*Zc**4
    a=omegaa*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
    omegab=0.026-0.181*Zc+0.061*Zc**2
    b=omegab*R_atml*compuesto.Tc/compuesto.Pc.atm
    omegac=0.578-1.904*Zc
    c=omegac*R_atml*compuesto.Tc/compuesto.Pc.atm
    return a*alfa, b, c

def HPW_lib(compuesto, T, P, alfa):
    """Hederer, H.; Peter, S., Wenzel, H. Calculation of Thermodynamic Properties from a Modified Redlich-Kwong Equation of State.  Chem. Eng. J. 1976, 11, 183."""
    #TODO: hallar la forma de calcular el valor de alfa a partir de los datos de la base de datos, se supone que representa la pendiente de la presión de vapor
    ac=0.42747*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
    a=ac*compuesto.tr(T)**alfa
    b=0.08664*R_atml*compuesto.Tc/compuesto.Pc.atm
    return  a, b


def TB_lib(compuesto, T):
    """Librería de cálculo de la ecuación de estado de Trebble-Bishnoi"""
    if compuesto.Tc<>0.0 and compuesto.Pc.atm<>0.0 and compuesto.vc<>0.0:
        Zc=compuesto.Pc.atm*compuesto.vc*compuesto.peso_molecular/R_atml/compuesto.Tc
    else:
        if compuesto.f_acent<0.225:
            TcPc=775.9+12009*compuesto.f_acent-57335*compuesto.f_acent**2+91393*compuesto.f_acent**3
        elif compuesto.f_acent<=1.:
            TcPc=1876-1160*compuesto.f_acent
        else:
            TcPc=compuesto.Tc*compuesto.Pc.MPa

        if compuesto.f_acent>=-0.14:
            Zc=0.29-0.0885*compuesto.f_acent-0.0005/((compuesto.Tc*compuesto.Pc.MPa)**0.5-TcPc**0.5)
        else:
            Zc=0.3024

    d=0.341*compuesto.vc*compuesto.peso_molecular-0.005
    Xc=1.075*Zc
    Dc=d*compuesto.Pc.atm/R_atml/compuesto.Tc
    Cc=1-3*Xc
    B=roots([1, 2-3*Xc, 3*Xc**2, -Dc**2-Xc**3])
    Bpositivos=[]
    for i in range(3):
        if B[i]>0:
            Bpositivos.append(B[i])
    Bc=min(Bpositivos)
    Ac=3*Xc**2+2*Bc*Cc+Bc+Cc+Bc**2+Dc**2

    if compuesto.indice==212 and compuesto.tr(T)<=1:
        q1=-0.31913
    elif compuesto.f_acent<-0.1:
        q1=0.66208+4.63961*compuesto.f_acent+7.45183*compuesto.f_acent**2
    elif compuesto.f_acent<=0.4:
        q1=0.35+0.7924*compuesto.f_acent+0.1875*compuesto.f_acent**2-28.93*(0.3-Zc)**2
    elif compuesto.f_acent>0.4:
        q1=0.32+0.9424*compuesto.f_acent-28.93*(0.3-Zc)**2

    if compuesto.f_acent<-0.0423:
        q2=0
    elif compuesto.f_acent<=0.3:
        q2=0.05246+1.15058*compuesto.f_acent-1.99348*compuesto.f_acent**2+1.5949*compuesto.f_acent**3-1.39267*compuesto.f_acent**4
    else:
        q2=0.17959+0.23471*compuesto.f_acent

    alfa=exp(q1*(1-compuesto.tr(T)))
    if T<=compuesto.Tc:
        beta=1.+q2*(1-compuesto.tr(T)+log(compuesto.tr(T)))
    else:
        beta=1
    ac=Ac*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
    bc=Bc*R_atml*compuesto.Tc/compuesto.Pc.atm
    c=Cc*R_atml*compuesto.Tc/compuesto.Pc.atm

    return ac*alfa, bc*beta, c, d

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
        Ho=self.Entalpia_ideal(T)
        Delta=self.TB_H_exceso(T, P)
        return unidades.Enthalpy(Delta+Ho)

    def TB_S_exceso(self, T, P):
        """Método de cálculo de la entropía de exceso mediante la ecuación de estado de Trebble-Bishnoi"""
        H=self.TB_H_exceso(T, P)
        f=self.TB_Fugacidad(T, P)
        return unidades.SpecificHeat(H.Jg/T-R/self.peso_molecular*log(f.atm/P), "JgK")

    def TB_Entropia(self, T, P):
        """Método de cálculo de la entropía mediante la ecuación de estado de Trebble-Bishnoi"""
        So=self.Entropia_ideal(T)
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


def DP_lib(compuesto, T):
    """Dohrn, R.; Prausnitz, J.M. A simple perturbation term for the Carnahan-Starling equation of state. Fluid Phase Eq. 1990, 61, 53."""
    Tr=T/compuesto.Tc
    if Tr >=1:
        m=1
    else: m=0

    a1=0.367845+0.055966*compuesto.f_acent
    a2=-1**m*(0.604709+0.008477*compuesto.f_acent)
    ac=0.550408*R_atml**2*compuesto.Tc**2/compuesto.Pc.atm
    a=ac*(a1*tanh(a2*abs(Tr-1)**0.7)+1)

    b1=0.356983-0.190003*compuesto.f_acent
    b2=-1**m*(1.37-1.898981*compuesto.f_acent)
    bc=0.187276*R_atml*compuesto.Tc/compuesto.Pc.atm
    b=bc*(b1*tanh(b2*abs(log(Tr))**0.8)+1)
    sigma=(3*b/2/pi/Avogadro)**(1./3)
    return a, b, sigma

    def DP_Z(self, T, P):
        """Factor de compresibilidad según la ecuación de estado de Dohrn-Prausnith"""
        V=self.DP_V(T, P)*self.peso_molecular
        return P*V/R_atml/T

    def DP_V(self, T, P):
        """Volumen según el modelo de Dohrn-Prausnith"""
        a, b, sigma=self.DP_lib(T, P)
        D=sigma
        E=sigma**2
        F=sigma**3

        if self.Fase(T, P)=="gas":
            V=25
        else: V=0.5

        def Vm(V):
            nu=b/4/V
            Zref=(1+(3*D*E/F-2)*nu+(3*E**3/F**2-3*D*E/F+1)*nu**2-E**3/F**2*nu**3)/(1-nu)**3
            Zpert=a/R_atml/T/V*(1-1.41*b/V/4+5.07*(b/V/4)**2)
            return V-(Zref+Zpert)*R_atml*T/P

        v=fsolve(Vm, V)

        return unidades.SpecificVolume(v/self.peso_molecular)

    def DP_RhoG(self, T, P):
        """Método para el cálculo de la densidad de gases haciendo uso de la ecuación de estado de Dohrn-Prausnith"""
        z=self.DP_Z(T, P)
        return unidades.SpecificVolume(P/z/R_atml/T*self.peso_molecular)


def ESD_lib(compuesto, T):
    """Elliot, J.R.; Suresh, S.J.; Donohue, M.D. A Simple Equation of State for Nonspherical and Associating Molecules. I&EC Res. 1990, 29, 1476."""
    #TODO: Añadir en la base de datos los parametros de esta ecuación, disponibles en la versión 6 de chemcad.
    c=1+3.535*compuesto.f_acent+0.533*compuesto.f_acent**2
    a=compuesto.Tc*(1+0.945*(c-1)+0.134*(c-1)**2)/(1.023+2.225*(c-1)+0.478*(c-1)**2)
    b=R_atml*compuesto.Tc/compuesto.Pc.atm*(0.0312+0.087*(c-1)+0.008*(c-1)**2)/(1+2.455*(c-1)+0.732*(c-1)**2)
    e=a*Bolzmann
    Y=exp(a/T)-1.0617
    return  a, b, c, e, Y

    def ESD_V(self, T, P):
        """Factor de compresibilidad según la ecuación de estado de Elliot-Suresh-Donohue"""
        a, b, c, e, Y=self.ESD_lib(T, P)

        if self.Fase(T, P)=="gas":
            V=25
        else:
            V=0.5

        def Vm(V):
            cnu=c*b/V
            qnuy=(1+1.90476*(c-1))*b*Y/V
            nu=b/V
            nuy=b*Y/V
            return P*V/R_atml/T-1-4*cnu/(1-1.9*nu)+9.5*qnuy/(1+1.7745*nuy)

        V=fsolve(Vm, V)
        return unidades.SpecificVolume(V/self.peso_molecular)

    def ESD_Z(self, T, P):
        """Volumen según el modelo de Elliot-Suresh-Donohue"""
        V=self.ESD_V(T, P)*self.peso_molecular
        return P*V/R_atml/T

    def ESD_RhoG(self, T, P):
        """Método para el cálculo de la densidad de gases haciendo uso de la ecuación de estado de Elliot-Suresh-Donohue"""
        z=self.ESD_Z(T, P)
        return unidades.SpecificVolume(P/z/R_atml/T*self.peso_molecular)


def Z_PT(self):
    """Factor de compresibilidad según la ecuación de estado de Patel_Teja"""
    ai=[]
    bi=[]
    ci=[]
    for componente in self.componente:
        a, b, c=componente.PT_lib(self.T)
        ai.append(a)
        bi.append(b)
        ci.append(c)
    a, b, c=self.Mixing_Rule([ai, bi], self.kij)
    Z=self.Z_Cubic_EoS(b, a, b+c, -c*b, b)
    self.titail, self.titaiv=self.Fugacidad_Cubic_EoS(Z, b, a, ai, bi, (b+c)/b, -c*b)
    return Z

def Z_PTV(self):
    """Factor de compresibilidad según la ecuación de estado de Patel_Teja Valderrama (1990)"""
    ai=[]
    bi=[]
    ci=[]
    for componente in self.componente:
        a, b, c=componente.PTV_lib(self.T)
        ai.append(a)
        bi.append(b)
        ci.append(c)
    a, b, c=self.Mixing_Rule([ai, bi], self.kij)
    Z=self.Z_Cubic_EoS(b, a, b+c, -c*b, b)
    self.titail, self.titaiv=self.Fugacidad_Cubic_EoS(Z, b, a, ai, bi, (b+c)/b, -c*b)
    return Z

def Z_PTVC(self):
    """Factor de compresibilidad según la ecuación de estado de Patel_Teja Valderrama-Cisternas (1986)"""
    ai=[]
    bi=[]
    ci=[]
    for componente in self.componente:
        a, b, c=componente.PTVC_lib(self.T)
        ai.append(a)
        bi.append(b)
        ci.append(c)
    a, b, c=self.Mixing_Rule([ai, bi], self.kij)
    Z=self.Z_Cubic_EoS(b, a, b+c, -c*b, b)
    self.titail, self.titaiv=self.Fugacidad_Cubic_EoS(Z, b, a, ai, bi, (b+c)/b, -c*b)
    return Z

def Z_TB(self):
    """Factor de compresibilidad según la ecuación de estado de Trebble and Bishnoi"""
    ai=[]
    bi=[]
    ci=[]
    di=[]
    for componente in self.componente:
        a, b, c, d=componente.TB_lib(self.T)
        ai.append(a)
        bi.append(b)
        ci.append(c)
        di.append(d)
    a, b, c, d=self.Mixing_Rule([ai, bi, ci, di], self.kij)
    Z=self.Z_Cubic_EoS(b, a, b+c, -(b*c+d**2), b)
    self.titail, self.titaiv=self.Fugacidad_Cubic_EoS(Z, b, a, ai, bi, (b+c)/b, -(b*c+d**2)/b**2)
    return Z



class BWRS(EoS):
    """Ecuación de estado de Benedict-Webb-Rubin-Starling"""
    __title__="Benedict-Webb-Rubin-Starling"
    __status__="BWRS"

    def __init__(self, T, P, mezcla):
        self.T=unidades.Temperature(T)
        self.P=unidades.Pressure(P, "atm")
        self.componente=mezcla.componente
        self.fraccion=mezcla.fraccion
        self.kij=mezcla.Kij(BWRS)

        Aoi=[]
        Boi=[]
        Coi=[]
        Doi=[]
        Eoi=[]
        ai=[]
        bi=[]
        ci=[]
        di=[]
        alfai=[]
        gammai=[]
        for compuesto in self.componente:
            Ao_, Bo_, Co_, Do_, Eo_, a_, b_, c_, d_, alfa_, gamma_=self._lib(compuesto)
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

        Ao=Co=Do=Eo=Bo=a=b=c=d=alfa=gamma=0
        for i in range(len(self.componente)):
            Bo+=self.fraccion[i]*Boi[i]
            a+=self.fraccion[i]*ai[i]**(1./3)
            b+=self.fraccion[i]*bi[i]**(1./3)
            c+=self.fraccion[i]*ci[i]**(1./3)
            d+=self.fraccion[i]*di[i]**(1./3)
            alfa+=self.fraccion[i]*alfai[i]**(1./3)
            gamma+=self.fraccion[i]*gammai[i]**0.5
        a=a**3
        b=b**3
        c=c**3
        d=d**3
        alfa=alfa**3
        gamma=gamma**2

        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
                Ao+=self.fraccion[i]*self.fraccion[j]*Aoi[i]**0.5*Aoi[j]**0.5*(1-self.kij[i][j])
                Co+=self.fraccion[i]*self.fraccion[j]*Coi[i]**0.5*Coi[j]**0.5*(1-self.kij[i][j])**3
                Do+=self.fraccion[i]*self.fraccion[j]*Doi[i]**0.5*Doi[j]**0.5*(1-self.kij[i][j])**4
                Eo+=self.fraccion[i]*self.fraccion[j]*Eoi[i]**0.5*Eoi[j]**0.5*(1-self.kij[i][j])**5

        self.Aoi=Aoi
        self.Boi=Boi
        self.Coi=Coi
        self.Doi=Doi
        self.Eoi=Eoi
        self.ai=ai
        self.bi=bi
        self.ci=ci
        self.di=di
        self.alfai=alfai
        self.gammai=gammai
        self.Ao=Ao
        self.Co=Co
        self.Do=Do
        self.Eo=Eo
        self.Bo=Bo
        self.a=a
        self.b=b
        self.c=c
        self.d=d
        self.alfa=alfa
        self.gamma=gamma


        Vm=lambda V: self.P.atm-R_atml*self.T/V - (Bo*R_atml*self.T-Ao-Co/self.T**2+Do/self.T**3-Eo/self.T**4)/V**2 - (b*R_atml*self.T-a-d/self.T)/V**3 - alfa*(a+d/self.T)/V**6 - c/self.T**2/V**3*(1+gamma/V**2)*exp(-gamma/V**2)

        #Usamos SRK para estimar los volumenes de ambas fases usados como valores iniciales en la iteeración
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
        """Librería de cálculo de la ecuación de estado de Benedict-Webb-Rubin-Starling"""
        a1, b1=0.443690, 0.115449
        a2, b2=1.28438, -0.920731
        a3, b3=0.356306, 1.70871
        a4, b4=0.544979, -0.270896
        a5, b5=0.528629, 0.349261
        a6, b6=0.484011, 0.754130
        a7, b7=0.0705233, -0.044448
        a8, b8=0.504087, 1.32245
        a9, b9=0.0307452, 0.179433
        a10, b10=0.0732828, 0.463492
        a11, b11=0.006450, -0.022143

        Bo=(a1+b1*cmp.f_acent)*cmp.Vc
        Ao=(a2+b2*cmp.f_acent)*R_atml*cmp.Tc*cmp.Vc
        Co=(a3+b3*cmp.f_acent)*R_atml*cmp.Tc**3*cmp.Vc
        gamma=(a4+b4*cmp.f_acent)*(cmp.Vc)**2
        b=(a5+b5*cmp.f_acent)*(cmp.Vc)**2
        a=(a6+b6*cmp.f_acent)*R_atml*cmp.Tc*(cmp.Vc)**2
        alfa=(a7+b7*cmp.f_acent)*(cmp.Vc)**3
        c=(a8+b8*cmp.f_acent)*R_atml*cmp.Tc**3*(cmp.Vc)**2
        Do=(a9+b9*cmp.f_acent)*R_atml*cmp.Tc**4*cmp.Vc
        d=(a10+b10*cmp.f_acent)*R_atml*cmp.Tc**2*(cmp.Vc)**2
        Eo=(a11+b11*cmp.f_acent*exp(-3.8*cmp.f_acent))*R_atml*cmp.Tc**5*cmp.Vc

        return Ao, Bo, Co, Do, Eo, a, b, c, d, alfa, gamma


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



class Grayson_Streed(EoS):
    """Ecuación de estado de Grayson Streed modificada por Chao-Seader
    Chao, K.C. and Seader, J.D.; A General Correlation of Vapor-Liquid Equilibria in Hydrocarbon Mixtures, AIChE Journal, 7, No 4 (December 1961"""
    __title__="Grayson Streed (1961)"
    __status__="GS"

    def __init__(self, T, P, mezcla):
        self.T=unidades.Temperature(T)
        self.P=unidades.Pressure(P, "atm")
        self.mezcla=mezcla
        self.componente=mezcla.componente
        self.fraccion=mezcla.fraccion
        p=unidades.Pressure(P, "atm")

        self.rk=RK(T, P, mezcla)
        self.Z=self.rk.Z
        self.V=self.rk.V
        self.x, self.xi, self.yi, self.Ki=self._Flash()

    def _k(self, xi, yi):
        fi=self.rk._fug(self.rk.Z[0], yi)
        suma1=0
        suma2=0
        Vi=[]
        for i in range(len(self.componente)):
            Vi.append(self.componente[i].Vc*self.componente[i].M)
            suma1+=xi[i]*Vi[i]*self.componente[i].parametro_solubilidad
            suma2+=xi[i]*Vi[i]
        dim=suma1/suma2

        gi=[]
        for i in range(len(self.componente)):
            gi.append(exp(Vi[i]/1e6*(self.componente[i].parametro_solubilidad-dim)**2/R/self.T))
        nio=[]
        for i in self.componente:
            tr=i.tr(self.T)
            pr=i.pr(self.P.atm)
            if i.indice==1:
                A=[1.50709, 2.74283, -0.02110, 0.00011, 0.0, 0.008585, 0., 0., 0., 0.]
            elif i.indice==2:
                A=[1.36822, -1.54831, 0., 0.02889, -0.01076, 0.10486, -0.02529, 0., 0., 0.]
            else:
                A=[2.05135, -2.10899, 0., -0.19396, 0.02282, 0.08852, 0., -0.00872, -0.00353, 0.00203]
            logn0=A[0]+A[1]/tr+A[2]*tr+A[3]*tr**2+A[4]*tr**3+(A[5]+A[6]*tr+A[7]*tr**2)*pr+(A[8]+A[9]*tr)*pr**2-log10(pr)
            logn1=-4.23893+8.65808*tr-1.2206/tr-3.15224*tr**3-0.025*(pr-0.6)
            nio.append(10**(logn0+i.f_acent*logn1))


        tital=[]
        for i in range(len(self.componente)):
            tital.append(nio[i]*gi[i])

        return tital, fi


    def _Flash(self):
        """Cálculo de los coeficientes de reparto entre fases, Ref Naji - Conventional and rapid flash claculations"""
        #Estimación inicial de K mediante correlación wilson Eq 19
        Ki=[]
        for i in self.componente:
            Ki.append(i.Pc/self.P*exp(5.37*(1.+i.f_acent)*(1.-i.Tc/self.T)))

        Rachford=lambda x: sum([zi*(ki-1.)/(1.-x+x*ki) for zi, ki in zip(self.fraccion, Ki)])

        if Rachford(0)>0 and Rachford(1)>0: #x>1, por tanto solo fase vapor
            xi=self.fraccion
            yi=self.fraccion
            x=1.
        elif Rachford(0)<0 and Rachford(1)<0: #x<0, por tanto solo fase líquida
            xi=self.fraccion
            yi=self.fraccion
            x=0.
        else:
            x=0.5
            while True:
                xo=x
                solucion=fsolve(Rachford, x, full_output=True)
                if solucion[2]!=1:
                    print solucion
                    break
                else:
                    x=solucion[0]
                    xi=[]
                    yi=[]
                    for zi, ki in zip(self.fraccion, Ki):
                        xi.append(zi/(1-x+x*ki))
                        yi.append(zi*ki/(1-x+x*ki))

                    titas=self._k(xi, yi)
                    tital=titas[0]
                    titav=titas[1]

                    fiv=[z*t*self.P for z, t in zip(yi, titav)]
                    fil=[z*t*self.P for z, t in zip(xi, tital)]
                    #criterio de convergencia Eq 21
                    if sum([abs(l/v-1) for l, v in zip(fil, fiv)])< 1e-12 and (x-xo)**2 < 1e-15:
                        break
                    else:
                        Ki=[l/v for l, v in zip(tital, titav)]

        return x, xi, yi, Ki




K=_cubic+[BWRS, Lee_Kesler, Grayson_Streed]
H=_cubic+[BWRS, Lee_Kesler]


if __name__ == "__main__":
    from corriente import Mezcla

#    eq=SRK_API(340, 1, Mezcla([10, 38, 22, 61], [0.3, 0.5, 0.05, 0.15]))
#    print eq.x
#    eq=SRK(340, 1., Mezcla([10, 38, 22, 61], [0.3, 0.5, 0.05, 0.15]))
#    print eq.x
#    print eq._Dew_T()

#    p=unidades.Pressure(2000, "psi")
#    t=unidades.Temperature(100, "F")
#    eq=van_Waals(t.K, p.atm, Mezcla([1, 4], [.805, 0.195]))
#    print eq.Z

#    mezcla=SRK(350, 1, Mezcla([4, 5, 6, 7, 8, 10, 11, 12, 13], [0.02361538, 0.2923077, 0.3638462, 0.02769231, 0.01153846, 0.01769231, 0.03007692, 0.2093846, 0.02384615]))
#    mezcla=Grayson_Streed(350, 1, Mezcla([4, 5, 6, 7, 8, 10, 11, 12, 13], [0.02361538, 0.2923077, 0.3638462, 0.02769231, 0.01153846, 0.01769231, 0.03007692, 0.2093846, 0.02384615]))

#    print mezcla._Bubble_T(), mezcla._Dew_T()
#    print mezcla.T, mezcla.x

    mezcla=Mezcla(ids=[10, 38, 22, 61], fraccionMolar=[0.3, 0.5, 0.05, 0.15])
#    for t in range(300, 345, 1):
#        eq=SRK(t, 1., mezcla )
#        print t, eq.x
    eq=SRK(340, 1., mezcla )
    print eq.x
