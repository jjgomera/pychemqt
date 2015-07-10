#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library to add EoS common functionality
###############################################################################

from scipy import exp, log, log10, tan, sinh, tanh, arctan, sqrt
from scipy import roots, r_
from scipy.constants import pi, Avogadro, R
from scipy.optimize import fsolve

from . import unidades
from . import config
from .physics import R_atml, factor_acentrico_octano

#from EoS import *


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
    def __init__(self, T, P, mezcla, **kwargs):
        self.T = unidades.Temperature(T)
        self.P = unidades.Pressure(P, "atm")
        self.mezcla = mezcla
        self.componente = mezcla.componente
        self.fraccion = mezcla.fraccion
        self.kwargs = kwargs

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
                    print(solucion)
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

        if x < 0:
            x = 0
        elif x > 1:
            x = 1
            
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








def PT_lib(compuesto, T):
    """Librería de cálculo de la ecuación de estado de Patel-Teja"""
    if compuesto.Tc!=0 and compuesto.Pc!=0 and compuesto.vc!=0:
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
    if compuesto.Tc!=0 and compuesto.Pc!=0 and compuesto.vc!=0:
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
    if compuesto.Tc!=0 and compuesto.Pc!=0 and compuesto.vc!=0:
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
    if compuesto.Tc!=0 and compuesto.Pc!=0 and compuesto.vc!=0:
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
    if compuesto.Tc!=0.0 and compuesto.Pc.atm!=0.0 and compuesto.vc!=0.0:
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









if __name__ == "__main__":
    from .corriente import Mezcla

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
    print(eq.x)
