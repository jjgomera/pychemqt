#!/usr/bin/python
# -*- coding: utf-8 -*-

import os

from scipy.optimize import fsolve, leastsq, brentq
from scipy.special import erf
from scipy.linalg import det
from scipy import roots, power, log, sqrt, log10, exp, sin, r_, real, imag, zeros, transpose
from pylab import triu, plot, grid, show, figtext, xlabel, ylabel
from PyQt4.QtGui import QApplication

from compuestos import Componente
from bip import srk
from physics import R_atml, R_cal, R, factor_acentrico_octano
from lib import unidades, config, EoS, mEoS, gerg, iapws, freeSteam, refProp, coolProp


class Mezcla(config.Entity):
    """Clase que define una mezcla definiendo, a partir de los componentes que la forman, cada una de sus propiedades físicas a partir de las reglas de mezcla"""
    kwargs={ "caudalMasico": 0.0,
                    "caudalMolar": 0.0,
                    "caudalUnitarioMolar": [],
                    "caudalUnitarioMasico": [],
                    "fraccionMolar": [],
                    "fraccionMasica": []}

    def __init__(self, tipo=0, **kwargs):
        """
        tipo: definición de la corriente
            0   -   No definida
            1   -   caudales unitarios másicos
            2   -   caudales unitarios molares
            3   -   caudal másico y fracción molar
            4   -   caudal másico y fracción másica
            5   -   caudal molar y fracción molar
            6   -   caudal molar y fracción másica
        kwargs:
            fraccionMolar
            fraccionMasica
            caudalMasico
            caudalMolar
            caudalUnitarioMasico
            caudalUnitarioMolar
        """
        self.kwargs.update(kwargs)
        if "ids" in self.kwargs:
            self.ids=self.kwargs.get("ids")
        else:
            Config=config.getMainWindowConfig()
            txt=Config.get("Components", "Components")
            if isinstance(txt, str):
                self.ids=eval(txt)
            else:
                self.ids=txt
        self.componente=[Componente(int(i)) for i in self.ids]
        fraccionMolar=self.kwargs.get("fraccionMolar", None)
        fraccionMasica=self.kwargs.get("fraccionMasica", None)
        caudalMasico=self.kwargs.get("caudalMasico", None)
        caudalMolar=self.kwargs.get("caudalMolar", None)
        caudalUnitarioMasico=self.kwargs.get("caudalUnitarioMasico", None)
        caudalUnitarioMolar=self.kwargs.get("caudalUnitarioMolar", None)

        # normalizar fracciones a la unidad
        if fraccionMolar:
            suma=float(sum(fraccionMolar))
            fraccionMolar=[x/suma for x in fraccionMolar]
        if fraccionMasica:
            suma=float(sum(fraccionMasica))
            fraccionMasica=[x/suma for x in fraccionMasica]

        if tipo==1:
            caudalMasico=sum(caudalUnitarioMasico)
            caudalUnitarioMolar=[mass/componente.M for mass, componente in zip(caudalUnitarioMasico, self.componente)]
            caudalMolar=sum(caudalUnitarioMolar)
            fraccionMolar=[mi/caudalMolar for mi in caudalUnitarioMolar]
            fraccionMasica=[mi/caudalMasico for mi in caudalUnitarioMasico]
        elif tipo==2:
            caudalMolar=sum(caudalUnitarioMolar)
            caudalUnitarioMasico=[mol*componente.M for mol, componente in zip(caudalUnitarioMolar, self.componente)]
            caudalMasico=sum(caudalUnitarioMasico)
            fraccionMolar=[mi/caudalMolar for mi in caudalUnitarioMolar]
            fraccionMasica=[mi/caudalMasico for mi in caudalUnitarioMasico]
        elif tipo==3:
            pesos=[x*componente.M for x, componente in zip(fraccionMolar, self.componente)]
            M=sum(pesos)
            caudalMolar=caudalMasico/M
            fraccionMasica=[peso/M for peso in pesos]
            caudalUnitarioMasico=[x*caudalMasico for x in fraccionMasica]
            caudalUnitarioMolar=[x*caudalMolar for x in fraccionMolar]
        elif tipo==4:
            caudalUnitarioMasico=[caudalMasico*x for x in fraccionMasica]
            caudalUnitarioMolar=[mass/componente.M for mass, componente in zip(caudalUnitarioMasico, self.componente)]
            caudalMolar=sum(caudalUnitarioMolar)
            fraccionMolar=[mi/caudalMolar for mi in caudalUnitarioMolar]
        elif tipo==5:
            caudalUnitarioMolar=[caudalMolar*x for x in fraccionMolar]
            caudalUnitarioMasico=[mol*componente.M for mol, componente in zip(caudalUnitarioMolar, self.componente)]
            caudalMasico=sum(caudalUnitarioMasico)
            fraccionMasica=[mi/caudalMasico for mi in caudalUnitarioMasico]
        elif tipo==6:
            moles=[x/componente.M for x, componente in zip(fraccionMasica, self.componente)]
            M=sum(moles)
            caudalMasico=caudalMolar*M
            fraccionMolar=[mol/caudalMolar for mol in moles]
            caudalUnitarioMasico=[x*caudalMasico for x in fraccionMasica]
            caudalUnitarioMolar=[x*caudalMolar for x in fraccionMolar]

        self.fraccion_normalizada=[unidades.Dimensionless(f) for f in fraccionMolar]
        self.componente_normalizado=self.componente[:]
        self.zeros=[]
        for i in range(fraccionMolar.count(0)):
            ind=fraccionMolar.index(0)
            self.zeros.insert(0, ind)
            del fraccionMolar[ind]
            del self.componente[ind]

        self.fraccion=fraccionMolar
        self.caudalmasico=unidades.MassFlow(caudalMasico)
        self.caudalmolar=unidades.MolarFlow(caudalMolar)
        self.fraccion_masica=[unidades.Dimensionless(f) for f in fraccionMasica]
        self.caudalunitariomasico=[unidades.MassFlow(q) for q in caudalUnitarioMasico]
        self.caudalunitariomolar=[unidades.MolarFlow(q) for q in caudalUnitarioMolar]

        self.Mixing_Rule=[self.Mix_van_der_Waals, self.Mix_Stryjek_Vera, self.Mix_Panagiotopoulos, self.Mix_Melhem][config.getMainWindowConfig().getint("Thermo","Mixing")]

        if tipo==0:
            self._bool=False
            self.status=0
            return
        else:
            self._bool=True
            self.status=1

        self.M=unidades.Dimensionless(caudalMasico/caudalMolar)


        """Método de cálculo de la temperatura crítica, API procedure 4B1.1 pag 304"""
        Vpc=sum([fraccion*componente.Vc for fraccion, componente in zip(self.fraccion, self.componente)])
        k=[fraccion*componente.Vc/Vpc for fraccion, componente in zip(self.fraccion, self.componente)]
        Tcm=sum([k*componente.Tc for k, componente in zip(k, self.componente)])
        self.Tc=unidades.Temperature(Tcm)

        """Cálculo de la temperatura pseudocrítica usada en muchos procedimientos"""
        tpc=sum([fraccion*componente.Tc for fraccion, componente in zip(self.fraccion, self.componente)])
        self.tpc=unidades.Temperature(tpc)

        """Cálculo de la presion pseudocrítica usada en muchos procedimientos"""
        ppc=sum([fraccion*componente.Pc for fraccion, componente in zip(self.fraccion, self.componente)])
        self.ppc=unidades.Pressure(ppc)

        """Método de cálculo de la presión crítica, API procedure 4B2.1 pag 307"""
        sumaw=0
        for i in range(len(self.componente)):
            sumaw+=self.fraccion[i]*self.componente[i].f_acent
        self.Pc=unidades.Pressure(self.ppc+self.ppc*(5.808+4.93*sumaw)*(self.Tc-self.tpc)/self.tpc)

        """Método de cálculo de factor acentrico #eq 6B2.2-6 pag 523"""
        self.f_acent=sum([fraccion*componente.f_acent for fraccion, componente in zip(self.fraccion, self.componente)])
        self.f_acent_mod=sum([fraccion*componente.f_acent_mod for fraccion, componente in zip(self.fraccion, self.componente)])

        """Método de cálculo del volumen crítico, API procedure 4B3.1 pag 314"""
        sumaxvc23=sum([fraccion*componente.Vc**(2./3) for fraccion, componente in zip(self.fraccion, self.componente)])
        k=[]
        k=[fraccion*componente.Vc**(2./3)/sumaxvc23 for fraccion, componente in zip(self.fraccion, self.componente)]
        #TODO: Generalizar valor de C en función de la naturaleza de los componentes. De momento se supondrá que se trata siempre de hidrocarburos (C=0)
        C=0
        V=[[-1.4684*abs((i.Vc-j.Vc)/(i.Vc+j.Vc))+C for j in self.componente] for i in self.componente]
        v=[[V[i][j]*(componentei.Vc+componentej.Vc)/2. for j, componentej in enumerate(self.componente)] for i, componentei in enumerate(self.componente)]
        suma1=sum([ki*componente.Vc for ki, componente in zip(k, self.componente)])
        suma2=sum([ki*kj*v[i][j] for j, kj in enumerate(k) for i, ki in enumerate(k)])
        self.Vc=unidades.SpecificVolume((suma1+suma2)*self.M)

        self.Tb=unidades.Temperature(sum([fraccion*componente.Tb for fraccion, componente in zip(self.fraccion, self.componente)]))
        self.SG=sum([fraccion*componente.SG for fraccion, componente in zip(self.fraccion, self.componente)])

    def recallZeros(self, fraccion, val=0):
        for indice in self.zeros:
            fraccion.insert(indice, val)


    def tr(self,T):
       return T/self.tpc
    def pr(self,P):
        return P/self.ppc

    def Kij(self, T=0, EOS=None):
        """Cálculo de los coeficientes de interacción binarios, usando los coeficientes BIP para la ecuación de estado SRK, y el API procedure 8D1.1 pag 819, ecuaciones pag 827"""
        if EOS:
            kij=[]
            for i in self.ids:
                kijk=[]
                for j in self.ids:
                    if i==j:
                        kijk.append(0)
                    else:
                        for indice in EOS:
                            if i in indice[0:2] and j in indice[0:2]:
                                kijk.append(indice[2])
                                break
                        else:
                            if i==1 or j==1:
                                kijk.append(1/(344.23*exp(-0.48586*T/databank.base_datos[1][3])+1))
                            elif i==2 or j==2:
                                kijk.append(0.014*(abs(self.componente[self.ids.index(i)].parametro_solubilidad-self.componente[self.ids.index(j)].parametro_solubilidad)))
                            elif i==46 or j==46:
                                kijk.append(0.0403*(abs(self.componente[self.ids.index(i)].parametro_solubilidad-self.componente[self.ids.index(j)].parametro_solubilidad)))
                            elif i==48 or j==48:
                                kijk.append(0)
                            elif i==49 or j==49:
                                kijk.append(0.1)
                            elif i==50 or j==50:
                                kijk.append(0.0316*(abs(self.componente[self.ids.index(i)].parametro_solubilidad-self.componente[self.ids.index(j)].parametro_solubilidad)))
                            else: kijk.append(0)
                kij.append(kijk)
        else:
            kij=zeros((len(self.ids), len(self.ids)))
        return kij


    def _Critical_API(self):
        """Método de cálculo de las propiedades críticas, haciendo uso de la ecuación de estado de Soave-Redlich-Kwong Modificada Thorwart-Daubert, API procedure 4B4.1 pag 317"""

        def parameters(T):
            ai=[]
            bi=[]
            for componente in self.componente:
                a, b=eos.SRK_Thorwart_lib(componente, T)
                ai.append(a)
                bi.append(b)
            b=sum([fraccion*b for fraccion, b in zip(self.fraccion, bi)])

            k=self.Kij(srk)

            aij=[[(ai[i]*ai[j])**0.5*(1-k[i][j]) for j in range(len(self.componente))] for i in range(len(self.componente))]
            a=sum([fraccioni*fraccionj*aij[i][j] for j, fraccionj in enumerate(self.fraccion) for i, fraccioni in enumerate(self.fraccion)])
            a_=[2*sum([fraccion*a for fraccion, a in zip(self.fraccion, aij[i])]) for i in range(len(self.fraccion))]

            return ai, bi, b, k, aij, a, a_


        def q(T, V):
            """Subrutina de cálculo del determinante de Q, eq 4B4.1-5"""
            ai, bi, b, k, aij, a, a_=parameters(T)
            B1=[[2*a*bi[i]*bi[j]-b*(a_[i]*bi[j]+a_[j]*bi[i]) for j in range(len(self.fraccion))] for i in range(len(self.fraccion))]
            B2=[[-B1[i][j]-2*aij[i][j]*b**2 for j in range(len(self.fraccion))] for i in range(len(self.fraccion))]
            d=[]
            for i in range(len(self.componente)):
                d_=[]
                for j in range(len(self.componente)):
                    if i==j:
                        d_.append(1)
                    else: d_.append(0)
                d.append(d_)

            Q=[[R_atml*T*(d[i][j]/self.fraccion[i]+(bi[i]+bi[j])/(V-b)+bi[i]*bi[j]/(V-b)**2)+a*bi[i]*bi[j]/b/(V+b)**2+B1[i][j]/b**2/(V+b)+B2[i][j]/b**3*log((V+b)/V)  for j in range(len(self.componente))] for i in range(len(self.componente))]
            q=triu(Q)
            return q


        def C(T, V):
            """Subrutina de cálculo del parámetro C, eq 4B4.1-19"""
            Q=q(T, V)

            ai, bi, b, k, aij, a, a_=parameters(T)

            deltaN=[1]
            for i in range(len(self.componente)-1, 0, -1):
                deltaN.insert(0, -deltaN[0]*Q[i-1][i]/Q[i-1][i-1])
            suma=sum([n**2 for n in deltaN])**0.5
            deltaN=[n/suma for n in deltaN]

            h=[]
            for i in range(len(self.componente)):
                hi=[]
                for j in range(len(self.componente)):
                    hij=[]
                    for k in range(len(self.componente)):
                        if i==j and j==k: hij.append(1)
                        elif i!=j and i!=k and j!=k: hij.append(6)
                        else: hij.append(3)
                    hi.append(hij)
                h.append(hi)

            F=[[[b_i*b_j*b_k for b_k in bi] for b_j in bi] for b_i in bi]
            D=[[[b*a_[i]*bi[j]*bi[k]+a_[j]*bi[i]*bi[k]+a_[k]*bi[i]*bi[j]-3*F[i][j][k]*a for k in range(len(self.componente))] for j in range(len(self.componente))] for i in range(len(self.componente))]
            E=[[[2*D[i][j][k]-2*b**2*(aij[i][j]*bi[k]+aij[i][k]*bi[j]+aij[j][k]*bi[i]) for k in range(len(self.componente))] for j in range(len(self.componente))] for i in range(len(self.componente))]

            d=[]
            for i in range(len(self.componente)):
                d_=[]
                for j in range(len(self.componente)):
                    d__=[]
                    for k in range(len(self.componente)):
                        if i==j and i==k: d__.append(1)
                        else: d__.append(0)
                    d_.append(d__)
                d.append(d_)

            delta=[[[R_atml*T*(d[i][j][k]/self.fraccion[i]**2+(bi[i]*bi[j]+bi[j]*bi[k]+bi[k]*bi[i])/(V-b)**2+2*F[i][j][k]/(V-b)**3)-2*a*F[i][j][k]/b/(V+b)**3+D[i][j][k]/b**2/(V+b)**2+E[i][j][k]/b**3/(V+b)-E[i][j][k]/b**4*log((V+b)/V)  for k in range(len(self.componente))] for j in range(len(self.componente))] for i in range(len(self.componente))]

            sumaijk=sum([sum([sum([h[i][j][k]*delta[i][j][k]*deltaN[i]*deltaN[j]*deltaN[k] for k in range(len(self.componente))]) for j in range(len(self.componente))]) for i in range(len(self.componente))])

            return ((V-b)/2/b)**2*sumaijk


        To=1.5*self.tpc
        Vo=sum([fraccion*R_atml*componente.Tc/3/componente.Pc.atm for fraccion, componente in zip(self.fraccion, self.componente)])
        funcion = lambda par: (det(q(par[0], par[1])), C(par[0], par[1]))
        Tc, Vc=fsolve(funcion, [To, Vo])

        ai, bi, b, k, aij, a, a_=parameters(Tc)
        Pc=R_atml*Tc/(Vc-b)-a/(Vc*(Vc+b))

        VcCorr=sum([fraccion/(R_atml*componente.Tc/3/componente.Pc.atm-componente.Vc*componente.M) for fraccion, componente in zip(self.fraccion, self.componente)])

        return unidades.Temperature(Tc), unidades.Pressure(Pc, "atm"), unidades.SpecificVolume(VcCorr/self.M, "lg")





#Reglas de mezcla
    def Mix_van_der_Waals(self, parameters, kij):
        """Reglas de mezcla de van der Waals"""
        ai=parameters[0]
        bi=parameters[1:]
        b=[0]*len(bi)
        a=0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j]+=self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                a+=self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j])
        return  tuple([a]+b)

    def Mix_Stryjek_Vera(self, parameters, kij):
        """Reglas de mezcla de Stryjek and Vera (1986)"""
        ai=parameters[0]
        bi=parameters[1:]
        b=[0]*len(bi)
        a=0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j]+=self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                if kij[i][j]==0 and kij[j][i]==0:
                    k=0.
                else:
                    k=kij[i][j]*kij[j][i]/(self.fraccion[i]*kij[i][j]+self.fraccion[j]*kij[j][i])
                a+=self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-k)
        return tuple([a]+b)

    def Mix_Panagiotopoulos(self, parameters, kij):
        """Reglas de mezcla de Panagiotopoulos (1985)"""
        ai=parameters[0]
        bi=parameters[1:]
        b=[0]*len(bi)
        a=0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j]+=self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                a+=self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j]+(kij[i][j]-kij[j][i])*self.fraccion[i])
        return  tuple([a]+b)

    def Mix_Melhem(self, parameters, kij):
        """Reglas de mezcla de Melhem (1991)"""
        ai=parameters[0]
        bi=parameters[1:]
        b=[0]*len(bi)
        a=0
        for i in range(len(self.componente)):
            for j in range(len(bi)):
                b[j]+=self.fraccion[i]*bi[j][i]
            for j in range(len(self.componente)):
                a+=self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j]+(kij[i][j]-kij[j][i])*self.fraccion[i]/(self.fraccion[i]+self.fraccion[j]))
        return  tuple([a]+b)


    def Lee_Kesler_Entalpia(self, T, P):
        """Método de cálculo de la entalpía haciendo uso de las propiedades críticas, método de Lee-Kesler
        Procedure API 7B3.7 Pag.643"""
        P=unidades.Pressure(P, "atm")
        H_adimensional=self.Lee_Kesler_Entalpia_lib(T, P.atm, self.Fase(T, P.atm))
        H_ideal=self.Entalpia_ideal(T)
        return unidades.Enthalpy(H_ideal.Jg-H_adimensional*R*self.Tc/self.M, "Jg")


    def Entalpia_ideal(self, T):
        """Regla de mezcla de la entalpia ideal"""
        entalpia=0
        for i in range(len(self.componente)):
            entalpia+=self.fraccion_masica[i]*self.componente[i].Entalpia_ideal(T)
        return unidades.Enthalpy(entalpia)

    def Hv_DIPPR(self, T):
        entalpia=0
        for i in range(len(self.componente)):
            if self.componente[i].calor_vaporizacion[-1]<T:
                Hi=self.componente[i].Hv_DIPPR(self.componente[i].calor_vaporizacion[-1])
            elif self.componente[i].calor_vaporizacion[-2]>T:
                Hi=self.componente[i].Hv_DIPPR(self.componente[i].calor_vaporizacion[-2])
            else:
                Hi=self.componente[i].Hv_DIPPR(T)
            entalpia+=self.fraccion[i]*Hi
        return unidades.Enthalpy(entalpia/self.M, "Jkg")


    def Entalpia(self, T, P):
        """Método de cálculo de la entalpía, API procedure 7B4.1, pag 645"""
        Ho=self.Entalpia_ideal(T)
        #TODO:


    def Calor_vaporizacion(self):
        """Método de cálculo del calor de vaporización, API procedure 7C2.1, pag 682"""
        pass


    def Cp_Liquido(self, T):
        """Método de cálculo de la capacidad calorífica del líquido a presión constante, API procedure 7D2.1, pag 700"""
        Cp=0
        for i, componente in enumerate(self.componente):
            Cp+=self.fraccion_masica[i]*componente.Cp_Liquido_DIPPR(T)
        return unidades.SpecificHeat(Cp)


    def Cp_Gas(self, T, P):
        """Método de cálculo de la capacidad calorífica del gas a presión constante, API procedure 7D4.1, pag 714"""
        Cp=0
        for i, componente in enumerate(self.componente):
            Cp+=self.fraccion_masica[i]*componente.Cp_Gas_DIPPR(T)
        return unidades.SpecificHeat(Cp)
        #TODO: Añadir correción a alta presión


    def Cp_v_ratio(self):
        """Calculo de la relación de capacidades caloríficas
        Procedure API 7E2.1 Pag.728"""
        pass


    def Entropia(self, T):
        """Método de cálculo de la entropia
        Procedure API 7F2.1, Pag.741"""
        pass


    def flash_water(self, T, P):
        """Método de cálculo del equilibrio líquido-vapor en sistemas hidrocarburo-agua, API procedure 9A6.1, pag 922"""
        #TODO:
        #Volumen molar
        a=1
        b=1
        V=roots([1, -R_atml*T/P, -(b**2+R_atml*T/P*b-a/P), a*b/P])

        return V


    def SOUR_water(self, T):
        """Método de cálculo del equilibrio líquido-vapor en sistemas H2O-NH3-H2S, API procedure 9A7.3 pag 930, en chemcad modelo k: sour water"""

        t=Temperature(T)

        if 63 in self.ids:
            amoniaco=True
            indice_amoniaco=self.ids.index(63)
        else: amoniaco=False
        if 50 in self.ids:
            sulfuro=True
            indice_sulfuro=self.ids.index(50)
        else: sulfuro=False
        if sulfuro and amoniaco:
            ternario=True
        else: ternario=False

        if ternario:
            ppmnh3=log10(self.fraccion_masica[indice_amoniaco]*1e6)
            ppmh2s=log10(self.fraccion_masica[indice_sulfuro]*1e6)
            ratio=self.fraccion[indice_amoniaco]/self.fraccion[indice_sulfuro]
            print 10**ppmh2s, ratio

            def Z(ind, M=ratio):
                """Subrutina que nos devuelve el valor de Zm"""
                if ind==1:
                   return 1/ratio**2
                elif ind==2:
                    return log10(ratio)
                else:
                    return sin(ratio)/ratio

            #Desglose de la tabla 9A7.4
            if t.R<537.7:
                if ratio>0.7:
                    parnh3=[Z(1), 0.02998, 0.838848, -1.34384, -2.8802]
                else:
                    parnh3=[Z(2), 0.13652, 0.074081, 1.12893, -3.8683]
                if ratio<1.:
                    parh2s=[Z(3), -0.00562, 1.028528, 8.33185, -8.9175]
                elif ratio<4.:
                    parh2s=[Z(2), 0.00227, 0.980646, -2.82752, -2.3105]
                else:
                    parh2s=[Z(2), 0.00138, 0.881842, -1.08346, -2.7605]
            elif t.R<559.7:
                if ratio>0.7:
                    parnh3=[Z(1), 0.03894, 0.775221, -1.17860, -2.5972]
                else:
                    parnh3=[Z(2), 0.06415, 0.659240, 1.24672, -4.5124]
                if ratio<1.:
                    parh2s=[Z(3), -0.00487, 1.025371, 7.46512, -7.9480]
                elif ratio<4.:
                    parh2s=[Z(2), 0.01018, 0.915687, -2.65953, -1.9435]
                else:
                    parh2s=[Z(2), 0.00668, 0.849756, -1.10875, -2.4133]
            elif t.R<579.7:
                if ratio>0.7:
                    parnh3=[Z(1), 0.02393, 0.881656, -1.08692, -2.5668]
                else:
                    parnh3=[Z(2), 0.09581, 0.420378, 1.30988, -3.7033]
                if ratio<1.:
                    parh2s=[Z(3), 0.00308, 0.976347, 6.81079, -7.1893]
                elif ratio<4.:
                    parh2s=[Z(2), 0.00562, 0.932481, -2.48535, -1.6691]
                else:
                    parh2s=[Z(2), 0.01577, 0.778382, -1.07217, -1.9946]
            elif t.R<609.7:
                if ratio>0.7:
                    parnh3=[Z(1), 0.02580, 0.856844, -0.90996, -2.2738]
                else:
                    parnh3=[Z(2), 0.07081, 0.625739, 1.23245, -3.5797]
                if ratio<0.9:
                    parh2s=[Z(2), -0.02064, 1.090006, -0.66045, -1.0928]
                elif ratio<2.9:
                    parh2s=[Z(2), 0.00144, 0.961073, -3.00981, -1.2917]
                else:
                    parh2s=[Z(2), 0.01306, 0.813499, -1.15162, -1.5845]
            elif t.R<659.7:
                if ratio>1.6:
                    parnh3=[Z(2), 0.06548, 0.621031, 0.24106, -1.9006]
                elif ratio>1.1:
                    parnh3=[Z(2), 0.05838, 0.683411, 1.95035, -2.4470]
                else:
                    parnh3=[Z(2), -0.07626, 1.462759, 2.57458, -3.5185]
                if ratio<1.5:
                    parh2s=[Z(2), 0.00398, 0.997000, -1.95892, -1.1374]
                else:
                    parh2s=[Z(2), -0.00651, 1.037056, -1.24395, -1.4799]
            elif t.R<679.7:
                if ratio>1.5:
                    parnh3=[Z(2), 0.05840, 0.647857, 0.28668, -1.8320]
                elif ratio>1.:
                    parnh3=[Z(2), 0.02447, 0.864904, 2.60250, -2.5186]
                else:
                    parnh3=[Z(2), -0.08280, 1.479719, 2.26608, -3.3007]
                if ratio<1.5:
                    parh2s=[Z(2), 0.00892, 0.957002, -1.84465, -0.9934]
                else:
                    parh2s=[Z(2), 0.00674, 0.910272, -1.16113, -1.0837]
            elif t.R<699.7:
                if ratio>1.6:
                    parnh3=[Z(2), 0.04542, 0.756870, 0.18305, -1.8041]
                elif ratio>1.:
                    parnh3=[Z(2), 0.03659, 0.790904, 2.00739, -2.2410]
                else:
                    parnh3=[Z(2), -0.07386, 1.445958, 1.79890, -3.1476]
                parh2s=[Z(2), 0.0, 1.002832, -1.40994, -0.9558]
            elif t.R<719.7:
                if ratio>1.6:
                    parnh3=[Z(2), 0.04603, 0.728930, 0.17864, -1.5652]
                elif ratio>1.:
                    parnh3=[Z(2), 0.03493, 0.796292, 1.69221, -2.0369]
                else:
                    parnh3=[Z(2), -0.07354, 1.437440, 1.73508, -2.8916]
                if ratio<1.5:
                    parh2s=[Z(2), 0.00871, 0.964568, -1.44142, -0.7861]
                else:
                    parh2s=[Z(2), -0.01035, 1.054256, -1.22074, -0.9966]
            else:
                if ratio>1.5:
                    parnh3=[Z(2), 0.04484, 0.728061, 0.25132, -1.5543]
                elif ratio>0.9:
                    parnh3=[Z(2), 0.00193, 0.998553, 1.87924, -2.1925]
                else:
                    parnh3=[Z(2), -0.06773, 1.391080, 1.68263, -2.6696]
                parh2s=[Z(2), 0.0, 1.000692, -1.23871, -0.7932]

            ppnh3=10**(parnh3[1]*ppmnh3**2+parnh3[2]*ppmnh3+parnh3[3]*parnh3[0]+parnh3[4])
            pph2s=10**(parh2s[1]*ppmh2s**2+parh2s[2]*ppmh2s+parh2s[3]*parh2s[0]+parh2s[4])
            return Pressure(ppnh3, "mmHg"), Pressure(pph2s, "mmHg")

        else:
            if amoniaco:
                ppm=log10(self.fraccion_masica[indice_amoniaco]*1e6)
                pp=10**(0.01129*ppm**2+0.9568*ppm+0.00719*t.R-6.83498)

            else:
                ppm=log10(self.fraccion_masica[indice_sulfuro]*1e6)
                pp=10**(-0.00322*ppm**2+1.0318*ppm+0.00248*t.R-1.95729)

            return Pressure(pp, "mmHg")


    def SOUR_water_ph(self, T):
        """Método de cálculo del ph de disoluciones de NH3-H2S en agua, API procedure 9A8.1 pag 934"""

        t=Temperature(T)

        if 63 in self.ids:
            amoniaco=True
            indice_amoniaco=self.ids.index(63)
        else: amoniaco=False
        if 50 in self.ids:
            sulfuro=True
            indice_sulfuro=self.ids.index(50)
        else: sulfuro=False
        if sulfuro and amoniaco:
            ternario=True
        else: ternario=False

        if ternario:
            M=self.fraccion[indice_amoniaco]/self.fraccion[indice_sulfuro]
            pH=10**(-0.0084*log10(log10(M))+0.1129*log10(M)-0.00032*T.R+1.0689)
        else:
            if amoniaco:
                ppm=log10(self.fraccion_masica[indice_amoniaco]*1e6)
                pH=10**(0.480*ppm-0.0106*t.R-15.236)

            else:
                ppm=log10(self.fraccion_masica[indice_sulfuro]*1e6)
                pH=10**(-0.505*ppm-0.0019*t.R+6.807)

        return pH

    def RhoL_Rackett(self, T):
        """Cálculo de de la densidad del líquido en su punto de burbuja mediante el método de Rackett, procedure API 6A3.1 pag.479
        Valor obtenido en mol/l"""

        # eq 6A3.1-2
        Zram=0
        for i in range(len(self.componente)):
            Zram+=self.fraccion[i]*self.componente[i].rackett

        # eq 6A3.1-5
        suma=0
        for i in range(len(self.componente)):
            suma+=self.fraccion[i]*self.componente[i].Vc
        fi=[]
        for i in range(len(self.componente)):
            fi.append(self.fraccion[i]*self.componente[i].Vc/suma)

        # eq 6A3.1-7
        k=[]
        for i in self.componente:
            ki=[]
            for j in self.componente:
                ki.append(1-(sqrt(i.Vc**(1./3)*j.Vc**(1./3))*2/(i.Vc**(1./3)+j.Vc**(1./3)))**3)
            k.append(ki)

        # eq 6A3.1-6
        Tc=[]
        for i in range(len(self.componente)):
            Tci=[]
            for j in range(len(self.componente)):
                Tci.append(sqrt(self.componente[i].Tc*self.componente[j].Tc)*(1-k[i][j]))
            Tc.append(Tci)

        # eq 6A3.1-4
        Tmc=0
        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
                Tmc+=fi[i]*fi[j]*Tc[i][j]

        Tr=T/Tmc        #eq 6A3.1-3

        #eq 6A3.1-1
        suma=0
        for i in range(len(self.componente)):
            suma+=self.fraccion[i]*self.componente[i].Tc/self.componente[i].Pc
        inv=R_atml*suma*Zram**(1+(1-Tr)**(2./7))

        return unidades.Density(1/inv*self.M, "gl")

    def _lib_Costald(self):
        #eq 6A3.1-2
        suma1=0
        suma2=0
        suma3=0
        for i in range(len(self.componente)):
            suma1+=self.fraccion[i]*self.componente[i].Vc
            suma2+=self.fraccion[i]*self.componente[i].Vc**(2./3)
            suma3+=self.fraccion[i]*self.componente[i].Vc**(1./3)
        Vm=(suma1+3*suma2*suma3)/4

        #eq 6A3.1-5
        suma=0
        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
                suma+=self.fraccion[i]*self.fraccion[j]*sqrt(self.componente[i].Vc*self.componente[j].Vc*self.componente[i].Tc*self.componente[j].Tc)
        Tmc=suma/Vm
        return Vm, Tmc


    def RhoL_Costald(self, T):
        """Cálculo de de la densidad del líquido en su punto de burbuja mediante el método de Costald, procedure API 6A3.2 pag.482
        Valor obtenido en mol/l"""

        Vm, Tmc=self._lib_Costald()
        Tr=T/Tmc    #eq 6A3.2-4
        Vr0=1-1.52816*(1-Tr)**(1./3)+1.43907*(1-Tr)**(2./3)-0.81446*(1-Tr)+0.190454*(1-Tr)**(4./3)        #eq 6A3.2-7
        Vrd=(-0.296123+0.386914*Tr-0.0427258*Tr**2-0.0480645*Tr**3)/(Tr-1.00001)    #eq 6A3.2-8

        #eq 6A3.2-1
        return unidades.Density(1/(Vm*Vr0*(1-self.f_acent_mod*Vrd))*self.M, "gl")


    def RhoL_Tait_Costald(self, T, P):
        """Cálculo de de la densidad del líquido comprimido mediante el método de Tait-Costald, procedure API 6A3.4 pag.489
        Presión dada en pascales
        Densidad obtenida en mol/l"""
        #FIXME. No sale
        densidad_bp=self.RhoL_Costald(T)
        Vm, Tmc=self._lib_Costald()
        Tr=T/Tmc
        zmc=0.291-0.08*self.f_acent_mod
        pmc=zmc*R_atml*Tmc/Vm

        #eq 6A3.4-12
        alfa=35.-36./Tr-96.736*log10(Tr)+Tr**6
        beta=log10(Tr)+0.03721754*alfa
        prm0=5.8031817*log10(Tr)+0.07608141*alfa
        prm1=4.86601*beta
        presion_bp=10**(prm0+self.f_acent_mod*prm1)*pmc

        #eq 6A3.4-2
        C=0.0861488+0.0344483*self.f_acent_mod
        e=exp(4.79594+0.250047*self.f_acent_mod+1.14188*self.f_acent_mod**2)
        B=pmc*(-1-9.070217*(1-Tr)**(1./3)+62.45326*(1-Tr)**(2./3)-135.1102*(1-Tr)+e*(1-Tr)**(4./3))
        return unidades.Density(densidad_bp/(1-C*log((B+P)/(B+presion_bp))))


    def RhoL_API(self, T, P):
        """Método de cálculo de la densidad de líquidos, API procedure 6A3.3, pag 485"""
        Tcm=Pcm=0
        for i in range(len(self.componente)):
            Tcm+=self.fraccion[i]*self.componente[i].Tc
            Pcm+=self.fraccion[i]*self.componente[i].Pc

#FIXME: Añadir corrección en componentes que son gases a 1 atm y 60ºF

        suma1=suma2=0
        for i in range(len(self.componente)):
            suma1+=self.fraccion[i]*self.componente[i].M
            suma2+=self.fraccion[i]*self.componente[i].M/self.componente[i].SG/1000
        rho1=suma1/suma2

        A02=1.6368-0.04615*self.pr(P, Pcm)+2.1138e-3*self.pr(P, Pcm)**2-0.7845e-5*self.pr(P, Pcm)**3-0.6923e-6*self.pr(P, Pcm)**4
        A12=-1.9693-0.21874*self.pr(P, Pcm)-8.0028e-3*self.pr(P, Pcm)**2-8.2328e-5*self.pr(P, Pcm)**3+5.2604e-6*self.pr(P, Pcm)**4
        A22=2.4638-0.36461*self.pr(P, Pcm)-12.8763e-3*self.pr(P, Pcm)**2+14.8059e-5*self.pr(P, Pcm)**3-8.6895e-6*self.pr(P, Pcm)**4
        A32=-1.5841-0.25136*self.pr(P, Pcm)-11.3805e-3*self.pr(P, Pcm)**2+9.5672e-5*self.pr(P, Pcm)**3+2.1812e-6*self.pr(P, Pcm)**4
        C2=A02+A12*self.tr(T, Tcm)+A22*self.tr(T, Tcm)**2+A32*self.tr(T, Tcm)**3
        A01=1.6368-0.04615*self.pr(1, Pcm)+2.1138e-3*self.pr(1, Pcm)**2-0.7845e-5*self.pr(1, Pcm)**3-0.6923e-6*self.pr(1, Pcm)**4
        A11=-1.9693-0.21874*self.pr(1, Pcm)-8.0028e-3*self.pr(1, Pcm)**2-8.2328e-5*self.pr(1, Pcm)**3+5.2604e-6*self.pr(1, Pcm)**4
        A21=2.4638-0.36461*self.pr(1, Pcm)-12.8763e-3*self.pr(1, Pcm)**2+14.8059e-5*self.pr(1, Pcm)**3-8.6895e-6*self.pr(1, Pcm)**4
        A31=-1.5841-0.25136*self.pr(1, Pcm)-11.3805e-3*self.pr(1, Pcm)**2+9.5672e-5*self.pr(1, Pcm)**3+2.1812e-6*self.pr(1, Pcm)**4
        C1=A01+A11*self.tr(288.71, Tcm)+A21*self.tr(288.71, Tcm)**2+A31*self.tr(288.71, Tcm)**3
        return unidades.Density(rho1*C2/C1)


    def Tension(self,T):
        """Método de cálculo de la tensión superficial a baja presión, API procedure 10A2.1, pag 991"""
        tension=sum([fraccion*componente.Tension(T) for fraccion, componente in zip(self.fraccion, self.componente)])
        return unidades.Tension(tension)


    def Tension_superficial_presion(self,T, parachor, fraccion_liquido, fraccion_vapor):
        """Método de cálculo de la tensión superficial a alta presión, API procedure 10A2.2, pag 993
        Este procedimiento tiene un uso interno, ya que necesita como parámetros además las fracciones molares de los componentes en ambas fases, datos provenientes de un cálculo flash previo"""
        #TODO: parachor tiene que ser indicado como parámetro mientras no sepa como calcularlo para cada elemento
        mv=ml=suma=0
        for i in range(len(self.componente)):
            mv+=self.componente[i].M*fraccion_vapor[i]
            ml+=self.componente[i].M*fraccion_liquido[i]
        rhoV=RhoG_Lee_Kesler(T, 1)
        suma=0
        for i in range(len(self.componente)):
            suma+=parachor[i]*(self.RhoL_Rackett(T)/ml*fraccion_liquido[i]-rhoV/mv*fraccion_vapor[i])

        return unidades.Tension(suma**4, "dyncm")


    def Tension_inferfacial_water(self, T):
        """Método de cálculo de la tensión interfacial entre agua e hidrocarburos, API procedure 10B1.3, pag 1007"""
        agua=Componente(62)
        sigma_w=agua.Tension_parametrica(T).dyncm
        sigma_h=self.Tension_superficial(T).dyncm
        return unidades.Tension(sigma_h+sigma_w-1.1*sqrt(sigma_h*sigma_w), "dyncm")


    def Mu_Liquido(self, T, P):
        """Método de cálculo de la viscosidad de líquidos, API procedure 11A3.1, pag 1051"""
        suma=0
        for i in range(len(self.componente)):
            suma+=self.fraccion[i]*self.componente[i].Mu_Liquido(T, P)**(1./3)
        return unidades.Viscosity(suma**3)


    def Mu_Gas_Wilke(self, T):
        """Método de cálculo de la viscosidad de gases, API procedure 11B2.1, pag 1102"""

        mui=[i.Mu_Gas(T, 1) for i in self.componente]

        kij=[]
        for i in range(len(self.componente)):
            kiji=[]
            for j in range(len(self.componente)):
                if i==j:
                    kiji.append(0)
                else:
                    kiji.append((1+(mui[i]/mui[j])**0.5*(self.componente[j].M/self.componente[i].M)**0.25)**2/8**0.5/(1+self.componente[i].M/self.componente[j].M)**0.5)
            kij.append(kiji)

        suma=[]
        for i in range(len(self.componente)):
            sumai=0
            if self.fraccion[i]!=0:
                for j in range(len(self.componente)):
                    sumai+=kij[i][j]*self.fraccion[j]/self.fraccion[i]
            suma.append(sumai)
        mu=0
        for i in range(len(self.componente)):
            mu+=mui[i]/(1.+suma[i])

        return unidades.Viscosity(mu)


    def Mu_Gas_Stiel(self, T, P, rhoG=0, muo=0):
        """Método de cálculo de la viscosidad de gases a alta presión, API procedure 11B4.1, pag 1107"""
        if muo == 0:
            muo=self.Mu_Gas_Wilke(T)
        if rhoG == 0:
            rhoG=P/self.Z/R_atml/T
        x=self.tpc**(1.0/6)/self.M**0.5/self.ppc**(2.0/3)
        rhor=rhoG*self.Vc/self.M
        mu=muo.cP+10.8e-5*(exp(1.439*rhor)-exp(-1.11*rhor**1.858))/x
        return unidades.Viscosity(mu, "cP")


    def Mu_Gas_Carr(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de no hidrocarburos gaseosos a alta presión, API procedure 11C1.2, pag 1113"""

        Tr=T/self.tpc
        Pr=P/self.ppc
        muo=self.Mu_Gas_Wilke(T)

        A1=83.8970*Tr**0.0105+0.6030*Tr**-0.0822+0.9017*Tr**-0.12-85.3080
        A2=1.514*Tr**-11.3036+0.3018*Tr**-0.6856+2.0636*Tr**-2.7611
        k=A1*1.5071*Pr**-0.4487+A2*(11.4789*Pr**0.2606-12.6843*Pr**0.1773+1.6953*Pr**-0.1052)
        return unidades.Viscosity(muo*k)

    def Mu_Gas(self, T, P):
        if P<2:
            return self.Mu_Gas_Wilke(T)
        else:
            return self.Mu_Gas_Stiel(T, P)
        #TODO: Añadir método de Carr cuando en la base de datos estén disponibles las clases químicas


    def ThCond_Liquido(self, T, P):
        """Método de cálculo de la conductividad térmica de líquidos, API procedure 12A2.1, pag 1145"""
        ki=[]
        for i in range(len(self.componente)):
            ki.append(self.componente[i].ThCond_Liquido(T, P))
        kij=[]
        for i in range(len(self.componente)):
            kiji=[]
            for j in range(len(self.componente)):
                kiji.append(2/(1/ki[i]+1/ki[j]))
            kij.append(kiji)
        xV=0
        Vi=[]
        for i in range(len(self.componente)):
            Vi.append(self.componente[i].M/self.componente[i].RhoL(T, P))
            xV+=self.fraccion[i]*Vi[i]
        fi_bruto=[]
        suma=0
        for i in range(len(self.componente)):
            fi_bruto.append(self.fraccion[i]*Vi[i]/xV)
            suma+=fi_bruto[i]
        fi=[]
        for i in range(len(self.componente)):
            fi.append(fi_bruto[i]/suma)
        k=0
        for i in range(len(self.componente)):
            for j in range(len(self.componente)):
               k+=fi[i]*fi[j]*kij[i][j]

        return unidades.ThermalConductivity(k)


    def ThCond_Gas(self, T, P):
        """Método de cálculo de la conductividad térmica de líquidos, API procedure 12A2.1, pag 1145"""
        ki=[]
        S=[]
        mu=[]
        for i in range(len(self.componente)):
            ki.append(self.componente[i].ThCond_Gas(T, P))
            if self.componente[i].indice==1:
                S.append(78.8888889)
            else:
                S.append(1.5*self.componente[i].Tb)
            mu.append(self.componente[i].Mu_Gas(T, P))

        Sij=[]
        for i in range(len(self.componente)):
            Siji=[]
            for j in range(len(self.componente)):
                Siji.append(sqrt(S[i]*S[j]))
            Sij.append(Siji)

        Aij=[]
        for i in range(len(self.componente)):
            Aiji=[]
            for j in range(len(self.componente)):
                Aiji.append(0.25*(1+sqrt(mu[i]/mu[j]*(self.componente[j].M/self.componente[i].M)**0.75*(1+S[i]/T)/(1+S[j]/T)))**2*(1+Sij[i][j]/T)/(1+S[i]/T))
            Aij.append(Aiji)

        sumaj=[]
        for i in range(len(self.componente)):
            suma=0
            for j in range(len(self.componente)):
                suma+=Aij[i][j]*self.fraccion[j]
            sumaj.append(suma)

        k=0
        for i in range(len(self.componente)):
            k+=ki[i]*self.fraccion[i]/sumaj[i]

        return unidades.ThermalConductivity(k)

    def Solubilidad_agua(self, T):
        """Método de cálculo de la solubilidad de agua en el componente, API procedure 9A1.1, 9A1.5, pag 897
        Solubilidad obtenida expresada en fracción molar"""
        #TODO: Rehacer para la última versión de API
        t=unidades.Temperature(T)
        modo1=[4, 6, 5, 8, 7, 10, 55, 11, 541, 12, 82, 14, 38, 60, 23, 27, 35, 28, 40, 41, 45, 42, 43, 71, 861, 885, 178]
        if self.indice in modo1:
            a1=[-6.158, -0.792, 10.289, -6.46, -6.361, -3.441, 11.511, -0.082, 12.831, 0.386, 1.934, 0.998, -3.615, 0.579, 16.123, 3.489, 1.202, 2.041, -0.595, -1.271, -0.587, 0.382, -1.024, 20.439, 0.879, 0.905, -21.475]
            a2=[-964, -2540, -6054, -572, -1007, -1306, -5918, -2428, -6281, -2536, -2939, -2689, -1263, -2569, -7437, -3270, -2661, -2489, -1591, -1386, -1631, -1896, -1705, -7147, -2192, -2209, 3902]
            a3=[0.8338e-2, 0.3597e-2, -0.464e-2, 0.7548e-2, 0.9174e-2, 0.4815e-2, -0.679e-2, 0.2379e-2, -0.789e-2, 0.1958e-2, 0.076e-2, 0.1416e-2, 0.4778e-2, 0.1654e-2, -1.052e-2, -0.073e-2, 0.1262e-2, -0.017e-2, 0.198e-2, 0.2451e-2, 0.1955e-2, 0.1123e-2, 0.2584e-2, -1.8040e-2, 0.1007e-2, 0.1010e-2, 2.0184e-2]
            indice=modo1.index(self.indice)
            return 10**(a1[indice]+a2[indice]/t.R+a3[indice]*t.R)
        else:
            H=self.composicion_molecular[1][self.composicion_molecular[0].index("H")]*Elemental.get_element(1).atomic_mass.value
            C=self.composicion_molecular[1][self.composicion_molecular[0].index("C")]*Elemental.get_element(6).atomic_mass.value
            return 10**(-(4200*H/C+1050)*(1.8/t.R-0.0016))


    def Solubilidad_en_agua(self, T):
        """Método de cálculo de la solubilidad en agua del componente en condiciones de equilibrio liquido-vapor, API procedure 9A2.1, 9A2.3, pag 899, método renovado en API v7
        Solubilidad obtenida expresada en fracción molar"""
        t=unidades.Temperature(T)
        modo1=[4, 6, 8, 10, 11, 12, 13, 36, 38, 39, 60, 389, 24, 35, 83, 40, 41, 45, 70, 377, 1486, 886, 71, 77, 862, 885, 178, 191, 1497]
        modo2=[185, 193, 194, 197, 409, 200, 203, 204, 206, ]
        if self.indice in modo1:
            a1=[-314.3825, -292.4658, -361.7959, -406.5873, -430.4197, -450.7616, -469.842, -235.6808, -326.816, -356.3782, -385.8427, -444.1288, -228.836, -291.413, -366.146, -195.5673, -218.6392, -241.8211, -330.4825, -375.6332, -420.7929, -465.8706, -262.5274, -248.0811, -311.3553, -308.4157, -288.0113, -220.5168, -247.3787]
            a2=[24225.156, 21452.23, 26167.45, 29388.83, 31018.136, 32355.695, 33782.08, 16843.48, 23264.01, 25331.92, 27399.83, 31535.66, 16995.95, 20436.66, 25673.53, 13544.694, 15119.661, 16694.626, 22994.478, 26144.406, 29294.352, 32444.28, 17838.556, 16349.596, 21748.495, 20014.725, 19278.918, 12603.309, 15072.32]
            a3=[41.50287, 38.58148, 47.97436, 53.89582, 56.95927, 59.55451, 61.94, 30.89572, 43.298, 47.1467, 50.9954, 58.6928, 30.0380, 38.4871, 48.3494, 25.8585, 28.8653, 31.8721, 43.8994, 49.913, 55.9266, 61.9402, 34.67748, 32.77122, 41.11214, 40.74593, 38.54914, 29.35657, 32.71111]
            indice=modo1.index(self.indice)
            return exp(a1[indice]+a2[indice]/t.R+a3[indice]*log(t.R))
        elif self.indice in modo2:
            a1=[-373.24613, -334.02054, -332.69462, -380.76757, -420.02594, -350.23574, -0.84491, -380.93125, -7.53314]
            a2=[21618.93, 17740.54, 17633.99, 20361.3, 21096.86, 17490.54, -9069.58, 19440.39, -8120.29]
            a3=[51.01156, 45.56672, 45.48078, 52.08392, 5753568, 47.9872, 0, 51.93974, 0]
            indice=modo2.index(self.indice)
            return exp(a1[indice]+a2[indice]/t.R+a3[indice]*log(t.R))


    def Solubilidad_en_agua_25(self):
        """Método de cálculo de la solubilidad en agua saturada a 25ºC, API procedure 9A2.6, pag 909"""
        #TODO: Díficil de implementar mientras no se añada a la base de datos alguna propiedad que indique la naturaleza química del compuesto

    def Solubilidad_Henry(self,T, P):
        """Solubilidad de gases en líquidos usando la ley de Henry
        Temperatura dada en ºR
        constante H obtenida en psia por unidad de fracción molar del gas
        lnH = A/T + B∗lnT + C∗T + D
        Solo disponible para algunos compuestos:
        Hydrogen, Helium, Argon, Neon, Krypton, Xenon, Oxygen, Nitrogen,
        Hydrogen sulfide, Carbon monoxide, Carbon dioxide, Sulfur dioxide,
        Nitrous oxide, Chlorine,Bromine, Iodine, Methane, Ethane, Propane,
        Ethylene, Ammonia.
        API procedure 9A7.1, pag 927
        Los parametros de la ecuación se encuentran en la base de datos
        en forma de lista en la posición décima
        Solubilidad obtenida en fracción molar"""
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        H=exp(self.henry[0]/T.R+self.henry[1]*log(T.R)+self.henry[2]*T.R+self.henry[3])
        return p.psi/H

class Solid(config.Entity):
    """Clase que define un sólido"""
    kwargs={"caudalSolido": [],
                    "diametroMedio": 0.0,
                    "distribucion_fraccion": [],
                    "distribucion_diametro": []
                    }

    def __init__(self, tipo=2, **kwargs):
        """
        tipo: definición de la corriente
            0   -   No definida
            1   -   diametro medio global
            2   -   distribución de tamaños
        kwargs:
            caudalSolido
            diametroMedio
            fraccionMasica
            diametros
        """
        self.Config=config.getMainWindowConfig()
        txt=self.Config.get("Components", "Solids")
        if isinstance(txt, str):
            self.ids=eval(txt)
        else:
            self.ids=txt
        self.componente=[Componente(int(i)) for i in self.ids]

        self.kwargs=kwargs
        caudal=self.kwargs.get("caudalSolido", [])
        diametro_medio=self.kwargs.get("diametroMedio", 0.0)
        fraccion=self.kwargs.get("distribucion_fraccion", [])
        diametros=self.kwargs.get("distribucion_diametro", [])

        if tipo==0:
            self._bool=False
            return
        else:
            self._bool=True

        self.caudalUnitario=[unidades.MassFlow(i) for i in caudal]
        self.caudal=unidades.MassFlow(sum(self.caudalUnitario))
        self.diametros=diametros
        self.fracciones=fraccion
        if tipo==2:
            self.diametros=[unidades.Length(i, magnitud="ParticleDiameter") for i in diametros]
            self.fracciones=fraccion
            diametro_medio=0
            self.fracciones_acumuladas=[0]
            for diametro, fraccion in zip(diametros, fraccion):
                diametro_medio+=diametro*fraccion
                self.fracciones_acumuladas.append(fraccion+self.fracciones_acumuladas[-1])
                del self.fracciones_acumuladas[0]
        self.diametro_medio=unidades.Length(diametro_medio, magnitud="ParticleDiameter")


    def RhoS(self, T):
        densidad=0
        for i in range(len(self.ids)):
            densidad+=self.caudalUnitario[i]/self.caudal*self.componente[i].RhoS(T)
        self.rho=unidades.Density(densidad)

    def __repr__(self):
        if self:
            return "%s con %s y dm %s" %(self.__class__, self.caudal.str, self.diametro_medio.str)
        else:
            return "%s empty" %(self.__class__)

    def ajustar_distribucion(self, eq=0):
        """
        Este método ajusta la distribución de tamaños disponible a alguna de las distribuciones estadisticas disponibles.
        eq: indice que indica el modelo al que ajustar la distribución de tamaños
            0   -   Rosin, Rammler, Sperling (distribución de Weibull)
            1   -   Gates, Gaudin, Shumann
            2   -   Gaudin, Meloy
            3   -   Broadbent, Callcott
            4   -   Distribución lognormal
            5   -   Harris

        Ref: Ahmed, Drzymala; Two-dimensional fractal linearization of distribution curves, pag 2
        """
        if eq==0:
            d=r_[self.diametros]
            y=r_[self.fracciones_acumuladas]
            funcion = lambda p, d: 1.-exp(-(d/p[0])**p[1]) # Función a ajustar
            residuo = lambda p, d, y: funcion(p, d) - y # Residuo
            inicio=r_[1, 1]
            ajuste, exito=leastsq(residuo,inicio,args=(d, y))
            return d,funcion(ajuste, d), "Rosin, Rammler, Sperling"

        elif eq==1:
            d=r_[self.diametros]
            y=r_[self.fracciones_acumuladas]
            funcion = lambda p, d: (d/p[0])**p[1] # Función a ajustar
            residuo = lambda p, d, y: funcion(p, d) - y # Residuo
            inicio=r_[1., 1.]
            ajuste, exito=leastsq(residuo,inicio,args=(d, y))
            return d,funcion(ajuste, d), "Gates, Gaudin, Shumann"

        elif eq==2:
            d=r_[self.diametros]
            y=r_[self.fracciones_acumuladas]
            funcion = lambda p, d: 1-(1-d/p[0])**p[1] # Función a ajustar
            residuo = lambda p, d, y: funcion(p, d) - y # Residuo
            inicio=r_[d[-1]*2, 0]
            ajuste, exito=leastsq(residuo,inicio,args=(d, y))
            return d,funcion(ajuste, d), "Gaudin, Meloy"

        elif eq==3:
            d=r_[self.diametros]
            y=r_[self.fracciones_acumuladas]
            funcion = lambda p, d: 1.-exp(-(d/p[0])**p[1])/(1-exp(-1.)) # Función a ajustar
            residuo = lambda p, d, y: funcion(p, d) - y # Residuo
            inicio=r_[1., 1.]
            ajuste, exito=leastsq(residuo,inicio,args=(d, y))
            return d,funcion(ajuste, d), "Broadbent, Callcott"

        elif eq==4:
            d=r_[self.diametros]
            y=r_[self.fracciones_acumuladas]
            funcion = lambda p, d: erf(log(d/p[0])/p[1]) # Función a ajustar
            residuo = lambda p, d, y: funcion(p, d) - y # Residuo
            inicio=r_[1., 1.]
            ajuste, exito=leastsq(residuo,inicio,args=(d, y))
            return d,funcion(ajuste, d), "lognormal"

        elif eq==5:
            d=r_[self.diametros]
            y=r_[self.fracciones_acumuladas]
            funcion = lambda p, d: 1-(1-d/p[0]**p[1])**p[2] # Función a ajustar
            residuo = lambda p, d, y: funcion(p, d) - y # Residuo
            inicio=r_[d[-1]*2, 1, 1]
            ajuste, exito=leastsq(residuo,inicio,args=(d, y))
            return d,funcion(ajuste, d), "Harris"


    def Separar(self, rendimiento_parcial):
        """Metodo que separa la corriente de sólido en dos a partir de un array con los rendimientos parciales
        Devuelve una tupla con los dos sólidos resultantes, el no recogido y el recogido"""
        rendimiento_global=0
        for i, fraccion in enumerate(self.fracciones):
            rendimiento_global+=rendimiento_parcial[i]*fraccion

        caudal_escapado=unidades.MassFlow(self.caudal*(1-rendimiento_global))
        caudal_separado=unidades.MassFlow(self.caudal*rendimiento_global)
        if rendimiento_global==1:
            return None, self
        elif rendimiento_global==0:
            return self, None
        else:
            fraccion_engas=[]
            fraccion_ensolido=[]
            for i in range(len(self.diametros)):
                fraccion_engas.append(self.caudal*self.fracciones[i]*(1-rendimiento_parcial[i])/caudal_escapado)
                fraccion_ensolido.append(self.caudal*self.fracciones[i]*rendimiento_parcial[i]/caudal_separado)
            Solido_NoCapturado=Solid(2, caudalSolido=[caudal_escapado], distribucion_diametro=self.diametros, distribucion_fraccion=fraccion_engas)
            Solido_Capturado=Solid(2, caudalSolido=[caudal_separado], distribucion_diametro=self.diametros, distribucion_fraccion=fraccion_ensolido)
            return Solido_NoCapturado, Solido_Capturado


class Corriente(config.Entity):
    """Clase que define una corriente"""
    kwargs={"T": 0.0,
                    "P": 0.0,
                    "x": None,

                    "caudalMasico": 0.0,
                    "caudalVolumetrico": 0.0,
                    "caudalMolar": 0.0,
                    "caudalUnitarioMolar": [],
                    "caudalUnitarioMasico": [],
                    "fraccionMolar": [],
                    "fraccionMasica": [],
                    "mezcla": None,

                    "solido": None,
                    "caudalSolido": [],
                    "diametroMedio": 0.0,
                    "distribucion_fraccion": [],
                    "distribucion_diametro": []}
    status=0
    msg="Unknown variables"
    kwargs_forbidden=["entrada", "mezcla", "solido"]
    solido=None

    def __init__(self, **kwargs):
        """Los parámetros necesarios para definirla son:
        -T: temperatura de la corriente en kelvin
        -P: presión de la corriente en atm
        -x: fracción de vapor

        -caudalMasico: caudal en k]/h (en el caso de corrientes con componentes sólidos, no se incluirán en el caudal total
        -caudalMolar: caudal en kmol/h (en el caso de corrientes con componentes sólidos, no se incluirán en el caudal total
        -caudalVolumetrico: caudal en m3/h medido en condiciones estandar de P y T
        -caudalUnitarioMolar: array con los caudales molares de cada componente
        -caudalUnitarioMasico: array con los caudales másicos de cada componente
        -fraccionMolar: array con las fracciones molares de los componentes
        -fraccionMasica: array con las fracciones másica de los componentes
        -mezcla: instancia de clase mezcla que define la mezcla

        -solido: instancia de clase solido que define el sólido
        -caudalSolido: caudal, en kg/h de sólidos, array con los valores para cada uno de los componentes sólidos
        -diametroMedio: diametro medio de las particulas, en micras
        -distribucion_fraccion: En el caso de indicarse una distribución de tamaños, fracción másica de los diferentes tamañanos
        -distribucion_diametro: En el caso de indicarse una distribución de tamaños, diámetros de las diferentes fracciones, en micras

        -notas: Texto explicativo sobre la corriente
        """
        self.Config=config.getMainWindowConfig()
        self.kwargs=Corriente.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        if kwargs.get("mezcla", None):
            kwargs.update(kwargs["mezcla"].kwargs)
        if kwargs.get("solido", None):
            kwargs.update(kwargs["solido"].kwargs)

        if kwargs.get("caudalUnitarioMasico", []):
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMolar"]=[]
            self.kwargs["fraccionMasica"]=[]
            self.kwargs["caudalMasico"]=None
            self.kwargs["caudalVolumetrico"]=None
            self.kwargs["caudalMolar"]=None

        elif kwargs.get("caudalUnitarioMolar", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["fraccionMolar"]=[]
            self.kwargs["fraccionMasica"]=[]
            self.kwargs["caudalMasico"]=None
            self.kwargs["caudalVolumetrico"]=None
            self.kwargs["caudalMolar"]=None

        elif kwargs.get("caudalMasico", None) and kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMasica"]=[]
            self.kwargs["caudalVolumetrico"]=None
            self.kwargs["caudalMolar"]=None

        elif kwargs.get("caudalMasico", None) and kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMolar"]=[]
            self.kwargs["caudalVolumetrico"]=None
            self.kwargs["caudalMolar"]=None

        elif kwargs.get("caudalMasico", None):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["caudalVolumetrico"]=None
            self.kwargs["caudalMolar"]=None

        elif kwargs.get("caudalMolar", None) and kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMasica"]=[]
            self.kwargs["caudalMasico"]=None
            self.kwargs["caudalVolumetrico"]=None

        elif kwargs.get("caudalMolar", None) and kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMolar"]=[]
            self.kwargs["caudalMasico"]=None
            self.kwargs["caudalVolumetrico"]=None

        elif kwargs.get("caudalMolar", None):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["caudalMasico"]=None
            self.kwargs["caudalVolumetrico"]=None

        elif kwargs.get("caudalVolumetrico", None) and kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMasica"]=[]
            self.kwargs["caudalMasico"]=None

        elif kwargs.get("caudalVolumetrico") and kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMolar"]=[]
            self.kwargs["caudalMasico"]=None

        elif kwargs.get("caudalVolumetrico", None):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["caudalMolar"]=[]
            self.kwargs["caudalMasico"]=None

        elif kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMolar"]=[]

        elif kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"]=[]
            self.kwargs["caudalUnitarioMolar"]=[]
            self.kwargs["fraccionMasica"]=[]

        self.kwargs.update(kwargs)

        self.kwargs.update()
        for key, value in self.kwargs.iteritems():
            if value:
#                self.__setattr__(key, value)
                self._bool=True

        if self.calculable:
            self.status=1
            self.calculo()
            self.msg=""
        elif self.tipoSolido:
            if self.kwargs["solido"]:
                self.solido=self.kwargs["solido"]
            else:
                self.solido=Solid(self.tipoSolido, **self.kwargs)
            if self.solido:
                self.solido.RhoS(self.kwargs["T"])

    @property
    def calculable(self):
        #Definición termodinámica
        self.tipoTermodinamica=""
        if self.kwargs["T"] and self.kwargs["P"]:
            self.tipoTermodinamica="TP"
        elif self.kwargs["T"] and self.kwargs["x"]:
            self.tipoTermodinamica="Tx"
        elif self.kwargs["P"] and self.kwargs["x"]:
            self.tipoTermodinamica="Px"
        #Definición de flujo
        self.tipoFlujo=0
        if self.kwargs["caudalUnitarioMasico"]:
            self.tipoFlujo=1

        elif self.kwargs["caudalUnitarioMolar"]:
            self.tipoFlujo=2

        elif self.kwargs["caudalMasico"] and self.kwargs["fraccionMolar"]:
            self.tipoFlujo=3

        elif self.kwargs["caudalMasico"] and self.kwargs["fraccionMasica"]:
            self.tipoFlujo=4

        elif self.kwargs["caudalMolar"] and self.kwargs["fraccionMolar"]:
            self.tipoFlujo=5

        elif self.kwargs["caudalMolar"] and self.kwargs["fraccionMasica"]:
            self.tipoFlujo=6

        elif self.kwargs["caudalVolumetrico"] and self.kwargs["fraccionMolar"]:
            self.kwargs["caudalMolar"]=1
            self.tipoFlujo=5

        elif self.kwargs["caudalVolumetrico"] and self.kwargs["fraccionMasica"]:
            self.kwargs["caudalMolar"]=1
            self.tipoFlujo=6

        elif self.kwargs["mezcla"]:
            self.tipoFlujo=7

        #Definición de sólido
        self.tipoSolido=0
        if sum(self.kwargs["caudalSolido"])>0:
            if self.kwargs["distribucion_fraccion"] and self.kwargs["distribucion_diametro"]:
                self.tipoSolido=2
            elif self.kwargs["diametroMedio"]:
                self.tipoSolido=1
        return self.tipoTermodinamica and self.tipoFlujo


    def calculo(self):
        if self.kwargs["mezcla"]:
            self.mezcla=self.kwargs["mezcla"]
        else:
            self.mezcla=Mezcla(self.tipoFlujo, **self.kwargs)

        self.ids=self.mezcla.ids
        self.componente=self.mezcla.componente_normalizado
        self.fraccion=self.mezcla.fraccion_normalizada
        self.caudalmasico=self.mezcla.caudalmasico
        self.caudalmolar=self.mezcla.caudalmolar
        self.fraccion_masica=self.mezcla.fraccion_masica
        self.caudalunitariomasico=self.mezcla.caudalunitariomasico
        self.caudalunitariomolar=self.mezcla.caudalunitariomolar

        T=unidades.Temperature(self.kwargs.get("T", None))
        P=unidades.Pressure(self.kwargs.get("P", None))
        x=self.kwargs.get("x", None)

        MEoS=self.Config.getboolean("Thermo","MEoS")
        COOLPROP=self.Config.getboolean("Thermo", "coolprop")
        REFPROP=self.Config.getboolean("Thermo", "refprop")
        IAPWS=self.Config.getboolean("Thermo", "iapws") and len(self.ids)==1 and self.ids[0]==62
        FREESTEAM=self.Config.getboolean("Thermo", "freesteam") and len(self.ids)==1 and self.ids[0]==62
        GERG=self.Config.getboolean("Thermo","GERG")
        setData=True

        mEoS_available=self.ids[0] in mEoS.id_mEoS
        GERG_available=True
        REFPROP_available=True
        COOLPROP_available=self.ids[0] in coolProp.__all__
        for id in self.ids:
            if id not in gerg.id_GERG:
                GERG_available=False
            if id not in refProp.__all__:
                REFPROP_available=False

        if IAPWS and FREESTEAM:
            compuesto=freeSteam.Freesteam(**self.kwargs)
        elif IAPWS:
            compuesto=iapws.IAPWS97(**self.kwargs)
        elif MEoS and REFPROP and REFPROP_available:
            fluido=[refProp.__all__[id] for id in self.ids]
            compuesto=refProp.RefProp(fluido=fluido, **self.kwargs)
        elif GERG and  GERG_available:
            ids=[]
            for id in self.ids:
                ids.append(gerg.id_GERG.index(id))
            kwargs=self.kwargs
            kwargs["mezcla"]=self.mezcla
            compuesto=gerg.GERG(componente=ids, fraccion=self.fraccion, **kwargs)
        elif MEoS and len(self.ids)==1 and COOLPROP and COOLPROP_available:
            compuesto=coolProp.CoolProp(fluido=self.ids[0], **self.kwargs)
        elif MEoS and len(self.ids)==1 and mEoS_available:
            if self.tipoTermodinamica=="TP":
                compuesto=mEoS.__all__[mEoS.id_mEoS.index(self.ids[0])](T=T, P=P.MPa)
            elif self.tipoTermodinamica=="Tx":
                compuesto=mEoS.__all__[mEoS.id_mEoS.index(self.ids[0])](T=T, x=x)
            elif self.tipoTermodinamica=="Px":
                compuesto=mEoS.__all__[mEoS.id_mEoS.index(self.ids[0])](P=P.MPa, x=x)

        else:
            setData=False
            self.M=unidades.Dimensionless(self.mezcla.M)
            self.Tc=self.mezcla.Tc
            self.Pc=self.mezcla.Pc
            self.SG=unidades.Dimensionless(self.mezcla.SG)

            if self.tipoTermodinamica=="TP":
                self.T=unidades.Temperature(T)
                self.P=unidades.Pressure(P)
                eos=EoS.K[self.Config.getint("Thermo","K")](self.T, self.P.atm, self.mezcla)
                self.eos=eos
                self.x=unidades.Dimensionless(eos.x)
            else:
                self.x=unidades.Dimensionless(x)

            self.mezcla.recallZeros(eos.xi)
            self.mezcla.recallZeros(eos.yi)
            self.mezcla.recallZeros(eos.Ki, 1.)

            if 0.<self.x<1.:
                self.Liquido=Mezcla(tipo=5, fraccionMolar=eos.xi, caudalMolar=self.caudalmolar*(1-self.x))
                self.Gas=Mezcla(tipo=5, fraccionMolar=eos.yi, caudalMolar=self.caudalmolar*self.x)
            elif self.x<=0:
                self.Liquido=self.mezcla
                self.Gas=Mezcla()
            else:
                self.Liquido=Mezcla()
                self.Gas=self.mezcla
            self.Gas.Z=unidades.Dimensionless(float(eos.Z[0]))
            self.Liquido.Z=unidades.Dimensionless(float(eos.Z[1]))

            if EoS.H[self.Config.getint("Thermo","H")].__title__ == EoS.K[self.Config.getint("Thermo","K")].__title__:
                eosH=eos
            else:
                eosH=EoS.H[self.Config.getint("Thermo","H")](self.T, self.P.atm, self.mezcla)
            self.H_exc=eosH.H_exc

            self.Liquido.Q=unidades.VolFlow(0)
            self.Gas.Q=unidades.VolFlow(0)
            self.Liquido.h=unidades.Power(0)
            self.Gas.h=unidades.Power(0)
            if self.x<1:      #hay corriente líquida
                Hl=(self.Liquido.Entalpia_ideal(self.T).Jg-self.Liquido.Hv_DIPPR(self.T).Jg)*self.Liquido.caudalmasico.gh
                self.Liquido.h=unidades.Power(Hl-R*self.T/self.M*self.H_exc[1]*(1-self.x)*self.Liquido.caudalmasico.gh, "Jh")
                self.Liquido.cp=self.Liquido.Cp_Liquido(T)
                self.Liquido.rho=self.Liquido.RhoL_Tait_Costald(T, self.P.atm)
                self.Liquido.mu=self.Liquido.Mu_Liquido(T, self.P.atm)
                self.Liquido.k=self.Liquido.ThCond_Liquido(T, self.P.atm)
                self.Liquido.epsilon=self.Liquido.Tension(T)
                self.Liquido.Q=unidades.VolFlow(self.Liquido.caudalmasico/self.Liquido.rho)
                self.Liquido.Pr=self.Liquido.cp*self.Liquido.mu/self.Liquido.k
            if self.x>0:            #hay corriente gaseosa
                Hg=self.Gas.Entalpia_ideal(self.T).Jg*self.Gas.caudalmasico.gh
                self.Gas.h=unidades.Power(Hg-R*self.T/self.M*self.H_exc[0]*self.x*self.Gas.caudalmasico.gh, "Jh")
                self.Gas.cp=self.Gas.Cp_Gas(T, self.P.atm)
                self.Gas.rho=unidades.Density(self.P.atm/self.Gas.Z/R_atml/self.T*self.M, "gl")
                self.Gas.rhoSd=unidades.Density(1./self.Gas.Z/R_atml/298.15*self.M, "gl")
                self.Gas.mu=self.Gas.Mu_Gas(T, self.P.atm)
                self.Gas.k=self.Gas.ThCond_Gas(T, self.P.atm)
                self.Gas.Q=unidades.VolFlow(self.Gas.caudalmasico/self.Gas.rho)
                self.Gas.Pr=self.Gas.cp*self.Gas.mu/self.Gas.k

            self.Q=unidades.VolFlow(self.Liquido.Q+self.Gas.Q)
            self.h=unidades.Power(self.Liquido.h+self.Gas.h)
            self.Molaridad=[caudal/self.Q.m3h for caudal in self.caudalunitariomolar]

            #TODO:
            self.cp_cv=0.5
            self.cp_cv_ideal=0.5


        if setData:
        #Asignación de valores comun
            self.T=compuesto.T
            self.P=compuesto.P
            self.x=compuesto.x
            self.M=unidades.Dimensionless(compuesto.M)
            self.Tc=compuesto.Tc
            self.Pc=compuesto.Pc
            self.h=unidades.Power(compuesto.h*self.caudalmasico)
            self.s=unidades.Entropy(compuesto.s*self.caudalmasico)
            self.rho=compuesto.rho
            self.Q=unidades.VolFlow(compuesto.v*self.caudalmasico)

#        if self.__class__!=mEoS.H2O:
#            agua=mEoS.H2O(T=self.T, P=self.P.MPa)
#            self.SG=unidades.Dimensionless(self.rho/agua.rho)
#        else:
            self.SG=unidades.Dimensionless(1.)

            self.Liquido=compuesto.Liquido
            self.Gas=compuesto.Gas

#            if self.x<1:      #Fase líquida
#                self.Liquido=compuesto.Liquido
#            if self.x>0:       #Fase gaseosa
#                self.Gas=compuesto
#            elif 0<self.x<1:    #Ambas fases
#                self.Liquido=compuesto.Liquido
#                self.Gas=compuesto.Gas
#                self.QLiquido=self.Q*(1-self.x)
#            elif self.x==1:            #Fase gaseosa
#                self.Gas=compuesto
#                self.Liquido=None
#                self.QLiquido=unidades.VolFlow(0)
#                self.Gas.rhoSd=unidades.Density(1./R_atml/298.15*self.M, "gl")

            if self.x>0:
                self.Gas.Q=unidades.VolFlow(self.Q*(1-self.x))
                self.Gas.caudalmasico=unidades.MassFlow(self.caudalmasico*self.x)
                self.Gas.caudalmolar=unidades.MolarFlow(self.caudalmolar*self.x)
                self.Gas.fraccion=[unidades.Dimensionless(1)]
                self.Gas.fraccion_masica=[unidades.Dimensionless(1)]
                self.Gas.caudalunitariomasico=[self.Gas.caudalmasico]
                self.Gas.caudalunitariomolar=[self.Gas.caudalmolar]

            if self.x<1:
                self.Liquido.Q=unidades.VolFlow(self.Q*(1-self.x))
                self.Liquido.caudalmasico=unidades.MassFlow(self.caudalmasico*(1-self.x))
                self.Liquido.caudalmolar=unidades.MolarFlow(self.caudalmolar*(1-self.x))
                self.Liquido.fraccion=[unidades.Dimensionless(1)]
                self.Liquido.fraccion_masica=[unidades.Dimensionless(1)]
                self.Liquido.caudalunitariomasico=[self.Liquido.caudalmasico]
                self.Liquido.caudalunitariomolar=[self.Liquido.caudalmolar]


        if self.Config.get("Components", "Solids"):
            if self.kwargs["solido"]:
                self.solido=self.kwargs["solido"]
            else:
                self.solido=Solid(self.tipoSolido, **self.kwargs)
            if self.solido:
                self.solido.RhoS(T)
        else:
            self.solido=None


        if self.kwargs["caudalVolumetrico"]:
            self.kwargs["caudalMolar"]*=self.kwargs["caudalVolumetrico"]/self.Q
            Q=self.kwargs["caudalVolumetrico"]
            self.kwargs["caudalVolumetrico"]=None
            self.calculo()
            self.kwargs["caudalVolumetrico"]=Q
            self.kwargs["caudalMolar"]=None


    def clone(self, **kwargs):
        old_kwargs=self.kwargs.copy()
        if kwargs.has_key("split"):
            split=kwargs["split"]
            del kwargs["split"]
            if self.kwargs["caudalUnitarioMasico"]:
                kwargs["caudalUnitarioMasico"]=[]
                for caudal in self.kwargs["caudalUnitarioMasico"]:
                    kwargs["caudalUnitarioMasico"].append(split*caudal)
            if self.kwargs["caudalUnitarioMolar"]:
                kwargs["caudalUnitarioMolar"]=[]
                for caudal in self.kwargs["caudalUnitarioMolar"]:
                    kwargs["caudalUnitarioMolar"].append(split*caudal)
            if self.kwargs["caudalMasico"]:
                kwargs["caudalMasico"]=split*self.kwargs["caudalMasico"]
            if self.kwargs["caudalMolar"]:
                kwargs["caudalMolar"]=split*self.kwargs["caudalMolar"]
        if kwargs.has_key("x"):
            del old_kwargs["T"]
        if kwargs.has_key("mezcla"):
            old_kwargs.update(kwargs["mezcla"].kwargs)
        old_kwargs.update(kwargs)
        return Corriente(**old_kwargs)


    def __repr__(self):
        if self:
            return "%s a %0.2f K y %0.2f atm" %(self.__class__, self.T, self.P.atm)
        else:
            return "%s empty" %(self.__class__)

    def txt(self):
        txt=str(self.notasPlain)+os.linesep+os.linesep
        txt+="#---------------"+QApplication.translate("pychemqt", "Input properties")+"-----------------#"+os.linesep
        for key, value in self.kwargs.iteritems():
            if value:
                txt+=key+": "+str(value)+os.linesep

        if self.calculable:
            txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Global stream")+"-------------------#"+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Temperature"), self.T.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure"), self.P.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Vapor Fraction"), self.x.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Molar Flow"), self.caudalmasico.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Mass Flow"), self.caudalmolar.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Volumetric Flow"), self.Q.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Enthalpy"), self.h.str)+os.linesep
            txt+="%-25s\t%s" %("Tc", self.Tc.str)+os.linesep
            txt+="%-25s\t%s" %("Pc", self.Pc.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "SG, water=1"), self.SG.str)+os.linesep
            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Molecular weight"), self.M.str)+os.linesep
            txt+="#"+QApplication.translate("pychemqt", "Molar Composition")+os.linesep
            for componente, fraccion in zip(self.componente, self.fraccion):
                txt+="%-25s\t %0.4f" %(componente.nombre, fraccion)+os.linesep

            if self.x>0:
                txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Vapor Only")+"--------------------#"+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Molar Flow"), self.Gas.caudalmasico.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Mass Flow"), self.Gas.caudalmolar.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Volumetric Flow"), self.Gas.Q.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Molecular weight"), self.Gas.M.str)+os.linesep
                txt+=os.linesep+"#"+QApplication.translate("pychemqt", "Molar Composition")+os.linesep
                for componente, fraccion in zip(self.componente, self.Gas.fraccion):
                    txt+="%-25s\t %0.4f" %(componente.nombre, fraccion)+os.linesep

                txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Density"), self.Gas.rho.str)+os.linesep
                txt+="%-25s\t%s" %("Compresibility", self.Gas.Z.str)+os.linesep

                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Enthalpy"), self.Gas.h.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Heat Capacity"), self.Gas.cp.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Viscosity"), self.Gas.mu.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thermal conductivity"), self.Gas.k.str)+os.linesep

            if self.x<1:
                txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Liquid Only")+"-------------------#"+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Molar Flow"), self.Liquido.caudalmasico.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Mass Flow"), self.Liquido.caudalmolar.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Volumetric Flow"), self.Liquido.Q.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Molecular weight"), self.Liquido.M.str)+os.linesep
                txt+=os.linesep+"#"+QApplication.translate("pychemqt", "Molar Composition")+os.linesep
                for componente, fraccion in zip(self.componente, self.Liquido.fraccion):
                    txt+="%-25s\t %0.4f" %(componente.nombre, fraccion)+os.linesep

                txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Density"), self.Liquido.rho.str)+os.linesep
                txt+="%-25s\t%s" %("Compresibility", self.Liquido.Z.str)+os.linesep

                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Enthalpy"), self.Liquido.h.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Heat Capacity"), self.Liquido.cp.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Viscosity"), self.Liquido.mu.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thermal Conductivity"), self.Liquido.k.str)+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Surface Tension"), self.Liquido.epsilon.str)+os.linesep

        else:
            txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "No Fluid Stream")+"-------------------#"+os.linesep

        if self.solido:
            txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Solid")+"-------------------#"+os.linesep
            for componente, caudal in zip(self.solido.componente, self.solido.caudalUnitario):
                txt+="%-25s\t%s" %(componente.nombre, caudal.str)+os.linesep
            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Density"), self.solido.rho.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Mean Diameter"), self.solido.diametro_medio.str)+os.linesep
            if self.solido.diametros:
                txt+=os.linesep+"#"+QApplication.translate("pychemqt", "Particle Size Distribution")+os.linesep
                txt+="%s, %s \t%s" %(QApplication.translate("pychemqt", "Diameter"), unidades.Length.text("ParticleDiameter"), QApplication.translate("pychemqt", "Fraction"))+os.linesep
                for diametro, fraccion in zip(self.solido.diametros, self.solido.fracciones):
                    txt+="%10.4f\t%0.4f\t" %(diametro.config("ParticleDiameter"), fraccion)+os.linesep
        return txt


    @classmethod
    def propertiesNames(cls):
        list=[(QApplication.translate("pychemqt", "Temperature"), "T", unidades.Temperature),
                (QApplication.translate("pychemqt", "Pressure"), "P", unidades.Pressure),
                (QApplication.translate("pychemqt", "Vapor Fraction"), "x", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Molar Flow"), "caudalmolar", unidades.MolarFlow),
                (QApplication.translate("pychemqt", "Mass Flow"), "caudalmasico", unidades.MassFlow),
                (QApplication.translate("pychemqt", "Volumetric Flow"), "Q", unidades.VolFlow),
                (QApplication.translate("pychemqt", "Enthalpy"), "h", unidades.Enthalpy),
                (QApplication.translate("pychemqt", "Critic Temperature"), "Tc", unidades.Temperature),
                (QApplication.translate("pychemqt", "Critic Pressure"), "Pc", unidades.Pressure),
                (QApplication.translate("pychemqt", "SG, water=1"), "SG", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Molecular weight"), "M", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Molar Composition"), "fraccion", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Mass Composition"), "fraccion_masica", unidades.Dimensionless),
                (QApplication.translate("pychemqt", "Molar Component Flow"), "caudalunitariomolar", unidades.MolarFlow),
                (QApplication.translate("pychemqt", "Mass Component Flow"),  "caudalunitariomasico", unidades.MassFlow),
                (QApplication.translate("pychemqt", "Notes"), "notasPlain", str)]
        return list

    def propertiesListTitle(self, index):
        """Define los titulos para los popup de listas"""
        lista=[comp.nombre for comp in self.componente]
        return lista


if __name__ == '__main__':

#    aire=Mezcla([46, 47, 49], [0.781, 0.209, 0.01])
#    print aire.RhoL_Rackett(20)
#    print aire.RhoL_Costald(20)[0]
#    print 1/aire.RhoG_Lee_Kesler(300, 1),  "l/mol"
#    print aire.RhoL_API(20, 1)
#    print aire.Viscosidad_liquido(20)
#    suma=0
#    for i in aire.fraccion_masica:
#        print i*0.9
#        suma+=i
#    print suma


#    ejemplo1=Mezcla([3, 11], [0.5871, 0.4129])
#    print Densidad(ejemplo1.densidad_Rackett(Temperature(91, "F"))).lbft3
#
#    ejemplo2=Mezcla([2, 14], [0.2, 0.8])
#    print Densidad(ejemplo2.densidad_Costald(Temperature(160, "F"))[0]).lbft3
#
#    ejemplo3=Mezcla([3, 14], [0.2, 0.8])
#    print Densidad(ejemplo3.densidad_Tait_Costald(Temperature(160, "F"), Pressure(3000, "psi").atm), "gl").lbft3

#    ejemplo5=Mezcla([1, 25, 49], [0.3, 0.415, 0.285], 300)
#    print Temperature(ejemplo5.Tc()).R
#    print ejemplo5.propiedades_criticas()
#    print ejemplo5.kij()


#    ejemplo6=Mezcla([45, 12], [0.3, 0.7])
#    print Pressure(ejemplo6.Pc(), "atm").psi

#    ejemplo7=Mezcla([6, 11], [0.63, 0.37])
#    print Volumen(ejemplo7.Vc()).ft3*lb

#    ejemplo8=Mezcla([2, 3], [0.1, 0.9])
#    print Temperature(ejemplo8.Tc()).R
#    print Pressure(ejemplo8.Pc(), "atm").psi
#    print Volumen(ejemplo8.Vc()).ft3*lb
#    print ejemplo8.kij()

#    ejemplo9=Mezcla([2, 3, 4, 6, 46], [0.868, 0.065, 0.025, 0.007, 0.034])
#    print ejemplo9.Tb().F
#    print ejemplo9.tc_gas().F

#    ejemplo10=Mezcla([133, 7], [7.56, 92.44])
#    print ejemplo10.Reid()

#    ejemplo11=Mezcla([3, 14], [0.2, 0.8])
#    t=Temperature(160, "F")
#    p=Pressure(3000, "psi")
#    print ejemplo11.RhoL_API(t, p.atm)

#    ejemplo12=Mezcla([3, 14], [0.2, 0.8])
#    print ejemplo12.flash_water(300, 1)

#    amoniaco=Mezcla([62, 63], [0.9894, 0.0106])
#    t=unidades.Temperature(248, "F")
#    print amoniaco.fraccion_masica[1]*1e6
#    print amoniaco.SOUR_water(t).psi

#    sour=Mezcla([62, 63, 50], [0.9796, 0.0173, 0.0031])
#    print sour.SOUR_water(t)[1].mmHg

#    sour=Mezcla([62, 63, 50], [0.97, 0.02, 0.01])
#    t=unidades.Temperature(80, "F")
#    print sour.SOUR_water_ph(t)


#    ejemplo13=Mezcla([40, 38], [0.379, 0.621])
#    print ejemplo13.Tension_superficial(unidades.Temperature(77, "F")).dyncm

#    ejemplo14=Mezcla([2, 4], [0, 1])
#    print ejemplo14.Tension_superficial_presion([72.6, 150.8], [0.418, 0.582], [0.788, 0.212])

#    ejemplo15=Mezcla([14], [1])
#    print ejemplo15.Tension_inferfacial_water(100+273.15).dyncm

#    ejemplo16=Mezcla([20, 40, 10], [0.2957, 0.3586, 0.3457])
#    print ejemplo16.Viscosidad_liquido(unidades.Temperature(77, "F")).cP

#    ejemplo17=Mezcla([4, 8, 38], [0.25, 0.5, 0.25])
#    print ejemplo17.Viscosidad_liquido(unidades.Temperature(160, "F")).cP
#    print ejemplo17.Viscosidad_gas(unidades.Temperature(160, "F"))

#    ejemplo18=Mezcla([1, 4], [0.5818, 0.4182])
#    t=unidades.Temperature(77, "F")
#    print ejemplo18.Mu_Gas_Wilke(t).cP

#    ejemplo=Mezcla([2, 3, 4, 46], [0.956, 0.036, 0.005, 0.003])
#    t=unidades.Temperature(85, "F")
#    print ejemplo.Mu_Gas_Wilke(t).cP

#    ejemplo19=Mezcla([2, 3, 4, 46], [95.6, 3.6, 0.5, 0.3])
#    print ejemplo19.Viscosidad_gas(unidades.Temperature(85, "F")).cP
#    print ejemplo19.Viscosidad_gas_presion(unidades.Temperature(85, "F"), 5).cP

#    ejemplo20=Mezcla([2, 4], [0.6, 0.4])
#    print ejemplo20.Viscosidad_gas(398.15).cP
#    print ejemplo20.Viscosidad_gas_presion(398.15, 102.07).cP
#    print ejemplo20.MuG_Carr(398.15, 102.07).cP

#    ejemplo21=Mezcla([11, 36], [0.68, 0.32])
#    t=unidades.Temperature(32, "F")
#    print ejemplo21.Conductividad_Termica_liquido(t, 1).BtuhftF

#    ejemplo22=Mezcla([8, 10], [0.2996, 0.7004])
#    t=unidades.Temperature(212, "F")
#    print ejemplo22.Conductividad_Termica_gas(t, 1).BtuhftF

#    ejemplo23=Mezcla([5, 10], [0.45, 0.55])
#    t=unidades.Temperature(200, "C")
#    print ejemplo23.RhoG_SRK(t, 1)
#    print ejemplo23.Z_SRK(t, 1)
#    print ejemplo23.Entalpia_SRK(t, 1)
#    print ejemplo23.RhoG_MSRK(t, 1)
#    print ejemplo23.Z_MSRK(t, 1)
#    print ejemplo23.Z_BWRS(t, 1)

#    ejemplo24=Mezcla([49, 4], [0.5, 0.5])
#    p=unidades.Pressure(140, "bar")
#    print ejemplo24.Z_RK(450, p.atm)
#    print ejemplo24.Z_SRK(450, p.atm)
#    print ejemplo24.Z_PR(450, p.atm)


#    ejemplo25=Mezcla([2, 3, 4, 50, 49, 46], [0.1, 0.3, 0.25, 0.05, 0.23, 0.07])
#    print ejemplo25.Z_PR(298.15, 1)
#    print ejemplo25.RhoG_PR(298.15, 1)
#    print ejemplo25.Fugacidad_PR(298.15, 1)
#    print ejemplo25.Mu_Gas_Wilke(298.15)
#    print ejemplo25.ThCond_Gas(298.15, 1)

#    ejemplo26=Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15])
#    t=unidades.Temperature(76, "C")
#    ejemplo26.Flash_SRK(320, 1)

#    ejemplo27=Mezcla( [46, 47], [0.781, 0.209])
#    t=unidades.Temperature(76, "C")
#    print ejemplo27.Flash_SRK(t, 1)

#    ejemplo28=Mezcla([8, 12], [39.2, 60.8]) #Ej pag 647
#    t=unidades.Temperature(590, "F")
#    p=unidades.Pressure(1400, "psi")
#    print ejemplo28.Lee_Kesler_Entalpia(t, p.atm)

#    agua=Corriente(298.15, 1, 3600, [[62], [1]])
#    print agua.RhoL



#    distribucion=[[17.5, 0.02],
#                                [22.4, 0.03],
#                                [26.2,  0.05],
#                                [31.8,  0.1],
#                                [37, 0.1],
#                                [42.4, 0.1],
#                                [48, 0.1],
#                                [54, 0.1],
#                                [60, 0.1],
#                                [69, 0.1],
#                                [81.3, 0.1],
#                                [96.5, 0.05],
#                                [109, 0.03],
#                                [127, 0.02]]
#    diametros=[]
#    fracciones=[]
#    for diametro, fraccion in distribucion:
#        diametros.append(diametro)
#        fracciones.append(fraccion)
#
#    solido=Solid(caudal=[5], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
#    solido.ajustar_distribucion(3)

#    t=unidades.Temperature(160, "F")
#    p=unidades.Pressure(3000, "psi")
#    mezcla=Mezcla([2, 14], [0.2, 0.8]) #Ej pag 482
#    print mezcla.RhoL_Costald(t).gml
#
#    mezcla=Mezcla([3, 14], [0.2, 0.8]) #Ej pag 490
#    print mezcla.RhoL_Tait_Costald(t, p.atm).gml

#    t=Temperature(27, "C")
#    P= Pressure(10, "atm")
#    agua=Corriente(t, P.atm, 3600, [[62], [1]])
#    print agua.RhoL

#    agua=Mezcla( [62], [1])
#    t=unidades.Temperature(76, "C")
#    print agua.Flash_SRK(t, 1)

#    ejemplo=Corriente(340, 1, 1, Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15]))
#    print  ejemplo.Entalpia.MJh

#    agua1=Corriente(298.15, 1, 1, Mezcla([62], [1.0]))
#    agua2=Corriente(350, 1, 1, Mezcla([62], [1.0]))
#    print agua2.Entalpia.MJkg-agua1.Entalpia.MJkg

#    ejemplo2=Mezcla([10, 38, 22, 61], [0.3045265648865772, 0.61740493041538613, 0.0010833278211912348, 0.076985176876845404])
#    print ejemplo2.RhoL_Tait_Costald(340, 1)

#    ejemplo=Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15])
#    print ejemplo.SRK_lib(300)

#    t=unidades.Temperature(100, "F")
#    p=unidades.Pressure(2000, "psi")
#    ejemplo=Corriente(t, p.atm, 1, Mezcla([1, 4], [.12, 0.88]))

#    ejemplo=Mezcla([2, 4], [0.6, 0.4])
#    t=unidades.Temperature(257, "F")
#    p=unidades.Pressure(1500, "psi")
#    print ejemplo.Mu_Gas_Stiel(t, p.atm)

#    ejemplo=Mezcla([11, 36], [0.68, 0.32])
#    t=unidades.Temperature(32, "F")
#    print ejemplo.ThCond_Liquido(t, 1).BtuhftF

#    ejemplo=Mezcla([10, 38, 22, 61], [.3, 0.5, 0.05, 0.15])
#    t=unidades.Temperature(32, "F")
#    print ejemplo.ThCond_Liquido(t, 1).BtuhftF

#    ejemplo=Mezcla([40, 38], [0.552, 0.448]) #Ej. pag 701
#    t=unidades.Temperature(72.1, "F")
#    print ejemplo.Cp_liquido(t).BtulbF
#
#
#    ejemplo=Mezcla([6, 11], [0.63, 0.37])
#    print ejemplo.Vc.ft3lb/ejemplo.M

#    ejemplo=Mezcla([2, 3], [0.1, 0.9])
#    print ejemplo.Tc.R, ejemplo.Pc.psi, ejemplo.Vc.ft3lb/ejemplo.M

#    agua1=Corriente(400, 1, 1, Mezcla([62], [1]))
#    agua2=Corriente(450, 1, 1, Mezcla([62], [1]))
#    print agua2.h.MJh-agua1.h.MJh
#    print agua1.Liquido.Z, agua1.Gas.Z

#    mezcla=Corriente(T=340, P=1, caudal=1, mezcla=Mezcla([10, 38, 22, 61], [0.3, 0.5, 0.05, 0.15]))
#    print mezcla.x, mezcla.H_exc

#    mezcla=Corriente()
#    if mezcla:
#        print True
#    else:
#        print False

#    mezcla=Corriente(T=300)
#    mezcla(P=101325)
#    mezcla(caudalMasico=2)
#    mezcla(fraccionMolar=[1, 0, 0, 0])
#    print mezcla.caudalmasico, mezcla.caudalmolar
#    mezcla(caudalMolar=3, fraccionMolar=[1, 0, 0, 0])
#    print mezcla.caudalmasico, mezcla.caudalmolar
#    print mezcla.x, mezcla.H_exc
#    print mezcla.T, mezcla.P.atm, mezcla.Q.ft3s
#    mezcla.clear()

#    agua=Corriente(T=300, P=1e5, caudalMasico=5, fraccionMolar=[1.])
#    agua2=agua.clone(P=101325, split=0.9)
#    print agua.P, agua.Liquido.caudalmasico

#    z=0.965
#    mez=Mezcla(tipo=3, fraccionMolar=[z, 1-z], caudalMasico=1.)
#    tb=mez.componente[0].Tb
#    print tb
#    corr=Corriente(T=tb, P=101325., mezcla=mez)
#    print corr.eos._Dew_T()

#    mezcla=Mezcla(caudalMasico=0.01, ids=[10, 38, 22, 61], fraccionMolar=[.3, 0.5, 0.05, 0.15], tipo=3)
#    corr=Corriente(T=300, P=101325., mezcla=mezcla)
#    print corr.x, corr.Q, corr.P, corr.caudalmasico

    corr=Corriente(T=300, P=101325.)
    if corr:
        print bool(corr)
