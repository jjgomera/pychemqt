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



from string import ascii_lowercase, digits
import tempfile
import time
import os

from scipy import exp, cosh, sinh, log, log10, roots, absolute, sqrt
from scipy.optimize import fsolve
from scipy.constants import R, Avogadro
from PyQt5.QtWidgets import QApplication

from lib.physics import R_atml, R_Btu, R_cal, factor_acentrico_octano
from lib import unidades, config, eos, sql


class Componente(object):
    """Clase que define los compuestos químicos con todas sus caracteristicas,
    algunas sacadas de la base de datos, otras calculadas a partir de estas
    Introduciendo el id del componente de la base de datos quedaría perfectamente definido"""

    def __init__(self, indice=None):
        if not indice:
            return
        self.indice=indice
        self.Config=config.getMainWindowConfig()
        componente=sql.getElement(indice)
        self.formula=componente[1]
        self.nombre=componente[2]
        self.M=componente[3]
        self.SG=componente[124]
        if componente[4]!=None:
            self.Tc=unidades.Temperature(componente[4])
        else:
            self.Tc=self._Tc_Nokay()
        self.Pc=unidades.Pressure(componente[5], "atm")
        self.Tb=unidades.Temperature(componente[131])
        self.Tf=unidades.Temperature(componente[132])
        if componente[125]!=0:
            self.f_acent=componente[125]
        elif self.Pc!=0 and self.Tc!=0:
            self.f_acent=self.factor_acentrico_Lee_Kesler()
        else:
            self.f_acent=0
        if componente[6]!=0:
            self.Vc=unidades.SpecificVolume(componente[6])
        elif self.f_acent!=0 and self.Tc!=0 and self.Pc!=0:
            self.Vc=self.vc_Riedel()
        else:
            self.Vc=0
        if self.Tc:
            self.Zc=self.Pc*self.Vc/R/self.Tc
        else:
            self.Zc=0
        if componente[7]!=0:
            self.API=componente[7]
        elif componente[124]!=0:
            self.API=141.5/componente[124]-131.5
        else:
            self.API=0
        self.cp=componente[8:14]
        self.antoine=componente[14:17]
        self.henry=componente[17:21]
        self.viscosidad_parametrica=componente[21:23]
        self.tension_superficial_parametrica=componente[23:25]
        self.densidad_solido=componente[25:33]
        self.densidad_liquido=componente[33:41]
        self.presion_vapor=componente[41:49]
        self.calor_vaporizacion=componente[49:57]
        self.capacidad_calorifica_solido=componente[57:65]
        self.capacidad_calorifica_liquido=componente[65:73]
        self.capacidad_calorifica_gas=componente[73:81]
        self.viscosidad_liquido=componente[81:89]
        self.viscosidad_gas=componente[89:97]
        self.conductividad_liquido=componente[97:105]
        self.conductividad_gas=componente[105:113]
        self.tension_superficial=componente[113:121]
        self.momento_dipolar=unidades.DipoleMoment(componente[121])
        if componente[123]!=0.0:
            self.rackett=componente[123]
        else:
            self.rackett=self.Rackett()
        if componente[122]!=0.0:
            self.Vliq=componente[122]
#        elif self.Pc!=0 and self.Tc>298.15 :
#            self.Vliq=self.Volumen_Liquido_Constante()
        else:
            self.Vliq=0

        if componente[126]!=0.0:
            self.parametro_solubilidad=unidades.SolubilityParameter(componente[126])
        else:
            self.parametro_solubilidad=self.Parametro_Solubilidad()
        self.Kw=componente[127]
        self.MSRK=componente[128:130]
        if componente[130]!=0.0:
            self.stiehl=componente[130]
        else:
            self.stiehl=0
        #FIXME: No esta bien
#            self.stiehl=self.Stiehl_Polar_factor()
        self.CASNumber=componente[133]
        self.formula_alternativa=componente[134]
        self.UNIFAC=eval(componente[135])
        self.diametro_molecular=componente[136]
        self.ek=componente[137]
        self.UNIQUAC_area=componente[138]
        self.UNIQUAC_volumen=componente[139]
        if componente[140]==0.0:
            self.f_acent_mod=componente[125]
        else:
            self.f_acent_mod=componente[140]
        self.calor_formacion=unidades.Enthalpy(componente[141]/self.M)
        self.energia_formacion=unidades.Enthalpy(componente[142]/self.M)
        self.wilson=componente[143]
        self.calor_combustion_neto=unidades.Enthalpy(componente[144]/self.M)
        self.calor_combustion_bruto=unidades.Enthalpy(componente[145]/self.M)
        self.nombre_alternativo=componente[146]
        self.V_char=componente[147]
        self.calor_formacion_solido=componente[148]
        self.energia_formacion_solido=componente[149]
        self.parametro_polar=componente[150]
        self.smile=componente[151]
        if self.smile!="" and os.environ["oasa"] == "True":
            import oasa
            # Install bkchem and create symbolic link for oasa
            # ln -s /usr/lib/bkchem/bkchem/oasa/oasa /usr/local/lib/python2.7/dist-packages
            molecula=oasa.smiles.text_to_mol(self.smile,calc_coords=40)
            self.archivo_imagen=tempfile.NamedTemporaryFile("w+r", suffix=".svg")
            oasa.svg_out.mol_to_svg(molecula, self.archivo_imagen.name)

        #TODO: Añadir las contribuciones de grupos de parachor a la base de datos
        self.parachor=[]
        #TODO: Añadir las parámetros de la ecuación de wagner de la presión de vapor, de momento usamos los datos del tetralin
        self.wagner=[-7.4138E+00, 9.5219E-01, -1.5111E+00, -4.9487E+00, 427, 1296]
        #TODO: Añadir los tipos de cada elemento
        self.hidrocarburo=True
        self.Van_Veltzen=[] #Quizá se pueda derivar de otro grupos, UNIFAC o similar
        #TODO: Añadir caraceristicas químicas del componente
        self.isCiclico=False
        self.isHidrocarburo=True
        self.isLineal=True
        self.isAlcohol=False

        #TODO: Añadir parametros S1,S2 a la base de datos, API databook, pag 823
        self.SRKGraboski=[0, 0]

        #TODO: Añadir parámetros, archivo /media/datos/Biblioteca/archivos/Melhem, Almeida - A data Bank of Parameters for the Attractive-Aznar Telles.pdf
        self.Melhem=[0, 0]          #Alcoholes en archivo de abajo
        self.Almeida=[0, 0]

        #TODO: Añadir parámetros, archivo /media/datos/Biblioteca/archivos/alfas.pdf
        self.Mathias=0
        self.MathiasCopeman=[0, 0, 0]
        self.Adachi=[0, 0]
        self.Andoulakis=[0, 0, 0]
        self.Yu_Lu=[0, 0, 0]


        #Desglosar formula en elementos y átomos de cada elemento
        formula=self.formula
        elementos=[]
        atomos=[]
        while len(formula)>0:
            letras=1
            numeros=0
            if len(formula)>1 and formula[1] in ascii_lowercase:
                letras+=1
            while len(formula)>letras+numeros and formula[letras+numeros] in digits:
                numeros+=1
            elementos.append(formula[:letras])
            if numeros==0:
                atomos.append(1)
            else:
                atomos.append(int(formula[letras:letras+numeros]))
            formula=formula[letras+numeros:]
        self.composicion_molecular=[elementos, atomos]
        if 'C' in elementos:
            self.C=atomos[elementos.index("C")]
        else:
            self.C=0
        if 'H' in elementos:
            self.H=atomos[elementos.index("H")]
        else:
            self.C=0
        if 'O' in elementos:
            self.O=atomos[elementos.index("O")]
        else:
            self.O=0
        if 'N' in elementos:
            self.N=atomos[elementos.index("N")]
        else:
            self.N=0
        if 'S' in elementos:
            self.S=atomos[elementos.index("S")]
        else:
            self.S=0

        if self.C and self.H:
            self.HC=self.H/self.C

    def tr(self,T):
       return T/self.Tc
    def pr(self,P):
        return P/self.Pc


#Metodos de estimacion de propiedades no disponibles en la base de datos
    def _Tc_Nokay(self, tipo):
        """Estimación de la temperatura crítica usando el método Nokay (Perry Cap II, Pag.342)
        Como parámetro opcional necesita el tipo de hidrocarburo de que se trata:
        0:Parafina
        1:Naftano
        2:Olefina
        3:Acetileno
        4:Diolefin
        5:Aromático
        Lo ideal sería no necesitar este parametro pero a ver como encuentro la forma, lease sin añadir ese parámetro en la base de datos de componentes"""
        #TODO: Puede que usando los grupos UNIFAC
        A=[1.35940, 0.65812, 1.09534, 0.74673, 0.1384, 1.0615]
        B=[0.43684, -0.07165, 0.27749, 0.30381, -0.39618, 0.22732]
        C=[0.56224, 0.81196, 0.65563, 0.79987, 0.99481, 0.66929]
        return unidades.Temperature(10**(A[tipo]+B[tipo]*log10(self.SG)+C[tipo]*log10(self.Tb)))

    def vc_Riedel(self):
        """Estimación del volumen crítico haciendo uso del método de Riedel. API procedure 4A3.1 pag. 302
        Volumen obtenido en l/mol"""
        riedel=5.811+4.919*self.f_acent
        return unidades.SpecificVolume(R_atml*self.Tc/self.Pc.atm/(3.72+0.26*(riedel-7))/self.M)

    def factor_acentrico_Lee_Kesler(self):
        """Estimación del factor acéntrico en componentes que no dispongan de ese valor, API,procedure 2A1.1 pag.210 """
        if self.Tb!=0:
            Tr=self.Tb/self.Tc
            Pr=1/self.Pc.atm
        else:
            Tr=0.7
            Pv=self.Pv_DIPPR(self.Tc*Tr)
            Pr=Pv.atm/self.Pc.atm
        return (log(Pr)-5.92714+6.09648/Tr+1.28862*log(Tr)-0.169347*Tr**6)/(15.2518-15.6875/Tr-13.4721*log(Tr)+0.43577*Tr**6)

    def factor_acentrico_Ambrose(self):
        """Estimación del factor acéntrico en componentes que no dispongan de ese valor,
       ref. propiedades físicas de líquidos y gases, pag.235 """
        if self.Tb!=0:
            Tr=self.Tb/self.Tc
            Pr=1/self.Pc.atm
        else:
            Tr=0.7
            pv=self.Pv_DIPPR(self.Tc*Tr)/self.Pc.atm
            t=1-Tr
            f0=(-5.97616*t+1.29874*t**1.5-0.60394*t**2.5-1.06841*t**5)/Tr
            f1=(-5.03365*t+1.11505*t**1.5-5.41217*t**2.5-7.46628*t**5)/Tr
            f2=(-0.64771*t+2.41539*t**1.5-4.26979*t**2.5+3.25259*t**5)/Tr
            coef=roots([f2, f1, f0-log(pv)])
            if absolute(coef[0])<absolute(coef[1]):
                self.f_acent=coef[0]
            else:
                self.f_acent=coef[1]

    def Rackett(self):
        """Método alternativa para el calculo de la constante de Rackett en aquellos componentes que no dispongan de ella a partir del factor acéntrico.
        Yamada, T., and R. Gunn. “Saturated Liquid Molar Volumes: The Rackett Equation,” Journal of Chemical Engineering Data 18, no. 2 (1973): 234–236."""
        return 0.29056-0.08775*self.f_acent

    def Volumen_Liquido_Constante(self):
        V=1/self.RhoL_Rackett(self.Tc)
        V=R_atml*1000*self.Tc/self.Pc.atm*self.rackett**(1+(1-self.tr(298.15))**(2.0/7)) #cm3/mol
        return V/(5.7+1611/self.Tc) #cm3/mol

    def Parametro_Solubilidad(self):
        """Método de cálculo del parametro de solubilidad, API databook pag 812"""
        if self.calor_vaporizacion==(0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0):
            Par=0
        else:
            DH=self.Hv_DIPPR(298.15).calg
            densidad=self.RhoL_DIPPR(298.15).gcc
            Par=sqrt(DH-R_cal*298.15/self.M*densidad)
        return unidades.SolubilityParameter(Par, "calcc")

    def Stiehl_Polar_factor(self):
        f0=5.92714-6.09648/0.6-1.28862*log(0.6)+0.169347*0.6**6
        f1=15.2518-15.6875/0.6-13.4721*log(0.6)+0.43577*0.6**6
        pv=exp(f0*0.6+self.f_acent*f1*0.6)
        return log10(self.pr(P)/pv)



#Propiedades ideales
    def Cp_ideal(self,T):
        """capacidad calorifica del gas ideal
        Temperatura dada en K
        cp obtenida en (cal/gmol-k)
        Los parámetros de la ecuación se encuentran en la base de datos
        en forma de lista en la posición octava [A,B,C,D,E,F]
        Cp = A + BT + CT^2 + DT^3 + ET^4 + FT^5"""
        cp=self.cp[0]+self.cp[1]*T+self.cp[2]*T**2+self.cp[3]*T**3+self.cp[4]*T**4+self.cp[5]*T**5
        return unidades.SpecificHeat(cp/self.M, "calgK")

    def Cv_ideal(self, T):
        """Capacidad calorífica isocórica (a volumen constante) ideal"""
        return unidades.SpecificHeat(self.Cp_ideal(T).JgK-R/self.M, "JgK")

    def Entalpia_ideal(self, T):
        """Entalpia del gas ideal usando los coeficientes de cp (su relación con la entalpia por integración) y suponiento H=0 a 0K
        API procedure 7A1.1 pag. 543"""
        H=self.cp[0]*T+self.cp[1]/2*T**2+self.cp[2]/3*T**3+self.cp[3]/4*T**4+self.cp[4]/5*T**5+self.cp[5]/6*T**6
        Ho=self.cp[0]*298.15+self.cp[1]/2*298.15**2+self.cp[2]/3*298.15**3+self.cp[3]/4*298.15**4+self.cp[4]/5*298.15**5+self.cp[5]/6*298.15**6
        return unidades.Enthalpy((H-Ho)/self. M, "calg")

    def Entropia_ideal(self, T):
        """Entropia del gas ideal usando los coeficientes de cp (su relación con la entropia por derivación) y suponiento S=0 a 1K
        API procedure 7A1.1 pag. 543"""
        s=self.cp[0]*log(T)+self.cp[1]*T+self.cp[2]/2*T**2+self.cp[3]/3*T**3+self.cp[4]/4*T**4+self.cp[5]/5*T**5
        return unidades.SpecificHeat(s/self.M, "calgK")


    def Entalpia_formacion(self, T):
        """Cálculo del calor de formación a una temperatura no estandart, API procedure 7A1.8, pag 581"""
        atomos=[64, 1, 47, 46, 992, 208, 105, 99, 213] #C,H,O,N,S,F,Cl,Br,I
        simbolos=["C", "H", "O", "N", "S", "F", "Cl", "Br", "I"]
        elementos=[]
        contribucion=[]
        for i in range(len(simbolos)):
            if simbolos[i] in self.composicion_molecular[0]:
                elementos.append(atomos[i])
                if i==0 or i==4:
                    contribucion.append(self.composicion_molecular[1][self.composicion_molecular[0].index(simbolos[i])])
                else:
                    contribucion.append(self.composicion_molecular[1][self.composicion_molecular[0].index(simbolos[i])]/2)

        Cpf=(self.cp[0]*T+self.cp[1]/2*T**2+self.cp[2]/3*T**3+self.cp[3]/4*T**4+self.cp[4]/5*T**5-(self.cp[0]*298.15+self.cp[1]/2*298.15**2+self.cp[2]/3*298.15**3+self.cp[3]/4*298.15**4+self.cp[4]/5*298.15**5))/self.M

        Cpi=0
        for i in range(len(elementos)):
            Cpi+=(databank.config.base_datos[elementos[i]][7][0]*contribucion[i]*T+databank.config.base_datos[elementos[i]][7][1]/2*contribucion[i]*T**2+databank.config.base_datos[elementos[i]][7][2]/3*contribucion[i]*T**3+databank.config.base_datos[elementos[i]][7][3]/4*contribucion[i]*T**4+databank.config.base_datos[elementos[i]][7][4]/5*contribucion[i]*T**5-(databank.config.base_datos[elementos[i]][7][0]*contribucion[i]*298.15+databank.config.base_datos[elementos[i]][7][1]/2*contribucion[i]*298.15**2+databank.config.base_datos[elementos[i]][7][2]/3*contribucion[i]*298.15**3+databank.config.base_datos[elementos[i]][7][3]/4*contribucion[i]*298.15**4+databank.config.base_datos[elementos[i]][7][4]/5*contribucion[i]*298.15**5))/databank.config.base_datos[elementos[i]][2]

        return unidades.Enthalpy(self.calor_formacion/self.M)+unidades.Enthalpy((Cpf-Cpi), "calg")


#Propieadades fisicas
    def DIPPR(self,T,parametros):
        """Función generica para calcular las propiedes fisicas paramétricas usando
        las ecuaciones DIPPR. Las propiedades disponibles y las unidades en que se obtienen los valores son:

            Solid density                       (kmol/m3)
            Liquid density                      (kmol/m3)
            Vapor pressure                      (Pa)
            Heat of vaporization                (J/kmol)
            Solid heat capacity                 (J/kmol-K)
            Liquid heat capacity                (J/kmol-K)
            Ideal gas heat capacity             (J/kmol-K)
            Liquid viscosity                    (Pa-sec)
            Vapor viscosity                     (Pa-sec)
            Liquid thermal conductivity         (W/m-K)
            Vapor thermal conductivity          (W/m-K)
            Surface Tension                     (N/m)

        A continuación las diferentes ecuaciones:

        Ecuación 1:     Y = A+B*T+C*T^2+D*T^3+E*T^4
        Ecuación 2:     Y = exp(A+B*T+C*ln(T)+D*T^E)
        Ecuación 3:     Y = A*T^B/(1+C*T+D*T^2)
        Ecuación 4:     Y = A+B*exp(-C/T^D)
        Ecuación 5:     Y = A + B*T + C*T^3 + D*T^8 + E*T^9
        Ecuación 6:     Y = A/(B^(1+(1-T/C)^D)
        Ecuación 7:     Y = A*(1-Tr)^(B+C*Tr+D*Tr^2+E*Tr^3)
        Ecuación 8:     Y = A+ B*((C/T)/sinh(C/T))^2 + D*((E/T)/cosh(E/T))^2
        Ecuación 9:     Y = A^2/Tr + B - 2ACTr - ADTr^2 - C^2Tr^3/3 - CDTr^4/2 - D^2Tr^5/5

        siendo: T la temperatura en kelvin
                Tr la temperatura reducida T/Tc
                A,B,C,D,E los parametros

        Estos parámetros vendrán dados en la base de datos, para cada propiedad física"""

        ecuacion=parametros[0]
        if ecuacion == 1:
            return parametros[1]+parametros[2]*T+parametros[3]*T**2+parametros[4]*T**3+parametros[5]*T**4
        elif ecuacion == 2:
            return exp(parametros[1]+parametros[2]/T+parametros[3]*log(T)+parametros[4]*T**parametros[5])
        elif ecuacion == 3:
            return parametros[1]*T**parametros[2]/(1+parametros[3]/T+parametros[4]/T**2)
        elif ecuacion == 4:
            return parametros[1]+parametros[2]*exp(-parametros[3]/T**parametros[4])
        elif ecuacion == 5:
            return parametros[1]+parametros[2]/T+parametros[3]/T**3+parametros[4]/T**8+parametros[5]/T**9
        elif ecuacion == 6:
            return parametros[1]/(parametros[2]**(1+((1-T/parametros[3])**parametros[4])))
        elif ecuacion == 7:
            return parametros[1]*(1-self.tr(T))**(parametros[2]+parametros[3]*self.tr(T)+parametros[4]*self.tr(T)**2+parametros[5]*self.tr(T)**3)
        elif ecuacion == 8:
            return parametros[1]+parametros[2]*(parametros[3]/T/sinh(parametros[3]/T))**2+parametros[4]*(parametros[5]/T/cosh(parametros[5]/T))**2
        elif ecuacion == 9:
            return parametros[1]**2/self.tr(T)+parametros[2]-2*parametros[1]*parametros[3]*self.tr(T)-parametros[1]*parametros[4]*self.tr(T)**2-parametros[3]**2*self.tr(T)**3/3-parametros[3]*parametros[4]*self.tr(T)**4/2-parametros[4]**2*self.tr(T)**5/5


    def RhoS(self,T):
        """Cálculo de la densidad del sólido usando las ecuaciones DIPPR"""
        return unidades.Density(self.DIPPR(T,self.densidad_solido)*self.M)


    def RhoL(self, T, P):
        """Procedimiento que define el método más apropiado para el calculo de la densidad del líquido"""
        rhoL=self.Config.getint("Transport","RhoL")
        corr=self.Config.getint("Transport","Corr_RhoL")
        if P<1013250:
            if rhoL==0 and self.densidad_liquido and self.densidad_liquido[6]<=T<=self.densidad_liquido[7]:
                return self.RhoL_DIPPR(T)
            elif rhoL==1 and self.rackett!=0 and T<self.Tc:
                return self.RhoL_Rackett(T)
            elif rhoL==2 and self.Vliq!=0:
                return self.RhoL_Cavett(T)
            elif rhoL==3:
                return self.RhoL_Costald(T)
            else:
                if self.densidad_liquido and self.densidad_liquido[6]<=T<=self.densidad_liquido[7]:
                    return self.RhoL_DIPPR(T)
                elif self.rackett!=0 and T<self.Tc:
                    return self.RhoL_Rackett(T)
                elif self.Vliq!=0:
                    return self.RhoL_Cavett(T)
                else:
                    return self.RhoL_Costald(T)
        else:
            if corr==0:
                return self.RhoL_Thomson_Brobst_Hankinson(T, P)
            else:
                return self.RhoL_API(T, P)

    def RhoL_DIPPR(self,T):
        """Cálculo de la densidad del líquido usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimocuarta posición
        densidad obtenida en kg/m3"""
        return unidades.Density(self.DIPPR(T,self.densidad_liquido)*self.M)

    def RhoL_Cavett(self, T):
        """Método alternativo para calcular la densidad de liquidos haciendo uso
        de la ecuación de Cavett:       V = Vol_Con * (5.7 + 3Tr)
        donde   V: volumen de liquido en cc/mol
                Vol_Con es una constante para cada compuesto, situada en la base de datos en el puesto veinticinco
                Tr es la temperatura reducida
                Densidad obtenida en g/l"""
        return unidades.Density(1/(self.Vliq*(5.7+3*self.tr(T)))*1000*self.M)

    def RhoL_Rackett(self, T):
        """Método alternativo para calcular la densidad de líquidos saturados haciendo uso
        de la ecuación modificada de Rackett:    V=R*Tc/Pc*(Zra)[1 + (1 − Tr)^2/7]
        API procedure 6A2.13 pag.454
        Spencer, F. F., and R. P. Danner. “Prediction of Bubble-Point Density of Mixtures,” Journal of Chemi-
   cal Engineering Data 18, no. 2 (1973): 230–234"""
        V=R_atml*self.Tc/self.Pc.atm*self.rackett**(1.+(1.-self.tr(T))**(2./7))
        return unidades.Density(1/V*self.M)

    def RhoL_Costald(self, T):
        """Método alternativo para el cálculo de la densidad de líquidos saturados
        API procedure 6A2.15 pag. 462"""
        if self.f_acent_mod!=0:
            w=self.f_acent_mod
        else: w=self.f_acent
        if self.V_char:
            V_=self.V_char
        else: V_=self.Vc

        Vr0=1-1.52816*(1-self.tr(T))**(1./3)+1.43907*(1-self.tr(T))**(2./3)-0.81446*(1-self.tr(T))+0.190454*(1-self.tr(T))**(4./3)
        Vr1=(-0.296123+0.386914*self.tr(T)-0.0417258*self.tr(T)**2-0.0480645*self.tr(T)**3)/(self.tr(T)-1.00001)
        #TODO: Añadiendo V* a la base de datos mejoraría la precisión de este método, en vez de usar el volumen critico, porque la constante de volumen de líquido no parece corresponder a esta constante
        return unidades.Density(1/(V_*Vr0*(1-w*Vr1))*self.M)

    def RhoL_Thomson_Brobst_Hankinson(self, T, P):
        """Método alternativo para el cálculo de la densidad de líquidos comprimidos
        API procedure 6A2.23 pag. 477"""
        if self.f_acent_mod:
            w=self.f_acent_mod
        else: w=self.f_acent
        rho_sat=self.RhoL(T, 101325)
        pv_sat=self.Pv_DIPPR(T)
        C=0.0861488+0.0344483*w
        e=exp(4.79594+0.250047*w+1.14188*w**2)
        B=self.Pc*(-1-9.070217*(1-self.tr(T))**(1./3)+62.45326*(1-self.tr(T))**(2./3)-135.1102*(1-self.tr(T))+e*(1-self.tr(T))**(4./3))
        return unidades.Density(rho_sat/(1-C*log((B+P)/(B+pv_sat))), "gl")

    def RhoL_API(self, T, P):
        """Método alternativo para calcular la densidad de líquidos haciendo uso del método API
        Dens2=dens1*C2/C1
        donde:  C2,C1: valor empirico que se calcula a partir de una tabla de valores, C2 a T y P dadas, C1 a 60F y 1 atm
                dens1: densidad a 60ºF y 1 atm, situada en la base de datos en el puesto veintisiete
                densidad obtenida en mol/l"""
        pr=P/self.Pc
        A02=1.6368-0.04615*pr+2.1138e-3*pr**2-0.7845e-5*pr**3-0.6923e-6*pr**4
        A12=-1.9693-0.21874*pr-8.0028e-3*pr**2-8.2328e-5*pr**3+5.2604e-6*pr**4
        A22=2.4638-0.36461*pr-12.8763e-3*pr**2+14.8059e-5*pr**3-8.6895e-6*pr**4
        A32=-1.5841-0.25136*pr-11.3805e-3*pr**2+9.5672e-5*pr**3+2.1812e-6*pr**4
        C2=A02+A12*self.tr(T)+A22*self.tr(T)**2+A32*self.tr(T)**3
        A01=1.6368-0.04615*self.pr(1)+2.1138e-3*self.pr(1)**2-0.7845e-5*self.pr(1)**3-0.6923e-6*self.pr(1)**4
        A11=-1.9693-0.21874*self.pr(1)-8.0028e-3*self.pr(1)**2-8.2328e-5*self.pr(1)**3+5.2604e-6*self.pr(1)**4
        A21=2.4638-0.36461*self.pr(1)-12.8763e-3*self.pr(1)**2+14.8059e-5*self.pr(1)**3-8.6895e-6*self.pr(1)**4
        A31=-1.5841-0.25136*self.pr(1)-11.3805e-3*self.pr(1)**2+9.5672e-5*self.pr(1)**3+2.1812e-6*self.pr(1)**4
        t1=unidades.Temperature(60, "F")
        C1=A01+A11*self.tr(t1)+A21*self.tr(t1)**2+A31*self.tr(t1)**3
        d2=self.SG*1000*C2/C1
        #FIXME: C2 no sale bien, sale muchas veces negativo, repaso, repaso, pero no se porqué

        return unidades.Density(d2, "gl")


    def Pv(self, T):
        """Procedimiento que define el método más apropiado para el cálculo de la presión de vapor"""
        Pv=self.Config.getint("Transport","Pv")
        T=unidades.Temperature(T)
        if Pv==0 and self.presion_vapor and self.presion_vapor[6]<=T<=self.presion_vapor[7]:
            return self.Pv_DIPPR(T)
        elif Pv==1 and self.antoine:
            return self.Pv_Antoine(T)
        elif Pv==2 and self.Pc and self.Tc and self.f_acent:
            return self.Pv_Lee_Kesler(T)
        elif Pv==3 and self.Kw and self.Tb:
            return self.Pv_Maxwell_Bonnel(T)
        elif Pv==4 and self.wagner and self.wagner[4]<=T.R<=self.wagner[5]:
            return self.Pv_Wagner(T)
        else:
            if self.presion_vapor and self.presion_vapor[6]<=T<=self.presion_vapor[7]:
                return self.Pv_DIPPR(T)
            elif self.antoine:
                return self.Pv_Antoine(T)
            elif self.Pc and self.Tc and self.f_acent:
                return self.Pv_Lee_Kesler(T)
            elif self.Kw and self.Tb:
                return self.Pv_Maxwell_Bonnel(T)
            elif self.wagner and self.wagner[4]<=T.R<=self.wagner[5]:
                return self.Pv_Wagner(T)
            else:
                print("Ningún método disponible")

    def Pv_DIPPR(self,T):
        """Cálculo de la presión de vapor usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimoquinta posición
        Presión de vapor obtenida en Pa"""
        return unidades.Pressure(self.DIPPR(T,self.presion_vapor))

    def Pv_Antoine(self,T,parameters=None):
        """Presión de vapor de Antoine
        Los parámetros de la ecuación se encuentran en la base de datos"""
        if parameters==None:
            parameters=self.antoine
        return unidades.Pressure(exp(parameters[0]-parameters[1]/(T+parameters[2])), "mmHg")

    def Pv_Lee_Kesler(self, T):
        """Denominado en Chemcad Curl Pitzer.
        Método alternativo para calcular la presión de vapor, usando las
        propiedades críticas, cuando no están disponibles los parametros DIPPR
        ni de Antoine pero si las propiedades críticas, API procedure 5A1.16, pag 390"""
        f0=5.92714-6.09648/self.tr(T)-1.28862*log(self.tr(T))+0.169347*self.tr(T)**6
        f1=15.2518-15.6875/self.tr(T)-13.4721*log(self.tr(T))+0.43577*self.tr(T)**6
        return unidades.Pressure(exp(f0+self.f_acent*f1)*self.Pc.atm, "atm")

    def Pv_Wagner(self, T):
        """Método alternativo para el cálculo de la presión de vapor, API procedure 5A1.3 pag 366"""
        Tr=self.tr(T)
        X1=(1-Tr)/Tr
        X2=(1-Tr)**1.5/Tr
        X3=(1-Tr)**2.6/Tr
        X4=(1-Tr)**5./Tr
        return unidades.Pressure(exp(self.wagner[0]*X1+self.wagner[1]*X2+self.wagner[2]*X3+self.wagner[3]*X4)*self.Pc)

    def Pv_Maxwell_Bonnel(self, T):
        """Método alternativo de cálculo de la presión de vapor, cuando no se dispone de los parametros DIPPR, ni de los valores de las propiedades críticas del elemento. Necesita el factor de Watson. API procedure 5A1.18  Pag. 394"""
        if self.Tb.F>400: f=1.0
        elif self.Tb.F <200: f=0.0
        else: f=(self.Tb.R-659.7)/200

        def P(Tb):
            X=(Tb/T/1.8-0.0002867*Tb)/(748.1-0.2145*Tb)
            if X > 0.0022: p=10**((3000.538*X-6.761560)/(43*X-0.987672))
            elif X < 0.0013: p=10**((2770.085*X-6.412631)/(36*X-0.989679))
            else: p=10**((2663.129*X-5.994296)/(95.76*X-0.972546))
            return p

        Tb=fsolve(lambda Tb: Tb-self.Tb.R+2.5*f*(self.Kw-12)*log10(P(Tb)/760), self.Tb)
        p=P(Tb)
        return unidades.Pressure(p, "mmHg")


    def ThCond_Liquido(self, T, P):
        """Procedimiento que define el método más apropiado para el cálculo de la conductividad térmica del líquido, pag 1135"""
        ThCondL=self.Config.getint("Transport","ThCondL")
        corr=self.Config.getint("Transport","Corr_ThCondL")
        p=unidades.Pressure(P, "atm")
        if p.psi<500:
            if ThCondL==0 and self.conductividad_liquido and self.conductividad_liquido[6]<=T<=self.conductividad_liquido[7]:
                return self.ThCond_Liquido_DIPPR(T)
            elif ThCondL==1 and T<self.Tc:
                return self.ThCond_Liquido_Pachaiyappan(T)
            else:
#                print "Warning: Thermal conductivity of %s out of range" % self.nombre
                return self.ThCond_Liquido_DIPPR(self.conductividad_liquido[7])

        else:
            if corr==0:
                return self.ThCond_Liquido_Lenoir(T, P)
            else:
                return self.ThCond_Liquido_Kanitkar_Thodos(T, P)

    def ThCond_Liquido_DIPPR(self,T):
        """Cálculo de la conductividad terminca del líquido usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en vigesimosegunda posición
        Conductividad termica obtenida en W/(m·K), API procedure 12A1.1, pag 1137"""
        return unidades.ThermalConductivity(self.DIPPR(T,self.conductividad_liquido))

    def ThCond_Liquido_Pachaiyappan(self, T):
        """Método alternativo para el cálculo de la conductividad de líquidos a baja presión, API procedure 12A1.2, pag 1141"""
        #TODO: Alguna forma mejor de definir los compuestos lineales hay que encontrar, quizá añadiendo un parámetro a la base de datos , branched True o False
        if self.indice in [1, 2, 3, 4, 6, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 90, 91, 92, 1763, 1765, 1767, 1738, 1769, 1770, 1844, 1741, 1842, 1771, 1742, 22, 23, 24, 25, 26, 28, 29, 30, 31, 35, 56, 57, 58]:
            n=1.001
            C=1.676e-3
        else:
            n=0.7717
            C=4.079e-3

        rho=unidades.Density(self.RhoL(293.15, 1)/self.M).lbft3
        Vm=1/rho
        k=C*self.M**n/Vm*(3+20*(1-self.tr(T))**(2./3))/(3+20*(1-293.15/self.Tc)**(2./3))
        return unidades.ThermalConductivity(k, "BtuhftF")

    def ThCond_Liquido_Riazi_Faghri(self, T):
        """Método de cálculo de la conductividad térmica de hidrocarburos pesados a baja presión,
        Riazi, M. R. and Faghri, A., "Thermal Conductivity of Liquid and Vapor Hydrocarbon Systems: Pentanes and Heavier at Low Pressures," Industrial and Engineering Chemistry, Process Design and Development, Vol. 24, No. 2, 1985, pp. 398-401."""
        t=(1.8*T-460)/100
        A=exp(-4.5093-0.6844*t-0.1305*t**2)
        B=0.3003+0.0918*t+0.01195*t**2
        C=0.1029+0.0894*t+0.0292*t**2
        k=1.7307*A*self.Tb.R**B*self.SG**C
        return unidades.Conductividad_termica(k)


    def ThCond_Liquido_Kanitkar_Thodos(self, T, P):
        """Método alternativo para el cálculo de la conductividad de líquidos por encima del punto de ebullición, API procedure 12A1.3, pag 1143"""
        Pr=self.pr(P)
        rho=self.RhoL(T, P)/self.M
        rhor=rho*self.Vc

        l=self.Tc.R**(1./6)*self.M**0.5/self.Pc.atm**(2./3)
        b=0.4+0.986/exp(0.58*l)
        a=7.137e-3/b**3.322
        k=(-1.884e-6*Pr**2+1.442e-3*Pr+a*exp(b*rhor))/l
        return unidades.ThermalConductivity(k, "BtuhftF")

    def ThCond_Liquido_Lenoir(self, T, P, ko=[]):
        """Método alternativo para el cálculo de la conductividad de líquidos a alta presión, API procedure 12A4.1, pag 1156
        Opcionalmente puede aceptar otro parametros ko que indica un valor experimental de la conductividad térmica (WmK, así como la temperatura y presión a la que se da, en un array [k,T,P]"""
        Tr=self.tr(T)
        if ko==[]:
            k1=self.ThCond_Liquido(T, 1)
            C1=17.77+0.065*self.pr(1)-7.764*Tr-2.065*Tr**2/exp(0.2*self.pr(1))
        else:
            k1=ko[0]
            C1=17.77+0.065*self.pr(ko[2])-7.764*self.tr(ko[1])-2.065*self.tr(ko[1])**2/exp(0.2*self.pr(ko[2]))
        C2=17.77+0.065*self.pr(P)-7.764*Tr-2.065*Tr**2/exp(0.2*self.pr(P))
        k=k1*C2/C1
        return unidades.ThermalConductivity(k)


    def ThCond_Gas(self, T, P):
        """Procedimiento que define el método más apropiado para el cálculo de la conductividad térmica del líquido, pag 1136"""
        ThCondG=self.Config.getint("Transport","ThCondG")
        p=unidades.Pressure(P)
        if p.psi<50:
            if ThCondG==0 and self.conductividad_gas and self.conductividad_gas[6]<=T<=self.conductividad_gas[7]:
                return self.ThCond_Gas_DIPPR(T)
            else:
                return self.ThCond_Gas_Misic_Thodos(T)
        elif self.indice in [1, 46, 47, 48, 50, 51, 111]:
            return self.ThCond_Gas_Nonhidrocarbon(T, P)
        else:
            # TODO: fix crooks methods with lost lee_kesler library
            return self.ThCond_Gas(T, 101325.)
#            return self.ThCond_Gas_Crooks(T, P)

    def ThCond_Gas_DIPPR(self,T):
        """Cálculo de la conductividad terminca del gas usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en vigesimotercera posición
        Conductividad termica obtenida en W/(m·K), API procedure 12B1.1 pag 1158"""
        return unidades.ThermalConductivity(self.DIPPR(T,self.conductividad_gas))

    def ThCond_Gas_Misic_Thodos(self, T):
        """Método alternativo para el cálculo de la conductividad térmica de gases a baja presión <4 atm, API Procedure 12B1.2 pag.1162"""
        l=(self.Tc.R)**(1./6)*self.M**0.5*(1/self.Pc.atm)**(2./3)
        cp=self.Cp_Gas_DIPPR(T)
        #TODO: Cuando se añada alguna propiedad en la base de datos que defina la naturaleza cíclica de los componentes se podrá generalizar este metodo a una temperatura reducida menor de 1 para compuestos cíclicos e hidrógeno, de momento todos se calculan por el metodo de tr mayor de 1.
#        if self.tr(T)<1:
#            k=1.188e-3*self.tr(T)*cp.BtulbF/l
#        else:
        k=2.67e-4*(14.52*self.tr(T)-5.14)**(2.0/3)*cp.BtulbF*self.M/l
        return unidades.ThermalConductivity(k, "BtuhftF")

    def ThCond_Gas_Crooks(self, T, P):
        """Método alternativo para el cálculo de la conductividad térmica de gases a alta presión >4 atm, API Procedure 12B4.1 pag.1170"""
        Tr=self.tr(T)
        Pr=self.pr(P)
        k=self.ThCond_Gas(T, 1)

        A=-0.0617*exp(1.91/Tr**9)
        B=2.29*exp(1.34/Tr**16)
        k_prima=1+(4.18/Tr**4+0.537*Pr/Tr**2)*(1-exp(A*Pr**B))+0.510*Pr/Tr**3*exp(A*Pr**B)
        # TODO: particularizar a compuestos cíclicos y no cíclicos
        k_prima2=1+1./Tr**5*Pr**4/(2.44*Tr**20+Pr**4)+0.012*Pr/Tr
        Cv=self.Cv_Lee_Kesler(T, P, 1)
        Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(Tr, Pr)
        Cv_prima=4.965-R_Btu*(1+Cv0)
        Cv_prima2=Cv.BtulbF-Cv_prima
        return unidades.ThermalConductivity(k*(Cv_prima/Cv.BtulbF*k_prima+Cv_prima2/Cv.BtulbF*k_prima2))

    def ThCond_Gas_Nonhidrocarbon(self, T, P):
        """Método de cálculo de la conductividad térmica de gases no hidrocarburos a alta presión, API procedure 12C1.1, pag 1174
        """
        t=unidades.Temperature(T)
        p=unidades.Pressure(P, "atm")
        if self.indice==1:
            par=[4.681e-3, 2.e-4, -3.6e-8, 0.0, 0.0, 0.0, 1.7e-3]
        elif self.indice==46:
            par=[4.561e-3, 1.61e-5, 0.0, 2.56e-9, 5.299e-3, 2.47e-3, 0.0]
        elif self.indice==47:
            par=[5.95e-4, 1.71e-5, 0.0, -2.10e-8, 5.869e-3, 6.995e-3, 0.0]
        elif self.indice==48:
            par=[1.757e-3, 1.55e-5, 0.0, 2.08e-8, 5.751e-3, 5.6e-3, 0.0]
        elif self.indice==50:
            par=[-1.51e-3, 2.25e-5, 3.32e-10, 0.0, 0.0, 0.0, 0.0]
        elif self.indice==51:
            par=[2.5826e-2, 1.35e-5, 0.0, -4.4e-7, 1.026e-2, -2.631e-2, 0.0]
        elif self.indice==111:
            par=[-1.02e-3, 1.35e-5, 4.17e-9, 0.0, 0.0, 0.0, 0.0]

        k=par[0]+par[1]*t.R+par[2]*t.R**2+par[3]*p.psi+par[4]*p.psi/t.R**1.2+par[5]/(0.4*p.psi-0.001*t.R)**0.015+par[6]*log(p.psi)
        return unidades.ThermalConductivity(k, "BtuhftF")


    def Mu_Gas(self, T, P):
        """Procedimiento que define el método más apropiado para el cálculo de la viscosidad del líquido, pag 1026"""
        MuG=self.Config.getint("Transport","MuG")
        if P/self.Pc<0.6:
           if MuG==0 and self.viscosidad_gas and self.viscosidad_gas[6]<=T<=self.viscosidad_gas[7]:
                return self.Mu_Gas_DIPPR(T)
           elif MuG==1:
               return self.Mu_Gas_Chapman_Enskog(T)
           else :
               return self.Mu_Gas_Thodos(T)
        else:
            if self.hidrocarburo:
                return self.Mu_Gas_Eakin_Ellingtong(T, P)
            else:
                return self.Mu_Gas_Carr(T, P)

    def Mu_Gas_DIPPR(self,T):
        """Cálculo de la viscosidad del gas usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en vigesimoprimera posición
        Viscosidad obtenida en Pa·s
        API procedure 11B1.1, pag 1091"""
        return unidades.Viscosity(self.DIPPR(T,self.viscosidad_gas))

    def Mu_Gas_Chapman_Enskog(self,T):
        """Método alternativo para calcular la viscosidad de gases (a baja presión):
        μg = 5/16(πMRT)^0.5/πO²Ω = 22.69*(MT)^0.5/O²Ω = Mp*P
        where:        M  = peso molecular
                      R  = constante del gas ideal
                      T  = temperatura
                      P  = presión
                      O  = diametro molécular
                      Ωv = colision integral
                      Mp = momento dipolar
        ref chemcad pag 81
        ref Properties gases and liquids pag. 470
        """
        #Metodo de Chung (pouling pag 473)
        if not self.diametro_molecular:
            diametro_molecular=0.809*self.Vc*self.M**(1./3)
        else: diametro_molecular=self.diametro_molecular
        if not self.ek:
            ek=self.Tc/1.2593
        else: ek=self.ek

        T_=T/ek
        if self.parametro_polar: #Polar, colisión integral de Brokaw (pouling pag 640)
            omega=1.03036/T_**0.15610+0.193/exp(0.47635*T_)+1.03587/exp(1.52996*T_)+1.76474/exp(3.89411*T_)+0.19*self.parametro_polar**2/T_
        else: #No polar, colisión integral de Neufeld
            omega=1.16145/T_**0.14874+0.52487/exp(0.7732*T_)+2.16178/exp(2.43787*T_)
        return unidades.Viscosity(26.69*(self.M*T)**0.5/diametro_molecular**2/omega, "microP")

    def Mu_Gas_Thodos(self, T):
        """Método alternativo para el cálculo de la viscosidad de gases a baja presión, solo necesita las propiedades críticas, API procedure 11B1.3, pag 1099"""
        t=unidades.Temperature(T)
        Tr=self.tr(T)
        if self.indice==1:
            if Tr<=1.5:
                mu=3.7e-5*t.R**0.94
            else:
                mu=9.071e-4*(7.639e-2*t.R-1.67)**0.625
        else:
            if Tr<=1.5:
                N=3.5e-4*Tr**0.94
            else:
                N=1.778e-4*(4.58*Tr-1.67)**0.625
            x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)
            mu=N/x
        return unidades.Viscosity(mu, "cP")

    def Mu_Gas_Jossi(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de hidrocarburos gaseosos pesados a alta presión,
        Jossi, J. A., Stiel, L. I., and Thodos, G., "The Viscosity of Pure Substances in the Dense Gaseous and Liquid Phases," American Institute of Chemical Engineers Journal, Vol. 8, 1962, pp. 59-63."""
        if muo==0:
            muo=self.Mu_Gas(T, 1)
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)
        rhor=self.RhoG_Lee_Kesler(T, P)*self.Vc*self.M
        mu=((0.1023+0.023367*rhor+0.058533*rhor**2-0.040758*rhor**3+0.0093324*rho**4)**4-1e-4)*x+muo
        return unidades.Viscosity(mu, "cP")

    def Mu_Gas_Eakin_Ellingtong(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de hidrocarburos gaseosos a alta presión, API procedure 11B4.1, pag 1107"""
        if muo==0:
            muo=self.Mu_Gas(T, 101325)
        else:
            muo=unidades.Viscosity(muo)
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)
        rhor=self.RhoG_Lee_Kesler(T, P)*self.Vc*self.M
        mu=muo.cP+10.8e-5*(exp(1.439*rhor)-exp(-1.11*rhor**1.858))/x
        return unidades.Viscosity(mu, "cP")

    def Mu_Gas_Carr(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de no hidrocarburos gaseosos a alta presión, API procedure 11C1.2, pag 1113"""
        Tr=self.tr(T)
        Pr=self.pr(P)
        if muo==0:
            muo=self.Mu_Gas(T, 101325)

        A1=83.8970*Tr**0.0105+0.6030*Tr**-0.0822+0.9017*Tr**-0.12-85.3080
        A2=1.514*Tr**-11.3036+0.3018*Tr**-0.6856+2.0636*Tr**-2.7611
        k=A1*1.5071*Pr**-0.4487+A2*(11.4789*Pr**0.2606-12.6843*Pr**0.1773+1.6953*Pr**-0.1052)
        return unidades.Viscosity(muo*k)

    def Mu_Gas_Stiel_Thodos(self, T, P, muo):
        """Método de cálculo de la viscosidad de gases polares a alta presión,
        Stiel, L. I. and Thodos, G., "The Viscosity of Polar Substances in the Dense Gaseous and Liquid Regions," American Institute of Chemical Engineers Journal, Vol. 10, No. 2, 1964, pp. 275-277."""
        if muo==0:
            muo=self.Mu_Gas(T, 1)
        rhor=self.RhoG_Lee_Kesler(T, P)*self.Vc*self.M
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)

        if rhor<=0.1:
            mu=1.626e-4*rhor**1.111/x+muo
        elif 0.1<rhor<=0.9:
            mu=6.07e-6*(9.045*rhor+0.63)**1.739/x+muo
        else:
            mu=10**(4-10**(0.6239-0.1005*rhor))/x/1e4+muo
        return unidades.Viscosity(mu, "cP")


    def Mu_Liquido(self, T, P):
        """Procedimiento que define el método más apropiado para el cálculo de la viscosidad del líquido, pag 1026"""
        #Comparacion de métodos: pag 405 Vismanath
        MuL=self.Config.getint("Transport","MuL")
        corr=self.Config.getint("Transport","Corr_MuL")
        if P<1013250:
            if MuL==0 and self.viscosidad_liquido and self.viscosidad_liquido[6]<=T<=self.viscosidad_liquido[7]:
                return self.Mu_Liquido_DIPPR(T)
            elif MuL==1 and self.viscosidad_parametrica[0]!=0 and self.viscosidad_parametrica[1]!=0:
                return self.Mu_Liquido_Parametrica(T)
            elif MuL==2 and self.Tc!=0 and self.Pc!=0 and self.f_acent!=0:
                return self.Mu_Liquido_Letsou_Steil(T)
            elif MuL==3 and self.Van_Veltzen:
                return self.Mu_Liquido_Van_Veltzen(T)
            else:
                if  self.viscosidad_liquido and self.viscosidad_liquido[6]<=T<=self.viscosidad_liquido[7]:
                    return self.Mu_Liquido_DIPPR(T)
                elif self.viscosidad_parametrica[0]!=0 and self.viscosidad_parametrica[1]!=0:
                    return self.Mu_Liquido_Parametrica(T)
                elif self.Tc!=0 and self.Pc!=0 and self.f_acent!=0:
                    return self.Mu_Liquido_Letsou_Steil(T)
                elif self.Van_Veltzen:
                    return self.Mu_Liquido_Van_Veltzen(T)
                else: print("Ningún método disponible")
        else:
            if corr==0 and self.Tb<650:
            #En realidad el criterio de corte es los hidrocarburos de menos de 20 átomos de carbono (hidrocarburos de bajo peso molecular), pero aprovechando que la temperatura de ebullición es proporcional al peso molecular podemos usar esta
                return self.Mu_Liquido_Graboski_Braun(T, P)
            elif corr==1:
                return self.Mu_Liquido_Kouzel(T, P)
            elif corr==2 and self.Pc and self.f_acent:
                return self.Mu_Liquido_Lucas(T, P)
            else:
                if self.Tb<650:
                    return self.Mu_Liquido_Graboski_Braun(T, P)
                elif self.Pc and self.f_acent:
                    return self.Mu_Liquido_Lucas(T, P)
                else:
                    return self.Mu_Liquido_Kouzel(T, P)

    def Mu_Liquido_DIPPR(self,T):
        """Cálculo de la viscosidad del líquido usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en vigésima posición
        Viscosidad obtenida en Pa·s
        API procedure 11A2.1, pag 1038"""
        return unidades.Viscosity(self.DIPPR(T,self.viscosidad_liquido))

    def Mu_Liquido_Parametrica(self, T, parameters=None):
        """Cálculo paramétrico de la viscosidad del líquido"""
        if parameters==None:
            parameters=self.viscosidad_parametrica
        return unidades.Viscosity(10**(parameters[0]*(1./T-1/parameters[1])), "cP")

    def Mu_Liquido_Van_Veltzen(self, T):
        """Método alternativo para calcular la viscosidad de líquidos haciendo uso de la contribución de los grupos moleculares, API procedure 11A2.3, pag 1048"""
        #TODO: Método dificil de implementar debido a que se tiene que calcular la contribución de grupos

    def Mu_Liquido_Letsou_Steil(self, T):
        """Método alternativo para el cálculo de la viscosidad en líquidos."""
        x= self.Tc**(1./6)/self.M**0.5/self.Pc.atm**(2./3)
        x0=0.015178-0.021351*self.tr(T)+0.007503*self.tr(T)**2
        x1=0.042559-0.07675*self.tr(T)+0.034007*self.tr(T)**2
        return unidades.Viscosity((x0+self.f_acent*x1)/x, "cP")

    def Mu_Liquido_Lucas(self, T, P, muo=0):
        """Método de cálculo de la viscosidad de líquidos a alta presión
        ref Viswanath, Viscosidad de líquidos, pag 145
        ref Poling, The Properties of Gases and Liquids, pag 521"""
        P=unidades.Pressure(P, "atm")
        Tr=self.tr(T)
        if muo==0:
            muo=self.Mu_Liquido(T, 1)
        else:
            muo=unidades.Viscosity(muo)
        Pvp=self.Pv(T)
        DPr=(P.atm-Pvp.atm)/self.Pc.atm
        A=0.9991-(4.674e-4/(1.0523*Tr**-0.03877-1.0513))
        C=-0.07921+2.1616*Tr-13.4040*Tr**2+44.1706*Tr**3-84.8291*Tr**4+96.1209*Tr**5-59.8127*Tr**6+15.6719*Tr**7
        D=0.3257/(1.0039-Tr**2.573)**0.2906-0.2086
        mu=muo*(1+D*(DPr/2.118)**A)/(1+C*self.f_acent*DPr)
        return unidades.Viscosity(mu)

    def Mu_Liquido_Graboski_Braun(self, T, P):
        """Método alternativo para el cálculo de la viscosidad en líquidos de bajo peso molécular a altas presiones, API procedure 11A5.1 (pag. 1074)"""
        Pr=self.pr(P)
        Tr=self.tr(T)

        A1=3.0294*Tr**9.0740+0.0032*Tr**10.9399-0.3689
        A2=-0.038*Tr**-7.2309+0.0229*Tr**11.7631+0.5781
        A3=-0.1415*Tr**27.2842+0.0778*Tr**-4.3406+0.0014
        A4=0.0028*Tr**69.4404-0.0042*Tr**3.3586+0.0062
        A5=0.0107*Tr**-7.4626-85.8276*Tr**0.1392+87.3164
        mur0=A1*log10(Pr)+A2*log10(Pr)**2+A3*Pr+A4*Pr**2+A5

        if Pr <= 0.75:
            B1=-0.2462*Tr**0.0484-0.7275*log(Tr)-0.0588*Tr+0.0079
            B2=-0.3199*Tr**17.0626-0.0695*log(Tr)+0.1267*Tr-0.0101
            B3=4.7217*Tr**-1.9831+19.2008*Tr**-1.7595+65.5728*log(Tr)+0.6110*Tr-19.1590
        else:
            B1=-0.0214*Tr**0.0484-0.1827*log(Tr)-0.0183*Tr+0.0090
            B2=-0.3588*Tr**5.0537-0.1321*log(Tr)+0.0204*Tr-0.0075
            B3=3.7266*Tr**-2.5689+52.1358*Tr**0.3514-13.0750*log(Tr)+0.6358*Tr-56.6687
        mur1=B1*Pr+B2*log(Pr)+B3

        mur=mur0+self.f_acent*mur1
        #TODO: Para implementar correctamente este método hay que añadir a la base de datos los datos de la viscosidad en el punto crítico. De momento usaremos el procedimiento DIPPR como alternativa para obtener un valor de viscosidad en las condiciones estandart y una minitabla para los elementos que aparecen en la tabla 11A5.4 del API Databook
        if 1<self.indice<=21 and self.indice!=7:
            muc=[0, 0, 0.014, 0.02, 0.0237, 0.027, 0.0245, 0, 0.0255, 0.0350, 0.0264, 0.0273, 0.0282, 0.0291, 0.0305, 0.0309, 0.0315, 0.0328, 0.0337, 0.0348, 0.0355, 0.0362][self.indice]
        elif self.indice==90:
            muc=0.0370
        elif self.indice==91:
            muc=0.0375
        elif self.indice==92:
            muc=0.0388
        else:
#            muc=self.Mu_critica().cP
            To=298.15
            Po=101325
            muo=self.Mu_Liquido(To, 101325)

            Pr=self.pr(Po)
            Tr=self.tr(To)

            A1=3.0294*Tr**9.0740+0.0032*Tr**10.9399-0.3689
            A2=-0.038*Tr**-7.2309+0.0229*Tr**11.7631+0.5781
            A3=-0.1415*Tr**27.2842+0.0778*Tr**-4.3406+0.0014
            A4=0.0028*Tr**69.4404-0.0042*Tr**3.3586+0.0062
            A5=0.0107*Tr**-7.4626-85.8276*Tr**0.1392+87.3164
            muor0=A1*log10(Pr)+A2*log10(Pr)**2+A3*Pr+A4*Pr**2+A5

            if Pr <= 0.75:
                B1=-0.2462*Tr**0.0484-0.7275*log(Tr)-0.0588*Tr+0.0079
                B2=-0.3199*Tr**17.0626-0.0695*log(Tr)+0.1267*Tr-0.0101
                B3=4.7217*Tr**-1.9831+19.2008*Tr**-1.7595+65.5728*log(Tr)+0.6110*Tr-19.1590
            else:
                B1=-0.0214*Tr**0.0484-0.1827*log(Tr)-0.0183*Tr+0.0090
                B2=-0.3588*Tr**5.0537-0.1321*log(Tr)+0.0204*Tr-0.0075
                B3=3.7266*Tr**-2.5689+52.1358*Tr**0.3514-13.0750*log(Tr)+0.6358*Tr-56.6687
            muor1=B1*Pr+B2*log(Pr)+B3
            muor=muor0+self.f_acent*muor1
            muc=muo.cP/muor

        return unidades.Viscosity(mur*muc, "cP")

    def Mu_Liquido_Kouzel(self, T, P, mua=0):
        """Método alternativo para el cálculo de la viscosidad en líquidos de alto peso molécular a altas presiones, API procedure 11A5.5 (pag. 1081)
        como parámetro opcional se puede indicar la viscosidad a presión atmosferica a tempratura T"""
        p=unidades.Pressure(P, "atm")
        if mua==0:
            mua=self.Mu_Liquido(T, 1)
        mup=mua*10**(p.psig/1000*(-0.0102+0.04042*mua.cP**0.181))
        return unidades.Viscosity(mup)

    def Mu_critica(self):
        """Procedimiento que define la viscosidad crítica si no viene en la base de datos
        Eq 8.9 Riazi - Characterization and Properties of Petroleum fractions, pag 348"""
        x=self.Tc**(1.0/6)/self.M**0.5/self.Pc.atm**(2.0/3)
        return unidades.Viscosity(7.7e-4/x, "cP")

    def Tension(self, T):
        """Procedimiento que define el método más apropiado para el cálculo de la tensión superficial"""
        tension=self.Config.getint("Transport","Tension")
        if tension==0 and self.tension_superficial and self.tension_superficial[6]<=T<=self.tension_superficial[7]:
            return self.Tension_DIPPR(T)
        elif tension==1 and self.tension_superficial_parametrica:
            return self.Tension_Parametrica(T)
        elif tension==2 and self.parachor:
            return self.Tension_Parachor(T)
        elif tension==3:
            return self.Tension_MIller(T)
        elif tension==4 and self.stiehl:
            return self.Tension_Hakim(T)
        elif tension==5 and self.Kw:
            return self.Tension_Hydrocarbon(T)
        else:
            if self.tension_superficial and self.tension_superficial[6]<=T<=self.tension_superficial[7]:
                return self.Tension_DIPPR(T)
            elif self.tension_superficial_parametrica:
                return self.Tension_Parametrica(T)
            elif self.parachor:
                return self.Tension_Parachor(T)
            elif self.stiehl:
                return self.Tension_Hakim(T)
            elif self.Kw:
                return self.Tension_Hydrocarbon(T)
            else:
                return self.Tension_MIller(T)

    def Tension_DIPPR(self,T):
        """Cálculo de la tensión superficial del líquido usando las ecuaciones DIPPR"""
        return unidades.Tension(self.DIPPR(T,self.tension_superficial))

    def Tension_Parametrica(self, T, parameters=None):
        """Cálculo paramétrico de la tensión superficial
        ST = A*(1-Tr)^B
        donde:  ST Tensión superficial in Newtons/metro
                Tr Temperatura reducida
        A,B: Los parámetros de la ecuación que se encuentran en la base de datos
        en forma de lista en la posición duodécima"""
        if parameters==None:
            parameters=self.tension_superficial_parametrica
        return unidades.Tension(parameters[0]*(1-self.tr(T))**parameters[1])

    def Tension_Hakim(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos.
        ref Chemcad pag 85
        ref Properties of gases and liquids pag 693 y sig."""
        Qp=0.1574+0.385*self.f_acent-1.769*self.stiehl-13.69*self.stiehl**2-0.510*self.f_acent**2+1.298*self.f_acent*self.stiehl
        m=1.21+0.5385*self.f_acent-14.61*self.stiehl-32.07*self.stiehl**2-1.656*self.f_acent**2+22.03*self.f_acent*self.stiehl
        return unidades.Tension(self.Pc.atm**0.67*self.Tc**0.33*Qp*((1-self.tr(T))/0.4)**m, "dyncm")

    def Tension_Block_Bird(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos.
        ref Eq.8.88 Riazi-Characterization of petroleum fractions pag 373"""
        Tbr=self.Tb/self.Tc
        Q=0.1196*(1+Tbr*log(self.Pc.atm)/(1-Tbr))-0.279
        return unidades.Tension(self.Pc.atm**0.67*self.Tc**0.33*Q*(1-self.tr(T))**(11./9), "dyncm")

    def Tension_Miqueu(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos.
        Miqueu, C., Broseta, D., Satherley, J., Mendiboure, B., Lachiase, J., and Graciaa, A., "An Extended Scaled Equation for the Temperature Dependence of the Surface Tension of Pure Compounds Inferred From an Analysis of Experimental Data," Fluid Phase Equilibria, Vol. 172, 2000, pp. 169-182."""
        t=1-self.tr(T)
        return unidades.Tension(Bolzmann*1e7*self.Tc*(Avogadro/self.Vc.ccg)**(2./3)*(4.35+4.14*self.f_acent)*(1+0.19*t**0.5-0.25*t)*t**1.26, "dyncm")

    def Tension_Hydrocarbon(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos"""
        t=unidades.Temperature(T)
        return unidades.Tension(673.7/self.Kw*((self.Tc.R-t.R)/self.Tc.R)**1.232, "dyncm")

    def Tension_MIller(self, T):
        """Método alternativo para el cálculo de la tensión superficial de líquidos"""
        Q=0.1207*(1+self.tr(self.Tb)*log(self.Pc.atm)/(1-self.tr(self.Tb)))-0.281
        return unidades.Tension(self.Pc.atm**0.67*self.Tc**0.33*Q*(1-self.tr(T))**(11.0/9), "dyncm")

    def Tension_Parachor(self, T, parachor):
        """Método alternativo para el cálculo de la tensión superficial de líquidos haciendo uso de las contribuciones de grupos
        API procedure 10A1.4, pag 987"""
        #TODO: Mientras no sepa como automatizar el cálculo de las contribuciones de grupo, habra que indicarlo como parámetro
        rhoL=self.RhoL_DIPPR(T)/1000
        rhoG=self.RhoG_Lee_Kesler(T, 1)*self.M/1000
        sigma=(parachor/self.M*(rhoL-rhoG))**4
        return unidades.Tension(sigma, "dyncm")





    def Hv_DIPPR(self,T):
        """Cálculo del calor de vaporización usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimosexta posición
        Calor de vaporización obtenido en (J/kmol)"""
        return unidades.Enthalpy(self.DIPPR(T,self.calor_vaporizacion)/self.M)

    def Hv_Lee_Kesler(self, T):
        """Método alternativo para el cálculo del calor de vaporización haciendo uso de las propiedades críticas
        Procedure API 7C1.16 Pag.680
        Valor en J/mol"""
        Pv=self.Pv_DIPPR(T)
        Tr=T/self.Tc
        Pr=Pv/self.Pc
        H_adimensional_vapor=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 1)
        H_adimensional_liquido=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, 0)
        return unidades.Enthalpy(R*self.Tc/self.M*(H_adimensional_vapor-H_adimensional_liquido), "Jg")

    def Cp_Solido_DIPPR(self,T):
        """Cálculo de la capacidad calorifica del solido usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimoseptima posición
        Capacidad calorifica obtenida en (J/kmol·K)"""
        return unidades.SpecificHeat(self.DIPPR(T,self.capacidad_calorifica_solido)/self.M)

    def Cp_Liquido_DIPPR(self,T):
        """Cálculo de la capacidad calorifica del liquido usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimoctava posición
        Capacidad calorifica obtenida en (J/kmol·K)"""
        return unidades.SpecificHeat(self.DIPPR(T,self.capacidad_calorifica_liquido)/self.M)

    def Cp_Hadden(self, T):
        """Método alternativo para el cálculo de la capacidad calorífica en líquidos por debajo de su punto de ebullición
        Procedure API 7D1.9 Pag.696"""
        #TODO: No fácil de implementar al depender del elemento

    def Cp_Gas_DIPPR(self,T):
        """Cálculo de la capacidad calorifica del gas ideal usando las ecuaciones DIPPR
        Los parámetros se encuentran en la base de datos en decimonovena posición
        Capacidad calorifica obtenida en (J/kmol·K)"""
        return unidades.SpecificHeat(self.DIPPR(T,self.capacidad_calorifica_gas)/self.M)

    def Cp_Lee_Kesler(self, T, P, fase=None):
        """Método alternativo para el cálculo de la capacidad calorífica
        Procedure API 7D3.6 Pag.711"""
        Tr=self.tr(T)
        if fase==None:
            fase=self.Fase(T, P)
        Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(T, P, fase)

        B=0.1181193-0.265728/Tr-0.154790/Tr**2-0.030323/Tr**3
        C=0.0236744-0.0186984/Tr
        D=0.155488e-4+0.623689e-4/Tr
        dpdt_0=1/vr0*(1+(0.1181193+0.154790/Tr**2+2*0.030323/Tr**3)/vr0+0.0236744/vr0**2+0.155488e-4/vr0**5-2*0.042724/Tr**3/vr0**2*((0.65392+0.060167/vr0**2)*exp(-0.060167/vr0**2)))
        dpdv_0=-Tr/vr0**2*(1+2*B/vr0+3*C/vr0**2+6*D/vr0**5+0.042724/Tr**3/vr0**2*(3*0.65392+(5-2*(0.65392+0.060167/vr0**2))*0.060167/vr0**2)*exp(-0.060167/vr0**2))
        Cp0=1+Tr*dpdt_0**2/dpdv_0+Cv0

        B=0.2026579-0.331511/Tr-0.027655/Tr**2-0.203488/Tr**3
        C=0.0313385-0.0503618/Tr+0.016901/Tr**3
        D=0.48736e-4+0.0740336e-4/Tr
        dpdt_h=1/vrh*(1+(0.2026579+0.027655/Tr**2+2*0.203488/Tr**3)/vrh+(0.0313385-2*0.016901/Tr**3)/vrh**2+0.48736e-4/vrh**5-2*0.041577/Tr**3/vrh**2*((1.226+0.03754/vrh**2)*exp(-0.03754/vrh**2)))
        dpdv_h=-Tr/vrh**2*(1+2*B/vrh+3*C/vrh**2+6*D/vrh**5+0.041577/Tr**3/vrh**2*(3*1.226+(5-2*(1.226+0.03754/vrh**2))*0.03754/vrh**2)*exp(-0.03754/vrh**2))
        Cph=1+Tr*dpdt_h**2/dpdv_h+Cvh

        Cp_adimensional=Cp0+self.f_acent/factor_acentrico_octano*(Cph-Cp0)
        return unidades.SpecificHeat(self.Cp_ideal(T).JgK-R/self.M*Cp_adimensional, "JgK")

    def Cv_Lee_Kesler(self, T, P, fase=None):
        """Método de cálculo de la capacidad calorífica a volumen constante
        Procedure API 7E1.6 Pag.726"""
        #FIXME: No sale, un factor de 100 tengo que añadir no sé de donde
        Pr=P/self.Pc
        Tr=T/self.Tc
        if fase==None:
            fase=self.Fase(T, P)
        Cpo=self.Cp_ideal(T)
        Cv0, Cvh, vr0, vrh=eos.Lee_Kesler_lib_Cp(Tr, Pr, fase)
        Cv_adimensional=Cv0+self.f_acent/factor_acentrico_octano*(Cvh-Cv0)
        return unidades.SpecificHeat(100*(Cpo.JgK-R/self.M*(1+Cv_adimensional)), "JgK")


    def Cp_Cv_Lee_Kesler(self, T, P):
        """Método de cálculo de la capacidad calorífica a volumen constante
        Procedure API 7E1.6 Pag.726"""
        Cv=self.Cv_Lee_Kesler(T, P.atm)
        Cp=self.Cp_Lee_Kesler(T, P.atm)
#        print Cp.BtulbF, Cv
        return Cp/Cv


    def Fugacidad_Lee_Kesler(self, T, P):
        """Método de cálculo de la fugacidad
        Procedure API 7G1.8 Pag.752"""
        Tr=T/self.Tc
        Pr=P/self.Pc
        f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P))
        return unidades.Pressure(P*exp(f), "atm")

    def Entropia_Lee_Kesler(self, T, P):
        """Método de cálculo de la entropia
        Procedure API 7F1.7 Pag.739"""
        Tr=T/self.Tc
        Pr=P/self.Pc
        S0=self.Entropia_ideal(T)
        H_adimensional=eos.Lee_Kesler_Entalpia_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        f=eos.Lee_Kesler_Fugacidad_lib(Tr, Pr, self.f_acent, self.Fase(T, P.atm))
        S=H_adimensional+f+log(P/101325)

        return unidades.SpecificHeat(S0.JgK-R*S/self.M, "JgK")

    def constante_Henry(self,T,parameters=None):
        """constante H obtenida en psia por unidad de fracción molar del gas
        lnH = A/T + B∗lnT + C∗T + D
        Solo disponible para algunos compuestos:
        Hydrogen, Helium, Argon, Neon, Krypton, Xenon, Oxygen, Nitrogen,
        Hydrogen sulfide, Carbon monoxide, Carbon dioxide, Sulfur dioxide,
        Nitrous oxide, Chlorine,Bromine, Iodine, Methane, Ethane, Propane,
        Ethylene, Ammonia.
        API procedure 9A7.1, pag 927
        Los parametros de la ecuación se encuentran en la base de datos
        en forma de lista en la posición décima
        """
        if parameters==None:
            parameters=self.henry
        t=unidades.Temperature(T)
        return exp(parameters[0]/t.R+parameters[1]*log(t.R)+parameters[2]*t.R+parameters[3])


    def Fase(self, T, P):
        """Método que calcula el estado en el que se encuentra la sustancia"""
        Pv=self.Pv(T).atm
        if Pv>P:
            return 1
        else:
            return 0







    def RhoG_Lee_Kesler(self, T, P):
        a, b=eos.SRK_lib(self, T)
        Z_srk=eos.Z_Cubic_EoS(T, P, b, a, b, 0, b)
        Vvo=Z_srk[0]*R_atml*T/P

        vr0v, vrhv, vr0l, vrhl=eos.Lee_Kesler_lib(T/self.Tc, P/self.Pc.atm, fase=1, Vvo=Vvo)
        z0v=P/self.Pc.atm*vr0v/T*self.Tc
        zhv=P/self.Pc.atm*vrhv/T*self.Tc
        z=z0v+self.f_acent/factor_acentrico_octano*(zhv-z0v)
        return P/z/R_atml/T


class newComponente(object):
    """Clase general que define la creaccion de nuevos componentes"""
    def export2Component(self):
        """Método que devuelve un array listo para añadirse a la base de datos"""
        elemento=[]
        elemento.append(self.formula)
        elemento.append(self.name)
        elemento.append(self.M)
        elemento.append(self.Tc)
        elemento.append(self.Pc.atm)
        elemento.append(self.Vc)
        elemento.append(self.API)
        elemento.append(self.cp)

        #Parametricas
        elemento.append([])
        elemento.append([])
        if self.mul:
            elemento.append(self.mul)
        else:
            elemento.append([])
        elemento.append([])

        #DIPPR
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])
        elemento.append([])

        #Otros
        elemento.append(0)
        elemento.append(self.Vliq)
        elemento.append(self.rackett)
        elemento.append(self.SG)
        elemento.append(self.f_acent)
        elemento.append(self.Parametro_solubilidad)
        elemento.append(self.watson)

        elemento.append([])

        elemento.append(0)
        elemento.append(self.Tb)
        elemento.append(self.Tf)
        elemento.append("")
        elemento.append("")

        elemento.append([])

        elemento.append(0)
        elemento.append(0)
        elemento.append(0)
        elemento.append(0)
        elemento.append(0)
        elemento.append(self.Hf)
        elemento.append(self.Gf)
        elemento.append(0)
        elemento.append(0)
        elemento.append(0)

        elemento.append("")
        elemento.append(0)
        elemento.append(0)
        elemento.append(0)

        elemento.append(0)
        elemento.append(0)
        elemento.append(0)
        elemento.append(0)

        elemento.append(0)
        elemento.append(0)
        elemento.append(0)
        elemento.append("")

        return elemento


class GroupContribution(newComponente):
    """Superclase que define un nuevo elemento según los diferentes métodos de contribución de grupos:
    -Joback: Preferido
    -Constantinou - Gani: Preferido
    -Wilson - Jasperson Contribuciones por átomos
    -Marrero - Pardillo: No predice propiedades termodinámicas
    -Elliot (UNIFAC) Preferido
    -Andrade: Hidrocarbonos sin heteroátomos"""

    kwargs={"group": [],
                    "contribution": [],
                    "M": 0.0,
                    "Tb": 0.0,
                    "SG": 0.0,
                    "name": "",

                    "ring": 0,
                    "atomos": 0,
                    "platt": 0}

    FirstOrder=0
    SecondOrder=0
    status=0
    _bool=False
    msg=""
    cp=[]
    Tb=0
    Tf=0
    Hf=0
    Gf=0
    mul=None
    Hm=0

    def __init__(self, **kwargs):
        """Lee la documentación de cada tipo para conocer los kwargs aceptados"""
        self.kwargs=self.kwargs.copy()
        if kwargs:
            self.__call__(**kwargs)

    def __call__(self, **kwargs):
        self.kwargs.update(kwargs)
        self._bool=True
        if self.isCalculable():
            self.calculo()

    def isCalculable(self):
        """Método que estima si el método es calculable en función de los datos disponibles, definido por cada método"""
        if not self.kwargs["group"] or not self.kwargs["contribution"]:
            self.msg=QApplication.translate("pychemqt", "undefined group")
            self.status=0
        else:
            self.status=1
            self.msg=""
            return True

    def calculo(self):
        """Procedimiento de cálculo propiamente dicho, definido por cada método"""
        if self.kwargs["name"]:
            self.name=str(self.kwargs["name"])
        else:
            self.name=self.__class__.__name__+"_"+time.strftime("%d/%m/%Y-%H:%M:%S")

        if "f_acent" not in self.__dict__:
            self.f_acent=self._factor_acentrico_Lee_Kesler()
        if "Hv" not in self.__dict__:
            self.Hv=self._Calor_vaporizacion()

        self.rackett=self._Rackett()
        self.Vliq=self._Volumen_Liquido_Constante()
        self.Parametro_solubilidad=self._Parametro_solubilidad()
        if self.kwargs["SG"]:
            self.SG=self.kwargs["SG"]
        else:
            self.SG=self._Gravedad_especifica()
        if "Vc" not in self.__dict__:
            self.Vc=self._Vc()
        self.watson=self._Watson()
        if "cp" not in self.__dict__:
            self.cp=self._cp()
        self.API=141.5/self.SG-131.5
        self.txt, self.formula=self.EmpiricFormula()

    def __bool__(self):
        return self._bool

    def clear(self):
        self.kwargs=self.__class__.kwargs
        self.__dict__.clear()
        self._bool=False


    def _Gravedad_especifica(self):
        #FIXME: No sale bién
        volumen=self.Vliq*(5.7+3*288.71/self.Tc)
        return 1/volumen*18

    def _Watson(self):
        return self.Tb.R**(1./3)/self.SG

    def _factor_acentrico_Lee_Kesler(self):
        Tr=self.Tb/self.Tc
        Pr=101325./self.Pc
        return (log(Pr)-5.92714+6.09648/Tr+1.28862*log(Tr)-0.169347*Tr**6)/(15.2518-15.6875/Tr-13.4721*log(Tr)+0.43577*Tr**6)

    def _cp(self):
        """Método de cálculo de la capacidad calorifica del gas ideal"""
        cp=[0, 0, 0, 0, 0, 0]
        cp[0]=unidades.Enthalpy((0.036863384*self.watson-0.4673722)*self.M, "Btulb").kcalkg/self.M*1.8
        cp[1]=unidades.Enthalpy((3.1865e-5*self.watson+0.001045186)*self.M, "Btulb").kcalkg/self.M*1.8**2
        cp[2]=unidades.Enthalpy(-4.9572e-7*self.M, "Btulb").kcalkg/self.M*1.8**3
        return cp


    def _Vc(self):
        """Método de cálculo del volumen crítico"""
        if self.Tc.R<536.67:
            D=8.75+1.987*log(self.Tb.R)+self.Tb.R/1.8
        elif self.Tc.R>593:
            #FIXME: No sale bién
            D=(0.398907*self.SG*(self.f_acent-592.4439)/self.M)**0.5
        else:
            f=((self.Tc.R-536.67)/(self.Tc.R-self.Tb.R))**0.38
            D=8.75+1.987*log(self.Tb.R)+self.Tb.R/1.8*f

        Zc=1/(3.43+6.7e-9*D**2)
        return Zc*self.Tc.R*10.73/self.Pc.psi

    def _Calor_vaporizacion(self):
        """Método de cálculo del calor de vaporización,
        ref. chemcad pag 60"""
        tbr=self.Tb/self.Tc
        return unidades.Enthalpy(1.093*R*1000*self.Tc*(tbr*(log(self.Pc.atm)-1)/(0.930-tbr))/self.M, "calg")


    def _Rackett(self):
        """ref 64"""
        return 0.29056-0.08775*self.f_acent

    def _Volumen_Liquido_Constante(self):
        V=R_atml*1000*self.Tc/self.Pc.atm*self.rackett**(1+(1-298.15/self.Tc)**(2.0/7)) #cm3/mol
        return V/(5.7+1611/self.Tc) #cm3/mol


    def _Parametro_solubilidad(self):
        V=R_atml/1000*self.Tc/self.Pc.atm*self.rackett**(1+(1-298.15/self.Tc)**(2.0/7)) #m3/mol
        return unidades.SolubilityParameter(((self.Hv-298*R)/V)**0.5)

    def EmpiricFormula(self):
        if "txt" not in self.coeff:
            return "", ""
        C=H=N=O=S=F=Cl=Br=I=0
        for i, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if i<self.FirstOrder:
                C+=self.coeff["txt"][i][1].get("C", 0)*contribucion
                H+=self.coeff["txt"][i][1].get("H", 0)*contribucion
                N+=self.coeff["txt"][i][1].get("N", 0)*contribucion
                O+=self.coeff["txt"][i][1].get("O", 0)*contribucion
                S+=self.coeff["txt"][i][1].get("S", 0)*contribucion
                F+=self.coeff["txt"][i][1].get("F", 0)*contribucion
                Cl+=self.coeff["txt"][i][1].get("Cl", 0)*contribucion
                Br+=self.coeff["txt"][i][1].get("Br", 0)*contribucion
                I+=self.coeff["txt"][i][1].get("I", 0)*contribucion
        string=""
        formula=""
        for elemento, txt in zip([C, H, N, O, S, F, Cl, Br, I], ["C", "H", "N", "O", "S", "F",  "Cl", "Br", "I"]):
            if elemento>1:
                string+="%s<sub>%i</sub>" % (txt, elemento)
                formula+="%s%i" % (txt, elemento)
            elif elemento==1:
                string+="%s" %txt
                formula+="%s" %txt

        return string, formula


class Joback(GroupContribution):
    """Joback, K. G.: ‘‘A Uniﬁed Approach to Physical Property Estimation Using Multivariate Statistical Techniques,’’ S.M. Thesis, Department of Chemical Engineering, Massachusetts Institute of Technology, Cambridge, MA, 1984.
    ref, properties of gases and liquids, pag 782
    grupos: grupos que forman la molécula
    contribuciones: contribuciones de cada grupo
    M: peso molecular, opcional
    Tb: Temperatura de ebullición, opcional
    SG: gravedad específica, opcional

    >>> desconocido=Joback(group=[0, 1, 13, 14, 20], contribution=[1, 1, 4, 2, 1])
    >>> print desconocido.Tb, desconocido.Tc
    489.74 715.745692022

    #http://en.wikipedia.org/wiki/Joback_method
    >>> acetona=Componente(140)
    >>> print acetona.Tc, acetona.Pc.bar, acetona.Tb, acetona.Tf
    508.2 47.0150127825 329.4 178.45
    >>> joback_acetona=Joback(group=[0, 23], contribution=[2, 1])
    >>> print joback_acetona.Tc, joback_acetona.Pc.bar, joback_acetona.Tb, joback_acetona.Tf
    500.248204912 48.0249960499 321.91 173.5
    """
    coeff={
        "M": [15.03452, 14.02658, 13.01864, 12.0107, 14.02658, 13.01864, 12.0107, 12.0107, 13.01864, 12.0107, 14.02658, 13.01864, 12.0107, 13.01864, 12.0107, 18.9984, 35.453, 79.904, 126.90447, 17.00734, 17.00734, 15.9994, 15.9994, 28.0101, 28.0101, 29.01804, 45.01744, 44.0095, 15.9994, 16.02258, 15.01464, 15.01464, 14.0067, 14.0067, 14.0067, 15.01464, 26.0174, 46.0055, 33.07294, 32.065, 32.065],
        "atomos": [4, 3, 2, 1, 3, 2, 1, 1, 2, 1, 3, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 3, 4, 3, 1, 3, 2, 2, 1, 1, 1, 2, 2, 3, 2, 1, 1],
        "tc": [0.0141, 0.0189, 0.0164, 0.0067, 0.0113, 0.0129, 0.0117, 0.0026, 0.0027, 0.002, 0.01, 0.0122, 0.0042, 0.0082, 0.0143, 0.0111, 0.0105, 0.0133, 0.0068, 0.0741, 0.024, 0.0168, 0.0098, 0.038, 0.0284, 0.0379, 0.0791, 0.0481, 0.0143, 0.0243, 0.0295, 0.0130, 0.0169, 0.0255, 0.0085, 0.0, 0.0496, 0.0437, 0.0031, 0.0119, 0.0019],
        "Pc": [-0.0012, 0.0, 0.002, 0.0043, -0.0028, -0.0006, 0.0011, 0.0028, -0.0008, 0.0016, 0.0025, 0.0004, 0.0061, 0.0011, 0.0008, -0.0057, -0.0049, 0.0057, -0.0034, 0.0112, 0.0184, 0.0015, 0.0048, 0.0031, 0.0028, 0.0030, 0.0077, 0.0005, 0.0101, 0.0109, 0.0077, 0.0114, 0.0074, -0.0099, 0.0076, 0.0, -0.0101, 0.0064, 0.0084, 0.0049, 0.0051],
        "vc": [65, 56, 41, 27, 56, 46, 38, 36, 46, 37, 48, 38, 27, 41, 32, 27, 58, 71, 97, 28, -25, 18, 13, 62, 55, 82, 89, 82, 36, 38, 35, 29, 9, 0, 34, 0, 91, 91, 63, 54, 38],
        "tb": [23.58, 22.88, 21.74, 18.25, 18.18, 24.96, 24.14, 26.15, 9.20, 27.38, 27.15, 21.78, 21.32, 26.73, 31.01, -0.03, 38.13, 66.86, 93.84, 92.88, 76.34, 22.42, 31.22, 76.75, 94.97, 72.24, 169.09, 81.10, -10.5, 73.23, 50.17, 52.82, 11.74, 74.6, 57.55, 0.0, 125.66, 152.54, 63.56, 68.78, 52.10],
        "tf": [-5.10, 11.27, 12.64, 46.43, -4.32, 8.73, 11.14, 17.78, -11.18, 64.32, 7.75, 19.88, 60.15, 8.13, 37.02, -15.78, 13.55, 43.43, 41.69, 44.45, 82.83, 22.23, 23.05, 61.20, 75.97, 36.9, 155.5, 53.6, 2.08, 66.89, 52.66, 101.51, 48.84, 0, 68.4, 0.0, 59.89, 127.24, 20.09, 34.4, 79.93],
        "hf": [-76.45, -20.64, 29.89, 82.23, -9.63, 37.97, 83.99, 142.14, 79.30, 115.51, -26.8, 8.67, 79.72, 2.09, 46.43, -251.92, -71.55, -29.48, 21.06, -208.04, -221.65, -132.22, -138.16, -133.22, -164.50, -162.03, -426.72, -337.92, -247.61, -22.02, 53.47, 31.65, 123.34, 23.61, 55.52, 93.7, 88.43, -66.57, -17.33, 41.87, 39.1],
        "gf": [-43.96, 8.42, 58.36, 116.02, 3.77, 48.53, 92.36, 136.7, 77.71, 109.82, -3.68, 40.99, 87.88, 11.30, 54.05, -247.19, -64.31, -38.06, 5.74, -189.20, -197.37, -105.0, -98.22, -120.50, -126.27, -143.48, -387.87, -301.95, -250.83, 14.07, 89.39, 75.61, 163.16, 0.0, 79.93, 119.66, 89.22, -16.83, -22.99, 33.12, 27.73],
        "hv": [567, 532, 404, 152, 412, 527, 511, 636, 276, 789, 573, 464, 154, 608, 731, -160, 1083, 1573, 2275, 4021, 2987, 576, 1119, 2144, 1588, 2173, 4669, 2302, 1412, 2578, 1538, 1656, 453, 797, 1560, 2908, 3071, 4000, 1645, 1629, 1430],
        "hm": [217, 619, 179, -349, -113, 643, 732, 1128, 555, 992, 117, 775, -328, 263, 572, 334, 601, 861, 651, 575, 1073, 284, 1405, 1001, 0, 764, 2641, 1663, 866, 840, 1197, 1790, 1124, 0, 872, 0, 577, 2313, 564, 987, 372],
        "cpa": [19.5, -0.909, -23.0, -66.2, -23.6, -8.0, -28.1, 27.4, 24.5, 7.87, -6.03, 8.67, -90.9, -2.14, -8.25, 26.5, 33.3, 28.6, 32.1, 25.7, -2.81, 25.5, 12.2, 6.45, 30.4, 30.9, 24.1, 24.5, 6.82, 26.9, -1.21, 11.8, -31.1, 0.0, 8.83, 5.69, 36.5, 25.9, 35.3, 19.6, 16.7],
        "cpb": [-8.08e-3, 9.5e-2, 2.04e-1, 4.27e-1, -3.81e-2, 1.05e-1, 2.08e-1, -5.57e-2, -2.71e-2, 2.01e-2, 8.54e-2, 1.62e-1, 5.57e-1, 5.74e-2, 1.01e-1, -9.13e-2, -9.63e-2, -6.49e-2, -6.41e-2, -6.91e-2, 1.11e-1, -6.32e-2, -1.26e-2, 6.7e-2, -8.29e-2, -3.36e-2, 4.27e-2, 4.02e-2, 1.96e-2, -4.12e-2, 7.62e-2, -2.3e-2, 2.27e-1, 0.0, -3.84e-3, -4.12e-3, -7.33e-2, -3.74e-3, -7.58e-2, -5.61e-3, 4.81e-3],
        "cpc": [1.53e-4, -5.44e-5, -2.65e-4, -6.41e-4, 1.72e-4, -9.63e-5, -3.06e-4, 1.01e-4, 1.11e-4, -8.33e-6, -8.0e-6, -1.6e-4, -9.0e-4, -1.64e-6, -1.42e-4, 1.91e-4, 1.87e-4, 1.36e-4, 1.26e-4, 1.77e-4, -1.16e-4, 1.11e-4, 6.03e-5, -3.57e-5, 2.36e-4, 1.6e-4, 8.04e-5, 4.02e-5, 1.27e-5, 1.64e-4, -4.86e-5, 1.07e-4, -3.2e-4, 0.0, 4.35e-5, 1.28e-4, 1.84e-4, 1.29e-4, 1.85e-4, 4.02e-5, 2.77e-5],
        "cpd": [-9.67e-8, 1.19e-8, 1.2e-7, 3.01e-7, -1.03e-7, 3.56e-8, 1.46e-7, -5.02e-8, -6.78e-8, 1.39e-9, -1.8e-8, 6.24e-8, 4.69e-7, -1.59e-8, 6.78e-8, -1.03e-7, -9.96e-8, -7.45e-8, -6.87e-8, -9.88e-8, 4.94e-8, -5.48e-8, -3.86e-8, 2.86e-9, -1.31e-7, -9.88e-8, -6.87e-8, -4.52e-8, -1.78e-8, -9.76e-8, 1.05e-8, -6.28e-8, 1.46e-7, 0.0, -2.6e-8, -8.88e-8, -1.03e-7, -8.88e-8, -1.03e-7, -2.76e-8, -2.11e-8],
        "mua": [548.29, 94.16, -322.15, -573.56, 495.01, 82.28, 0, 0, 0, 0, 307.53, -394.29, 0, 259.65,-245.74, 0, 625.45, 738.91, 809.55, 2173.72, 3018.17, 122.09, 440.24, 340.35, 0, 740.92, 1317.23, 483.88, 675.24,0,0,0,0,0,0,0,0,0,0,0,0],
        "mub": [-1.719,-0.199,1.187,2.307,-1.539,-0.242,0,0,0,0,-0.798,1.251,0,-0.702,0.912,0,-1.814,-2.038,-2.224,-5.057,-7.314,-0.386,-0.953,-0.35,0,-1.713,-2.578,-0.966,-1.34,0,0,0,0,0,0,0,0,0,0,0,0],
        "txt": [("CH3", {"C": 1, "H": 3}),
                        ("CH2", {"C": 1, "H": 2}),
                        ("CH", {"C": 1, "H": 1}),
                        ("C", {"C": 1}),
                        ("=CH2", {"C": 1, "H": 2}),
                        ("=CH", {"C": 1, "H": 1}),
                        ("=C", {"C": 1}),
                        ("=C=", {"C": 1}),
                        ("≡CH", {"C": 1, "H": 1}),
                        ("≡C", {"C": 1}),
                        ("CH2 (cyclic)", {"C": 1, "H": 2}),
                        ("CH (cyclic)", {"C": 1, "H": 1}),
                        ("C (cyclic)", {"C": 1,}),
                        ("-CH (Aromatic)", {"C": 1, "H": 1}),
                        ("=C (Aromatic)", {"C": 1}),
                        ("F", {"F": 1}),
                        ("Cl", {"Cl": 1}),
                        ("Br", {"Br": 1}),
                        ("I", {"I": 1}),
                        ("-OH", {"O": 1, "H": 1}),
                        ("-OH (Aromatic)", {"O": 1, "H": 1}),
                        ("-O-", {"O": 1}),
                        ("-O- (cyclic)", {"O": 1}),
                        ("C=O", {"C": 1, "O": 1}),
                        ("C=O (cyclic)", {"C": 1, "O": 1}),
                        ("CH=O", {"C": 1, "H": 1, "O": 1}),
                        ("COOH", {"C": 1, "H": 1, "O": 2}),
                        ("COO", {"C": 1, "O": 2}),
                        ("=O", {"O": 1}),
                        ("NH2", {"N": 1, "H": 2}),
                        ("NH", {"N": 1, "H": 1}),
                        ("NH (cyclic)", {"N": 1, "H": 1}),
                        ("N", {"N": 1}),
                        ("=N-", {"N": 1}),
                        ("=N- (cyclic)", {"N": 1}),
                        ("=NH", {"N": 1, "H": 1}),
                        ("CN", {"C": 1, "N": 1}),
                        ("NO2", {"N": 1, "O": 2}),
                        ("SH", {"S": 1, "H": 1}),
                        ("S", {"S": 1}),
                        ("S (cyclic)",  {"S": 1})]}

    FirstOrder=41

    def calculo(self):
        if self.kwargs["M"]:
            self.M=self.kwargs["M"]
        else:
            self.M=sum([self.coeff["M"][grupo]*self.kwargs["contribution"][i] for i, grupo in enumerate(self.kwargs["group"])])

        if self.kwargs["Tb"] :
            self.Tb=unidades.Temperature(self.kwargs["Tb"])
        else:
            self.Tb=unidades.Temperature(198+sum([self.coeff["tb"][grupo]*self.kwargs["contribution"][i] for i, grupo in enumerate(self.kwargs["group"])]))

        atomos=tcsuma=pcsuma=vcsuma=0
        Tf=122.5
        Hf=68.29
        Gf=53.88
        Hv=15.3
        Hm=-0.88
        cpa=-37.93
        cpb=0.21
        cpc=-3.91e-4
        cpd=2.06e-7
        mua, mub=0, 0
        for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            Tf+=contribucion*self.coeff["tf"][grupo]
            atomos+=contribucion*self.coeff["atomos"][grupo]
            tcsuma+=contribucion*self.coeff["tc"][grupo]
            pcsuma+=contribucion*self.coeff["Pc"][grupo]
            vcsuma+=contribucion*self.coeff["vc"][grupo]
            Hf+=contribucion*self.coeff["hf"][grupo]
            Gf+=contribucion*self.coeff["gf"][grupo]
            Hv+=contribucion*self.coeff["hv"][grupo]
            Hm+=contribucion*self.coeff["hm"][grupo]
            cpa+=contribucion*self.coeff["cpa"][grupo]
            cpb+=contribucion*self.coeff["cpb"][grupo]
            cpc+=contribucion*self.coeff["cpc"][grupo]
            cpd+=contribucion*self.coeff["cpd"][grupo]
            mua+=contribucion*self.coeff["mua"][grupo]
            mub+=contribucion*self.coeff["mub"][grupo]
        self.Tf=unidades.Temperature(Tf)
        self.Tc=unidades.Temperature(self.Tb/(0.584+0.965*tcsuma-tcsuma**2))
        self.Pc=unidades.Pressure((0.113+0.0032*atomos-pcsuma)**-2, "bar")
        self.Vc=unidades.SpecificVolume((vcsuma+17.5)/1000/self.M)
        self.Hf=unidades.Enthalpy(Hf/self.M, "kJg")
        self.Gf=unidades.Enthalpy(Gf/self.M, "calg")
        self.Hv=unidades.Enthalpy(Hv/self.M, "calg")
        self.Hm=unidades.Enthalpy(Hm/self.M, "calg")
        self.cp=[cpa, cpb, cpc, cpd, 0, 0]
        #Ajuste de la expresion de la viscosidad a la viscosidad parametrica
        self.mul=[mua-597.82, -(mua-597.82)/(log(self.M)+mub-11.202)]

        GroupContribution.calculo(self)


class Constantinou_Gani(GroupContribution):
    """Constantinou, L., R. Gani, and J. P. O’Connell: Fluid Phase Equil., 104: 11 (1995).
    ref, properties of gases and liquids, pag 782
    grupos: grupos que forman la molécula
    contribuciones: contribuciones de cada grupo
    M: peso molecular, opcional
    Tb: Temperatura de ebullición, opcional
    SG: gravedad específica, opcional

    >>> unknown=Constantinou_Gani(group=[0, 1, 15], contribution=[1, 3, 1])
    >>> print unknown.f_acent, unknown.Tc
    0.380647791898 558.910978427
    """
    coeff={
            "M": [15.03452, 14.02658, 13.01864, 12.0107, 27.04522, 26.03728, 26.03728, 25.02934, 24.0214, 39.05592, 13.01864, 12.0107, 27.04522, 26.03728, 25.02934, 17.00734, 29.01804, 43.04462, 42.03668, 29.01804, 59.04402, 58.03608, 45.01744, 31.03392, 30.02598, 29.01804, 49.0243832, 30.04916, 29.04122, 28.03328, 29.04122, 28.03328,  28.03328, 78.09196, 77.08402, 40.04398, 45.01744, 49.47958, 48.47164, 47.4637, 83.92464, 118.3697, 82.9167, 47.4637, 60.03208, 59.02414, 58.0162, 47.09952, 126.90447, 79.904, 25.02934, 24.0214, 59.4744, 31.0091032, 71.0779, 69.0059096, 50.0075064, 31.0091032, 44.0095, 101.9151032, 67.4700432, 85.4605064, 18.9984032, 44.03268, 58.05926, 57.05132, 72.08584, 70.06996, 70.06996, 64.08372, 63.07578, 47.09952, 46.09158, 45.08364, 83.13162, 82.12368],
            "tc": [1.6781, 3.4920, 4.0330, 4.8823, 5.0146, 7.3691, 6.5081, 8.9582, 11.3764, 9.9318, 3.7337, 14.6409, 8.2130, 10.3239, 10.4664, 9.7292, 25.9145, 13.2896, 14.6273, 10.1986, 12.5965, 13.8116, 11.6057, 6.4737, 6.0723, 5.0663, 9.5059, 12.1726, 10.2075, 9.8544, 10.4677, 7.2121, 7.6924, 5.5172, 28.7570, 29.1528, 27.9464, 20.3781, 23.7593, 11.0752, 10.8632, 11.3959, 16.3945, 0, 18.5875, 14.1565, 24.7369, 23.2050, 34.5870, 13.8058, 17.3947, 10.5371, 7.5433, 11.4501, 5.4334, 2.8977, 0, 2.4778, 1.7399, 3.5192, 12.1084, 9.8408, 0, 4.8923, 1.5974, 65.1053, 0, 0, 36.1403, 0, 0, 17.9668, 0, 14.3969, 17.7916, 0, 0, 0, -0.5334, -0.5143, 1.0699, 1.9886, 5.8254, -2.3305, -1.2978, -0.6785, 0.8479, 3.6714, 0.4402, 0.0167, -0.5231, -0.3850, 2.1160, 2.0427, -1.5826, 0.2996, 0.5018, 2.9571, 1.1696, -1.7493, 6.1279, -1.3406, 2.5413, -2.7617, -3.4235, -2.8035, -3.5442, 5.4941, 0.3233, 5.4864, 2.0699, 2.1345, 1.0159, -5.3307, 4.4847, -0.4996, -1.9334, 0, -2.2974, 2.8907, 0],
            "Pc": [0.0199, 0.0106, 0.0013, -0.0104, 0.0250, 0.0179, 0.0223, 0.0126, 0.0020, 0.0313, 0.0075, 0.0021, 0.0194, 0.0122, 0.0028, 0.0051, -0.0074, 0.0251, 0.0178, 0.0141, 0.0290, 0.0218, 0.0138, 0.0204, 0.0151, 0.0099, 0.0090, 0.0126, 0.0107, 0.0126, 0.0104, -0.0005, 0.0159, 0.0049, 0.0011, 0.0296, 0.0257, 0.0361, 0.0115, 0.0198, 0.0114, 0.0031, 0.0268, 0, 0.0349, 0.0131, 0.0210, 0.0122, 0.0150, 0.0136, 0.0028, -0.0018, 0.0148, 0.0041, 0.0160, 0.0130, 0, 0.0442, 0.0129, 0.0047, 0.0113, 0.0354, 0, 0.0390, 0.0144, 0.0043, 0, 0, 0.0401, 0, 0, 0.0254, 0, 0.0160, 0.0111, 0, 0, 0, 0.000488, 0.001410, -0.001850, -0.005200, -0.013230, 0.003714, 0.001171, 0.000424, 0.002257, -0.009800, 0.004186, -0.000180, 0.003538, 0.005675, -0.002550, 0.005175, 0.003659, 0.001474, -0.002300, 0.003818, -0.002480, 0.004920, 0.000344, 0.000659, 0.001067, -0.004880, -0.000540, -0.004390, 0.000178, 0.005052, 0.006917, 0.001408, 0.002148, -0.005950, -0.000880, -0.002250, 0, 0.000319, 0, 0, 0.009027, 0.008247, 0],
            "vc": [0.0750, 0.0558, 0.0315, -0.0003, 0.1165, 0.0954, 0.0918, 0.0733, 0.0762, 0.1483, 0.0422, 0.0398, 0.1036, 0.1010, 0.0712, 0.0390, 0.0316, 0.1340, 0.1119, 0.0863, 0.1589, 0.1365, 0.1056, 0.0875, 0.0729, 0.0587, 0.0686, 0.1313, 0.0753, 0.1215, 0.0996, 0.0916, 0.1260, 0.0670, 0.0636, 0.2483, 0.1703, 0.1583, 0.1019, 0.1156, 0.1035, 0.0792, 0.1695, 0, 0.2103, 0.1016, 0.1653, 0.1423, 0.1426, 0.1025, 0.1081, 0.0828, 0.0933, 0.0763, 0.0569, 0.0567, 0, 0.1148, 0.0952, 0, 0.0859, 0.1821, 0, 0.1475, 0.0378, 0.1443, 0, 0, 0.2503, 0, 0, 0.1675, 0, 0.1302, 0.1165, 0, 0, 0, 0.00400, 0.00572, -0.00398, -0.01081, -0.02300, -0.00014, -0.00851, -0.00866, 0.01636, -0.02700, -0.00781, -0.00098, 0.00281, 0.00826, -0.01755, 0.00227, -0.00664, -0.00510, -0.00122, -0.01966, 0.00664, 0.00559, -0.00415, -0.00293, -0.00591, -0.00144, 0.02605, -0.00777, 0.01511, 0.00397, -0.02297, 0.00433, 0.00580, -0.01380, 0.00297, -0.00045, 0, -0.00596, 0.00510, 0, -0.00832, -0.00341, 0],
            "w": [0.296, 0.147, -0.071, -0.351, 0.408, 0.252, 0.223, 0.235, -0.210, 0.152, 0.027, 0.334, 0.146, -0.088, 1.524, 0.737, 1.015, 0.633, 0.963, 1.133, 0.756, 0.765, 0.526, 0.442, 0.218, 0.509, 0.800, 0, 0.953, 0.550, 0.386, 0.384, 0.075, 0.793, 0, 0, 0, 1.670, 0.570, 0, 0, 0.716, 0, 0.617, 0, 0.296, 0, 0, 0, 0, 0.233, 0.278, 0.618, 0, 0, 0.263, 0.500, 0, 0, 0, 0, 0.503, 0, 0.547, 0, 0, 0, 0, 0, 0, 0, 0.428, 0, 0, 0.438, 0.739, 0, 0, 0.01740, 0.01922, -0.00475, -0.02883, -0.08632, 0.17563, 0.22216, 0.16284, -0.03065, -0.02094, 0.01648, 0.00619, -0.01150, 0.02778, -0.11024, -0.11240, 0, -0.20789, -0.16571, 0, 0, 0.08774, 0, -0.26623, 0, 0.91939, 0, 0.03654, 0.21106, 0, 0, 0, 0, -0.13106, 0, 0, -0.01509, 0, 0, 0, -0.03078, 0.00001, 0],
            "tb": [0.8894, 0.9225, 0.6033, 0.2878, 1.7827, 1.8433, 1.7117, 1.7957, 1.8881, 3.1243, 0.9297, 1.6254, 1.9669, 1.9478, 1.7444, 3.2152, 4.4014, 3.5668, 3.8967, 2.8526, 3.6360, 3.3953, 3.1459, 2.2536, 1.6249, 1.1557, 2.5892, 3.1656, 2.5983, 3.1376, 2.6127, 1.5780, 2.1647, 1.2171, 5.4736, 6.2800, 5.9234, 5.0525, 5.8337, 2.9637, 2.6948, 2.2073, 3.9300, 3.5600, 4.5797, 2.6293, 5.7619, 5.0767, 6.0837, 3.2914, 3.6650, 2.6495, 2.3678, 2.5645, 1.7824, 0.9442, 7.2644, 1.2880, 0.6115, 1.1739, 2.6446, 2.8881, 2.3086, 1.9163, 1.0081, 10.3428, 0, 0, 7.6904, 0, 6.7822, 5.5566, 5.4248, 3.6796, 3.6763, 2.6812, 5.7093, 5.8260, -0.1157, -0.0489, 0.1798, 0.3189, 0.7273, 0.4745, 0.3563, 0.1919, 0.1957, 0.3489, 0.1589, 0.0668, -0.1406, -0.0900, 0.0511, 0.6884, -0.1074, 0.0224, 0.0920, 0.5580, 0.0735, -0.1552, 0.7801, -0.2383, 0.4456, -0.1977, 0.0835, -0.5385, -0.6331, 1.4108, -0.0690, 1.0682, 0.4247, 0.2499, 0.1134, -0.2596, 0.4408, -0.1168, -0.3201, -0.4453, -0.6776, -0.3678, 0],
            "tf": [0.4640, 0.9246, 0.3557, 1.6479, 1.6472, 1.6322, 1.7899, 2.0018, 5.1175, 3.3439, 1.4669, 0.2098, 1.8635, 0.4177, -1.7567, 3.5979, 13.7349, 4.8776, 5.6622, 4.2927, 4.0823, 3.5572, 4.2250, 2.9248, 2.0695, 4.0352, 4.5047, 6.7684, 4.1187, 4.5341, 6.0609, 3.4100, 4.0580, 0.9544, 10.1031, 0, 12.6275, 4.1859, 11.5630, 3.3376, 2.9933, 9.8409, 5.1638, 0, 10.2337, 2.7336, 5.5424, 4.9738, 8.4724, 3.0044, 4.6089, 3.7442, 3.9106, 9.5793, 1.5598, 2.5015, 0, 3.2411, 0, 0, 3.4448, 7.4756, 0, 2.7523, 1.9623, 31.2786, 0, 0, 11.3770, 0, 0, 0, 0, 5.0506, 3.1468, 0, 0, 0, 0.0381, -0.2355, 0.4401, -0.4923, 6.0650, 1.3772, 0, 0.6824, 1.5656, 6.9709, 1.9913, 0.2476, -0.5870, -0.2361, -2.8298, 1.4880, 2.0547, -0.2951, -0.2986, 0.7143, -0.6697, -3.1034, 28.4324, 0.4838, 0.0127, -2.3598, -2.0198, -0.5480, 0.3189, 0.9124, 9.5209, 2.7826, 2.5114, 1.0729, 0.2476, 0.1175, -0.2914, -0.0514, -1.6425, 0, 2.5832, -1.5511, 0],
            "hf": [-45.947, -20.763, -3.766, 17.119, 53.712, 69.939, 64.145, 82.528, 104.293, 197.322, 11.189, 27.016, -19.243, 9.404, 27.671, -181.422, -164.609, -182.329, -164.410, -129.2, -389.737, -359.258, -332.822, -163.569, -151.143, -129.488, -140.313, -15.505, 3.320, 5.432, 23.101, 26.718, 54.929, 69.885, 20.079, 134.062, 139.758, 88.298, -396.242, -73.568, -63.795, -57.795, -82.921, 0, -107.188, -16.752, -66.138, -59.142, -7.365, -8.253, 57.546, 1.834, 220.803, 227.368, -36.097, -161.740, 0, -679.195, 0, 0, -313.545, -258.960, 0, -446.835, -223.398, -203.188, -67.778, -182.096, -189.888, -46.562, 0, -344.125, 0, -2.084, 18.022, 0, 0, 0, -0.860, -1.338, 6.771, 7.205, 14.271, 104.800, 99.455, 13.782, -9.660, 15.465, -8.392, 0.474, 1.472, 4.504, 1.252, -2.792, -2.092, 0.975, 4.753, 14.145, -3.173, 1.279, 12.245, -7.807, 37.462, -16.097, -9.874, -3.887, -24.125, 0.366, -16.333, -2.992, 2.855, 0.351, -8.644, 1.532, -0.329, 0, 11.989, 0, 12.285, 11.207, 11.740],
            "gf": [-8.030, 8.231, 19.848, 37.977, 84.926, 92.900, 88.402, 93.745, 116.613, 221.308, 22.533, 30.485, 22.505, 41.228, 52.948, -158.589, -132.097, -131.366, -132.386, -107.858, -318.616, -291.188, -288.902, -105.767, -101.563, -92.099, -90.883, 58.085, 63.051, 82.471, 95.888, 85.001, 128.602, 132.756, 68.861, 199.958, 199.288, 121.544, -349.439, -33.373, -31.502, -25.261, -35.814, 0, -53.332, -0.596, 17.963, 18.088, 60.161, 16.731, 46.945, -1.721, 217.003, 216.328, -28.148, -144.549, 0, -626.580, 0, 0, -281.495, -209.337, 0, -392.975, -212.718, -136.742, 0, 0, -65.642, 0, 0, -241.373, 0, 30.222, 38.346, 0, 0, 0, 0.297, -0.399, 6.342, 7.466, 16.224, 94.564, 92.573, 5.733, -8.180, 20.597, -5.505, 0.950, 0.699, 1.013, 1.041, -1.062, -1.359, 0.075, 0, 23.539, -2.602, 2.149, 10.715, -6.208, 29.181, -11.809, -7.415, -6.770, -20.770, 3.805, -5.487, -1.600, 1.858, 8.846, -13.167, -0.654, -2.091, 0, 12,373, 0, 14.161, 12.530, 0],
            "hv": [4.116, 4.650, 2.771, 1.284, 6.714, 7.370, 6.797, 8.178, 9.342, 12.318, 4.098, 12.552, 9.776, 10.185, 8.834, 24.529, 40.246, 18.999, 20.041, 12.909, 22.709, 17.759, 0, 10.919, 7.478, 5.708, 11.227, 14.599, 11.876, 14.452, 14.481, 0, 6.947, 6.918, 28.453, 31.523, 31.005, 23.340, 43.046, 13.780, 11.985, 9.818, 19.208, 17.574, 0, 11.883, 30.644, 26.277, 0, 14.931, 14.364, 11.423, 7.751, 11.549, 0, 4.877, 0, 8.901, 1.860, 8.901, 0, 13.322, 0, 8.301, 0, 0, 0, 51.787, 0, 0, 0, 0, 0, 16.921, 17.117, 13.265, 27.966, 0, 0.292, -0.720, 0.868, 1.027, 2.426, 0, 0, -0.568, -0.905, -0.847, 2.057, -0.073, -0.369, 0.345, -0.114, 0, 0.207, -0.668, 0.071, 0.744, -3.410, 0, 8.502, -3.345, 0, 1.517, 0, -1.398, 0.320, -3.661, 4.626, 0, 0, 2.311, 0, 0, 0.972, 0, 0, 0, -7.488, -4.864, 0],
            "vliq": [0.0261, 0.0164, 0.0071, -0.0038, 0.0373, 0.0269, 0.0270, 0.0161, 0.0030, 0.0434, 0.0132, 0.0044, 0.0289, 0.0192, 0.0099, 0.0055, 0.0113, 0.0365, 0.0282, 0.0200, 0.0450, 0.0357, 0.0267, 0.0327, 0.0231, 0.0180, 0.0206, 0.0265, 0.0195, 0.0267, 0.0232, 0.0181, 0.0191, 0.0168, 0.0137, 0.0608, 0.0524, 0.0331, 0.0223, 0.0337, 0.0266, 0.0202, 0.0468, 0.0620, 0, 0.0241, 0.0338, 0.0262, 0.0250, 0.0345, 0.0279, 0.0214, 0, 0.0145, 0.0153, 0.0173, 0, 0, 0, 0, 0.0192, 0.0538, 0, 0.0538, 0, 0, 0, 0, 0.0548, 0, 0, 0.0410, 0, 0.0348, 0.0273, 0, 0, 0, 0.00133, 0.00179, -0.00203, -0.00243, -0.00744, 0, 0, 0.00213, 0.00063, -0.00519, -0.00188, 0.00009, 0.00012, 0.00142, -0.00107, 0, -0.00009, -0.00030, -0.00108, -0.00111, -0.00036, -0.00050, 0.00777, 0.00083, 0.00036, 0.00198, 0.00001, -0.00092, 0.00175, 0.00235, -0.00250, 0.00046, 0, -0.00179, -0.00206, 0.01203, -0.00023, 0, 0, 0, 0.00178, 0.00171, 0],
            "cpa": [35.1152, 22.6346, 8.9272, 0.3456, 49.2506, 35.2248, 37.6299, 21.3528, 10.2797, 66.0574, 16.3794, 10.4283, 42.8569, 32.8206, 19.9504, 27.2107, 39.7712, 59.3032, 0, 40.7501, 66.8423, 0, 51.5048, 50.5604, 39.5784, 25.6750, 0, 57.6861, 44.1122, 53.7012, 44.6388, 0, 41.4064, 30.1561, 47.1311, 84.7602, 0, 58.2837, 46.5577, 48.4648, 36.5885, 29.1848, 60.8262, 56.1685, 78.6054, 33.6450, 63.7851, 51.1442, 0, 58.2445, 29.1815, 28.0260, 45.9768, 26.7371, 25.8094, 30.1696, 0, 63.2024, 44.3567, 0, 0, 0, 0, 0, 22.2082, 0, 0, 0, 0, 0, 0, 0, 0, 57.7670, 45.0314, 40.5275, 80.3010, 0, 0.5830, 0.3226, 0.9668, -0.3082, -0.1201, 8.5546, 3.1721, -5.9060, -3.9682, -3.2746, 2.6142, -1.3913, 0.2630, 6.5145, 4.1707, 0, 0, 3.7978, 0, 0, 0, 0, -15.7667, 0, 0, -6.4072, 0, 2.4484, -1.5252, 0, 0, 0, 0, 0, 0, 0, -2.7407, 0, -1.6978, 0, -2.2923, -0.3162, 0],
            "cpb": [39.5923, 45.0933, 59.9786, 74.0368, 59.3840, 62.1924, 62.1285, 66.3947, 65.5372, 69.3936, 32.7433, 25.3634, 65.6464, 70.4153, 81.8764, 2.7609, 35.5676, 67.8149, 0, 19.6990, 102.4553, 0, 44.4133, 38.9681, 41.8177, 24.7281, 0, 64.0768, 77.2155, 71.7948, 68.5041, 0, 85.0996, 81.6814, 51.3326, 177.2513, 0, 49.6388, 48.2322, 37.2370, 47.6004, 52.3817, 41.9908, 46.9337, 32.1318, 23.2759, 83.4744, 94.2934, 0, 46.9958, -9.7846, -7.1651, 20.6417, 21.7676, -5.2241, 26.9738, 0, 51.9366, 44.5875, 0, 0, 0, 0, 0, -2.8385, 0, 0, 0, 0, 0, 0, 0, 0, 44.1238, 55.1432, 55.0141, 132.7786, 0, -1.2002, 2.1309, -2.0762, 1.8969, 4.2846, -22.9771, -10.0834, -1.8710, 17.7889, 32.1670, 4.4511, -1.5496, -2.3428, -17.5541, -3.1964, 0, 0, -7.3251, 0, 0, 0, 0, -0.1174, 0, 0, 15.2583, 0, -0.0765, -7.6380, 0, 0, 0, 0, 0, 0, 0, 11.1033, 0, 1.0477, 0, 3.1142, 2.3711, 0],
            "cpc": [-9.9232, -15.7033, -29.5143, -45.7878, -21.7908, -24.8156, -26.0637, -29.3703, -30.6057, -25.1081, -13.1692, -12.7283, -21.0670, -28.9361, -40.2864, 1.3060, -15.5875, -20.9948, 0, -5.4360, -43.3306, 0, -19.6155, -4.7799, -11.0837, 4.2419, 0, -21.0480, -33.5086, -22.9685, -26.7106, 0, -35.6318, -36.1441, -25.0276, -72.3213, 0, -15.6291, -20.4868, -13.0635, -22.8148, -30.8526, -20.4091, -31.3325, -19.4033, -12.2406, -35.1171, -45.2029, 0, -10.5106, 3.4554, 2.4332, -8.3297, -6.4481, 1.4542, -13.3722, 0, -28.6308, -23.2820, 0, 0, 0, 0, 0, 1.2679, 0, 0, 0, 0, 0, 0, 0, 0, -9.5565, -18.7776, -31.7190, -58.3241, 0, -0.0584, -1.5728, 0.3148, -1.6454, -2.0262, 10.7278, 4.9674, 4.2945, -3.3639, -17.8246, -5.9808, 2.5899, 0.8975, 10.6977, -1.1997, 0, 0, 2.5312, 0, 0, 0, 0, 6.1191, 0, 0, -8.3149, 0, 0.1460, 8.1795, 0, 0, 0, 0, 0, 0, 0, -11.0878, 0, 0.2002, 0, -1.4995, -1.4825, -0.0584],
        "txt": [("CH3", {"C": 1, "H": 3}),
                        ("CH2", {"C": 1, "H": 2}),
                        ("CH", {"C": 1, "H": 1}),
                        ("C", {"C": 1}),
                        ("CH2=CH", {"C": 2, "H": 3}),
                        ("CH=CH", {"C": 2, "H": 2}),
                        ("CH2=C", {"C": 2, "H": 2}),
                        ("CH=C", {"C": 2, "H": 1}),
                        ("C=C", {"C": 2}),
                        ("CH2=C=CH", {"C": 3, "H": 3}),
                        ("-CH (Aromatic)", {"C": 1, "H": 1}),
                        ("=C (Aromatic)", {"C": 1}),
                        ("-CCH3 (Aromatic)", {"C": 2, "H": 3}),
                        ("-CCH2 (Aromatic)", {"C": 2, "H": 2}),
                        ("-CCH (Aromatic)", {"C": 2, "H": 1}),
                        ("-OH", {"O": 1, "H": 1}),
                        ("-OH (Aromatic)", {"O": 1, "H": 1}),
                        ("CH3COO", {"C": 2, "H": 3, "O": 2}),
                        ("CH2COO", {"C": 2, "H": 2, "O": 2}),
                        ("HCOO", {"C": 1, "H": 1, "O": 2}),
                        ("CH3O", {"C": 1, "H": 3, "O": 1}),
                        ("CH2O", {"C": 1, "H": 2, "O": 1}),
                        ("CH-O", {"C": 1, "H": 1, "O": 1}),
                        ("FCH2O", {"C": 1, "H": 2, "O": 1, "F": 1}),
                        ("CH2NH2", {"C": 1, "H": 4, "N": 1}),
                        ("CHNH2", {"C": 1, "H": 3, "N": 1}),
                        ("CH3NH", {"C": 1, "H": 4, "N": 1}),
                        ("CH2NH", {"C": 1, "H": 3, "N": 1}),
                        ("CHNH", {"C": 1, "H": 2, "N": 1}),
                        ("CH3N", {"C": 1, "H": 3, "N": 1}),
                        ("CH2N", {"C": 1, "H": 2, "N": 1}),
                        ("=CNH2 (Aromatic)", {"C": 1, "H": 2, "N": 1}),
                        ("C5H4N", {"C": 5, "H": 4, "N": 1}),
                        ("C5H3N", {"C": 5, "H": 3, "N": 1}),
                        ("CH2CN", {"C": 2, "H": 2, "N": 1}),
                        ("COOH", {"C": 1, "H": 1, "O": 2}),
                        ("CH2Cl", {"C": 1, "H": 2, "Cl": 1}),
                        ("CHCl", {"C": 1, "H": 1, "Cl": 1}),
                        ("CCl", {"C": 1, "Cl": 1}),
                        ("CHCl2", {"C": 1, "H": 1, "Cl": 2}),
                        ("CCl3", {"C": 1, "Cl": 3}),
                        ("CCl2", {"C": 1, "Cl": 2}),
                        ("=CCl (Aromatic)", {"C": 1, "Cl": 1}),
                        ("CH2NO2", {"C": 1, "H": 2, "O": 2, "N": 1}),
                        ("CHNO2", {"C": 1, "H": 1, "O": 2, "N": 1}),
                        ("=CNO2 (Aromatic)", {"C": 1, "O": 2, "N": 1}),
                        ("CH2SH", {"C": 1, "H": 3, "S": 1}),
                        ("I", {"I": 1}),
                        ("Br", {"Br": 1}),
                        ("CH≡C", {"C": 2, "H": 1}),
                        ("C≡C", {"C": 2}),
                        ("Cl-C=C", {"C": 2, "Cl": 1}),
                        ("=CF (Aromatic)", {"C": 1, "F": 1}),
                        ("HCON(CH2)2", {"C": 3, "H": 5, "O": 1, "N": 1}),
                        ("CF3", {"C": 1, "F": 3}),
                        ("CF2", {"C": 1, "F": 2}),
                        ("CF", {"C": 1, "F": 1}),
                        ("COO", {"C": 1, "O": 2}),
                        ("CCl2F", {"C": 1, "F": 1,  "Cl": 2}),
                        ("HCClF", {"C": 1, "H": 1, "F": 1, "Cl": 1}),
                        ("CClF2", {"C": 1, "F": 2,  "Cl": 1}),
                        ("F (others)", {"F": 1}),
                        ("CONH2", {"C": 1, "H": 2, "O": 1, "N": 1}),
                        ("CONHCH3", {"C": 2, "H": 4, "O": 1, "N": 1}),
                        ("CONHCH2", {"C": 2, "H": 3, "O": 1, "N": 1}),
                        ("CON(CH3)2", {"C": 3, "H": 6, "O": 1, "N": 1}),
                        ("CONCH2CH2", {"C": 3, "H": 4, "O": 1, "N": 1}),
                        ("CON(CH2)2", {"C": 3, "H": 4, "O": 1, "N": 1}),
                        ("C2H5O2", {"C": 2, "H": 5, "O": 2}),
                        ("C2H4O2", {"C": 2, "H": 4, "O": 2}),
                        ("CH3S", {"C": 1, "H": 3, "S": 1}),
                        ("CH2S", {"C": 1, "H": 2, "S": 1}),
                        ("CHS", {"C": 1, "H": 1, "S": 1}),
                        ("C4H3S", {"C": 4, "H": 3, "S": 1}),
                        ("C4H2S", {"C": 4, "H": 2, "S": 1}),

                        ("CH(CH3)2", ),
                        ("C(CH3)3", ),
                        ("CHCH3CHCH3", ),
                        ("CH(CH3)C(CH3)2", ),
                        ("C(CH3)2C(CH3)2", ),
                        ("3 membered ring", ),
                        ("4 membered ring", ),
                        ("5 membered ring", ),
                        ("6 membered ring", ),
                        ("7 membered ring", ),
                        ("CHn=CHm-CHp=CHk (m, p (0,1); k, n (0,2)", ),
                        ("CH3-CHm=CHn (m (0,1); n (0,2))", ),
                        ("CH2-CHm=CHn (m (0,1); n (0,2))", ),
                        ("CH-CHm=CHn (m (0,1); n (0,2))", ),
                        ("Alicyclic side-chain CcyclicCm", ),
                        ("CH3CH3", ),
                        ("CHCHO or CCHO", ),
                        ("CH3COCH2", ),
                        ("CH3COCH or CH3COC", ),
                        ("Ccyclic=O", ),
                        ("ACCHO", ),
                        ("CHCOOH or CCOOH", ),
                        ("ACCOOH", ),
                        ("CH3COOCH or CH3COOC", ),
                        ("COCH2COO or COCHCOO or COCCOO", ),
                        ("CO-O-CO", ),
                        ("ACCOO", ),
                        ("CHOH", ),
                        ("COH", ),
                        ("CHm(OH)CHn(OH) (0<m,n<2)", ),
                        ("CHm cyclic-OH (0<m<1)", ),
                        ("CHn(OH)CHm(NHp) (0<m<1); (0<n,p<2)", ),
                        ("CHm(NH2)CHn(NH2) (0<m,n<2)", ),
                        ("CHm cyclic-NHp-CHn cyclic (0<n,m,p<1)", ),
                        ("CHn-O-CHm=CHp (0<m<1); (0<n,p<2)", ),
                        ("AC-O-CHm (0<m<3)", ),
                        ("CHm cyclic-S-CHn cyclic (0<n,m<1)", ),
                        ("CHn=CHm-F (0<m<1); (0<n<2)", ),
                        ("CHn=CHm-Br (0<m<1); (0<n<2)", ),
                        ("CHn=CHm-I (0<m<1); (0<n<2)", ),
                        ("ACBr", ),
                        ("ACI", ),
                        ("CHm(NH2)-COOH (0<m<2)", )]
            }

    FirstOrder=75
    SecondOrder=118

    def calculo(self):
        if self.kwargs["M"]:
            self.M=self.kwargs["M"]
        else:
            self.M=sum([self.coeff["M"][grupo]*contribucion for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]) if grupo<=self.FirstOrder])

        tc=Pc=tf=tb=vc=w=gf=hf=hv=vliq=cpa=cpb=cpc=0
        for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tf+=self.coeff["tf"][grupo]*contribucion
            tb+=self.coeff["tb"][grupo]*contribucion
            tc+=self.coeff["tc"][grupo]*contribucion
            Pc+=contribucion*self.coeff["Pc"][grupo]
            vc+=contribucion*self.coeff["vc"][grupo]
            w+=contribucion*self.coeff["w"][grupo]
            hf+=contribucion*self.coeff["hf"][grupo]
            gf+=contribucion*self.coeff["gf"][grupo]
            hv+=contribucion*self.coeff["hv"][grupo]
            vliq+=contribucion*self.coeff["vliq"][grupo]
            cpa+=contribucion*self.coeff["cpa"][grupo]
            cpb+=contribucion*self.coeff["cpb"][grupo]
            cpc+=contribucion*self.coeff["cpc"][grupo]

        self.Tf=unidades.Temperature(102.425*log(tf))
        self.Tb=unidades.Temperature(204.359*log(tb))
        self.Tc=unidades.Temperature(181.128*log(tc))
        self.Pc=unidades.Pressure((Pc+0.10022)**-2+1.3705, "bar")
        self.Vc=unidades.SpecificVolume((vc-0.00435)/self.M)
        self.f_acent=0.4085*log(w+1.1507)**(1/0.5050)
        self.Hf=unidades.Enthalpy((hf+10.835)/self.M, "kJg")
        self.Gf=unidades.Enthalpy((gf-14.83)/self.M, "kJg")
        self.Hv=unidades.Enthalpy((hv+6.829)/self.M, "kJg")
        self.Vliq=unidades.SpecificVolume((vliq-0.00435)/self.M, "ccg")
        self.cp=[cpa-19.7779, cpb+22.5981, cpc-10.7983, 0, 0, 0]

        GroupContribution.calculo(self)


class Wilson_Jasperson(GroupContribution):
    """Des. Dev., 20: 94 (1981). Wilson. G. M., and L. V. Jasperson: ‘‘Critical Constants Tc, Pc, Estimation Based on Zero, First and Second Order Methods,’’ AIChE Spring Meeting, New Orleans, LA, 1996.
    ref, properties of gases and liquids, pag 2.9
    grupos: grupos que forman la molécula
    contribuciones: contribuciones de cada grupo
    Tb: Temperatura de ebullición
    Ring: anillos que forman la molecula
    M: peso molecular, opcional
    SG: gravedad específica, opcional

    >>> cresol_wilson=Wilson_Jasperson(group=[3, 0, 5, 41], contribution=[8, 10, 1, 1], M=108.14, Tb=464.15, ring=1)
    >>> print cresol_wilson.Tc, cresol_wilson.f_acent
    654.074200805 0.620089916983
    """
    coeff={
        "M": [1.00794, 4.002602, 10.811, 12.0107, 14.0067, 15.9994, 18.9984032, 20.1797, 26.9815386, 28.0855, 30.973762, 32.065, 35.453, 39.948, 47.867, 50.9415, 69.723, 72.64, 74.9216, 78.96, 79.904, 83.798, 85.4678, 91.224, 92.90638, 95.94, 118.71, 121.76, 127.6, 126.90447, 131.293, 132.9054519, 178.49, 180.94788, 183.84, 186.207, 190.23, 200.59, 208.9804, 222, 238.02891],
        "tc": [0.002793, 0.320000, 0.019000, 0.008532, 0.019181, 0.020341, 0.008810, 0.036400, 0.088000, 0.020000, 0.012000, 0.007271, 0.011151, 0.016800, 0.014000, 0.018600, 0.059000, 0.031000, 0.007000, 0.010300, 0.012447, 0.013300, -0.027000, 0.175000, 0.017600, 0.007000, 0.020000, 0.010000, 0.000000, 0.005900, 0.017000, -0.027500, 0.219000, 0.013000, 0.011000, 0.014000, -0.050000, 0.000000, 0.000000, 0.007000, 0.015000, 0.0350, 0.0100, -0.0075, -0.0040, 0.0000, -0.0550, 0.0170, -0.0150, 0.0170, -0.0200, 0.0020, 0.0000, -0.0250],
        "Pc": [0.12660, 0.43400, 0.91000, 0.72983, 0.44805, 0.43360, 0.32868, 0.12600, 6.05000, 1.34000, 1.22000, 1.04713, 0.97711, 0.79600, 1.19000, 0.0, 0.0, 1.42000, 2.68000, 1.20000, 0.97151, 1.11000, 0.0, 1.11000, 2.71000, 1.69000, 1.95000, 0.0, 0.43000, 1.315930, 1.66000, 6.33000, 1.07000, 0.0, 1.08000, 0.0, 0.0, -0.08000, 0.69000, 2.05000, 2.04000, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.50, 0.00, 1.50, 1.00, 0.00, 0.00, -0.50],
        "txt": [("H", ),
                    ("He",),
                    ("B",),
                    ("C",),
                    ("N",),
                    ("O",),
                    ("F",),
                    ("Ne",),
                    ("Al",),
                    ("Si",),
                    ("P",),
                    ("S",),
                    ("Cl",),
                    ("Ar",),
                    ("Ti",),
                    ("V",),
                    ("Ga",),
                    ("Ge",),
                    ("As",),
                    ("Se",),
                    ("Br",),
                    ("Kr",),
                    ("Rb",),
                    ("Zr",),
                    ("Nb",),
                    ("Mo",),
                    ("Sn",),
                    ("Sb",),
                    ("Te",),
                    ("I",),
                    ("Xe",),
                    ("Cs",),
                    ("Hf",),
                    ("Ta",),
                    ("W",),
                    ("Re",),
                    ("Os",),
                    ("Hg",),
                    ("Bi",),
                    ("Rn",),
                    ("U",),
                    ("-OH, C4 or less",),
                    ("-OH, C5 or more",),
                    ("-O-",),
                    ("-NH2, >NH, >N-",),
                    ("-CHO",),
                    (">CO",),
                    ("-COOH",),
                    ("-COO-",),
                    ("-CN",),
                    ("-NO2",),
                    ("Organic Halides (once / molecule)",),
                    ("-SH, -S-, -SS-",),
                    ("Siloxane bond",)]}

    FirstOrder=41
    SecondOrder=54

    def isCalculable(self):
        if not self.kwargs["Tb"]:
            self.msg=QApplication.translate("pychemqt", "undefined boiling point")
            self.status=0
        else:
            return GroupContribution.isCalculable(self)


    def calculo(self):
        self.Tb=unidades.Temperature(self.kwargs["Tb"])
        if self.kwargs["M"]:
            self.M=self.kwargs["M"]
        else:
            self.M=sum([self.coeff["M"][grupo]*contribucion for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]) if grupo<=self.FirstOrder])

        tc=Pc=0
        for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc+=self.coeff["tc"][grupo]*contribucion
            Pc+=contribucion*self.coeff["Pc"][grupo]

        self.Tc=unidades.Temperature(self.Tb/(0.048271-0.019846*self.kwargs["ring"]+tc)**0.2)
        self.Pc=unidades.Pressure(0.0186233*self.Tc/(-0.96601+exp(-0.00922295-0.0290403*self.kwargs["ring"]+0.041*Pc)), "bar")

        GroupContribution.calculo(self)

    def EmpiricFormula(self):
        string=""
        formula=""
        for i, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            if i<self.FirstOrder:
                if contribucion>1:
                    string+="%s<sub>%i</sub>" % (self.coeff["txt"][i][0], contribucion)
                    formula+="%s%i" % (self.coeff["txt"][i][0], contribucion)
                elif contribucion==1:
                    string+="%s" %self.coeff["txt"][i][0]
                    formula+="%s" %self.coeff["txt"][i][0]
        return string, formula


class Marrero_Pardillo(GroupContribution):
    """Marrero-Marejon, J., and E. Pardillo-Fontdevila: AIChE J., 45: 615 (1999).
    ref, properties of gases and liquids, pag 2.9
    grupos: grupos que forman la molécula
    contribuciones: contribuciones de cada grupo
    Atomos: átomos que forman la molecula
    M: peso molecular
    Tb: Temperatura de ebullición, opcional
    SG: gravedad específica, opcional

    >>> cresol_marrero=Marrero_Pardillo(group=[1, 36, 129, 130, 132, 140, 148], contribution=[1, 1, 1, 2, 2, 1, 1], M=122.17, atomos=19)
    >>> print cresol_marrero.Tc, cresol_marrero.Pc.bar
    704.321804655 42.2204405196
    """
    coeff={
        "tc": [-0.0213, -0.0227, -0.0223, -0.0189, 0.8526, 0.1792, 0.3818, -0.0214, 0.1117, 0.0987, -0.0370, -0.9141, -0.9166, -0.9146, -0.0876, -0.0205, -0.0362, -0.0606, -0.0890, 0.0267, -0.0974, -0.0397, -0.0313, -0.0199, -0.0766, -0.0591, -0.9192, -0.0181, -0.0206, -0.0134, -0.0098, 0.8636, 0.1874, 0.4160, -0.0149, 0.1193, 0.1012, -0.0255, -0.0162, -0.0205, -0.0210, -0.0786, -0.0205, -0.0256, -0.0267, -0.0932, 0.0276, -0.0993, -0.0301, -0.0248, -0.0161, -0.0654, -0.0137, -0.0192, -0.0039, 0.0025, 0.8547, 0.1969, 0.0025, 0.1187, -0.0200, -0.0142, -0.0757, -0.0162, -0.0194, -0.0406, -0.0918, -0.1054, -0.0286, -0.0158, 0.0084, 0.8767, 0.2061, 0.0207, 0.0049, 0.1249, -0.0176, -0.0133, -0.0084, -0.0780, -0.0156, -0.0114, -0.1008, -0.9129, -0.8933, -0.4158, -0.0123, -1.7660, -1.2909, -0.8945, 1.7377, 1.0731, 1.2865, 0.9929, 0.8623, 0.8613, 0.8565, 0.8246, 0.7862, 0.8818, 0.7780, 0.8122, -0.8155, -0.4009, 0.3043, 0.1868, 0.1886, -0.0159, -0.0288, -0.4222, -0.7958, -0.0098, -0.0093, -0.1386, 0.0976, 0.1089, -0.0092, -0.0148, -0.0139, -0.0071, -0.0055, -0.1341, 0.0, 0.0, -0.0218, -0.0737, 0.0329, 0.0, -0.0314, -0.2246, -0.3586, 0.3913, 0.2089, 0.2190, 0.1000, 0.0947, -0.4067, 0.1027, -0.4848, 0.2541, 0.2318, 0.2424, 0.1104, -0.3972, 0.1069, 0.1028, 0.1060, 0.1075, 0.0931, 0.0997, 0.1112, 0.0919, 0.0313, 0.0241, 0.0830, 0.0978, 0.0938, 0.0768, -0.0191, -0.1926, -0.5728, -0.3553, -0.0422, -0.0690, -0.0781, -0.0301, -0.0124],
        "Pc": [-0.0618, -0.0430, -0.0376, -0.0354, 0.0654, 0.0851, -0.2320, -0.0396, -0.0597, -0.0746, -0.0345, -0.0231, -0.0239, -0.0241, -0.0180, -0.0321, -0.0363, -0.0466, -0.0499, 0.1462, -0.2290, -0.0288, -0.0317, -0.0348, -0.0507, -0.0385, -0.0244, -0.0305, -0.0272, -0.0219, -0.0162, 0.0818, 0.1010, -0.2199, -0.0265, -0.0423, -0.0626, -0.0161, -0.0150, -0.0140, -0.0214, -0.0119, -0.0184, -0.0204, -0.0210, -0.0253, 0.1561, -0.2150, -0.0214, -0.0203, -0.0170, -0.0329, -0.0163, -0.0173, -0.0137, -0.0085, 0.0816, 0.1080, -0.0168, -0.0556, -0.0147, -0.0131, -0.0093, -0.0155, -0.0112, -0.0280, -0.2098, -0.0358, -0.0212, -0.0162, 0.0002, 0.0953, 0.1109, 0.0213, -0.0111, -0.0510, -0.0161, -0.0129, -0.0121, -0.0094, -0.0103, -0.0085, -0.0455, -0.0476, -0.1378, -0.2709, -0.0239, -0.2291, -0.3613, -0.1202, 0.1944, 0.2146, -0.1087, 0.0533, 0.0929, 0.0919, 0.0947, 0.0801, 0.0806, 0.2743, -0.1007, 0.0771, -0.4920, -0.2502, 0.0705, 0.1064, 0.1102, -0.0010, -0.0226, 0.1860, 0.3933, -0.0221, -0.0181, 0.0081, -0.1034, -0.0527, -0.0119, -0.0177, -0.0127, 0.0, -0.0088, 0.0162, 0.0, 0.0, -0.0091, -0.0220, -0.0071, 0.0, -0.0119, 0.1542, 0.1490, 0.1356, -0.1822, -0.1324, -0.0900, 0.0, -0.1491, -0.0916, 0.1432, 0.0, -0.0809, -0.0792, -0.0374, -0.0971, -0.0504, -0.0512, -0.0548, -0.0514, -0.0388, -0.0523, -0.0528, -0.0597, -0.0684, -0.2573, -0.0579, -0.0471, -0.0462, -0.0625, -0.0125, -0.0878, 0.0, -0.0176, -0.0123, 0.0, -0.1878, 0.0, 0.0],
        "vc": [123.2, 88.6, 78.4, 69.8, 81.5, 57.7, 65.8, 58.3, 49.0, 71.7, 88.1, 113.8, 0.0, 0.0, 92.9, 66.0, 88.9, 128.9, 145.9, 93.3, 108.2, 0.0, 0.0, 76.3, 147.9, 148.1, 119.7, 87.9, 56.6, 40.2, 32.0, 50.7, 24.0, 33.9, 31.9, 0.0, 52.1, 49.3, 80.8, 101.3, 0.0, 45.2, 34.5, 62.3, 106.1, 114.0, 69.9, 79.1, 63.3, 49.4, 32.7, 113.5, 93.3, 57.9, 18.3, 8.6, 48.9, 4.3, 0.0, 0.0, 37.7, 68.6, 45.6, 23.7, 39.3, 92.2, 72.3, 110.2, 39.2, 0.0, 22.7, 23.4, 8.8, 0.0, 0.0, 0.0, 30.0, 63.7, 85.7, 40.6, 40.8, 62.1, 89.0, 105.3, 77.4, 99.2, 68.4, 47.8, 73.6, 43.6, 42.1, 16.6, 0.0, 0.0, 41.4, 68.7, 36.4, 0.0, 107.4, 55.2, 64.1, 107.4, 93.7, 58.1, 0.0, 14.6, 43.3, 51.4, 87.6, 73.1, 64.3, 47.2, 47.5, 49.9, 42.5, 0.0, 29.2, 50.7, 38.8, 0.0, 33.9, 0.0, 0.0, 0.0, 19.2, 0.0, 36.2, 0.0, 18.4, 36.5, 34.4, 8.3, 39.3, 29.8, 40.3, 0.0, 65.9, 40.8, 37.8, 0.0, 20.6, 51.7, -0.3, 35.6, 23.7, 60.3, 83.2, 110.2, 8.5, 0.0, 46.3, 0.0, 100.2, 55.2, 33.2, 0.0, 0.0, 0.0, 84.0, 0.0, 0.0, 0.0, 0.0, 0.0, 51.2, 0.0, 0.0],
        "tb": [113.12, 194.25, 194.27, 186.41, 137.18, 182.20, 194.40, 176.16, 180.60, 145.56, 160.83, 453.70, 758.44, 1181.44, 736.93, 228.01, 445.61, 636.49, 1228.84, 456.92, 510.65, 443.76, 293.86, 207.75, 891.15, 1148.58, 588.31, 409.85, 244.88, 244.14, 273.26, 201.80, 242.47, 207.49, 238.81, 260.00, 167.85, 166.59, 517.62, 875.85, 1262.80, 673.24, 243.37, 451.27, 648.70, 1280.39, 475.65, 541.29, 452.30, 314.71, 240.08, 869.18, 612.31, 451.03, 291.41, 344.06, 179.96, 249.10, 295.33, 132.66, 68.80, 438.47, 585.99, 215.94, 434.45, 630.07, 497.58, 1270.16, 388.44, 260.32, 411.56, 286.30, 286.42, 456.90, 340.00, 188.99, -16.64, 360.79, 610.26, 540.38, 267.26, 373.71, 1336.54, 51.13, 205.73, 245.27, 183.55, 334.64, 354.41, 316.46, 174.18, 228.38, 174.39, 184.20, 5.57, 370.60, 204.81, 658.53, 1245.86, 423.86, 525.35, 761.36, 399.58, 321.02, 250.88, -37.99, 367.05, 160.42, 120.85, 222.40, 333.26, 201.89, 209.40, 182.74, 218.07, 106.21, 225.52, 451.74, 283.55, 424.13, 210.66, 220.24, 254.50, 184.36, 169.17, 597.82, 348.23, 111.51, -41.35, 112.00, 291.15, 221.55, 285.07, 237.22, 171.59, 420.54, 321.44, 348.00, 477.77, 334.09, 180.07, 123.05, 134.23, 174.31, -48.79, 347.33, 716.23, 1294.98, 456.25, 199.70, 437.51, 700.06, 1232.55, 437.78, 517.75, 411.29, 422.51, 682.19, 532.24, 1012.51, 382.25, 385.36, 387.17, 1022.45, 298.12, 673.59, 597.59],
        "txt": [("CH3- & CH3-",),
                    ("CH3- & -CH2-",),
                    ("CH3- & >CH-",),
                    ("CH3- & >C<",),
                    ("CH3- & =CH-",),
                    ("CH3- & =C<",),
                    ("CH3- & ≡C-",),
                    ("CH3- & >CH- [r]",),
                    ("CH3- & >C< [r]",),
                    ("CH3- & =C< [r]",),
                    ("CH3- & F-",),
                    ("CH3- & Cl-",),
                    ("CH3- & Br-",),
                    ("CH3- & I-",),
                    ("CH3- & -OH",),
                    ("CH3- & -O-",),
                    ("CH3- & >CO",),
                    ("CH3- & -CHO",),
                    ("CH3- & -COOH",),
                    ("CH3- & -COO[-]",),
                    ("CH3- & [-]COO-",),
                    ("CH3- & -NH2",),
                    ("CH3- & -NH-",),
                    ("CH3- & >N-",),
                    ("CH3- & -CN",),
                    ("CH3- & -NO2",),
                    ("CH3- & -SH",),
                    ("CH3- & -S-",),
                    ("-CH2- & -CH2-",),
                    ("-CH2- & >CH-",),
                    ("-CH2- & >C<",),
                    ("-CH2- & =CH-",),
                    ("-CH2- & =C<",),
                    ("-CH2- & ≡C-",),
                    ("-CH2- & >CH- [r]",),
                    ("-CH2- & >C< [r]",),
                    ("-CH2- & =C< [r]",),
                    ("-CH2- & F-",),
                    ("-CH2- & Cl-",),
                    ("-CH2- & Br-",),
                    ("-CH2- & I-",),
                    ("-CH2- & -OH",),
                    ("-CH2- & -O-",),
                    ("-CH2- & >CO",),
                    ("-CH2- & -CHO",),
                    ("-CH2- & -COOH",),
                    ("-CH2- & -COO[-]",),
                    ("-CH2- & [-]COO-",),
                    ("-CH2- & -NH2",),
                    ("-CH2- & -NH-",),
                    ("-CH2- & >N-",),
                    ("-CH2- & -CN",),
                    ("-CH2- & -SH",),
                    ("-CH2- & -S-",),
                    (">CH- & CH-",),
                    (">CH- & >C<",),
                    (">CH- & =CH-",),
                    (">CH- & =C<",),
                    (">CH- & >CH- [r]",),
                    (">CH- & =C< [r]",),
                    (">CH- & F-",),
                    (">CH- & Cl-",),
                    (">CH- & -OH",),
                    (">CH- & -O-",),
                    (">CH- & >CO",),
                    (">CH- & -CHO",),
                    (">CH- & [-]COO-",),
                    (">CH- & -COOH",),
                    (">CH- & -NH2",),
                    (">CH- & -NH-",),
                    (">C< & >C<",),
                    (">C< & =CH-",),
                    (">C< & =C<",),
                    (">C< & >C< [r]",),
                    (">C< & >CH- [r]",),
                    (">C< & =C< [r]",),
                    (">C< & F-",),
                    (">C< & Cl-",),
                    (">C< & Br-",),
                    (">C< & -OH",),
                    (">C< & -O-",),
                    (">C< & >CO",),
                    (">C< & -COOH",),
                    ("[=]CH2 & [=]CH2",),
                    ("[=]CH2 & -CH[=]",),
                    ("[=]CH2 & >C[=]",),
                    ("[=]CH2 & =C[=]",),
                    ("-CH[=] & -CH[=]",),
                    ("-CH[=] & >C[=]",),
                    ("-CH[=] & =C[=]",),
                    ("=CH- & =CH-",),
                    ("=CH- & =C<",),
                    ("=CH- & ≡C-",),
                    ("=CH- & =C< [r]",),
                    ("=CH- & F-",),
                    ("=CH- & Cl-",),
                    ("=CH- & -O-",),
                    ("=CH- & -CHO",),
                    ("=CH- & -COOH",),
                    ("=CH- & -COO[-]",),
                    ("=CH- & [-]COO-",),
                    ("=CH- & -CN",),
                    (">C[=] & >C[=]",),
                    (">C[=] & =C[=]",),
                    ("=C< & =C< [r]",),
                    ("=C< & F-",),
                    ("=C< & Cl-",),
                    ("=C[=] & O[=]",),
                    ("CH[≡] & CH[≡]",),
                    ("CH[≡] & -C[≡]",),
                    ("-C[≡] & -C[≡]",),
                    ("-CH2- [r] & -CH2- [r]",),
                    ("-CH2- [r] & >CH- [r]",),
                    ("-CH2- [r] & >C< [r]",),
                    ("-CH2- [r] & =CH- [r]",),
                    ("-CH2- [r] & =C< [r]",),
                    ("-CH2- [r] & -O- [r]",),
                    ("-CH2- [r] & >CO [r]",),
                    ("-CH2- [r] & -NH- [r]",),
                    ("-CH2- [r] & -S- [r]",),
                    (">CH- [r] & >CH- [r]",),
                    (">CH- [r] & >C< [r]",),
                    (">CH- [r] & >CH- [rr]",),
                    (">CH- [r] & >C[=] [rr]",),
                    (">CH- [r] & -O- [r]",),
                    (">CH- [r] & -OH",),
                    (">C< [r] & >C< [r]",),
                    (">C< [r] & =C< [r]",),
                    (">C< [r] & F-",),
                    ("-CH[=] [r] & -CH[=] [r]",),
                    ("-CH[=] [r] & >C[=] [r]",),
                    ("-CH[=] [r] & -N[=] [r]",),
                    ("=CH- [r] & =CH- [r]",),
                    ("=CH- [r] & =C< [r]",),
                    ("=CH- [r] & -O- [r]",),
                    ("=CH- [r] & -NH- [r]",),
                    ("=CH- [r] & =N- [r]",),
                    ("=CH- [r] & -S- [r]",),
                    (">C[=] [r] & >C[=] [r]",),
                    (">C[=] [r] & -N[=] [r]",),
                    ("=C< [r] & =C< [r]",),
                    ("=C< [r] & =C< [rr]",),
                    ("=C< [r] & -O- [r]",),
                    ("=C< [r] & =N- [r]",),
                    ("=C< [r] & F-",),
                    ("=C< [r] & Cl-",),
                    ("=C< [r] & Br-",),
                    ("=C< [r] & I-",),
                    ("=C< [r] & -OH",),
                    ("=C< [r] & -O-",),
                    ("=C<[r] & >CO",),
                    ("=C<[r] & -CHO",),
                    ("=C<[r] & -COOH",),
                    ("=C<[r] & [-]COO-",),
                    ("=C<[r] & -NH2",),
                    ("=C<[r] & -NH-",),
                    ("=C<[r] & >N-",),
                    ("=C<[r] & -CN",),
                    ("Cl- & >CO",),
                    ("[-]COO- & [-]COO-",),
                    ("-O- [r] & =N- [r]",),
                    (">CO & -O-",),
                    ("-H & -CHO",),
                    ("-H & -COOH",),
                    ("-H & [-]COO-",),
                    ("-NH- & -NH2",),
                    ("-S- & -S-",)]}

    FirstOrder=29

    def isCalculable(self):
        if not self.kwargs["atomos"]:
            self.msg=QApplication.translate("pychemqt", "undefined atoms number")
            self.status=0
        elif not self.kwargs["M"]:
            self.msg=QApplication.translate("pychemqt", "undefined molecular weight")
            self.status=0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        if self.kwargs["M"]:
            self.M=self.kwargs["M"]
        else:
            self.M=sum([self.coeff["M"][grupo]*contribucion for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"])])

        tc=Pc=vc=tb=0
        for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tb+=self.coeff["tb"][grupo]*contribucion
            tc+=self.coeff["tc"][grupo]*contribucion
            Pc+=contribucion*self.coeff["Pc"][grupo]
            vc+=contribucion*self.coeff["vc"][grupo]

        if self.kwargs["Tb"]:
            self.Tb=unidades.Temperature(self.kwargs["Tb"])
        else:
            self.Tb=unidades.Temperature(self.M**-0.404*tb+156.)
        self.Tc=unidades.Temperature(self.Tb/(0.5851-0.9286*tc-tc**2))
        self.Pc=unidades.Pressure((0.1285-0.0059*self.kwargs["atomos"]-Pc)**-2, "bar")
        self.Vc=unidades.SpecificVolume((25.1+vc)/self.M, "ccg")

        GroupContribution.calculo(self)

    def EmpiricFormula(self):
        return "", ""


class Elliott(GroupContribution):
    """Zuppo and Elliott, Ind. Eng. Chem. Res. Submitted (1999).
    ref, chemcad propiedades fisicas pag 62
    grupos: grupos que forman la molécula
    contribuciones: contribuciones de cada grupo
    M: peso molecular
    Tb: Temperatura de ebullición, opcional
    SG: gravedad específica, opcional

    >>> elliot=Elliott(group=[0, 5], contribution=[4, 1], M=72)
    >>> print elliot.Tb, elliot.Tc
    333.268576405 829.20395796
    """
    coeff={
        "tc": [0.135, 0.131, 0.077, 0.073, 0.070, -0.015, 0.070, 0.169, 0.169, 0.169, 0.169, 0.169, 0.338, 0.069, 0.099, 0.221, 0.207, 0.136, 0.554, 0.0, 0.0, 0.278, 0.387, 0.383, 0.299, 0.457, 0.453, 0.305, 0.234, 0.230, 0.175, 0.140, 0.0, 0.301, 0.247, 0.306, 0.301, 0.247, 0.148, 0.144, 0.270, 0.0, 0.433, 0.433, 0.0, 0.512, 0.615, 0.0, 0.236, 0.178, 0.090, 0.0, 0.283, 0.196, 0.0, 0.326, 0.0, 0.165, 0.0, 0.440, 0.440, 0.440, 0.0, 0.0, 0.203, 0.0, 0.0, 0.056, 0.056, 0.125, 0.125, 0.0, 0.0, 0.082, 0.147, 0.0, 0.0, 0.340, 0.222, 0.103, 0.327, 0.209, 0.205, 0.151, 0.144, 0.245, 0.245, 0.215, 0.148, 0.0, 0.314, 0.0, 0.209, 0.327, 0.0, 0.0, 0.0, 0.0, 0.422, 0.557, 0.553, 0.670, 0.666, 0.662, 0.839, 0.609, 0.207, 0.203, 0.149, 0.0, 0.0, 0.379, 0.372, 0.0],
        "Pc": [0.232, 0.224, 0.177, 0.186, 0.195, 0.143, 0.204, 0.360, 0.360, 0.360, 0.360, 0.360, 0.720, 0.153, 0.173, 0.375, 0.370, 0.356, 0.075, 0.0, 0.0, 0.126, 0.513, 0.504, 0.324, 0.712, 0.704, 0.455, 0.367, 0.358, 0.311, 0.249, 0.0, 0.316, 0.269, 0.324, 0.316, 0.269, 0.313, 0.304, 0.211, 0.0, 0.869, 0.869, 0.0, 0.564, 0.511, 0.0, 0.542, 0.504, 0.461, 0.0, 0.822, 0.779, 0.0, 1.161, 0.0, 0.460, 0.0, 0.617, 0.617, 0.617, 0.0, 0.0, 0.476, 0.0, 0.0, 0.816, 0.522, 0.274, 0.274, 0.0, 0.0, 0.318, 0.340, 0.0, 0.0, 0.886, 0.638, 0.391, 0.485, 0.398, 0.298, 0.251, 0.269, 0.675, 0.675, 0.645, 0.200, 0.0, 1.027, 0.0, 0.709, 0.956, 0.0, 0.0, 0.0, 0.0, 0.372, 0.605, 0.596, 0.946, 0.937, 0.929, 0.658, 0.761, 0.485, 0.476, 0.429, 0.0, 0.0, 0.960, 0.978, 0.0],
        "vc": [40, 41, 25, 30, 37, 5, 55, 32, 32, 32, 32, 32, 64, 16, 87, 68, 95, 107, -25, 0.0, 0.0, -20, 77, 78, -8, 102, 103, -6, 41, 42, 27, -57, 0.0, 78, 62, 77, 78, 62, 111, 112, 24, 0.0, 107, 107, 0.0, 27, -31, 0.0, 79, 68, 43, 0.0, 107, 82, 0.0, 124, 0.0, 47, 0.0, 34, 34, 34, 0.0, 0.0, 65, 0.0, 0.0, -7, 6, -12, -12, 0.0, 0.0, 23, 27, 0.0, 0.0, 188, 127, 66, 47, -6, 41, 25, 37, 108, 108, 108, -15, 0.0, 143, 0.0, 104, 165, 0.0, 0.0, 0.0, 0.0, 73, 114, 115, 101, 102, 103, 55, 109, 64, 65, 49, 0.0, 0.0, 125, 137, 0.0],
        "tb": [123, 121, 138,  97, 107,  74,  20, 257, 257, 257, 257, 257, 514, 124, 247, 282, 303, 191, 474,  0.0,  0.0, 525, 514, 512, 396, 451, 573, 426, 288, 286, 262, 323,  0.0, 437, 412, 444, 442, 418, 293, 291, 655,  0.0, 942, 942,  0.0, 794, 858,  0.0, 360, 336, 313,  0.0, 575, 552,  0.0, 598,  0.0, 358,  0.0, 692, 668, 818,  0.0,  0.0, 515,  0.0,  0.0, 525, 353, 288, 288,  0.0,  0.0, 190, 135,  0.0,  0.0, 141, 108,  91, 338, 164, 164, 164, 164,  44,  44,  61, 225,  0.0, 569,  0.0, 477, 348,  0.0,  0.0,  0.0,  17, 707, 835, 833, 862, 860, 858, 830, 495, 473, 471, 447,  0.0,  0.0,   0, 0, 0.0],
        "hf": [-45.947, -20.763, -20.763, -3.766, -3.766, 17.119, 17.119, 53.712, 69.939, 64.145, 82.528, 104.293, 197.322, 11.189, 27.016, -19.243, 9.404, 27.671, -181.422, 0.0, 0.0, -164.609, -182.329, -164.41, -129.158, -389.737, -359.258, -332.822, -163.569, -151.143, -129.488, -140.313, 0.0, -15.505, 3.32, 5.432, 23.101, 26.718, 54.929, 69.885, 20.079, 0.0, 134.062, 139.758, 0.0, 88.298, -396.242, 0.0, -73.568, -63.795, -57.795, 0.0, -82.921, 0.0, 0.0, -107.188, 0.0, -16.752, 0.0, -66.138, -59.142, -7.365, 0.0, 0.0, -8.253, 0.0, 0.0, 57.546, 1.834, 220.803, 227.368, 0.0, 0.0, -36.097, -161.74, 0.0, 0.0, -679.195, 0.0, 0.0, -313.545, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -258.96, 0.0, 0.0, -446.835, 0.0, 0.0, 0.0, -223.398, -203.188, -67.778, -182.005, -189.888, -46.562, 0.0, -344.125, 0.0, -2.084, 18.022, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "gf": [-8.03, 8.231, 8.231, 19.848, 19.848, 37.977, 37.977, 84.926, 92.9, 88.402, 93.745, 116.613, 221.308, 22.533, 30.485, 22.505, 41.228, 52.948, -158.589, 0.0, 0.0, -132.097, -131.366, -132.386, -107.858, -318.616, -291.188, -288.902, -105.767, -101.563, -92.099, -90.883, 0.0, 58.085, 63.051, 82.471, 95.888, 85.001, 128.602, 132.756, 68.861, 0.0, 199.958, 199.288, 0.0, 121.544, -349.439, 0.0, -33.373, -31.502, -25.261, 0.0, -35.814, 0.0, 0.0, -53.332, 0.0, -0.50, 0.0, 17.963, 18.088, 60.161, 0.0, 0.0, 16.731, 0.0, 0.0, 46.945, -1.721, 217.003, 216.328, 0.0, 0.0, -28.148, -144.549, 0.0, 0.0, -626.58, 0.0, 0.0, -281.495, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -209.337, 0.0, 0.0, -392.975, 0.0, 0.0, 0.0, 212.718, 136.742, 0.0, 0.0, -65.642, 0.0, 0.0, 241.373, 0.0, 30.222, 38.346, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "hv": [4.116, 4.65, 4.65, 2.771, 2.771, 1.284, 1.284, 6.714, 7.37, 6.797, 8.178, 9.342, 12.318, 4.098, 12.552, 9.776, 10.185, 8.834, 24.529, 0.0, 0.0, 40.246, 18.999, 20.041, 12.909, 22.709, 17.759, 0.0, 10.919, 7.478, 5.708, 11.227, 0.0, 14.599, 11.876, 14.452, 14.481, 0.0, 6.947, 6.918, 28.453, 0.0, 31.523, 31.005, 0.0, 23.34, 43.046, 0.0, 13.78, 11.985, 9.818, 0.0, 19.208, 17.574, 0.0, 0.0, 0.0, 11.883, 0.0, 30.644, 26.277, 0.0, 0.0, 0.0, 14.931, 0.0, 0.0, 14.364, 11.423, 7.751, 11.549, 0.0, 0.0, 4.877, 0.0, 0.0, 8.901, 1.86, 8.901, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.322, 0.0, 0.0, 8.301, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 51.787, 0.0, 0.0, 0.0, 0.0, 0.0, 16.921, 17.117, 13.265, 0.0, 0.0, 27.966, 0.0, 0.0],
        "txt": [("CH3-",),
                        ("CH2<",),
                        ("RCH2<",),
                        ("CH",),
                        (">RCH-",),
                        (">C<",),
                        (">RC<",),
                        ("CH2=CH",),
                        ("CH=CH",),
                        ("CH2=C",),
                        ("CH=C",),
                        ("C=C",),
                        ("CH2=C=CH",),
                        ("ACH",),
                        ("AC-",),
                        ("ACCH3",),
                        ("ACCH2",),
                        ("ACCH",),
                        ("OH",),
                        ("CH3OH",),
                        ("H2O",),
                        ("ACOH",),
                        ("CH3CO",),
                        ("CH2CO",),
                        ("CHO",),
                        ("CH3COO",),
                        ("CH2COO",),
                        ("HCOO",),
                        ("CH3O",),
                        ("CH2O",),
                        ("CH-O",),
                        ("FCH2O",),
                        ("CH3NH2",),
                        ("CH2NH2",),
                        ("CHNH2",),
                        ("CH3NH",),
                        ("CH2NH",),
                        ("CHNH",),
                        ("CH3-RN",),
                        ("CH2-RN",),
                        ("ACNH2",),
                        ("C5H5N",),
                        ("C5H4N",),
                        ("C5H3N",),
                        ("CH3CN",),
                        ("CH2CN",),
                        ("COOH",),
                        ("HCOOH",),
                        ("CH2CL",),
                        ("CHCL",),
                        ("CCL",),
                        ("CH2CL2",),
                        ("CHCL2",),
                        ("CCL2",),
                        ("CHCL3",),
                        ("CCL3",),
                        ("CCL4",),
                        ("ACCL",),
                        ("CH3NO2",),
                        ("CH2NO2",),
                        ("CHNO2",),
                        ("ACNO2",),
                        ("CS2",),
                        ("CH3SH",),
                        ("CH2SH",),
                        ("FURFURAL",),
                        ("<CH2OH>2",),
                        ("I",),
                        ("Br",),
                        ("CH===C",),
                        ("C===C",),
                        ("ME2SO",),
                        ("ACRY",),
                        ("CL<C=C>",),
                        ("ACF",),
                        ("DMF-1",),
                        ("DMF-2",),
                        ("CF3",),
                        ("CF2",),
                        ("CF",),
                        ("COO",),
                        ("SiH3",),
                        ("SiH2",),
                        ("SiH",),
                        ("Si",),
                        ("SiH2O",),
                        ("SiHO",),
                        ("SiO",),
                        ("TERT-N",),
                        ("CCL3F",),
                        ("CCL2F",),
                        ("HCCL2F",),
                        ("HCCLF",),
                        ("CCLF2",),
                        ("HCCLF2",),
                        ("CCLF3",),
                        ("CCL2F2",),
                        ("F (exceptions)",),
                        ("CONH2",),
                        ("CONHCH3",),
                        ("CONHCH2",),
                        ("CON<CH3>2",),
                        ("CONCH3CH2",),
                        ("CON<CH2>2",),
                        ("C2H5O2",),
                        ("C2H4O2",),
                        ("CH3S",),
                        ("CH2S",),
                        ("CHS",),
                        ("MORPH",),
                        ("C4H4S",),
                        ("C4H3S",),
                        ("C4H2S",),
                        ("NMP",)]
        }

    def isCalculable(self):
        """Método que estima si el método es calculable en función de los datos disponibles, definido por cada método"""
        if not self.kwargs["M"]:
            self.msg=QApplication.translate("pychemqt", "undefined molecular weight")
            self.status=0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        self.M=self.kwargs["M"]
        tc=Pc=vc=tb=hv=gf=hf=0
        for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tb+=contribucion*self.coeff["tb"][grupo]
            tc+=contribucion*self.coeff["tc"][grupo]
            Pc+=contribucion*self.coeff["Pc"][grupo]
            vc+=contribucion*self.coeff["vc"][grupo]
            hv+=contribucion*self.coeff["hv"][grupo]
            gf+=contribucion*self.coeff["gf"][grupo]
            hf+=contribucion*self.coeff["hf"][grupo]

        if self.kwargs["Tb"] :
            self.Tb=unidades.Temperature(self.kwargs["Tb"])
        else:
            self.Tb=unidades.Temperature(1000/(0.5+35.7/tb**0.5+1000/(142+tb)))
        self.Tc=unidades.Temperature(self.Tb*(1+(1.28*tc)**-1))
        self.Pc=unidades.Pressure(self.M/(0.346+Pc)**2, "bar")
        self.Vc=unidades.SpecificVolume((172+vc)/self.M, "ccg")
        self.Hv=unidades.Enthalpy((hv+6.829)/self.M, "kJg")
        self.Hf=unidades.Enthalpy((hf+10.835)/self.M, "kJg")
        self.Gf=unidades.Enthalpy((gf-14.828)/self.M, "kJg")

        GroupContribution.calculo(self)


class Ambrose(GroupContribution):
    """Ambrose, D., “Correlation and Estimation of Vapor-Liquid Critical Properties, II. Critical Pressures and Critical Volumes of OrganicCompounds,”National Physical Laboratory, Teddington, NPL  Report 98 (May 1979).
    ref, API procedure 4A1.1, pag.294
    grupos: grupos que forman la molécula
    contribuciones: contribuciones de cada grupo
    platt: número de Platt, The Platt number is the number of pairs of carbon atoms which are separated by three carbon-carbon bonds and is an indicator of the degree of branching in the molecule. The Platt number of an n-alkane is equal to the number of carbons minus three. Further discussion of the Platt number is given by Wiener, J . Am. Chem. Soc., 69, 17(1947).
    Tb: Temperatura de ebullición
    M: peso molecular
    SG: gravedad específica, opcional

    >>> desconocido=Ambrose(group=[0, 1, 2, 3], contribution=[5, 1, 1, 1], Tb=unidades.Temperature(229.72, "F"), M=114.23, platt=3)
    >>> print desconocido.Tc.F, desconocido.Pc.psi, desconocido.Vc.ft3lb
    555.826906339 400.749899167 0.0639229274271
    """
    coeff={
            "Pc": [0.2260, 0.2260, 0.22, 0.1960, 0.1935, 0.1935, 0.1875, 0.1610, 0.1410, 0.1410, 0.1820, 0.1820, 0.1820, 0.1820, 0.1495, 0.1495, 0.1170, 0.9240, 0.8940, 0.9440, 0.9440, 0.8640, 0.9140, 0.8340, 0.8840, 0.8840, 0.8040, 0.7240, 0.5150],
            "tc": [0.138, 0.138, 0.095, 0.018, 0.113, 0.113, 0.070, 0.088, 0.038, 0.038, 0.09, 0.09, 0.03, 0.09, 0.075, 0.075, 0.06, 0.458, 0.448, 0.488, 0.488, 0.438, 0.478, 0.428, 0.468, 0.468, 0.418, 0.368, 0.22],
            "vc": [55.1, 55.1, 47.1, 38.1, 45.1, 45.1, 37.1, 35.1, 35.1, 35.1, 44.5, 44.5, 44.5, 44.5, 37, 37, 29.5, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 148],
            "txt": [("CH3", {"C": 1, "H": 3}),
                    ("CH2", {"C": 1, "H": 2}),
                    ("CH", {"C": 1, "H": 1}),
                    ("C", {"C": 1}),
                    ("=CH2", {"C": 1, "H": 2}),
                    ("=CH-", {"C": 1, "H": 1}),
                    ("=C", {"C": 1}),
                    ("=C=", {"C": 1}),
                    ("≡CH", {"C": 1, "H": 1}),
                    ("≡C-", {"C": 1}),
                    ("-CH2- (Cyclic)", {"C": 1, "H": 2}),
                    ("-CH< (Ciclic)", {"C": 1, "H": 1}),
                    ("-CH< (in fused ring)", {"C": 1, "H": 1}),
                    (">C< (Ciclic)", {"C": 1}),
                    ("=CH- (Cyclic)", {"C": 1, "H": 1}),
                    ("=C< (Cyclic)", {"C": 1}),
                    ("=C= (Cyclic)", {"C": 1}),
                    ("Phenyl- ", {"C": 6, "H": 5}),
                    ("o-Phenyl- ", {"C": 6, "H": 4}),
                    ("m-Phenyl- ", {"C": 6, "H": 4}),
                    ("p-Phenyl- ", {"C": 6, "H": 4}),
                    ("1,2,3-Phenyl- ", {"C": 6, "H": 3}),
                    ("1,2,4-Phenyl- ", {"C": 6, "H": 3}),
                    ("1,2,3,4-Phenyl- ", {"C": 6, "H": 2}),
                    ("1,2,3,5-Phenyl- ", {"C": 6, "H": 2}),
                    ("1,2,4,5-Phenyl- ", {"C": 6, "H": 2}),
                    ("1,2,3,4,5-Phenyl- ", {"C": 6, "H": 1}),
                    ("1,2,4,5,6-Phenyl- ", {"C": 6}),
                    ("=CH-CH= (in fused Aromatic ring)", {"C": 2, "H": 2})] }

    FirstOrder=29

    def isCalculable(self):
        if not self.kwargs["M"]:
            self.msg=QApplication.translate("pychemqt", "undefined molecular weight")
            self.status=0
        elif not self.kwargs["Tb"]:
            self.msg=QApplication.translate("pychemqt", "undefined boiling point")
            self.status=0
        else:
            return GroupContribution.isCalculable(self)

    def calculo(self):
        self.Tb=unidades.Temperature(self.kwargs["Tb"])
        self.M=self.kwargs["M"]

        Pc=tc=vc=0
        for grupo, contribucion in zip(self.kwargs["group"], self.kwargs["contribution"]):
            tc+=contribucion*self.coeff["tc"][grupo]
            Pc+=contribucion*self.coeff["Pc"][grupo]
            vc+=contribucion*self.coeff["vc"][grupo]

        self.Tc=unidades.Temperature(self.Tb.R*(1+1/(1.242+tc-0.023*self.kwargs["platt"])), "R")
        self.Pc=unidades.Pressure(14.5*self.M/(0.339+Pc-0.026*self.kwargs["platt"])**2, "psi")
        self.Vc=unidades.SpecificVolume(0.01602*(40+vc)/self.M, "ft3lb")

        GroupContribution.calculo(self)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#    cresol=Componente(177)
#    print cresol.Tb, cresol.Tc
#    joback=Joback(group=[0, 1, 13, 14, 20], contribution=[1, 1, 4, 2, 1])
#    print joback.Tb, joback.Tc

#    butanol_1=Componente(160)
#    print butanol_1.f_acent, butanol_1.Tc
#    unknown=Constantinou_Gani(group=[0, 1, 15], contribution=[1, 3, 1])
#    print unknown.f_acent, unknown.Tc

#    ic5=Componente(7)
#    print ic5.Tb, ic5.Tc
#    elliot=Elliott(group=[0, 5], contribution=[4, 1], M=72)
#    print elliot.Tb, elliot.Tc

#    cresol=Componente(177)
#    print cresol.Tc, cresol.f_acent
#    cresol_wilson=Wilson_Jasperson(group=[3, 0, 5, 41], contribution=[8, 10, 1, 1], M=108.14, Tb=464.15, ring=1)
#    print cresol_wilson.Tc, cresol_wilson.f_acent

#    cresol=Componente(177)
#    print cresol.Tc, cresol.Pc.bar
#    cresol_marrero=Marrero_Pardillo(group=[1, 36, 129, 130, 132, 140, 148], contribution=[1, 1, 1, 2, 2, 1, 1], M=122.17, atomos=19)
#    print cresol_marrero.Tc, cresol_marrero.Pc.bar

#    trimetilpentano=Componente(541)
#    print trimetilpentano.Tc, trimetilpentano.Pc.psi, trimetilpentano.Vc.ft3lb
#    desconocido=Ambrose(group=[0, 1, 2, 3], contribution=[5, 1, 1, 1], Tb=unidades.Temperature(229.72, "F"), M=114.23, platt=3)
#    print desconocido.Tc.F, desconocido.Pc.psi, desconocido.Vc.ft3lb

#http://en.wikipedia.org/wiki/Joback_method
#    acetona=Componente(140)
#    print acetona.Tc, acetona.Pc.bar
#    joback_acetona=Joback(group=[0, 23], contribution=[2, 1])
#    print joback_acetona.Tc, joback_acetona.Pc.bar, joback_acetona.Tb, joback_acetona.Tf



#    etilbenceno=Componente(45)
#    t=unidades.Temperature(180, "F")
#    print "DIPPR: ", etilbenceno.Tension_DIPPR(t).dyncm
#    print "Paramétrica: ", etilbenceno.Tension_Parametrica(t).dyncm
#    print "Hakim: ", etilbenceno.Tension_Hakim(t).dyncm
#    print "Miller: ", etilbenceno.Tension_MIller(t).dyncm
#    print "Hydrocarbon: ", etilbenceno.Tension_Hydrocarbon(t).dyncm
#    print "Parachor: ", etilbenceno.Tension_Parachor(t, 285.1).dyncm
#    print "Miqueu: ", etilbenceno.Tension_Miqueu(t).dyncm
#    print "Block Bird: ", etilbenceno.Tension_Block_Bird(t).dyncm

#    ipentano=Componente(7)
#    t=unidades.Temperature(212, "F")
#    print "DIPPR: ", ipentano.ThCond_Gas_DIPPR(t).BtuhftF
#    print "Misic-Thodos: ", ipentano.ThCond_Gas_Misic_Thodos(t).BtuhftF


#    heptano=Componente(11)
#    t=unidades.Temperature(572, "F")
#    p=unidades.Pressure(1450, "psi")
#    print "Crooks: ", heptano.ThCond_Gas_Crooks(t, p.atm).BtuhftF

#    oxigeno=Componente(47)
#    t=unidades.Temperature(984.6, "R")
#    p=unidades.Pressure(6075, "psi")
#    print "Nonhidrocarbon: ", oxigeno.ThCond_Gas_Nonhidrocarbon(t, p.atm).BtuhftF

#    butilbenceno=Componente(78)
#    t=unidades.Temperature(140, "F")
#    print "DIPPR: ", butilbenceno.ThCond_Liquido_DIPPR(t).BtuhftF
#    print "Pachaiyappan: ", butilbenceno.ThCond_Liquido_Pachaiyappan(t).BtuhftF

#    heptano=Componente(11)
#    t=unidades.Temperature(320, "F")
#    print "Kanitkar Thodos: ", heptano.ThCond_Liquido_Kanitkar_Thodos(t, 197.4).BtuhftF
#    print "Lenoir: ", heptano.ThCond_Liquido_Lenoir(t, 197.4).BtuhftF

#    decano=Componente(14)
#    t=unidades.Temperature(104, "F")
#    print "DIPPR: ", decano.Mu_Liquido_DIPPR(t).cP
#    print "Paramétrico: ", decano.Mu_Liquido_Parametrica(t).cP
#    print "Letsou Steil: ", decano.Mu_Liquido_Letsou_Steil(t).cP
#
#    pentano=Componente(50)
#    t=unidades.Temperature(200, "F")
#    p=unidades.Pressure(3000, "psi")
#    print "Graboski Broun: ", pentano.Mu_Liquido_Graboski_Braun(t, p.atm).cP
#    print "Lucas: ", pentano.Mu_Liquido_Lucas(t, p.atm).cP

#    tetralin=Componente(376)
#    t=unidades.Temperature(302, "F")
#    print "DIPPR:  %0.5f" % tetralin.Pv_DIPPR(t).psi
#    print "Antoine:   %0.5f" % tetralin.Pv_Antoine(t).psi
#    print "Lee-Kesler:  %0.5f" % tetralin.Pv_Lee_Kesler(t).psi
#    print "Maxwell-Bonnel: %0.5f" % tetralin.Pv_Maxwell_Bonnel(t).psi
#    print "Wagner: %0.5f" % tetralin.Pv_Wagner(t).psi

#    propano=Componente(4)
#    t=unidades.Temperature(30, "F")
#    print "DIPPR: ", propano.RhoL_DIPPR(t).gml
#    print "Rackett: ", propano.RhoL_Rackett(t).gml
#    print "Cavett: ", propano.RhoL_Cavett(t).gml
#    print "Costald: ", propano.RhoL_Costald(t).gml

#    octano=Componente(12)
#    t=unidades.Temperature(212, "F")
#    p=unidades.Pressure(4410, "psi")
#    print "Thomson Brobst Hankinson: ", octano.RhoL_Thomson_Brobst_Hankinson(t, p.atm).kgl
#    print "API: ", octano.RhoL_API(t, p.atm).kgl


#    ciclohexano=Componente(38)      #ej pag 637
#    t=unidades.Temperature(300, "F")
#    p=unidades.Pressure(1000, "psi")
#    print ciclohexano.Z_SRK(t, p.atm)
#    print ciclohexano.Lee_Kesler_Entalpia(t, p.atm).Btulb
#    print ciclohexano.Entropia(t, p.atm).BtulbF*1.8

#    print ciclohexano.Hv_Lee_Kesler(422.04), ciclohexano.Calor_vaporizacion(422.04)
#    print ciclohexano.Cp_Lee_Kesler(422.04, 68.046), ciclohexano.Cv_Lee_kesler(422.04, 68.046)


#    isobutano=Componente(5)
#    t=unidades.Temperature(370, "F")
#    p=unidades.Pressure(4000, "psi")
#    print isobutano.Lee_Kesler_Fugacidad(t, p.atm).psi #Ej pag 745
#    t=unidades.Temperature(475, "F")
#    print isobutano.Lee_Kesler_Entropia(t, p.atm).BtulbF #Ej pag 733

#    print "     SRK    Lee_Kesler    BWRS"
#    print "Z  %5.4f   %7.4f   %5.4f" % (isobutano.Z_SRK(t, p.atm), isobutano.Z_Lee_Kesler(t, p.atm), isobutano.Z_BWRS(t, p.atm))
#    print isobutano.RhoG_Lee_Kesler(t, p.atm)
#    print isobutano.RhoG_SRK(t, p.atm)
#    print isobutano.RhoG_BWRS(t, p.atm)
#    print isobutano.Entalpia_SRK(t, p.atm)

#    buteno=Componente(24)
#    print buteno.f_acent
#    print buteno.factor_acentrico()


#    butano=Componente(6)
#    T=unidades.Temperature(200, "F")
#    print unidades.Enthalpy(butano.Entalpia_formacion(T)).Btulb



#    decano=Componente(14)
#    t=unidades.Temperature(104, "F")
#    print "DIPPR: ", decano.Mu_Liquido_DIPPR(t).cP
#    print "Paramétrico: ", decano.Mu_Liquido_Parametrica(t).cP
#    print "Letsou Steil: ", decano.Mu_Liquido_Letsou_Steil(t).cP

#    pentano=Componente(8)
#    t=unidades.Temperature(200, "F")
#    p=unidades.Pressure(3000, "psi")
#    print "Graboski Broun: ", pentano.Mu_Liquido_Graboski_Braun(t, p.atm).cP
#    print "Lucas: ", pentano.Mu_Liquido_Lucas(t, p.atm).cP

#    metilciclohexano=Componente(39)
#    t=unidades.Temperature(300, "K")
#    p=unidades.Pressure(500, "bar")
#    muo=unidades.Viscosity(0.68, "cP")
#    print "Graboski Broun: ", metilciclohexano.Mu_Liquido_Graboski_Braun(t, p.atm).cP
#    print "Lucas: ", metilciclohexano.Mu_Liquido_Lucas(t, p.atm).cP


#    decano=Componente(14)
#    print decano.MuL_Kouzel(unidades.Temperature(120, "F"), unidades.Pressure(9940, "psi").atm, unidades.Viscosity(52.7, "cP")).cP

#    propano=Componente(4)
#    t=unidades.Temperature(176, "F")
#    print propano.MuG_Thodos(t).cP

#    metano=Componente(2)
#    t=unidades.Temperature(543, "F")
#    print metano.Mu_Gas(t, 1).cP
#    print metano.Mu_Gas_Thodos(t).cP
#    print metano.Mu_Gas_Eakin_Ellingtong(t, 50).cP

#    nitrogeno=Componente(46)
#    t=unidades.Temperature(-58, "F")
#    p=unidades.Pressure(1677, "psi")
#    print nitrogeno.Mu_Gas(t, p).cP
#    print nitrogeno.Mu_Gas_Carr(t, p.atm).cP


#    from pylab import arange, plot, show
#    nonano=Componente(13)
#    t=linspace(0.3,1, 10)
#    C1=[]
#    C2=[]
#    C3=[]
#    C5=[]
#    C10=[]
#    C30=[]
#    for i in t:
#        C1.append(nonano.C_API(i, 1))
#        C2.append(nonano.C_API(i, 2))
#        C3.append(nonano.C_API(i, 3))
#        C5.append(nonano.C_API(i, 5))
#        C10.append(nonano.C_API(i, 10))
#        C30.append(nonano.C_API(i, 30))
#        #        C3.append(nonano.RhoL_API(i, nonano.Pc*3))
##        C5.append(nonano.RhoL_API(i, nonano.Pc*5))
##        C10.append(nonano.RhoL_API(i, nonano.Pc*10))
##        C30.append(nonano.RhoL_API(i, nonano.Pc*30))
#
#    plot(t, C1, t, C2)
#    show()


#    agua=Componente(62)
#    print agua.composicion_molecular
#    oxigeno=Componente(47)
#    print oxigeno.composicion_molecular
#    benceno=Componente(40)
#    print benceno.composicion_molecular
#    cfc=Componente(241)
#    print cfc.composicion_molecular

#    i_pentano=Componente(7)
#    t=unidades.Temperature(68, "F")
#    print i_pentano.Solubilidad_agua(t)
#    benceno=Componente(40)
#    t=unidades.Temperature(104, "F")
#    print benceno.Solubilidad_agua(t)

#    hexano=Componente(10)
#    t=unidades.Temperature(212, "F")
#    print hexano.Solubilidad_en_agua(t)

#    fluorene=Componente(197)
#    t=unidades.Temperature(122, "F")
#    print fluorene.Solubilidad_en_agua(t)

#    sulfuro=Componente(50)
#    t=unidades.Temperature(77, "F")
#    print sulfuro.Solubilidad_Henry(t, 1)

#    agua=Componente(62)
#    T=unidades.Temperature(100, "C")
#
#    print "Liquid Thermal Conductivity: ", agua.ThCond_Liquido(T, 1), "W/mK"
#    print "Liquid viscosity: ", agua.Mu_Liquido(T, 1), "Pa·s"
#    print "Liquid surface tension: ", agua.Tension(T), "N/m"
#    print "Gas Thermal Conductivity: ", agua.ThCond_Gas(T, 1), "W/mK"
#    print "Gas viscosity: ", agua.Mu_Gas(T, 1), "Pa·s"
#
#    print "Vapor pressure: ", agua.Pv(T).atm, "atm"

#    propeno=Componente(23)
#    t=unidades.Temperature(302, "F")
#    p=unidades.Pressure(2290, "psi")
#    print propeno.Cp_Lee_Kesler(t, p.atm).BtulbF
#    print propeno.Cp_Cv_ratio(t, p.atm)

#    SO2=Componente(51)
#    t=unidades.Temperature(300, "C")
#    print SO2.Mu_Gas_Chapman_Enskog(t, 1).microP

#    from pylab import arange, plot, show
#    nonano=Componente(13)
#    p=arange(0.2*nonano.Pc,5*nonano.Pc,1)
#    C1=[]
#    C2=[]
#    C3=[]
#    C5=[]
#    C10=[]
#    C30=[]#    for i in p:
#        C1.append(nonano.pr(i)*nonano.Lee_Kesler(nonano.Tc*1, i)[0]/1)
#        C11.append(nonano.Lee_Kesler(nonano.Tc*1.1, i))
#        C12.append(nonano.Lee_Kesler(nonano.Tc*1.2, i))
#        C13.append(nonano.Lee_Kesler(nonano.Tc*1.3, i))
#        C15.append(nonano.Lee_Kesler(nonano.Tc*1.5, i))
#        C17.append(nonano.Lee_Kesler(nonano.Tc*1.7, i))
#        C2.append(nonano.pr(i)*nonano.Lee_Kesler(nonano.Tc*2, i)[0]/2)
#        C25.append(nonano.Lee_Kesler(nonano.Tc*2.5, i))
#        C3.append(nonano.Lee_Kesler(nonano.Tc*3, i))
#        C4.append(nonano.Lee_Kesler(nonano.Tc*4, i))
#    plot(p/nonano.Pc, C1)
#    show()

#    Hidrogeno=Componente(1)
#    print unidades.Temperature(Hidrogeno.Tc).R
#    print unidades.Pressure(Hidrogeno.Pc, "atm").psi
#    print Hidrogeno.f_acent

#    agua=Componente(62)
#    print agua.SRK_Z(298.15, 1)
#    print agua.SRK_RhoG(298.15, 1).kgm3
#    print agua.SRK_Entalpia(298.15, 1).MJkg

#    agua=Componente(62)
#    t=400
#    print agua.BWRS_Z(298.15, 1)
#    print agua.van_Waals_Z(t, 1), agua.PR_Z(t, 1), agua.RK_Z(t, 1), agua.HPW_Z(t, 1, -0.5)
#    print agua.RK_Z(t, 1), agua.Wilson_Z(t, 1), agua.SRK_Z(t, 1)
#    print agua.BWRS_RhoG(298.15, 1).kgm3
#    print agua.BWRS_Entalpia(298.15, 1).MJkg

#    print agua.PR_V(298.15, 1)
#    print agua.PR_RhoG(298.15, 1)
#    print agua.PR_Entalpia(t, 1).MJkg, agua.Lee_Kesler_Entalpia(t, 1).MJkg, agua.iapws_Entalpia(t, 1).MJkg
#    print agua.Lee_Kesler_Z(t, 1), agua.SRK_Z(t, 1)

#    print agua.Cp_Gas_DIPPR(400), iapws_Cp(400, 1), agua.Cp_ideal(400)
#    print agua.Hv_Lee_Kesler(t).MJkg, agua.Hv_DIPPR(t).MJkg
#    print agua.Lee_Kesler_Entalpia(t, 1).MJkg, agua.iapws_Entalpia(t, 1).MJkg, agua.Entalpia_ideal(t).MJkg
#    print agua.Cp_Lee_Kesler(t, 1).JkgK*agua.M, agua.iapws_Cp(t, 1).JkgK*agua.M
#    print agua.Lee_Kesler_Entropia(t, 1).JkgK, agua.iapws_Entropia(t, 1).JkgK

#    agua=Componente(62)
#    from scipy import arange
#    from pylab import plot, grid, show
#    d=arange(270, 500, 10.)
#    y=[]
#    y2=[]
#    y3=[]
#    delta=[]
#    for i in d:
#        y.append(agua.Lee_Kesler_Entalpia(i, 1))
#        y2.append(agua.TB_Entalpia(i, 1))
#        y3.append(agua.iapws_Entalpia(i, 1))
##        delta.append(y3[-1]-y2[-1])
#    plot(d, y, d, y2, d, y3)
#    grid(True)
#    show()
#


#    sulfuro=Componente(50)
#    t=300
#    p=1
#    print sulfuro.H2S_V(t, p).ccg*sulfuro.M
#    print sulfuro.H2S_RhoG(t, p).gcc
#    print sulfuro.H2S_Z(t, p), sulfuro.TB_Z(t, p)
#    print sulfuro.H2S_Fugacidad(t, p)
#    print sulfuro.H2S_Entalpia(t, p).Jg*sulfuro.M

#    agua=Componente(62)
#    t=273
#    p=1
#    print agua.TB_Fugacidad(t, p), agua.Lee_Kesler_Fugacidad(t, p)
#    print agua.TB_U_exceso(t, p), agua.TB_H_exceso(t, p), agua.TB_S_exceso(t, p), agua.TB_Cv_exceso(t, p)
#    print agua.TB_Entalpia(t, p).MJkg, agua.Lee_Kesler_Entalpia(t, p).MJkg, agua.iapws_Entalpia(t, p).MJkg
#    print agua.TB_Joule_Thomson(t, p)

#    solido=Componente(533)
#    print solido.PT_lib(300)
#
#    Hexano=Componente(10)
#    print Hexano.Mu_Liquido(340, 1)




#    agua=Componente(62)
#    print [agua.Tension_Parametrica(t) for t in range(300, 350, 10)]
#    print agua.RhoL_Tait_Costald(300, 1)
#    print agua.Tc, agua.Pc.bar, agua.f_acent


