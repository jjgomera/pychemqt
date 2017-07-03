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
import os

from scipy import exp, cosh, sinh, log, log10, roots, absolute, sqrt
from scipy.optimize import fsolve
from scipy.constants import R, Avogadro

from lib.physics import R_atml, R_Btu, R_cal, factor_acentrico_octano
from lib import unidades, config, eos, sql


__doi__ = {
    "1":
        {"autor": "Lee, B. I. and Kesler, M. G.",
         "title": "A Generalized Thermodynamic Correlation Based on"
                  "Three-Parameter Corresponding States",
         "ref": "American Institute of Chemical Engineers Journal, 21, 1975",
         "doi": "10.1002/aic.690210313"},
    "2":
        {"autor": "Tarek Ahmed",
         "title": "Equations of State and PVT Analysis: Applications for"
                  "Improved Reservoir Modeling, 2nd Edition",
         "ref": "Gulf Professional Publishing, 2016, ISBN 9780128015704,",
         "doi": "10.1016/B978-0-12-801570-4.00002-7"},

}


def Pv_Lee_Kesler(T, Tc, Pc, w):
    """Calculates vapor pressure of a fluid using the Lee-Kesler correlation

    The vapor pressure is given by:

    .. math::
        \ln P_r = f^{(0)} + \omega f^{(1)}

        f^{(0)} = 5.92714-\frac{6.09648}{T_r}-1.28862\ln T_r + 0.169347T_r^6

        f^{(1)} = 15.2518-\frac{15.6875}{T_r} - 13.4721 \ln T_r + 0.43577T_r^6

    Parameters
    ----------
    T : float
        Temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure [Pa]
    w : float
        Acentric factor [-]

    Returns
    -------
    Pv : float
        Vapor pressure at T [Pa]

    Examples
    --------
    Example 1.2 from [2]_; propane at 80ºF

    >>> T = unidades.Temperature(80, "F")
    >>> Tc = unidades.Temperature(666.01, "R")
    >>> Pc = unidades.Pressure(616.3, "psi")
    >>> "%0.0f" % Pv_Lee_Kesler(T, Tc, Pc, 0.1522).psi
    '144'

    References
    ----------
    .. [1] Lee, B. I. and Kesler, M. G., A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States. American
       Institute of Chemical Engineers Journal, Vot. 21, 1975
    .. [2] Tarek Ahmed. Equations of State and PVT Analysis: Applications for
       Improved Reservoir Modeling, 2nd Edition. Gulf Professional Publishing,
       2016, ISBN 9780128015704
    """
    # Eq 17, pag 525
    Tr = T/Tc
    f0 = 5.92714 - 6.09648/Tr - 1.28862*log(Tr) + 0.169347*Tr**6
    f1 = 15.2518 - 15.6875/Tr - 13.4721*log(Tr) + 0.43577*Tr**6
    return unidades.Pressure(exp(f0 + w*f1)*Pc)


def f_acent_Lee_Kesler(Tb, Tc, Pc):
    """Calculates acentric factor of a fluid using the Lee-Kesler correlation

    Parameters
    ----------
    Tb : float
        Boiling temperature [K]
    Tc : float
        Critical temperature [K]
    Pc : float
        Critical pressure [Pa]

    Returns
    -------
    w : float
        Acentric factor [-]

    References
    ----------
    .. [1] Lee, B. I. and Kesler, M. G., A Generalized Thermodynamic
       Correlation Based on Three-Parameter Corresponding States. American
       Institute of Chemical Engineers Journal, Vot. 21, 1975
    """
    Tr = Tb/Tc
    Pr = 101325/Pc
    w = (log(Pr) - 5.92714 + 6.09648/Tr + 1.28862*log(Tr) - 0.169347*Tr**6)/(
        15.2518 - 15.6875/Tr - 13.4721*log(Tr) + 0.43577*Tr**6)

    return unidades.Dimensionless(w)


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
        return Pv_Lee_Kesler(T, self.Tc, self.Pc, self.f_acent)

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


if __name__ == '__main__':
    import doctest
    doctest.testmod()

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


