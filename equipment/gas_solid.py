#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Librería de definición de equipos de separación gas-sólido:
#     -Ciclones
#     -Cámaras de sedimentación
#     -Filtros de mangas (fabric filters, baghouse)
#     -Precipitadores electrostáticos
###############################################################################

import os

from PyQt5.QtWidgets import  QApplication
from scipy import exp, sqrt, ceil, roots
from scipy.constants import pi, g, e, epsilon_0
from scipy.optimize import fsolve

from lib.unidades import (Length, Pressure, DeltaP, Speed, Time, Area,
                          PotencialElectric, Currency, Dimensionless, MassFlow)
# from lib.datasheet import pdf
from lib.corriente import Corriente, Solid
from .parents import equipment, UI_equip


class Separador_SolidGas(equipment):
    """Clase generica con la funcionalidad comun de los equipos de separacion gas-solido"""

    def calcularRendimiento(self, rendimientos):
        rendimiento_global = 0
        for i, fraccion in enumerate(self.entrada.solido.fracciones):
            rendimiento_global += rendimientos[i]*fraccion
        return Dimensionless(rendimiento_global)

    def CalcularSalidas(self):
        Solido_NoCapturado, Solido_Capturado=self.entrada.solido.Separar(self.rendimiento_parcial)
        self.salida=[]
        if self.rendimiento==0:
            self.salida.append(self.entrada.clone(P=self.entrada.P-self.deltaP))
            self.salida.append(self.entrada.clone(split=0, P=self.entrada.P-self.deltaP))
        else:
            self.salida.append(self.entrada.clone(solido=Solido_NoCapturado, P=self.entrada.P-self.deltaP))
            self.salida.append(Corriente(P=self.entrada.P-self.deltaP, T=self.entrada.T, solido=Solido_Capturado))

        self.Pin=self.entrada.P
        self.Pout=self.salida[0].P
        self.Min=self.entrada.solido.caudal
        self.Dmin=self.entrada.solido.diametro_medio
        self.Mr=self.salida[0].solido.caudal
        self.Dmr=self.salida[0].solido.diametro_medio
        self.Ms=self.salida[1].solido.caudal
        self.Dms=self.salida[1].solido.diametro_medio


class Ciclon(Separador_SolidGas):
    """Clase que modela un ciclón

    Parámetros:
        Entrada: instancia de la clase corriente que define la corriente a tratar
        tipo_calculo:
            0   -   Evaluación del funcionamiento de un ciclón de dimensiones conocidas (diametro, forma, número de ciclones)
            1   -   Diseño de un ciclón a partir del rendimiento o perdida de carga admisibles
        modelo_rendimiento: Modelo de simulación del ciclón:
            0   -   Modelo Rosin-Rammler-Intelmann
            1   -   Modelo Leith-Licht
        modelo_DeltaP: Método de cálculo de la perdida de presión:
            0   -   Cálculo simplificado teniendo en cuenta únicamente las pérdidas cinéticas
            1   -   Método de Casal-Martinez-Benet
            2   -   Método de Leith-Licht
            3   -   Método de Sheferd, Lapple y Ter Linden
        modelo_ciclón: Define la forma del ciclón a partir del diametro del ciclón y una serie de modelo habituales
            0   -   Stairmand (Alta η)
            1   -   Swift (Alta η)
            2   -   Lapple (Baja η)
            3   -   Swift (Baja η)
            4   -   Peterson/Whitby (Baja η)
            5   -   Lorenz I
            6   -   Lorenz II
            7   -   Lorenz III
            8   -   Personalizado
        diametro: diametro del ciclón
        num_ciclones: Número de ciclones dispuestos en paralelo para tratar la corriente
        dimensiones: En el caso de que se trate de un modelo de ciclón personalizado, se pueden indicar las dimensiones principales en un array de 7 elementos
        rendimientoAdmisible: Rendimiento admisible del ciclón (en diseño)
        DeltaPAdmisible: Perdida de presión admisible del ciclon
        velocidadAdmisible: Velocidad de admisión del ciclón

    Coste
        f_install: factor instalación
        base: indice base
        actual: Indice actual
        tipo_costo:
            0   -   Heavy duty
            1   -   Standart duty
            2   -   Multiciclone

    >>> diametros=[17.5, 22.4, 26.2, 31.8, 37, 42.4, 48, 54, 60, 69, 81.3, 96.5, 109, 127]
    >>> fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    >>> corriente=Corriente(T=300, P=101325, caudalMasico=0.3,  fraccionMolar=[1.], caudalSolido=[0.001], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    >>> ciclon=Ciclon(entrada=corriente, tipo_calculo=1, rendimientoAdmisible=0.95, velocidadAdmisible=5)
    >>> print ciclon.C_instTotal, ciclon.C_adqTotal
    2260.1946585 1614.42475607
    """
    title = QApplication.translate("pychemqt", "Cyclone")
    help = os.environ["pychemqt"] + "doc/Ciclones.htm"
    kwargs = {"entrada": None,
              "tipo_calculo": 0,
              "modelo_rendimiento": 0,
              "modelo_DeltaP": 0,
              "modelo_ciclon": 0,
              "Dc": 0.0,
              "num_ciclones": 0,
              "dimensiones": [],
              "rendimientoAdmisible": 0.0,
              "DeltaPAdmisible": 0.0,
              "velocidadAdmisible": 0.0,

              "f_install": 1.4,
              "Base_index": 0.0,
              "Current_index": 0.0,
              "tipo_costo": 0}
    kwargsInput = ("entrada", )
    kwargsValue = ("Dc", "num_ciclones", "rendimientoAdmisible", "velocidadAdmisible", "DeltaPAdmisible")
    kwargsList = ("tipo_calculo", "modelo_rendimiento", "modelo_DeltaP", "modelo_ciclon", "tipo_costo")
    calculateValue = ("deltaP", "V", "rendimientoCalc", "NCalc", "Dcc", "Hc", "Bc", "Jc", "Lc", "Zc", "De", "Sc")
    calculateCostos = ("C_adq", "C_inst", "num_ciclonesCoste", "Q")
    indiceCostos = 2

    TEXT_TIPO = [QApplication.translate("pychemqt", "Rating: Calculate ΔP and efficiency", None),
                 QApplication.translate("pychemqt", "Design: Calculate cyclone dimensions")]
    TEXT_MODEL = ["Rossin, Rammler & Intelmann", "Leith & Licht"]
    TEXT_MODEL_DELTAP = [QApplication.translate("pychemqt", "Standart"),
                         "Casal & Martinez-Benet", "Leith & Licht", "Sheferd, Lapple & Ter Linden"]
    TEXT_MODEL_CICLON = ["Stairmand ("+QApplication.translate("pychemqt", "High η", None)+")",
                         "Swift ("+QApplication.translate("pychemqt", "High η", None)+")",
                         "Lapple ("+QApplication.translate("pychemqt", "Low η", None)+")",
                         "Swift ("+QApplication.translate("pychemqt", "Low η", None)+")",
                         "Peterson/Whitby ("+QApplication.translate("pychemqt", "Low η", None)+")",
                         "Lorenz I", "Lorenz II", "Lorenz III",
                         QApplication.translate("pychemqt", "Custom")]
    TEXT_COST = [QApplication.translate("pychemqt", "Heavy duty"),
                 QApplication.translate("pychemqt", "Standard dury"),
                 QApplication.translate("pychemqt", "Multicyclone")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and self.kwargs["Current_index"]:
            self.statusCoste = True
        else:
            self.statusCoste = False

        if not self.kwargs["entrada"]:
            self.msg = QApplication.translate("pychemqt", "undefined input")
            self.status = 0
            return

        if self.kwargs["tipo_calculo"]:
            if self.kwargs["modelo_ciclon"] == 8:
                if self.kwargs["Dc"] and self.kwargs["Hc"] and self.kwargs["Bc"] and self.kwargs["De"]:
                    self.msg = ""
                    self.status = 1
                    return True
                else:
                    self.msg = QApplication.translate("pychemqt", "undefined cyclone dimension")
                    self.status = 0
            else:
                if (self.kwargs["DeltaPAdmisible"] or self.kwargs["velocidadAdmisible"]) and self.kwargs["rendimientoAdmisible"]:
                    self.msg = ""
                    self.status = 1
                    return True
                elif self.kwargs["rendimientoAdmisible"]:
                    self.msg = QApplication.translate("pychemqt", "undefined efficiency")
                    self.status = 0
                else:
                    self.msg = QApplication.translate("pychemqt", "undefined loss pressure specification")
                    self.status = 0

        else:
            if self.kwargs["Dc"] and self.kwargs["num_ciclones"]:
                self.msg = ""
                self.status = 1
                return True
            elif self.kwargs["Dc"]:
                self.msg = QApplication.translate("pychemqt", "undefined cyclone number")
                self.status = 0
            else:
                self.msg = QApplication.translate("pychemqt", "undefined cyclone diameter")
                self.status = 0

    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.Dc=Length(self.kwargs["Dc"])
        self.num_ciclones=Dimensionless(self.kwargs["num_ciclones"])
        self.dimensiones=self.kwargs["dimensiones"]
        self.rendimientoAdmisible=self.kwargs["rendimientoAdmisible"]
        self.DeltaPAdmisible=Pressure(self.kwargs["DeltaPAdmisible"])

        if self.kwargs["tipo_calculo"]==0: #Calculo
            if self.kwargs["modelo_ciclon"]!=8:
                self.dimensiones=self.dimensionado(Dc=self.Dc)  #Hc, Bc, Jc, Lc, Zc, De, Sc,kf, G
            else:
                self.dimensiones=self.dimensionado(dimensiones=self.dimensiones)
            self.Hc=Length(self.dimensiones[0])
            self.Bc=Length(self.dimensiones[1])
            self.Jc=Length(self.dimensiones[2])
            self.Lc=Length(self.dimensiones[3])
            self.Zc=Length(self.dimensiones[4])
            self.De=Length(self.dimensiones[5])
            self.Sc=Length(self.dimensiones[6])
            self.kf=Length(self.dimensiones[7])
            self.G=Length(self.dimensiones[8])

            self.V=Speed(self.entrada.Q/(self.Bc*self.Hc*self.num_ciclones))
            if self.kwargs["modelo_rendimiento"]==0:
                self.N=11.3*(self.Hc*self.Bc/self.Sc**2)**2+3.33
                self.dc=Length(sqrt(9*self.Bc*self.entrada.Gas.mu/(2*pi*self.N*self.V*(self.entrada.solido.rho-self.entrada.Gas.rho))), "ParticleDiameter")
            else:
                self.N=self.V*(0.1079-0.00077*self.V+1.924e-6*self.V**2)

            self.rendimiento_parcial=self.calcularRendimientos_parciales()
            self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)
            self.DeltaP=self.PerdidaPresion()
        else:
            if self.DeltaPAdmisible and self.kwargs["velocidadAdmisible"]:
                V_fpresion=self.velocidad_f_presion()
                if self.kwargs["velocidadAdmisible"]>V_fpresion:
                    self.kwargs["velocidadAdmisible"]=V_fpresion
            elif self.DeltaPAdmisible:
                self.kwargs["velocidadAdmisible"]=self.velocidad_f_presion()

            def f(diametro):
                self.dimensiones=self.dimensionado(Dc=diametro)
                self.Dc=Length(diametro)
                self.Hc=Length(self.dimensiones[0])
                self.Bc=Length(self.dimensiones[1])
                self.Jc=Length(self.dimensiones[2])
                self.Lc=Length(self.dimensiones[3])
                self.Zc=Length(self.dimensiones[4])
                self.De=Length(self.dimensiones[5])
                self.Sc=Length(self.dimensiones[6])
                self.kf=Length(self.dimensiones[7])
                self.G=Length(self.dimensiones[8])
                self.num_ciclones=Dimensionless(ceil(self.entrada.Q/self.Bc/self.Hc/self.kwargs["velocidadAdmisible"]))
                self.V=Speed(self.entrada.Q/(self.Bc*self.Hc*self.num_ciclones))
                if self.kwargs["modelo_rendimiento"]==0:
                    self.N=Dimensionless(11.3*(self.Hc*self.Bc/self.Sc**2)**2+3.33)
                    self.dc=Length(sqrt(9*self.Bc*self.entrada.Gas.mu/(2*pi*self.N*self.V*(self.entrada.solido.rho-self.entrada.Gas.rho))), magnitud="ParticleDiameter")
                else:
                    self.N=Dimensionless(self.V*(0.1079-0.00077*self.V+1.924e-6*self.V**2))
                rendimiento_parcial=self.calcularRendimientos_parciales()
                rendimiento=self.calcularRendimiento(rendimiento_parcial)
                return rendimiento-self.kwargs["rendimientoAdmisible"]

            Dc0=1
            diametro=fsolve(f, Dc0)
            self.Dc=Length(diametro[0])
            self.dimensiones=self.dimensionado(self.Dc)
            self.Hc=Length(self.dimensiones[0])
            self.Bc=Length(self.dimensiones[1])
            self.Jc=Length(self.dimensiones[2])
            self.Lc=Length(self.dimensiones[3])
            self.Zc=Length(self.dimensiones[4])
            self.De=Length(self.dimensiones[5])
            self.Sc=Length(self.dimensiones[6])
            self.kf=Length(self.dimensiones[7])
            self.G=Length(self.dimensiones[8])
            self.rendimiento_parcial=self.calcularRendimientos_parciales()
            self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)
            self.deltaP=self.PerdidaPresion()

        self.num_ciclonesCoste=self.num_ciclones
        self.NCalc=self.num_ciclones
        self.Q=self.entrada.Q
        self.rendimientoCalc=self.rendimiento
        self.Dcc=self.Dc

        self.CalcularSalidas()


    def dimensionado(self, Dc=0, dimensiones=[]):
        coef=[[0.5, 0.2, 0.375, 1.5, 2.5, 0.5, 0.5, 6.4, 551.3], #Stairmand (Alta η)
                    [0.44, 0.21, 0.4, 1.4, 2.5, 0.4, 0.5, 9.2, 699.2], #Swift (Alta η)
                    [0.5, 0.25, 0.25, 2, 2, 0.5, 0.625, 8, 402.9], #Lapple (Baja η)
                    [0.5, 0.25, 0.4, 1.75, 2, 0.5, 0.6, 7.6, 381.8], #Swift (Baja η)
                    [0.583, 0.208, 0.5, 1.333, 1.837, 0.5, 0.583, 16*0.583*0.208/0.5/0.5, 0], #Peterson/Whitby (Baja η)
                    [0.533, 0.133, 0.333, 0.693, 1.887, 0.333, 0.733, 16*0.583*0.133/0.333/0.333, 0], #Lorenz I
                    [0.533, 0.133, 0.333, 0.693, 1.887, 0.233, 0.733, 16*0.583*0.133/0.233/0.233, 0], #Lorenz II
                    [0.4, 0.1, 0.333, 0.693, 1.887, 0.233, 0.733, 16*0.4*0.1/0.233/0.233, 0]] #Lorenz III

        if self.kwargs["modelo_ciclon"]==8:
            coef=dimensiones
            dimensiones.append(16*dimensiones[0]*dimensiones[1]/dimensiones[5]**2)
            dimensiones.append(0)
        else:
            coef=coef[self.kwargs["modelo_ciclon"]]

        Hc=Dc*coef[0]
        Bc=Dc*coef[1]
        Jc=Dc*coef[2]
        Lc=Dc*coef[3]
        Zc=Dc*coef[4]
        De=Dc*coef[5]
        Sc=Dc*coef[6]
        kf=coef[7]
        G=coef[8]

        return Hc, Bc, Jc, Lc, Zc, De, Sc, kf, G


    def calcularRendimientos_parciales(self):
        rendimiento_parcial=[]
        if self.kwargs["modelo_rendimiento"]:     #modelo Leith-Licht
            Vo=self.entrada.Q/self.num_ciclones
            n=1-(1-(12*self.Dc.ft)**0.14/2.5)*((self.entrada.T.F+460)/530)**0.3
            if self.G==0:
                Vs=self.entrada.solido.caudal/self.entrada.solido.rho/self.num_ciclones
                self.G=4*self.Dc*(2*Vs+Vo)/self.num_ciclones/self.Hc**2/self.Bc**2
            for i in range(len(self.entrada.solido.diametros)):
                t=self.entrada.solido.rho*(self.entrada.solido.diametros[i])**2/18/self.entrada.Gas.mu
                rendimiento_parcial.append(1-exp(-2*(self.G*t*Vo/self.Dc**3*(n+1))**(0.5/(n+1))))
        else:
            #modelo Rosin-Rammler-Intelmann
            for i in range(len(self.entrada.solido.diametros)):
                rendimiento_parcial.append((self.entrada.solido.diametros[i]/self.dc)**2/(1+(self.entrada.solido.diametros[i]/self.dc)**2))
        return rendimiento_parcial


    def PerdidaPresion(self):
        if self.kwargs["modelo_DeltaP"]==0:   # Cálculo simplificado teniendo en cuenta unicamente las pérdidas cinéticas
            DeltaP=Pressure(self.kf/2.*self.entrada.Gas.rho*self.V**2)
        elif self.kwargs["modelo_DeltaP"]==1:     # Casal-Martinez-Benet
            DeltaP=Pressure(0.003*self.N*self.entrada.Gas.rho*self.V**2, "inH2O")
        elif self.kwargs["modelo_DeltaP"]==2:     # Leith-Licht
            DeltaP=Pressure(0.003*self.entrada.Gas.rho*16*(self.entrada.Q/self.num_ciclones)**2/self.De**2/self.Bc/self.Hc, "mmHg")
        else:       # Sheferd, Lapple y Ter Linden
            Ae=self.Bc*self.Hc
            As=pi/4*self.De**2
            densidad_media=self.entrada.Gas.rho+self.entrada.solido.caudal/self.entrada.Gas.rho/(self.entrada.caudalmasico+self.entrada.solido.caudal/self.entrada.Gas.rho)*(self.entrada.solido.rho-self.entrada.Gas.rho)
            DeltaP=Pressure(1.078*(Ae/As)**1.21*self.entrada.Gas.rho*self.V**2, "mmH2O")
        return DeltaP

    def velocidad_f_presion(self):
        if self.kwargs["modelo_DeltaP"]==0:
            velocidad=sqrt(self.DeltaPAdmisible*2/self.kf/self.entrada.densidad_gas)
        elif self.kwargs["modelo_DeltaP"]==1:
            velocidad=sqrt(self.DeltaPAdmisible.inH2O/self.vueltas_supuestas/0.003/self.entrada.densidad_gas)
        elif self.kwargs["modelo_DeltaP"]==2:
            velocidad=sqrt(self.DeltaPAdmisible*self.De**2*self.Bc*self.Hc/0.003/self.entrada.densidad_gas/16)
        else:
            Ae=self.Bc*self.Hc
            As=pi/4*self.De**2
            densidad_media=self.entrada.Gas.rho+self.entrada.solido.caudal/self.entrada.Gas.rho/(self.entrada.caudal+self.entrada.solido.caudal/self.entrada.Gas.rho)*(self.entrada.solido.rho-self.entrada.Gas.rho)
            velocidad=sqrt(self.DeltaPAdmisible.mmH2O/1.078/(Ae/As)**1.21/densidad_media)
        return velocidad

    def coste(self):
        Q=self.entrada.Q.kft3min*self.entrada.Gas.rho/self.entrada.Gas.rho#Sd
        if self.num_ciclones!=1:
            Q=Q/self.num_ciclones

        if self.kwargs["tipo_costo"]==0:
            C=1.39*Q**0.98*1000     #Heavy duty
        elif self.kwargs["tipo_costo"]==1:
            C=0.65*Q**0.91*1000     #Standart duty
        else:
            C=1.56*Q**0.68*1000     #multiciclone

        self.C_adq=Currency(C * self.kwargs["Current_index"] / self.kwargs["Base_index"])
        self.C_inst=Currency(self.C_adq * self.kwargs["f_install"])
        self.C_adqTotal=Currency(self.C_adq*self.num_ciclones)
        self.C_instTotal=Currency(self.C_inst*self.num_ciclones)


    def Pdf(self):
        archivo=pdf("CICLÓN")
        archivo.ciclon("Ciclon limpieza polvo", self)
        archivo.dibujar()
        #Si queremos abrirlo con algún visor de pdf
        os.system("xpdf datasheet.pdf")


    def propTxt(self):
        txt=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_TIPO[self.kwargs["tipo_calculo"]].split(":")[0])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Model"), self.TEXT_MODEL[self.kwargs["modelo_rendimiento"]])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Pressure Loss Model"), self.TEXT_MODEL_DELTAP[self.kwargs["modelo_DeltaP"]])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Ciclon Model"), self.TEXT_MODEL_CICLON[self.kwargs["modelo_ciclon"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure Loss"), self.deltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Critic Particle Diameter"), self.dc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Internal Cycles"), self.N.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Speed"), self.V.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimiento.str)+os.linesep
        txt+="%-25s\t %i" %(QApplication.translate("pychemqt", "Cyclone number"), self.num_ciclones)+os.linesep

        txt+=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Cyclone Dimensions")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Ciclon Diameter"), self.Dc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Inlet Height"), self.Hc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Inlet Width"), self.Bc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Diameter"), self.Jc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Cylinder Cyclone Section Length"), self.Lc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Conical Cyclone Section Length"), self.Zc.str)+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Clean Gas Output Diameter"), self.De.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Clean Gas Inlet Orifice Length"), self.Sc.str)+os.linesep
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mass Flow"), self.entrada.solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mean Diameter"), self.entrada.solido.diametro_medio.str)+os.linesep
        if len(self.entrada.solido.diametros)>=1:
            txt+="%10s %10s %10s %10s %10s" %("Dp(µm)", "Xi", "ηi", "Xis", "Xig")+os.linesep
            for i in range(len(self.rendimiento_parcial)):
                txt+="%10.1f %10.2f %10.3f %10.3f %10.3f" % (self.entrada.solido.diametros[i].config("ParticleDiameter"), self.entrada.solido.fracciones[i],  self.rendimiento_parcial[i], self.salida[1].solido.fracciones[i], self.salida[0].solido.fracciones[i])+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), self.salida[0].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), self.salida[0].solido.diametro_medio.str)+os.linesep
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mass Flow"), self.salida[1].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mean Diameter"), self.salida[1].solido.diametro_medio.str)+os.linesep

        if self.statusCoste:
            txt+=os.linesep
            txt+="#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Model"), self.TEXT_COST[self.kwargs["tipo_costo"]])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt


    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Mode"), ("TEXT_TIPO", "tipo_calculo"), str),
                (QApplication.translate("pychemqt", "Model"), ("TEXT_MODEL", "modelo_rendimiento"), str),
                (QApplication.translate("pychemqt", "Pressure Loss Model"), ("TEXT_MODEL_DELTAP", "modelo_DeltaP"), str),
                (QApplication.translate("pychemqt", "Ciclon Model"), ("TEXT_MODEL_CICLON", "modelo_ciclon"), str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Critic Particle Diameter"), "dc", Length),
                (QApplication.translate("pychemqt", "Gas Internal Cycles"), "N", Dimensionless),
                (QApplication.translate("pychemqt", "Gas Speed"), "V", Speed),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimiento", Dimensionless),
                (QApplication.translate("pychemqt", "Cyclone number"), "num_ciclones", Dimensionless),
                (QApplication.translate("pychemqt", "Ciclon Diameter"), "Dc", Length),
                (QApplication.translate("pychemqt", "Inlet Height"), "Hc", Length),
                (QApplication.translate("pychemqt", "Inlet Width"), "Bc", Length),
                (QApplication.translate("pychemqt", "Solid Output Diameter"), "Jc", Length),
                (QApplication.translate("pychemqt", "Cylinder Cyclone Section Length"), "Lc", Length),
                (QApplication.translate("pychemqt", "Conical Cyclone Section Length"), "Zc", Length),
                (QApplication.translate("pychemqt", "Clean Gas Output Diameter"), "De", Length),
                (QApplication.translate("pychemqt", "Clean Gas Inlet Orifice Length"), "Sc", Length),
                (QApplication.translate("pychemqt", "Input Solid Mass Flow"), "Min", MassFlow),
                (QApplication.translate("pychemqt", "Input Solid Mean Diameter"), "Dmin", Length),
                (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), "Mr", MassFlow),
                (QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
                (QApplication.translate("pychemqt", "Solid Output Mass Flow"), "Ms", MassFlow),
                (QApplication.translate("pychemqt", "Solid Output Mean Diameter"), "Dms", Length),
                (QApplication.translate("pychemqt", "Cost Mode"), ("TEXT_COST", "tipo_costo"), str),
                (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", Currency),
                (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", Currency)]
        return list


class GravityChamber(Separador_SolidGas):
    """Clase que define la cámara de sedimentación por gravedad

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la cámara
        metodo: Integer que indica el tipo de cálculo a realizar:
            0   -   Calculo del funcionamiento conociendo la geometria y corriente de entrada, calcula el rendimiento
            1   -   Diseño del equipo conocida la entrada y el rendimiento admisible para obtener la geometría
        modelo: Integer que indica el tipo de modelo de cálculo a aplicar:
            0   -   Modelo de flujo piston sin mezcla vertical
            1   -   Modelo de mezcla perfecta
        W, H, L: dimensiones de la cámara de sedimentación
        rendimientoAdmisible: rendimiento requerido para el equipo, necesario en el caso de que estemos diseñando el equipo
        velocidadAdmisible: velocidad admisible del gas en la cámara, si no se indica se considera el valor de 1 m/s
        deltaP: Perdida de presion

    >>> diametros=[17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, 54e-4, 60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
    >>> fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    >>> solido=Solid(caudalSolido=[12426.28], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    >>> entrada=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    >>> camara=GravityChamber(entrada=entrada, modelo=1, H=1, W=3, L=10)
    >>> print camara.Vgas, camara.rendimiento
    0.283284753885 0.124577139186
    """
    title=QApplication.translate("pychemqt", "Gravity settling chamber")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "modelo": 0,
                    "W": 0.0,
                    "H": 0.0,
                    "L": 0.0,
                    "rendimientoAdmisible": 0.0,
                    "velocidadAdmisible": 0.0,
                    "deltaP": 0.0}

    kwargsInput=("entrada", )
    kwargsValue=("W", "H", "L", "rendimientoAdmisible", "velocidadAdmisible", "deltaP")
    kwargsList=("metodo", "modelo")
    calculateValue=("Q", "LCalc", "WCalc", "HCalc", "Vgas",  "rendimiento")

    TEXT_TIPO=[QApplication.translate("pychemqt", "Rating: Known dimensions, calculate efficiency"),
                        QApplication.translate("pychemqt", "Design: Calculate dimensions to fit a requerid efficiency")]
    TEXT_MODEL=[QApplication.translate("pychemqt", "Plug flow without vertical mix"),
                                    QApplication.translate("pychemqt", "Vertical mix")]

    @property
    def isCalculable(self):
        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
            return

        if self.kwargs["metodo"]:
            if self.kwargs["rendimientoAdmisible"]:
                self.msg=""
                self.status=1
                return True
            else:
                self.msg=QApplication.translate("pychemqt", "undefined efficiency")
                self.status=0
        else:
            if not self.kwargs["W"]:
                self.msg=QApplication.translate("pychemqt", "undefined width")
                self.status=0
            elif self.kwargs["W"]  and self.kwargs["H"] and self.kwargs["L"]:
                self.msg=""
                self.status=1
                return True
            elif not self.kwargs["H"]:
                self.msg=QApplication.translate("pychemqt", "height undefined, using default")
                self.status=3
                return True
            else:
                self.msg=""
                self.status=1
                return True

    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.W=Length(self.kwargs["W"])
        self.L=Length(self.kwargs["L"])
        self.deltaP=DeltaP(self.kwargs["deltaP"])
        if self.kwargs["H"]:
            self.H=Length(self.kwargs["H"])
        else:
            self.H=Length(1)
        if self.kwargs["velocidadAdmisible"]:
            velocidadAdmisible=self.kwargs["velocidadAdmisible"]
        else:
            velocidadAdmisible=1.

        if self.kwargs["metodo"]==0: #Calculo
            self.Vgas=Speed(self.entrada.Q/self.H/self.W)
            self.rendimiento_parcial=self.calcularRendimientos_parciales(self.L)
            self.rendimiento=Dimensionless(self.calcularRendimiento(self.rendimiento_parcial))
        else: #Diseño
            self.Vgas=Speed(velocidadAdmisible)
            self.W=Length(self.entrada.Q/velocidadAdmisible/self.H)

            def f(longitud):
                rendimientos=self.calcularRendimientos_parciales(longitud)
                return self.calcularRendimiento(self.calcularRendimientos_parciales(longitud))-self.kwargs["rendimientoAdmisible"]

            if f(0.1)>self.kwargs["rendimientoAdmisible"]:
                longitud=1
            else:
                longitud_=0.1
                longitud=fsolve(f, longitud_)
            self.L=Length(float(longitud))
            self.rendimiento_parcial=self.calcularRendimientos_parciales(self.L)
            self.rendimiento=Dimensionless(self.calcularRendimiento(self.rendimiento_parcial))

        self.LCalc=self.L
        self.WCalc=self.W
        self.HCalc=self.H
        self.Q=self.entrada.Q

        self.CalcularSalidas()

    def calcularRendimientos_parciales(self, L):
        """Método que calcula el rendimiento de separación, tanto global, como desglosado según la distribución de sólidos"""
        rendimientos=[]
        for diametro in self.entrada.solido.diametros:
            Ar=diametro**3*self.entrada.Gas.rho*(self.entrada.solido.rho-self.entrada.Gas.rho)*g/self.entrada.Gas.mu**2
            Vt=self.entrada.Gas.mu/diametro*self.entrada.Gas.rho*((14.42+1.827*Ar**0.5)**0.5-3.798)**2
            if self.kwargs["modelo"]==0:
                r=Vt*L/self.Vgas/self.H
            else:
                r=1-exp(-Vt*L/self.Vgas/self.H)
            if r>1:
                rendimientos.append(1)
            else:
                rendimientos.append(r)
        return rendimientos


    def propTxt(self):
        txt=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_TIPO[self.kwargs["metodo"]].split(":")[0])+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Model"), self.TEXT_MODEL[self.kwargs["modelo"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure Loss"), self.deltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Height"), self.HCalc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Width"), self.WCalc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Length"), self.LCalc.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Speed"), self.Vgas.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimiento.str)+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mass Flow"), self.entrada.solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mean Diameter"), self.entrada.solido.diametro_medio.str)+os.linesep
        if len(self.entrada.solido.diametros)>=1:
            txt+="%10s %10s %10s %10s %10s" %("Dp(µm)", "Xi", "ηi", "Xis", "Xig")+os.linesep
            for i in range(len(self.rendimiento_parcial)):
                txt+="%10.1f %10.2f %10.3f %10.3f %10.3f" % (self.entrada.solido.diametros[i].config("ParticleDiameter"), self.entrada.solido.fracciones[i],  self.rendimiento_parcial[i], self.salida[1].solido.fracciones[i], self.salida[0].solido.fracciones[i])+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), self.salida[0].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), self.salida[0].solido.diametro_medio.str)+os.linesep
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mass Flow"), self.salida[1].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mean Diameter"), self.salida[1].solido.diametro_medio.str)+os.linesep

        return txt


    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Mode"), ("TEXT_TIPO", "metodo"), str),
                (QApplication.translate("pychemqt", "Model"), ("TEXT_MODEL", "modelo"), str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Height"), "HCalc", Length),
                (QApplication.translate("pychemqt", "Width"), "WCalc", Length),
                (QApplication.translate("pychemqt", "Length"), "LCalc", Length),
                (QApplication.translate("pychemqt", "Gas Speed"), "Vgas", Speed),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimiento", Dimensionless),
                (QApplication.translate("pychemqt", "Input Solid Mass Flow"), "Min", MassFlow),
                (QApplication.translate("pychemqt", "Input Solid Mean Diameter"), "Dmin", Length),
                (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), "Mr", MassFlow),
                (QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
                (QApplication.translate("pychemqt", "Solid Output Mass Flow"), "Ms", MassFlow),
                (QApplication.translate("pychemqt", "Solid Output Mean Diameter"), "Dms", Length)]
        return list



class Baghouse(Separador_SolidGas):
    """Clase que define los filtros de mangas

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la cámara
        metodo: Integer que indica el tipo de cálculo a realizar:
            0   -   Calcular caida de presión
            1   -   Calcular el tiempo de filtración
            2   -   Calcular el número de filtros necesarios
        num_filtros: número de filtros
        tiempo: tiempo de filtración que lleva operando
        deltaP: perdida de presión a traves de la membrana
        resistenciaFiltro: coeficiente de perdida de presión del filtros, (in water)/(cP)(ft/min)
        resistenciaTorta: coeficiente de perdida de presión de la torta,  (in water)/(cP)(gr/ft2)(ft/min)
        limpieza: indica el número de filtros que están en limpieza
        membranasFiltro: número de membranas que tiene cada filtro
        diametroMembrana: diametro  de cada membrana
        areaMembrana: area de filtración que tiene cada membrana
        rendimientos: array con los rendimientos de fabrica de la membrana

    >>> diametros=[17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, 54e-4, 60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
    >>> fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    >>> solido=Solid(caudalSolido=[12426.28], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    >>> corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    >>> filtro=Baghouse(entrada=corriente, metodo=1, num_filtros=4, tiempo=3600, deltaP=0.1)
    >>> print filtro.floorArea, filtro.Vgas.ftmin
    7.24643712 0.480966666862
    """
    title=QApplication.translate("pychemqt", "Baghouse")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "num_filtros": 0,
                    "tiempo": 0.0,
                    "deltaP": 0.0,
                    "resistenciaFiltro": 0.0,
                    "resistenciaTorta": 0.0,
                    "limpieza": 0,
                    "membranasFiltro": 0,
                    "diametroMembrana": 0.0,
                    "areaMembrana": 0.0,
                    "rendimientos": []}
    kwargsInput=("entrada", )
    kwargsValue=("num_filtros", "tiempo", "deltaP", "resistenciaFiltro", "resistenciaTorta", "limpieza", "membranasFiltro", "diametroMembrana", "areaMembrana")
    kwargsList=("metodo", )
    calculateValue=("floorArea", "rendimiento", "Vgas", "num_filtrosCalc", "tiempoCalc", "deltaPCalc")

    TEXT_TIPO=[QApplication.translate("pychemqt", "Calculate Pressure drop"),
                        QApplication.translate("pychemqt", "Calculate time of filtration"),
                        QApplication.translate("pychemqt", "Calculate number of cells")]


    @property
    def isCalculable(self):
        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
            return

        if self.kwargs["metodo"]==0 and (not self.kwargs["num_filtros"] or not self.kwargs["tiempo"]):
            self.msg=QApplication.translate("pychemqt", "undefined values")
            self.status=0
            return
        elif self.kwargs["metodo"]==1 and (not self.kwargs["num_filtros"] or not self.kwargs["deltaP"]):
            self.msg=QApplication.translate("pychemqt", "undefined values")
            self.status=0
            return
        elif self.kwargs["metodo"]==2 and (not self.kwargs["tiempo"] or not self.kwargs["deltaP"]):
            self.msg=QApplication.translate("pychemqt", "undefined values")
            self.status=0
            return

        if self.kwargs["metodo"]==2 and self.kwargs["limpieza"]>=self.kwargs["num_filtros"]:
            self.msg=QApplication.translate("pychemqt", "All filters cleaned")
            self.status=0
            return

        if not self.kwargs["rendimientos"]:
            self.msg=QApplication.translate("pychemqt", "using default efficiency")
            self.status=3
            return True
        if len(self.kwargs["rendimientos"])!=self.kwargs["entrada"].solido.diametros:
            self.msg=QApplication.translate("pychemqt", "using default efficiency")
            self.status=3
            return True

        self.msg=""
        self.status=1
        return True


    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.metodo=self.kwargs["metodo"]
        self.num_filtros=self.kwargs["num_filtros"]
        self.tiempo=Time(self.kwargs["tiempo"])
        self.deltaP=Pressure(self.kwargs["deltaP"])

        if self.kwargs["resistenciaFiltro"]:
            self.resistenciaFiltro=Dimensionless(self.kwargs["resistenciaFiltro"])
        else:
            self.resistenciaFiltro=Dimensionless(0.84)
        if self.kwargs["resistenciaTorta"]:
            self.resistenciaTorta=Dimensionless(self.kwargs["resistenciaTorta"])
        else:
            self.resistenciaTorta=Dimensionless(0.1)
        if self.kwargs["limpieza"]:
            self.limpieza=self.kwargs["limpieza"]
        else:
            self.limpieza=1
        if self.kwargs["membranasFiltro"]:
            self.membranasFiltro=self.kwargs["membranasFiltro"]
        else:
            self.membranasFiltro=78

        if self.kwargs["diametroMembrana"]:
            self.diametroMembrana=Length(self.kwargs["diametroMembrana"])
        else:
            self.diametroMembrana=Length(0.5, "ft")

        if self.kwargs["areaMembrana"]:
            self.areaMembrana=Area(self.kwargs["areaMembrana"])
        else:
            self.areaMembrana=Area(16, "ft2")

        if any(self.kwargs["rendimientos"]):
            self.rendimiento_parcial=self.kwargs["rendimientos"]
        else:
            self.rendimiento_parcial=self.defaultRendimiento()
        self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)

        if self.kwargs["metodo"]==0:
            self.area=Area(self.areaMembrana*self.membranasFiltro*(self.num_filtros-self.limpieza))
            ms=self.entrada.solido.caudal.gmin*self.tiempo.min/self.area.ft2
            self.Vgas=Speed(self.entrada.Q/self.area)
            self.deltaP=Pressure(self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin+self.resistenciaTorta*self.entrada.Gas.mu.cP*ms*self.Vgas.ftmin, "inH2O")
        elif self.kwargs["metodo"]==1:
            self.area=Area(self.areaMembrana*self.membranasFiltro*(self.num_filtros-self.limpieza))
            self.Vgas=Speed(self.entrada.Q/self.area)
            deltaPMinimo=Pressure(self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin, "inH2O")
            if deltaPMinimo>self.deltaP:
                self.deltaP=deltaPMinimo
                self.tiempo=Time(0)
            else:
                ms=(self.deltaP.inH2O-self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin)/(self.resistenciaTorta*self.entrada.Gas.mu.cP*self.Vgas.ftmin)
                self.tiempo=Time(ms*self.area.ft2/self.entrada.solido.caudal.gmin, "min")
        else:
            N=roots([self.deltaP.inH2O, -self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.entrada.Q.ft3min/self.areaMembrana.ft2/self.membranasFiltro,-self.resistenciaTorta*self.entrada.Gas.mu.cP*self.entrada.Q.ft3min*self.entrada.solido.caudal.gmin*self.tiempo.min/self.areaMembrana.ft2**2/self.membranasFiltro**2])
            if N[0]>0:
                self.num_filtros=ceil(N[0])+self.limpieza
            else:
                self.num_filtros=ceil(N[1])+self.limpieza
            self.area=Area(self.areaMembrana*self.membranasFiltro*(self.num_filtros-self.limpieza))
            ms=self.entrada.solido.caudal.gmin*self.tiempo.min/self.area.ft2
            self.Vgas=Speed(self.entrada.Q/self.area)
            self.deltaP=Pressure(self.resistenciaFiltro*self.entrada.Gas.mu.cP*self.Vgas.ftmin+self.resistenciaTorta*self.entrada.Gas.mu.cP*ms*self.Vgas.ftmin, "inH2O")

        self.floorArea=Area(self.num_filtros*self.membranasFiltro*self.diametroMembrana**2)
        self.num_filtrosCalc=self.num_filtros
        self.tiempoCalc=self.tiempo
        self.deltaPCalc=self.deltaP

        self.CalcularSalidas()


    def defaultRendimiento(self):
        """Modelo de rendimiento sylvan, se usara si no se indica el rendimiento"""
        rendimiento=[(d/0.3177e-6)/(1+(d/0.3177e-6)) for d in self.entrada.solido.diametros]
        return rendimiento


    def propTxt(self):
        txt=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_TIPO[self.kwargs["metodo"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure Loss"), self.deltaP.str)+os.linesep
        txt+="%-25s\t %i" %(QApplication.translate("pychemqt", "Filter Number"), self.num_filtros)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Operation Time"), self.tiempo.str)

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Cloth resistence"), self.resistenciaFiltro.str)
        if not self.kwargs["resistenciaFiltro"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Cake resistence"), self.resistenciaTorta.str)
        if not self.kwargs["resistenciaTorta"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t %i" %(QApplication.translate("pychemqt", "Cells cleaned"), self.limpieza)
        if not self.kwargs["limpieza"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t %i" %(QApplication.translate("pychemqt", "Bags per cell"), self.membranasFiltro)
        if not self.kwargs["membranasFiltro"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Bag diameter"), self.diametroMembrana.str)
        if not self.kwargs["diametroMembrana"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Area per bag"), self.areaMembrana.str)
        if not self.kwargs["areaMembrana"]:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Speed"), self.Vgas.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Surface"), self.floorArea.str)+os.linesep

        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimiento.str)+os.linesep
        if self.entrada.solido.diametros:
            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mass Flow"), self.entrada.solido.caudal.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mean Diameter"), self.entrada.solido.diametro_medio.str)+os.linesep
            if len(self.entrada.solido.diametros)>=1:
                txt+=os.linesep+"#"+QApplication.translate("pychemqt", "Partial Efficiency")
                if not any(self.kwargs["rendimientos"]):
                    txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")

                txt+=os.linesep+"%10s %8s %10s %10s %10s" %("Dp("+Length.text("ParticleDiameter")+")", "Xi", "ŋi", "Xis", "Xig")+os.linesep
                for i in range(len(self.rendimiento_parcial)):
                    txt+="%10.4f %10.4f %10.4f %10.4f %10.4f" % (self.entrada.solido.diametros[i].config("ParticleDiameter"), self.entrada.solido.fracciones[i],  self.rendimiento_parcial[i], self.salida[1].solido.fracciones[i], self.salida[0].solido.fracciones[i])+os.linesep

            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), self.salida[0].solido.caudal.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), self.salida[0].solido.diametro_medio.str)+os.linesep
            txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mass Flow"), self.salida[1].solido.caudal.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mean Diameter"), self.salida[1].solido.diametro_medio.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Mode"), ("TEXT_TIPO", "metodo"), str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Filter Number"), "num_filtros", int),
                (QApplication.translate("pychemqt", "Operation Time"), "tiempo", Time),
                (QApplication.translate("pychemqt", "Cloth resistence"), "resistenciaFiltro", Dimensionless),
                (QApplication.translate("pychemqt", "Cake resistence"), "resistenciaTorta", Dimensionless),
                (QApplication.translate("pychemqt", "Cells cleaned"), "limpieza", int),
                (QApplication.translate("pychemqt", "Bags per cell"), "membranasFiltro", int),
                (QApplication.translate("pychemqt", "Bag diameter"), "diametroMembrana", Length),
                (QApplication.translate("pychemqt", "Area per bag"), "areaMembrana", int),
                (QApplication.translate("pychemqt", "Speed"), "Vgas", Speed),
                (QApplication.translate("pychemqt", "Surface"), "floorArea", Area),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimiento", Dimensionless),
                (QApplication.translate("pychemqt", "Input Solid Mass Flow"), "Min", MassFlow),
                (QApplication.translate("pychemqt", "Input Solid Mean Diameter"), "Dmin", Length),
                (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), "Mr", MassFlow),
                (QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
                (QApplication.translate("pychemqt", "Solid Output Mass Flow"), "Ms", MassFlow),
                (QApplication.translate("pychemqt", "Solid Output Mean Diameter"), "Dms", Length)]
        return list


class ElectricPrecipitator(Separador_SolidGas):
    """ Clase que define el precipitador electrostático

    Parámetros:
        entrada: Instancia de clase corriente que define el gas que fluye por la cámara
        metodo: indice que indica el tipo de cálculo a realizar
            0   -   Evaluación, conocidas las dimensiones, calcular el rendimiento del equipo
            1   -   Diseño, imponiendo el rendimiento a obtener calcular el area y potencia del carga
        potencialCarga: potencial de carga del equipo en V/m
        potencialDescarga: potencial de descarga del equipo en V/m
        area: area de recolección del equipo, m2
        epsilon: constante dieléctrica relativa del material
        rendimientoAdmisible: rendimiento exigido al equipo (modo diseño)
        deltaP: parámetro opcional que indica la perdida de presión del equipo

    >>> diametros=[17.5e-4, 22.4e-4, 26.2e-4, 31.8e-4, 37e-4, 42.4e-4, 48e-4, 54e-4, 60e-4, 69e-4, 81.3e-4, 96.5e-4, 109e-4, 127e-4]
    >>> fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    >>> solido=Solid(caudalSolido=[12426.28], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    >>> corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)
    >>> precipitador=ElectricPrecipitator(entrada=corriente, metodo=1, rendimientoAdmisible=0.9)
    >>> print precipitador.area, precipitador.rendimiento
    2389.41624164 0.9
    """
    title=QApplication.translate("pychemqt", "Electrostatic precipitator")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "potencialCarga": 0.0,
                    "potencialDescarga": 0.0,
                    "area": 0.0,
                    "epsilon": 0.0,
                    "rendimientoAdmisible": 0.0,
                    "deltaP": 0.0}
    kwargsInput=("entrada", )
    kwargsValue=("potencialCarga", "potencialDescarga", "area", "rendimientoAdmisible", "epsilon", "deltaP")
    kwargsList=("metodo", )
    calculateValue=("areaCalculada", "rendimiento")

    TEXT_TIPO=[QApplication.translate("pychemqt", "Rating: Calculate efficiency"),
                        QApplication.translate("pychemqt", "Design: Calculate dimensions to fit a requerid efficiency")]


    @property
    def isCalculable(self):
        if not self.kwargs["entrada"]:
            self.msg=QApplication.translate("pychemqt", "undefined input")
            self.status=0
            return

        if self.kwargs["metodo"]==0 and not self.kwargs["area"]:
            self.msg=QApplication.translate("pychemqt", "undefined area")
            self.status=0
            return
        elif self.kwargs["metodo"]==1 and not self.kwargs["rendimientoAdmisible"]:
            self.msg=QApplication.translate("pychemqt", "undefined efficiency")
            self.status=0
            return


        if not self.kwargs["potencialCarga"]:
            self.msg=QApplication.translate("pychemqt", "using default charging field")
            self.status=3
            return True
        if not self.kwargs["potencialDescarga"]:
            self.msg=QApplication.translate("pychemqt", "using default collecting field")
            self.status=3
            return True
        if not self.kwargs["epsilon"]:
            self.msg=QApplication.translate("pychemqt", "using default dielectric constant")
            self.status=3
            return True

        self.msg=""
        self.status=1
        return True


    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.areaCalculada=Area(self.kwargs["area"])
        self.deltaP=Pressure(self.kwargs["deltaP"])

        if self.kwargs["potencialCarga"]:
            self.potencialCarga=PotencialElectric(self.kwargs["potencialCarga"])
        else:
            self.potencialCarga=PotencialElectric(24000)
        if self.kwargs["potencialDescarga"]:
            self.potencialDescarga=PotencialElectric(self.kwargs["potencialDescarga"])
        else:
            self.potencialDescarga=PotencialElectric(24000)
        if self.kwargs["epsilon"]:
            self.epsilon=Dimensionless(self.kwargs["epsilon"])
        else:
            self.epsilon=Dimensionless(4.)

        if self.kwargs["metodo"]==0:
            self.rendimiento_parcial=self.calcularRendimientos_parciales(self.areaCalculada)
            self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)
        else:
            def f(area):
                rendimientos=self.calcularRendimientos_parciales(area)
                return self.calcularRendimiento(rendimientos)-self.kwargs["rendimientoAdmisible"]

            self.areaCalculada=Area(fsolve(f, 100)[0])
            self.rendimiento_parcial=self.calcularRendimientos_parciales(self.areaCalculada)
            self.rendimiento=self.calcularRendimiento(self.rendimiento_parcial)

        self.CalcularSalidas()


    def calcularRendimientos_parciales(self, A):
        """Método que calcula el rendimiento de separación, desglosado según la distribución de sólidos"""
        rendimiento_parcial=[]
        for dp in self.entrada.solido.diametros:
            if dp <= 1e-6:
                q=dp*e*1e8
            else:
                q=pi*epsilon_0*self.potencialCarga*dp**2*(1+2*(self.epsilon-1.)/(self.epsilon+2.))
            U=q*self.potencialDescarga/(3*pi*dp*self.entrada.Gas.mu)
            rendimiento_parcial.append(1-exp(-U*A/self.entrada.Q))
        return rendimiento_parcial


    def propTxt(self):
        txt=os.linesep+"#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Mode"), self.TEXT_TIPO[self.kwargs["metodo"]].split(":")[0])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure Loss"), self.deltaP.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Charging field"), self.potencialCarga.str)
        if self.kwargs["potencialCarga"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Collecting field"), self.potencialDescarga.str)
        if self.kwargs["potencialDescarga"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Dielectric constant"), self.epsilon.str)
        if self.kwargs["epsilon"]==0.0:
            txt+=" (%s)" % QApplication.translate("pychemqt", "stimated")+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Area"), self.areaCalculada.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimiento.str)+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mass Flow"), self.entrada.solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Solid Mean Diameter"), self.entrada.solido.diametro_medio.str)+os.linesep
        if len(self.entrada.solido.diametros)>=1:
            txt+="%10s %10s %10s %10s %10s" %("Dp(µm)", "Xi", "ηi", "Xis", "Xig")+os.linesep
            for i in range(len(self.rendimiento_parcial)):
                txt+="%10.1f %10.2f %10.3f %10.3f %10.3f" % (self.entrada.solido.diametros[i].config("ParticleDiameter"), self.entrada.solido.fracciones[i],  self.rendimiento_parcial[i], self.salida[1].solido.fracciones[i], self.salida[0].solido.fracciones[i])+os.linesep

        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), self.salida[0].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), self.salida[0].solido.diametro_medio.str)+os.linesep
        txt+=os.linesep+"%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mass Flow"), self.salida[1].solido.caudal.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Solid Output Mean Diameter"), self.salida[1].solido.diametro_medio.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Mode"), ("TEXT_TIPO", "metodo"), str),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Pressure Loss"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Charging field"), "potencialCarga", PotencialElectric),
                (QApplication.translate("pychemqt", "Collecting field"), "potencialDescarga", PotencialElectric),
                (QApplication.translate("pychemqt", "Dielectric constant"), "epsilon", Dimensionless),
                (QApplication.translate("pychemqt", "Area"), "areaCalculada", Area),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimiento", Dimensionless),
                (QApplication.translate("pychemqt", "Input Solid Mass Flow"), "Min", MassFlow),
                (QApplication.translate("pychemqt", "Input Solid Mean Diameter"), "Dmin", Length),
                (QApplication.translate("pychemqt", "Gas Output Solid Mass Flow"), "Mr", MassFlow),
                (QApplication.translate("pychemqt", "Gas Output Solid Mean Diameter"), "Dmr", Length),
                (QApplication.translate("pychemqt", "Solid Output Mass Flow"), "Ms", MassFlow),
                (QApplication.translate("pychemqt", "Solid Output Mean Diameter"), "Dms", Length)]
        return list


if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    solido=Solid(caudalSolido=[0.1], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    corriente=Corriente(T=300, P=101325, caudalMasico=1.,  fraccionMolar=[1.], solido=solido)

#    ciclon=Ciclon(entrada=corriente, tipo_calculo=1, rendimientoAdmisible=0.95, velocidadAdmisible=5)
#    print ciclon.msg, ciclon.status
#    print ciclon.C_instTotal, ciclon.C_adqTotal

#    camara=GravityChamber(entrada=corriente, metodo=1, modelo=1, H=1, rendimientoAdmisible=0.9)
#    camara.show()

#
    filtro=Baghouse(entrada=corriente, metodo=1, num_filtros=4, tiempo=3600, deltaP=0.1)
    print((filtro.floorArea, filtro.Vgas.ftmin))


#    precipitador=ElectricPrecipitator(entrada=corriente, metodo=1, rendimientoAdmisible=0.9)
#    print precipitador.msg, precipitador.status
#    print precipitador.areaCalculada, precipitador.rendimiento
