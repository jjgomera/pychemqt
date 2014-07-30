#!/usr/bin/python
# -*- coding: utf-8 -*-

# Modulo que define los equipos de cambio de presión de gases: compresores, turbinas, sopladores...

import os

from PyQt4.QtGui import QApplication
from scipy import log, exp
from scipy.constants import R
from scipy.optimize import fsolve

from lib.unidades import DeltaT, DeltaP, Temperature, Pressure, MassFlow, Power, Currency, Dimensionless
from lib.corriente import Corriente
from parents import equipment


class Compressor(equipment):
    """Clase que define un compresor

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la tubería
        metodo: valor que indica el típo de cálculo necesario y las magnitudes a especificar
            0   -   Presión de salida y eficiencia
            1   -   Razón de compresión y eficiencia
            2   -   Trabajo real y eficiencia
            3   -   Presión de salida y trabajo real
            4   -   Razón de compresión y trabajo real
            5   -   Calcular caudal comprimido conocidos presión a la salida, trabajo y rendimiento
            6   -   Curva de funcionamiento (No implementado)
        termodinamica: indice con el que se indica el típo de proceso termodinámico:
            0   -   Adiabático
            1   -   Politrópico
            2   -   Isotérmico
        Pout: presión a la salida
        razon: razón de compresión
        rendimiento: rendimiento del compresor
        etapas: numero de etapas del compresor
        trabajo: trabajo conocido que realiza el compresor

    Coste
        f_install: factor instalación
        base: indice base
        actual: Indice actual
        compresor: Tipo de compresor
            0   -   Centrifugal compressor
            1   -   Reciprocating compressor
            2   -    Screw compressor
        transmision: Tipo de trasmisión
            0   -   Belt drive coupling
            1   -   Chain drive coupling
            2   -   Variable speed drive coupling
        motor: Tipo de motor
            0   -   Open drip-proof
            1   -   Totally enclosed, fan-cooled
            2   -   Explosion-proof
        rpm
            0   -   3600 rpm
            1   -   1800 rpm
            2   -   1200 rpm

    >>> corriente=Corriente(T=400, P=101325, caudalMasico=0.1, fraccionMasica=[1., 0, 0, 0])
    >>> compresor=Compressor(entrada=corriente, metodo=1, termodinamica=0, razon=3, rendimiento=0.75, etapas=1)
    >>> print compresor.power.kW
    31.1007420914
    >>> compresor(compresor=2, transmision=1, motor=0, rpm=1)
    >>> print compresor.C_inst
    60464.196881
    """
    title=QApplication.translate("pychemqt", "Compressor")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "termodinamica": 0,
                    "Pout": 0.0,
                    "razon": 0.0,
                    "rendimiento": 0.0,
                    "etapas": 0,
                    "trabajo": 0.0,

                    "f_install": 1.3,
                    "Base_index": 0.0,
                    "Current_index": 0.0,
                    "compresor": 0,
                    "transmision": 0,
                    "motor": 0,
                    "rpm": 0}
    kwargsInput=("entrada", )
    kwargsValue=("Pout", "razon", "rendimiento", "etapas", "trabajo")
    kwargsList=("metodo", "termodinamica", "compresor", "transmision", "motor", "rpm")
    calculateValue=("power", "cp_cv", "razonCalculada", "rendimientoCalculado")
    calculateCostos=("C_comp", "C_motor", "C_trans", "C_adq", "C_inst")
    indiceCostos=7

    TEXT_METODO=[QApplication.translate("pychemqt", "Specify out pressure and efficiency"),
                                QApplication.translate("pychemqt", "Specify pressure ratio and efficiency"),
                                QApplication.translate("pychemqt", "Specify actual power and efficiency"),
                                QApplication.translate("pychemqt", "Specify out pressure and actual power"),
                                QApplication.translate("pychemqt", "Specify pressure ratio and actual power"),
                                QApplication.translate("pychemqt", "Calculate input flowrate")]
    TEXT_TERMODINAMICA=[QApplication.translate("pychemqt", "Adiabatic"),
                                            QApplication.translate("pychemqt", "Polytropic"),
                                            QApplication.translate("pychemqt", "Isotermic")]
    TEXT_COMPRESOR=[QApplication.translate("pychemqt", "Centrifugal compressor"),
                                    QApplication.translate("pychemqt", "Reciprocating compressor"),
                                    QApplication.translate("pychemqt", "Screw compressor")]
    TEXT_TRANSMISION=[QApplication.translate("pychemqt", "Belt drive coupling"),
                                        QApplication.translate("pychemqt", "Chain drive coupling"),
                                        QApplication.translate("pychemqt", "Variable speed drive coupling")]
    TEXT_MOTOR=[QApplication.translate("pychemqt", "Open drip-proof"),
                            QApplication.translate("pychemqt", "Totally enclosed, fan-cooled"),
                            QApplication.translate("pychemqt", "Explosion-proof")]
    TEXT_RPM=["3600 RPM", "1800 RPM", "1200 RPM"]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and self.kwargs["Current_index"]:
            self.statusCoste=True
        else:
            self.statusCoste=False

        if not self.kwargs["etapas"]:
            self.msg=QApplication.translate("pychemqt", "undefined variables")
            self.status=0
        if not self.kwargs["entrada"]:
            if self.kwargs["metodo"]==5:
                if (self.kwargs["razon"] or self.kwargs["Pout"]) and self.kwargs["trabajo"] and self.kwargs["rendimiento"]:
                    self.status=1
                    self.msg=""
                    return True
                else:
                    self.msg=QApplication.translate("pychemqt", "undefined variables")
                    self.status=0
            else:
                self.msg=QApplication.translate("pychemqt", "undefined input")
                self.status=0
        else:
            if self.kwargs["metodo"]==0:
                valores=self.kwargs["Pout"] and self.kwargs["rendimiento"]
            elif self.kwargs["metodo"]==1:
                valores=self.kwargs["razon"] and self.kwargs["rendimiento"]
            elif self.kwargs["metodo"]==2:
                valores=self.kwargs["trabajo"] and self.kwargs["rendimiento"]
            elif self.kwargs["metodo"]==3:
                valores=self.kwargs["Pout"] and self.kwargs["trabajo"]
            elif self.kwargs["metodo"]==4:
                valores=self.kwargs["razon"] and self.kwargs["trabajo"]
            elif self.kwargs["metodo"]==5:
                valores=self.kwargs["razon"] and self.kwargs["rendimiento"]

            if valores:
                self.status=1
                self.msg=""
                return True
            else:
                self.msg=QApplication.translate("pychemqt", "undefined variables")
                self.status=0

    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        metodo=self.kwargs["metodo"]
        self.Pout=Pressure(self.kwargs["Pout"])
        razon=self.kwargs["razon"]
        self.rendimientoCalculado=Dimensionless(self.kwargs["rendimiento"])
        if self.kwargs["etapas"]:
            self.etapas=self.kwargs["etapas"]
        else:
            self.etapas=1.
        self.power=Power(self.kwargs["trabajo"])

        def f(Pout, rendimiento):
            W_ideal=self.__Wideal(Pout)
            power=W_ideal*self.entrada.caudalmasico.gs/rendimiento
            return power

        if metodo in [0, 3] or (metodo==5 and self.Pout):
            if self.etapas==1:
                razon=self.Pout.atm/self.entrada.P.atm
            else:
                razon=(self.Pout.atm/self.entrada.P.atm)**(1./self.etapas)
        elif metodo in [1, 4]  or (metodo==5 and razon):
            if self.etapas==1:
                self.Pout=Pressure(self.entrada.P*razon)
            else:
                self.Pout=Pressure(razon**self.etapas*self.entrada.P)

        if metodo in [0, 1]:
            Wideal=self.__Wideal(self.Pout.atm)
            power=Wideal*self.entrada.caudalmasico.gs/self.rendimientoCalculado
            self.power=Power(power*self.etapas)
        elif metodo == 2:
            funcion=lambda P: f(P, self.rendimientoCalculado)-self.power
            self.Pout=Pressure(fsolve(funcion, self.entrada.P+1))
            if self.etapas==1:
                razon=self.Pout/self.entrada.P
            else:
                razon=(self.Pout/self.entrada.P)**(1./self.etapas)
        elif metodo in [3, 4]:
            funcion=lambda rendimiento: f(self.Pout.atm, rendimiento)-self.power
            self.rendimientoCalculado=Dimensionless(fsolve(funcion, 0.5))
        elif metodo ==5:
            Wideal=self.__Wideal(self.Pout.atm)
            caudalmasico=MassFlow(self.power*self.rendimientoCalculado/Wideal, "gs")
            self.Tout=self.__Tout(Wideal)
            self.entrada=self.entrada.clone(caudalMasico=caudalmasico)

        Wideal=self.__Wideal(self.Pout.atm)
        self.Tout=Temperature(self.__Tout(Wideal))
        self.salida=[self.entrada.clone(T=self.Tout, P=self.Pout)]
        self.razonCalculada=Dimensionless(self.Pout/self.entrada.P)
        self.deltaT=DeltaT(self.salida[0].T-self.entrada.T)
        self.deltaP=DeltaP(self.salida[0].P-self.entrada.P)
        self.cp_cv=self.entrada.Gas.cp_cv
        self.Pin=self.entrada.P
        self.Tin=self.entrada.T


    def __Wideal(self, Pout):
        """Calculo del trabajo ideal"""
        if self.kwargs["termodinamica"]==0:
            cpv=self.entrada.Gas.cp_cv
        else:
            cpv=1.40388
        if self.kwargs["termodinamica"]==2:
            W_ideal=R*self.entrada.T/self.entrada.M*log(Pout/self.entrada.P)
        else:
            W_ideal=R*self.entrada.T*cpv/self.entrada.M/(cpv-1)*((Pout/self.entrada.P.atm)**((cpv-1)/cpv)-1)
        return W_ideal

    def __Tout(self, W_ideal):
        """Cálculo de la temperatura de salida"""
        if self.kwargs["termodinamica"]==0:
            cpv=self.entrada.Gas.cp_cv
        else:
            cpv=1.40388
        if self.kwargs["termodinamica"]==2:
            Tout=self.entrada.T
        else:
            Tout=(self.Pout.atm/self.entrada.P.atm)**((cpv-1)/cpv)*self.entrada.T+(1-self.rendimientoCalculado)/self.rendimientoCalculado*(W_ideal/self.entrada.Gas.cv)
        return Tout


    def coste(self):
        HP=self.power.hp/self.etapas

        if self.kwargs["compresor"]==0:
            #Centrifugal compressor
            C=6.49*HP**0.62*1000
        elif self.kwargs["compresor"]==1:
            #Reciprocating compressor
            C=5.96*HP**0.61*1000
        elif self.kwargs["compresor"]==2:
            #Screw compressor
            C=1.49*HP**0.71*1000

        C_comp=self.etapas*C*self.kwargs["Current_index"] / self.kwargs["Base_index"]

        if self.kwargs["motor"] == 0:                                     #Open, drip-proof
            if self.kwargs["rpm"] == 0 and HP <= 7.5:
                a1, a2, a3 = 4.8314, 0.0966, 0.10960
            elif self.kwargs["rpm"] == 0 and HP <= 250. and HP > 7.5:
                a1, a2, a3 = 4.1514, 0.5347, 0.05252
            elif self.kwargs["rpm"] == 0 and HP > 250.:
                a1, a2, a3 = 4.2432, 1.03251, -0.03595
            elif self.kwargs["rpm"] == 1 and HP <= 7.5:
                a1, a2, a3 = 4.7075, -0.01511, 0.22888
            elif self.kwargs["rpm"] == 1 and HP <= 250. and HP > 7.5:
                a1, a2, a3 = 4.5212, 0.47242, 0.04820
            elif self.kwargs["rpm"] == 1 and HP > 250.:
                a1, a2, a3 = 7.4044, -0.06464, 0.05448
            elif self.kwargs["rpm"] == 2 and HP <= 7.5:
                a1, a2, a3 = 4.9298, 0.30118, 0.12630
            elif self.kwargs["rpm"] == 2 and HP <= 250. and HP > 7.5:
                a1, a2, a3 = 5.0999, 0.35861, 0.06052
            elif self.kwargs["rpm"] == 2 and HP > 250.:
                a1, a2, a3 = 4.6163, 0.88531, -0.02188
        elif self.kwargs["motor"] == 1:                                   #Totally enclosed, fan-cooled
            if self.kwargs["rpm"] == 0 and HP <= 7.5:
                a1, a2, a3 = 5.1058, 0.03316, 0.15374
            elif self.kwargs["rpm"] == 0 and HP <= 250. and HP > 7.5:
                a1, a2, a3 = 3.8544, 0.83311, 0.02399
            elif self.kwargs["rpm"] == 0 and HP > 250.:
                a1, a2, a3 = 5.3182, 1.08470, -0.05695
            elif self.kwargs["rpm"] == 1 and HP <= 7.5:
                a1, a2, a3 = 4.9687, -0.00930, 0.22616
            elif self.kwargs["rpm"] == 1 and HP > 7.5:
                a1, a2, a3 = 4.5347, 0.57065, 0.04609
            elif self.kwargs["rpm"] == 2 and HP <= 7.5:
                a1, a2, a3 = 5.1532, 0.28931, 0.14357
            elif self.kwargs["rpm"] == 2 and HP > 7.5:
                a1, a2, a3 = 5.3858, 0.31004, 0.07406
        elif self.kwargs["motor"] == 2:                                    #Explosion-proof
            if self.kwargs["rpm"] == 0 and HP <= 7.5:
                a1, a2, a3 = 5.3934, -0.00333, 0.15475
            elif self.kwargs["rpm"] == 0 and HP > 7.5:
                a1, a2, a3 = 4.4442, 0.60820, 0.05202
            elif self.kwargs["rpm"] == 1 and HP <= 7.5:
                a1, a2, a3 = 5.2851, 0.00048, 0.19949
            elif self.kwargs["rpm"] == 1 and HP > 7.5:
                a1, a2, a3 = 4.8178, 0.51086, 0.05293
            elif self.kwargs["rpm"] == 2 and HP <= 7.5:
                a1, a2, a3 = 5.4166, 0.31216, 0.10573
            elif self.kwargs["rpm"] == 2 and HP > 7.5:
                a1, a2, a3 = 5.5655, 0.31284, 0.07212

        LnHP = log(HP)
        C_motor = self.etapas*1.2 * exp(a1 + a2 * LnHP + a3 * LnHP**2)

        if self.kwargs["compresor"]==2:
            C_trans=0
        elif self.kwargs["transmision"]==0:
            C_trans=1.2*exp(3.689+0.8917*log(HP))
        elif self.kwargs["transmision"]==1:
            C_trans=1.2*exp(5.329+0.5048*log(HP))
        elif self.kwargs["transmision"]==2:
            C_trans=12000/(1.562+7.877/HP)
        C_trans=self.etapas*C_trans * self.kwargs["Current_index"] / self.kwargs["Base_index"]

        C_adq=C_comp+C_motor+C_trans
        C_inst=C_adq * self.kwargs["f_install"]

        self.C_comp=Currency(C_comp)
        self.C_motor=Currency(C_motor)
        self.C_trans=Currency(C_trans)
        self.C_adq=Currency(C_adq)
        self.C_inst=Currency(C_inst)


    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Temperature"), self.entrada.T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Temperature"), self.salida[0].T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Temperature increase"), self.deltaT.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure increase"), self.deltaP.str)+os.linesep
        txt+="%-25s\t %0.4f" %(QApplication.translate("pychemqt", "Pressure ratio"), self.razonCalculada)+os.linesep
        txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Thermo Condition"), self.TEXT_TERMODINAMICA[self.kwargs["termodinamica"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Power"), self.power.str)+os.linesep
        txt+="%-25s\t %0.4f" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimientoCalculado)+os.linesep
        txt+="%-25s\t %0.4f" %("Cp/Cv", self.entrada.Gas.cp_cv)+os.linesep

        if self.statusCoste:
            txt+=os.linesep
            txt+="#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep

            txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Compressor type"), self.TEXT_COMPRESOR[self.kwargs["compresor"]])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Compressor Cost"), self.C_comp.str)+os.linesep
            if self.kwargs["compresor"]!=2:
                txt+="%-25s\t %s" %(QApplication.translate("pychemqt", "Transmission type"), self.TEXT_TRANSMISION[self.kwargs["transmision"]])+os.linesep
                txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Transmission Cost"), self.C_trans.str)+os.linesep
            txt+="%-25s\t %s, %s" %(QApplication.translate("pychemqt", "Motor type"), self.TEXT_MOTOR[self.kwargs["motor"]], self.TEXT_RPM[self.kwargs["rpm"]])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Motor Cost"), self.C_motor.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Input Temperature"), "Tin", Temperature),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Temperature"), "Tout", Temperature),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Temperature increase"), "deltaT", DeltaT),
                (QApplication.translate("pychemqt", "Pressure increase"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Pressure ratio"), "razonCalculada", Dimensionless),
                (QApplication.translate("pychemqt", "Thermodinamic mode"), ("TEXT_TERMODINAMICA", "termodinamica"),  str),
                (QApplication.translate("pychemqt", "Power"), "power", Power),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimientoCalculado", Dimensionless),
                (QApplication.translate("pychemqt", "Especific capacities ratio"), "cp_cv", Dimensionless),
                (QApplication.translate("pychemqt", "Compressor Type"), ("TEXT_COMPRESOR", "compresor"),  str),
                (QApplication.translate("pychemqt", "Transmission Type"), ("TEXT_TRANSMISION", "transmision"),  str),
                (QApplication.translate("pychemqt", "Motor Type"), ("TEXT_MOTOR", "motor"),  str),
                (QApplication.translate("pychemqt", "Motor RPM"), ("TEXT_RPM", "rpm"),  str),
                (QApplication.translate("pychemqt", "Cost compressor"), "C_comp", Currency),
                (QApplication.translate("pychemqt", "Cost Transmission"), "C_trans", Currency),
                (QApplication.translate("pychemqt", "Cost motor"), "C_motor", Currency),
                (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", Currency),
                (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", Currency)]
        return list


class Turbine(equipment):
    """Clase que define una turbina

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la tubería
        metodo: valor que indica el típo de cálculo necesario y las magnitudes a especificar
            0   -   Presión de salida y eficiencia
            1   -   Razón de presiones y eficiencia
            2   -   Trabajo real y eficiencia
            3   -   Presión de salida y trabajo real
            4   -   Razón de compresión y trabajo real
            5   -   Calcular caudal comprimido conocidos presión a la salida, trabajo y rendimiento
        termodinamica: indice con el que se indica el típo de proceso termodinámico:
            0   -   Adiabático
            1   -   Politrópico
            2   -   Isotérmico
        Pout: presión a la salida
        razon: razón de presiones
        rendimiento: rendimiento del compresor
        trabajo: trabajo conocido que realiza el compresor

    Coste
        f_install: factor instalación
        base: indice base
        actual: Indice actual

    >>> corriente=Corriente(T=400, P=101325, caudalMasico=0.1, fraccionMasica=[1., 0, 0, 0])
    >>> turbina=Turbine(entrada=corriente, metodo=1, razon=0.3, rendimiento=1.)
    >>> print turbina.power.MJh
    -69.138664576
    >>> print turbina.C_inst
    30713.7301133
    """

    title=QApplication.translate("pychemqt", "Turbine")
    help=""
    kwargs={"entrada": None,
                    "metodo": 0,
                    "termodinamica": 0,
                    "Pout": 0.0,
                    "razon": 0.0,
                    "rendimiento": 0.0,
                    "trabajo": 0.0,

                    "f_install": 1.5,
                    "Base_index": None,
                    "Current_index": None}
    kwargsInput=("entrada", )
    kwargsValue=("Pout", "razon", "rendimiento", "trabajo")
    kwargsList=("metodo", "termodinamica")
    calculateValue=("power", "cp_cv", "razonCalculada", "rendimientoCalculado")
    calculateCostos=("C_adq", "C_inst")
    indiceCostos=2

    TEXT_METODO=[QApplication.translate("pychemqt", "Specify out pressure and efficiency"),
                                QApplication.translate("pychemqt", "Specify pressure ratio and efficiency"),
                                QApplication.translate("pychemqt", "Specify actual power and efficiency"),
                                QApplication.translate("pychemqt", "Specify out pressure and actual power"),
                                QApplication.translate("pychemqt", "Specify pressure ratio and actual power"),
                                QApplication.translate("pychemqt", "Calculate input flowrate")]
    TEXT_TERMODINAMICA=[QApplication.translate("pychemqt", "Adiabatic"),
                                            QApplication.translate("pychemqt", "Polytropic"),
                                            QApplication.translate("pychemqt", "Isotermic")]

    @property
    def isCalculable(self):
        if self.kwargs["f_install"] and self.kwargs["Base_index"] and self.kwargs["Current_index"]:
            self.statusCoste=True
        else:
            self.statusCoste=False

        if not self.kwargs["entrada"]:
            if self.kwargs["metodo"]==5:
                if (self.kwargs["razon"] or self.kwargs["Pout"]) and self.kwargs["trabajo"] and self.kwargs["rendimiento"]:
                    self.status=1
                    self.msg=""
                    return True
                else:
                    self.msg=QApplication.translate("pychemqt", "undefined variables")
                    self.status=0
            else:
                self.msg=QApplication.translate("pychemqt", "undefined input")
                self.status=0
        else:
            if self.kwargs["metodo"]==0:
                valores=self.kwargs["Pout"] and self.kwargs["rendimiento"]
            elif self.kwargs["metodo"]==1:
                valores=self.kwargs["razon"] and self.kwargs["rendimiento"]
            elif self.kwargs["metodo"]==2:
                valores=self.kwargs["trabajo"] and self.kwargs["rendimiento"]
            elif self.kwargs["metodo"]==3:
                valores=self.kwargs["Pout"] and self.kwargs["trabajo"]
            elif self.kwargs["metodo"]==4:
                valores=self.kwargs["razon"] and self.kwargs["trabajo"]
            elif self.kwargs["metodo"]==5:
                valores=self.kwargs["razon"] and self.kwargs["rendimiento"]

            if valores:
                self.status=1
                self.msg=""
                return True
            else:
                self.msg=QApplication.translate("pychemqt", "undefined variables")
                self.status=0

    def calculo(self):
        self.entrada=self.kwargs["entrada"]
        self.Pout=Pressure(self.kwargs["Pout"])
        self.razon=Dimensionless(self.kwargs["razon"])
        self.rendimientoCalculado=Dimensionless(self.kwargs["rendimiento"])
        self.power=Power(-abs(self.kwargs["trabajo"]))

        def f(Pout, rendimiento):
            W_ideal=self.__Wideal(Pout)
            power=W_ideal*self.entrada.caudalmasico.gs*rendimiento
            return power

        self.cp_cv=self.entrada.Gas.cp_cv

        if self.kwargs["metodo"] in [0, 3] or (self.kwargs["metodo"]==5 and self.Pout):
            self.razon=Dimensionless(self.Pout/self.entrada.P)
        elif self.kwargs["metodo"] in [1, 4]  or (self.kwargs["metodo"]==5 and self.razon):
                self.Pout=Pressure(self.entrada.P*self.razon)

        if self.kwargs["metodo"] in [0, 1]:
            Wideal=self.__Wideal(self.Pout)
            self.power=Power(Wideal*self.entrada.caudalmasico.gs*self.rendimientoCalculado)
        elif self.kwargs["metodo"] == 2:
            funcion=lambda P: f(P, self.rendimientoCalculado)-self.power
            self.Pout=Pressure(fsolve(funcion, self.entrada.P.atm+1), "atm")
            self.razon=Dimensionless(self.Pout/self.entrada.P)
        elif self.kwargs["metodo"] in [3, 4]:
            funcion=lambda rendimiento: f(self.Pout.atm, rendimiento)-self.power
            self.rendimientoCalculado=fsolve(funcion, 0.5)
        elif self.kwargs["metodo"] ==5:
            Wideal=self.__Wideal(self.Pout)
            caudalmasico=MassFlow(self.power/self.rendimientoCalculado/Wideal, "gs")
            self.Tout=self.__Tout(Wideal)
            self.entrada=self.entrada.clone(caudalMasico=caudalmasico)

        Wideal=self.__Wideal(self.Pout)
        self.Tout=Temperature(self.__Tout(Wideal))
        self.razonCalculada=Dimensionless(self.Pout/self.entrada.P)
        self.salida=[self.entrada.clone(T=self.Tout, P=self.Pout)]
        self.deltaT=DeltaT(self.salida[0].T-self.entrada.T)
        self.deltaP=DeltaP(self.salida[0].P-self.entrada.P)
        self.Pin=self.entrada.P
        self.Tin=self.entrada.T

    def __Wideal(self, Pout):
        """Calculo del trabajo ideal"""
        if self.kwargs["termodinamica"]==0:
            cpv=self.cp_cv
        else:
            cpv=1.40388
        if self.kwargs["termodinamica"]==2:
            W_ideal=R*self.entrada.T/self.entrada.M*log(Pout/self.entrada.P)
        else:
            W_ideal=R*self.entrada.T*cpv/self.entrada.M/(cpv-1)*((Pout/self.entrada.P)**((cpv-1)/cpv)-1)
        return W_ideal

    def __Tout(self, W_ideal):
        """Cálculo de la temperatura de salida"""
        if self.kwargs["termodinamica"]==0:
            cpv=self.cp_cv
        else:
            cpv=1.40388
        if self.kwargs["termodinamica"]==2:
            Tout=self.entrada.T
        else:
            Tout=(self.Pout/self.entrada.P)**((cpv-1)/cpv)*self.entrada.T+(1-self.rendimientoCalculado)/self.rendimientoCalculado*(W_ideal/self.entrada.Gas.cv)
        return Tout

    def coste(self):
        HP=abs(self.power.hp)

        if self.salida[0].P.psi>=14.696:
            C=0.31*HP**0.81*1000
        else:
            C=0.69*HP**0.81*1000

        self.C_adq=Currency(C * self.kwargs["Current_index"] / self.kwargs["Base_index"])
        self.C_inst=Currency(self.C_adq * self.kwargs["f_install"])


    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Calculate properties")+"-----------------#"+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Temperature"), self.entrada.T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Temperature"), self.salida[0].T.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Temperature increase"), self.deltaT.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Input Pressure"), self.entrada.P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Output Pressure"), self.salida[0].P.str)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Pressure increase"), self.deltaP.str)+os.linesep
        txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Pressure ratio"), self.razon)+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Thermo Condition"), self.TEXT_TERMODINAMICA[self.kwargs["termodinamica"]])+os.linesep
        txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Power"), self.power.str)+os.linesep
        txt+="%-25s\t %0.4f" %(QApplication.translate("pychemqt", "Efficiency"), self.rendimientoCalculado)+os.linesep
        txt+="%-25s\t %0.4f" %("Cp/Cv", self.cp_cv)+os.linesep

        if self.statusCoste:
            txt+=os.linesep
            txt+="#---------------"+QApplication.translate("pychemqt", "Preliminary Cost Estimation")+"-----------------#"+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Base index"), self.kwargs["Base_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Current index"), self.kwargs["Current_index"])+os.linesep
            txt+="%-25s\t %0.2f" %(QApplication.translate("pychemqt", "Install factor"), self.kwargs["f_install"])+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Purchase Cost"), self.C_adq.str)+os.linesep
            txt+="%-25s\t%s" %(QApplication.translate("pychemqt", "Installed Cost"), self.C_inst.str)+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Input Temperature"), "Tin", Temperature),
                (QApplication.translate("pychemqt", "Input Pressure"), "Pin", Pressure),
                (QApplication.translate("pychemqt", "Output Temperature"), "Tout", Temperature),
                (QApplication.translate("pychemqt", "Output Pressure"), "Pout", Pressure),
                (QApplication.translate("pychemqt", "Temperature increase"), "deltaT", DeltaT),
                (QApplication.translate("pychemqt", "Pressure increase"), "deltaP", DeltaP),
                (QApplication.translate("pychemqt", "Pressure ratio"), "razonCalculada", Dimensionless),
                (QApplication.translate("pychemqt", "Thermodinamic mode"), ("TEXT_TERMODINAMICA", "termodinamica"),  str),
                (QApplication.translate("pychemqt", "Power"), "power", Power),
                (QApplication.translate("pychemqt", "Efficiency"), "rendimientoCalculado", Dimensionless),
                (QApplication.translate("pychemqt", "Especific capacities ratio"), "cp_cv", Dimensionless),
                (QApplication.translate("pychemqt", "Purchase Cost"), "C_adq", Currency),
                (QApplication.translate("pychemqt", "Installed Cost"), "C_inst", Currency)]
        return list


if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

#    corriente=Corriente(T=400, P=101325, caudalMasico=0.1, fraccionMasica=[1., 0, 0, 0])
#    compresor=Compressor(entrada=corriente, metodo=1, termodinamica=0, razon=3, rendimiento=0.75, etapas=1)
#    print compresor.power.kW
#    compresor(compresor=2, transmision=1, motor=0, rpm=1)
#    print compresor.C_inst

    corriente=Corriente(T=873.15, P=1013250, caudalMasico=0.1, fraccionMasica=[1.])
    turbina=Turbine(entrada=corriente, metodo=1, razon=0.3, rendimiento=0.5)
    print turbina.power.MJh
    print turbina.C_inst
    print turbina.salida[0].T

#    agua=Corriente(T=500, P=101325., caudalMasico=1, fraccionMasica=[1.])
#    print agua.P*2
#    compresor=Compressor(entrada=agua, Pout=5*101325, rendimiento=0.75)
