#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################################
# librería de definición de equipos de separación gas-sólido-liquido:
#     -Secadores de sólidos
#     -Lavadores de gases
#######################################################################

from PyQt4.QtGui import QApplication
from scipy import pi

from lib import unidades
from lib.corriente import Corriente
from lib.psycrometry import Punto_Psicrometrico
from parents import equipment


class Dryer(equipment):
    """Clase que define un equipo de secado de sólidos

    Parámetros:
        entradaSolido: Instancia de clase Corriente que define la entrada de sólidos
        entradaAire: Instancia de clase Psychrometry que define le aire entrante
        entrada: array con ambas entradas al equipo
        tipoCalculo: integer que indica el tipo de cálculo
            0   -   Cálculo, se calculan la humedad del aire a la salida
            1   -   Diseño, se imponen las condiciones de salida del gas y se calcula el caudal necesario
        HR: Humedad relativa del aire a la salida, por defecto 100%
        TemperaturaSolid: Temperatura del sólido a la salida, por defecto se considerará la temperatura de salida del vapor
        HumedadResidual: Humedad residual del sólido, en kg/kg por defecto se considerará 0
        Heat: Calor intercambiado por el equipo, por defecto se considerará funcionamiento adiabático, heat=0
        deltaP: Perdida de presión del equipo

    """
    title=QApplication.translate("pychemqt", "Solid dryer")
    help=""
    kwargs={"entradaSolido": None,
                    "entradaAire": None,
                    "entrada": [],
                    "tipoCalculo": 0,
                    "HR": 0.0,
                    "TemperaturaSolid": 0.0,
                    "HumedadResidual": 0.0,
                    "Heat": 0.0,
                    "deltaP": 0.0}

    @property
    def isCalculable(self):
        if not self.kwargs["entradaSolido"]:
            self.msg=QApplication.translate("pychemqt", "undefined solid stream input")
            self.status=0
        elif not self.kwargs["entradaAire"]:
            self.msg=QApplication.translate("pychemqt", "undefined air stream input")
            self.status=0
        elif not self.kwargs["HR"]:
            self.msg=QApplication.translate("pychemqt", "using default air output relative humid 100%")
            self.status=3
            return True
        else:
            self.msg=""
            self.status=1
            return True

    def cleanOldValues(self, **kwargs):
        """Actualización de los kwargs con los nuevos introducidos si es necesario para cada equipo"""
        if kwargs.has_key("entrada"):
            kwargs["entradaSolido"]=kwargs["entrada"][0]
            kwargs["entradaSolido"]=kwargs["entrada"][1]
            del kwargs["entrada"]
        self.kwargs.update(kwargs)


    def calculo(self):
        #TODO: De momento, no se implementan las cuestiones de cinetica de intercambio de calor y materia que definirían las dimensiones necesarias del equipo
        HR=self.kwargs.get("HR", 100)
        self.Heat=unidades.Power(self.kwargs["Heat"])
        self.deltaP=unidades.Pressure(self.kwargs["deltaP"])
        self.entradaAire=self.kwargs["entradaAire"]

        Pout=min(self.kwargs["entradaSolido"].P.atm, self.kwargs["entradaAire"].P.atm)-self.deltaP.atm

        aguaSolidoSalida=self.kwargs["HumedadResidual"]*self.kwargs["entradaSolido"].solido.caudal.kgh
        aguaSolidoEntrada=self.kwargs["entradaSolido"].caudalmasico.kgh
        if self.kwargs["tipoCalculo"]==0:
            Caudal_aguaenAireSalida=aguaSolidoEntrada-aguaSolidoSalida+self.entradaAire.caudal.kgh*self.entradaAire.Xw
            Caudal_airesalida=self.entradaAire.caudal.kgh*self.entradaAire.Xa
            if self.entradaAire.Hs>Caudal_aguaenAireSalida/Caudal_airesalida:
                H=Caudal_aguaenAireSalida/Caudal_airesalida
            else:
                H=self.entradaAire.Hs
                aguaSolidoSalida+=Caudal_aguaenAireSalida/Caudal_airesalida-self.entradaAire.Hs
            self.SalidaAire=Punto_Psicrometrico(caudal=Caudal_aguaenAireSalida+Caudal_airesalida, tdb=self.entradaAire.Tdb, H=H)
            self.SalidaSolido=self.kwargs["entradaSolido"].clone(T=self.SalidaAire.Tdb, P=Pout, split=aguaSolidoSalida/aguaSolidoEntrada)
        else:
            pass


#        if self.HumedadResidual==0:
#            self.SalidaSolido=Corriente(self.entradaAire.Tdb, self.entradaAire.P.atm-self.deltaP.atm, 0, self.entradaAire.mezcla, self.entradaAire.solido)
#        else:
#            self.SalidaSolido=Corriente(self.entradaAire.Tdb, self.entradaAire.P.atm-self.deltaP.atm, 0, self.entradaAire.mezcla, self.entradaAire.solido)
#        self.SalidaAire=Corriente(self.entradaAire.T, self.entradaAire.P.atm-self.deltaP.atm, self.entradaAire.caudalmasico.kgh, self.entradaAire.mezcla)


class Scrubber(equipment):
    """ Clase que define el scrubber

    Parámetros:
        entradaGas: Instancia de clase corriente que define el gas que fluye por el equipo
        entradaLiquido: Instancia de clase corriente que define el líquido lavador
        tipo_calculo:
            0   -   Evaluación del funcionamiento de un venturi de dimensiones conocidas
            1   -   Diseño de un venturi a partir del rendimiento o perdida de carga admisibles
        diametro: diametro del venturi
        rendimiento: Rendimiento admisible
        k: contante empirica del venturi, valores entre 500 y 100
        modelo_rendimiento: Modelo de cálculo del venturi:
            0   -   Johnstone (1954)
            1   -   Calvert (1972)
        modelo_DeltaP:
            0   -   Young (1977)
    """
    title = QApplication.translate("pychemqt", "Scrubber")
    help = ""
    kwargs = {"entradaGas": None,
              "entradaLiquido": None,
              "tipo_calculo": 0,
              "diametro": 0.0,
              "rendimiento": 0.0,
              "k": 0.0,
              "modelo_rendimiento": 0,
              "modelo_DeltaP": 0, 
              "Lt": 0.0}
    kwargsInput = ("entradaGas", "entradaLiquido")
    kwargsValue = ("diametro", "rendimiento", "k", "Lt")
    kwargsList = ("tipo_calculo", "modelo_rendimiento", "modelo_DeltaP")
    calculateValue = ("deltaP", "rendimientoCalc")

    TEXT_TIPO = [QApplication.translate("pychemqt", "Rating: Calculate efficiency"),
                 QApplication.translate("pychemqt", "Design: Calculate diameter")]
    TEXT_MODEL = ["Johnstone (1954)", "Calvert (1972)"]
    TEXT_MODEL_DELTAP = ["Young (1977)", 
                                            "Hesketh (1974)",
                                            "Gleason (1971)",
                                            "Volgin (1968)"]

    @property
    def isCalculable(self):
        self.status = 1
        self.msg = ""
        
        self.statusDeltaP=1
        if self.kwargs["modelo_DeltaP"] == 3 and not self.kwargs["Lt"]:
            self.statusDeltaP=0
            


        if not self.kwargs["entradaGas"]:
            self.msg = QApplication.translate("pychemqt", "undefined gas stream")
            self.status = 0
            return
        if not self.kwargs["entradaLiquido"]:
            self.msg = QApplication.translate("pychemqt", "undefined liquid stream")
            self.status = 0
            return

        if self.kwargs["tipo_calculo"] == 0 and not self.kwargs["diametro"]:
            self.msg = QApplication.translate("pychemqt", "undefined diameter")
            self.status = 0
            return
        elif self.kwargs["tipo_calculo"] == 1 and not self.kwargs["rendimiento"]:
            self.msg = QApplication.translate("pychemqt", "undefined efficiency")
            self.status = 0
            return

        if not self.kwargs["k"]:
            self.msg = QApplication.translate("pychemqt", "undefined venturi constant")
            self.status = 3

        return True

    def calculo(self):
        self.entradaGas=self.kwargs["entradaGas"]
        self.entradaLiquido=self.kwargs["entradaLiquido"]
        self.Dt=unidades.Length(self.kwargs["diametro"])
        self.k=self.kwargs.get("k", 1000)
        
        self.At=unidades.Area(pi/4*self.Dt**2)
        self.Vg=unidades.Speed(self.entradaGas.Q/self.At)
        self.R=self.entradaLiquido.Q/self.entradaGas.Q
        
        self.rendimientoCalc=0
        
        if self.statusDeltaP:
            self.deltaP=self._deltaP()
        else:
            self.deltaP=unidades.DeltaP(0)
        
        self.salida=[self.entradaGas, self.entradaLiquido]
#        rendimiento_fraccional=[]
#        rendimiento_global=0
#        for i in self.entradaGas.distribucion_solidos:
#            l=1e-7
#            d0=0.585*(self.entradaLiquido.tension_superficial/self.entradaLiquido.RhoL/V**2)**0.5+53.*(self.entradaLiquido.viscosidad**2/self.entradaLiquido.tension_superficial/liquido.RhoL)**0.225*R**1.5
#            A=1.257+0.4*exp(-1.1*i[0]/2/l)
#            C=1+2*A*l/i[0]
#            fi=C*self.entradaGas.densidad_solido*V*i[0]**2/(9*d0*self.entradaGas.viscosidad)
#            rendimiento=1-exp(-k*R*fi**0.5)
#            rendimiento_fraccional.append(rendimiento)
#            print i[1]
#            rendimiento_global+=rendimiento*i[1]
#
#        self.rendimiento=rendimiento_global
#
#        DeltaP=Pressure(1.002*V**2*R, "kPa")
#
#        self.DeltaP=DeltaP
#
#        print V, DeltaP
        """        1. Calvert: 'dp = 5.4e-04 * (v^2) * rho_gas * (L/G)'
        2a. Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * (0.56 + 0.125*L/G + 0.0023*(L/G)^2) / 507'
        2b. Simplified Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * ((L/G)^0.78) / 1270'
                """
    def _deltaP(self):
        if self.kwargs["modelo_DeltaP"] == 0: #Young (1977)
            deltaP=0
        elif self.kwargs["modelo_DeltaP"] == 1: #Hesketh (1974)
            deltaP=1.36e-4*self.Vg.cms**2*self.entradaGas.Gas.rho.gcc*self.At.cm2**0.133*(0.56+935*self.R+1.29e-2*self.R**2)
        elif self.kwargs["modelo_DeltaP"] == 2: #Gleason (1971)
            deltaP=2.08e-5*self.Vg.cms**2*(0.264*self.entradaLiquido.Liquido.Q.ccs+73.8)
        elif self.kwargs["modelo_DeltaP"] == 3: #Volgin (1968)
            deltaP=3.32e-6*self.Vg.cms**2*self.R*0.26*l**1.43

#        elif self.kwargs["modelo_rendimiento"] == 0: #Matrozov (1953)
#            deltaP=dPd+1.38e-3*self.Vg.cms**1.08*self.R**0.63
#        elif self.kwargs["modelo_rendimiento"] == 1: #Yoshida (1960)
#            pass
#        elif self.kwargs["modelo_rendimiento"] == 3: #Tohata (1964)
#            pass
#        elif self.kwargs["modelo_rendimiento"] == 4: #Geiseke (1968)
#            pass
#        elif self.kwargs["modelo_rendimiento"] == 6: #Calvert (1968)
#            deltaP=1.03e-3*self.Vg.cms**2*self.R
#        elif self.kwargs["modelo_rendimiento"] == 8: #Boll (1973)
#            pass
#        elif self.kwargs["modelo_rendimiento"] == 9: #Behie & Beeckman (1973)
#            pass
        
        return unidades.DeltaP(deltaP, "cmH2O")




if __name__ == '__main__':
#    import doctest
#    doctest.testmod()

    from lib.corriente import Corriente, Solid
    diametros=[17.5e-6, 22.4e-6, 26.2e-6, 31.8e-6, 37e-6, 42.4e-6, 48e-6, 54e-6, 60e-6, 69e-6, 81.3e-6, 96.5e-6, 109e-6, 127e-6]
    fracciones=[0.02, 0.03, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.03, 0.02]
    solido=Solid(caudalSolido=[0.1], distribucion_diametro=diametros, distribucion_fraccion=fracciones)
    aire=Corriente(T=300, P=101325, caudalMasico=1., ids=[475], fraccionMolar=[1.], solido=solido)
    agua=Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.])
    secador=Scrubber(entradaGas=aire, entradaLiquido=agua, diametro=0.05, modelo_rendimiento=10)
    print secador.status, secador.msg
    print secador.DeltaP.bar
