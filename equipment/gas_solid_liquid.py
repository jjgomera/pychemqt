#!/usr/bin/python
# -*- coding: utf-8 -*-

#######################################################################
# librería de definición de equipos de separación gas-sólido-liquido:
#     -Secadores de sólidos
#     -Lavadores de gases
#######################################################################

from PyQt4.QtGui import QApplication

from lib.unidades import Pressure, Power
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
        self.Heat=Power(self.kwargs["Heat"])
        self.deltaP=Pressure(self.kwargs["deltaP"])
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
            0   -   Evaluación del funcionamiento de un venturi de dimensiones conocidas (diametro, forma, número de ciclones)
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
              "modelo_DeltaP": 0}
    kwargsInput = ("entradaGas", "entradaLiquido")
    kwargsValue = ("diametro", "rendimiento", "k")
    kwargsList = ("tipo_calculo", "modelo_rendimiento", "modelo_DeltaP")
    #    calculateValue = ("deltaP", "V", "rendimientoCalc", "NCalc", "Dcc", "Hc", "Bc", "Jc", "Lc", "Zc", "De", "Sc")

    TEXT_TIPO = [QApplication.translate("pychemqt", "Rating: Calculate efficiency"),
                 QApplication.translate("pychemqt", "Design: Calculate diameter")]
    TEXT_MODEL = ["Johnstone (1954)", "Calvert (1972)"]
    TEXT_MODEL_DELTAP = ["Young (1977)"]

    @property
    def isCalculable(self):
        pass

    def calculo(self, **kwargs):
        self.entradaGas=kwargs.get("entrada", None)
        self.entradaLiquido=kwargs.get("entradaLiquido", 0)
        self.D=unidades.Length(kwargs.get("diametro", 0))
        self.k=kwargs.get("k", 707)

        R=self.entradaLiquido.caudal_volumetrico/self.entradaGas.caudal_volumetrico
        V=entradaGas.caudal_volumetrico/(pi/4*self.D**2)

        rendimiento_fraccional=[]
        rendimiento_global=0
        for i in self.entradaGas.distribucion_solidos:
            l=1e-7
            d0=0.585*(self.entradaLiquido.tension_superficial/self.entradaLiquido.RhoL/V**2)**0.5+53.*(self.entradaLiquido.viscosidad**2/self.entradaLiquido.tension_superficial/liquido.RhoL)**0.225*R**1.5
            A=1.257+0.4*exp(-1.1*i[0]/2/l)
            C=1+2*A*l/i[0]
            fi=C*self.entradaGas.densidad_solido*V*i[0]**2/(9*d0*self.entradaGas.viscosidad)
            rendimiento=1-exp(-k*R*fi**0.5)
            rendimiento_fraccional.append(rendimiento)
            print i[1]
            rendimiento_global+=rendimiento*i[1]

        self.rendimiento=rendimiento_global

        DeltaP=Pressure(1.002*V**2*R, "kPa")

        self.DeltaP=DeltaP

        print V, DeltaP
        """        1. Calvert: 'dp = 5.4e-04 * (v^2) * rho_gas * (L/G)'
        2a. Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * (0.56 + 0.125*L/G + 0.0023*(L/G)^2) / 507'
        2b. Simplified Hesketh: '(v^2) * rho_gas * (Throat_Area^0.133) * ((L/G)^0.78) / 1270'
                """



if __name__ == '__main__':
#    import doctest
#    doctest.testmod()


    from lib.corriente import Solid
    from lib.psycrometry import Punto_Psicrometrico
    diametros=[96.5, 105, 110, 118, 125, 130, 140, 150, 170]
    fraccion=[0.02, 0.05, 0.1, 0.15, 0.25, 0.2, 0.15, 0.05, 0.03]
    solido=Solid(caudalSolido=[5000], distribucion_fraccion=fraccion, distribucion_diametro=diametros)
    Solido=Corriente(T=300, P=101325., caudalMasico=50, fraccionMolar=[1, 0], solido=solido)
    aire=Punto_Psicrometrico(caudal=100, tdb=300, HR=50)
    secador=Dryer(entradaSolido=Solido, entradaAire=aire)

