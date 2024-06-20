#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2023, Juan José Gómez Romera <jjgomera@gmail.com>

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



############################
###  librería de definición de reactores   ###
############################

from scipy.optimize import fsolve
from tools.qt import QtWidgets

from lib import unidades
from lib.corriente import Corriente
from lib.reaction import Reaction
from .parents import equipment


class Reactor(equipment):
    """Clase que define un reactor

    Parámetros:
        entrada: Instancia de clase corriente que define la corriente que fluye por la tubería
        thermal: Comportamiento térmico
            0   -   Adiabático
            1   -   Isotérmico
            2   -   Flujo de calor
            3   -   Calcular el intercambio de calor conocidas U y Tª externa
        reaccion: lista con instanacias de clase Reaction
        P: Presión en el reactor
        T: temperatura en el reactor
        Q: Calor transmitido a través de las paredes de la tubería
        Text: Temperatura en el exterior de la tubería
        U: Coeficiente global de transimisión de calor entre la tubería y el exterior
    """
    title = QtWidgets.QApplication.translate("equipment", "Reactor")
    help=""
    kwargs = {"entrada": None,
                    "thermal": 0,
                    "deltaP": 0.0,
                    "Pout": 0.0,
                    "Hmax": 0.0,
                    "eficiencia": 0.0,
                    "poderCalorifico": 0.0,

                    "f_install": 1.3,
                    "Base_index": 0.0,
                    "Current_index": 0.0,
                    "tipo": 0,
                    "subtipo": 0,
                    "material": 0,
                    "P_dis": 0.0}

    kwargsInput=("entrada", )

    @property
    def isCalculable(self):
        return True

    def calculo(self):
        self.entrada=self.kwargs.get("entrada", None)
        self.thermal=self.kwargs.get("thermal", 0)

        #TODO: implementar para más de una reacción
        self.reaccion=self.kwargs.get("reaccion", 0)[0]


#        P=self.kwargs.get("P", 0)
#        if P:
#            self.entrada=Corriente(self.entrada.T, P, self.entrada.caudalmasico.kgh, self.entrada.mezcla, self.entrada.solido)
        T=self.kwargs.get("T", 0)
        if T:
            self.T=unidades.Temperature(T)
        else:
            self.T=self.entrada.T

        self.Q=unidades.Power(self.kwargs.get("Q", 0))

        self.Text=self.kwargs.get("Text", 0)
        self.U=unidades.HeatTransfCoef(self.kwargs.get("U", 0))

        if self.thermal in [0, 2]:

            def f(T):
                fracciones, h=self.reaccion.conversion(self.entrada, T)
                corriente=Corriente(T, self.entrada.P.atm, self.entrada.caudalmasico.kgh, Mezcla(self.entrada.ids, fracciones), self.entrada.solido)
                return corriente.h-self.Q-self.entrada.h-h
            T=fsolve(f, self.entrada.T)
            fracciones, h=self.reaccion.conversion(self.entrada, T)

        elif self.thermal==1:
            T=self.T
            fracciones, h=self.reaccion.conversion(self.entrada, T)

        elif self.thermal==3:
            pass

        print(fracciones)
        self.Salida=Corriente(T=T, P=self.entrada.P, caudalMasico=self.entrada.caudalmasico, fraccionMolar=fracciones, solido=self.entrada.solido)
        self.Heat=unidades.Power(self.Salida.h-self.entrada.h-h)



if __name__ == '__main__':
    from math import exp, log
    mezcla=Corriente(T=300, P=101325., caudalMasico=1.0, ids=[1, 46, 47, 62], fraccionMolar=[0.03, 0.96, 0.01, 0.])
    reaccion=Reaction(comp=[1, 46, 47, 62], coef=[-2, 0, -1, 2], tipo=0, base=0, conversion=0.9)
    print(reaccion)
    reactor=Reactor(entrada=mezcla, reaccion=[reaccion], thermal=1)
    print(reactor.status)
    print(reactor.Salida.fraccion, reactor.Salida.T, reactor.Heat.MJh)

