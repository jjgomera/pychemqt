#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=attribute-defined-outside-init

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2025, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Module with solid definition

    * :class:`Solid`: Solid entity

"""

from ast import literal_eval

from numpy import log, exp, r_
from scipy.optimize import leastsq
from scipy.special import erf
from tools.qt import translate

from lib.compuestos import Componente
from lib.config import Entity, getMainWindowConfig
from lib.unidades import Density, MassFlow, Length, Temperature, VolFlow


class Solid(Entity):
    """Class to define a solid entity
    Parameters:
        tipo: definition of solid
            0 - Undefined
            1 - Mean diameter
            2 - Particle solid distribution
        kwargs:
            caudalSolido
            diametroMedio
            fraccionMasica
            diametros
            T

            solidos: To override the config value
        """
    kwargs = {"caudalSolido": [],
              "diametroMedio": 0.0,
              "distribucion_fraccion": [],
              "distribucion_diametro": [],

              "solids": None}

    status = 0
    msg = translate("Solid", "undefined")

    def __call__(self, **kwargs):
        """All equipment are callables, so we can instance or add/change
        input value with flexibility"""
        Entity.__call__(self, **kwargs)
        if self._oldkwargs != self.kwargs and self.isCalculable:
            self.calculo()

    @property
    def isCalculable(self):
        """Procedure to check complete definition of instance"""
        self.status = 0
        if sum(self.kwargs["caudalSolido"]) > 0:
            if self.kwargs["distribucion_fraccion"] and \
                    self.kwargs["distribucion_diametro"]:
                self.status = 2
            elif self.kwargs["diametroMedio"]:
                self.status = 1
        return self.status

    def calculo(self):
        """Procedure to calculate instance properties"""
        if self.kwargs["solids"] is not None:
            self.ids = self.kwargs["solids"]
        else:
            Config = getMainWindowConfig()
            txt = Config.get("Components", "Solids")
            if isinstance(txt, str):
                self.ids = literal_eval(txt)
            else:
                self.ids = txt
        self.componente = [Componente(int(i)) for i in self.ids]

        caudal = self.kwargs.get("caudalSolido", [])
        diametro_medio = self.kwargs.get("diametroMedio", 0.0)
        fraccion = self.kwargs.get("distribucion_fraccion", [])
        dms = self.kwargs.get("distribucion_diametro", [])

        if self.status == 0:
            self._bool = False
            return

        self._bool = True

        self.caudalUnitario = [MassFlow(i) for i in caudal]
        self.caudal = MassFlow(sum(self.caudalUnitario))
        self.diametros = dms
        self.fracciones = fraccion
        if self.status == 2:
            self.diametros = [Length(i, "m", magnitud="ParticleDiameter")
                              for i in dms]
            self.fracciones = fraccion
            diametro_medio = 0
            self.fracciones_acumuladas = [0]
            for di, xi in zip(dms, fraccion):
                diametro_medio += di*xi
                self.fracciones_acumuladas.append(
                    xi+self.fracciones_acumuladas[-1])
            del self.fracciones_acumuladas[0]
        self.diametro_medio = Length(diametro_medio,
                                     magnitud="ParticleDiameter")
        self.RhoS(self.kwargs.get("T", 300))
        self.Q = VolFlow(self.caudal/self.rho)

    def RhoS(self, T):
        """Calculate solid density as a posible function of temperature"""
        densidad = 0
        for i in range(len(self.ids)):
            densidad += self.caudalUnitario[i]/self.caudal * \
                self.componente[i].RhoS(T)
        self.rho = Density(densidad)
        self.T = Temperature(T)

    def __repr__(self):
        if self.status:
            return f"Solid: Q:{self.caudal.str}, Dm:{self.diametro_medio.str}"

        return f"{self.__class__} empty"

    def ajustar_distribucion(self, eq=0):
        """
        Fit current distribution with any of standard.
        eq: index with distribuciont to fit
            0 - Rosin, Rammler, Sperling (Weibull distribution)
            1 - Gates, Gaudin, Shumann
            2 - Gaudin, Meloy
            3 - Broadbent, Callcott
            4 - Distribución lognormal
            5 - Harris
        Ref: Ahmed, Drzymala; Two-dimensional fractal linearization of
        distribution curves, pag 2
        """
        d = r_[self.diametros]
        y = r_[self.fracciones_acumuladas]
        if eq == 0:
            model = "Rosin, Rammler, Sperling"
            inicio = r_[1, 1]

            def function(p, d):
                return 1.-exp(-(d/p[0])**p[1])

        elif eq == 1:
            model = "Gates, Gaudin, Shumann"
            inicio = r_[1, 1]

            def function(p, d):
                return (d/p[0])**p[1]

        elif eq == 2:
            model = "Gaudin, Meloy"
            inicio = r_[d[-1]*2, 0]

            def function(p, d):
                return 1-(1-d/p[0])**p[1]

        elif eq == 3:
            model = "Broadbent, Callcott"
            inicio = r_[1, 1]

            def function(p, d):
                return 1.-exp(-(d/p[0])**p[1])/(1-exp(-1.))

        elif eq == 4:
            model = "lognormal"
            inicio = r_[1, 1]

            def function(p, d):
                return erf(log(d/p[0])/p[1])

        else:
            model = "Harris"
            inicio = r_[d[-1]*2, 1, 1]

            def function(p, d):
                return 1-(1-d/p[0]**p[1])**p[2]

        def residuo(p, d, y):
            return function(p, d) - y
        ajuste = leastsq(residuo, inicio, args=(d, y))
        return d, function(ajuste[0], d), model

    def Separar(self, etas):
        """Split solid with efficiency array input
        return two array with solids filtered and no filtered"""
        rendimiento_global = 0
        for i, fraccion in enumerate(self.fracciones):
            rendimiento_global += etas[i]*fraccion

        G_skip = MassFlow(self.caudal*(1-rendimiento_global))
        G_sep = MassFlow(self.caudal*rendimiento_global)
        if rendimiento_global == 1:
            return None, self

        if rendimiento_global == 0:
            return self, None

        f_gas = []
        f_solid = []
        for i in range(len(self.diametros)):
            f_gas.append(self.caudal*self.fracciones[i]*(1-etas[i])/G_skip)
            f_solid.append(self.caudal*self.fracciones[i]*etas[i]/G_sep)
        S_skip = Solid(caudalSolido=[G_skip],
                       distribucion_diametro=self.diametros,
                       distribucion_fraccion=f_gas)
        S_sep = Solid(caudalSolido=[G_sep],
                      distribucion_diametro=self.diametros,
                      distribucion_fraccion=f_solid)
        return S_skip, S_sep

    def writeStatetoJSON(self, data):
        """Write entity state to JSON file"""
        if self.status:
            data["status"] = self.status
            data["ids"] = self.ids
            data["unitFlow"] = self.caudalUnitario
            data["caudal"] = self.caudal
            data["diametros"] = self.diametros
            data["fracciones"] = self.fracciones
            data["fracciones_acumuladas"] = self.fracciones_acumuladas
            data["diametro_medio"] = self.diametro_medio
            data["rho"] = self.rho
            data["T"] = self.T

    def readStatefromJSON(self, data):
        """Read entity state from JSON file"""
        if data:
            self._bool = True
            self.status = data["status"]
            self.ids = data["ids"]
            self.componente = [Componente(int(i)) for i in self.ids]
            self.caudalUnitario = [MassFlow(q) for q in data["unitFlow"]]
            self.caudal = MassFlow(data["caudal"])
            self.diametros = [Length(d, "m", "ParticleDiameter")
                              for d in data["diametros"]]
            self.fracciones = data["fracciones"]
            self.fracciones_acumuladas = data["fracciones_acumuladas"]
            self.diametro_medio = Length(data["diametro_medio"])
            self.rho = Density(data["rho"])
            self.T = Temperature(data["T"])
        else:
            self._bool = False
            self.status = False


if __name__ == '__main__':
    distribucion = [[17.5, 0.02],
                    [22.4, 0.03],
                    [26.2, 0.05],
                    [31.8, 0.1],
                    [37, 0.1],
                    [42.4, 0.1],
                    [48, 0.1],
                    [54, 0.1],
                    [60, 0.1],
                    [69, 0.1],
                    [81.3, 0.1],
                    [96.5, 0.05],
                    [109, 0.03],
                    [127, 0.02]]
    diametros = []
    fracciones = []
    for diametro, frac in distribucion:
        diametros.append(diametro)
        fracciones.append(frac)

    solido = Solid(caudalSolido=[5], distribucion_diametro=diametros,
                   distribucion_fraccion=fracciones)
