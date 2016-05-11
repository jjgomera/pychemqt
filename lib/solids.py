#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2016, Juan José Gómez Romera <jjgomera@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>."""


###############################################################################
# Module with solid  definition
#   -Solid: Solid entity
###############################################################################


from scipy import log, exp, r_
from scipy.optimize import leastsq
from scipy.special import erf
from PyQt5.QtWidgets import QApplication

from lib.compuestos import Componente
from lib.config import Entity, getMainWindowConfig
from lib.unidades import Density, MassFlow, Length


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
        """
    kwargs = {"caudalSolido": [],
              "diametroMedio": 0.0,
              "distribucion_fraccion": [],
              "distribucion_diametro": []}
    status = 0
    msg = QApplication.trnaslate("pychemqt", "undefined")

    def __call__(self, **kwargs):
        """All equipment are callables, so we can instance or add/change
        input value with flexibility"""
        Entity.__call__(self, **kwargs)
        if self._oldkwargs != self.kwargs and self.isCalculable:
            self.calculo()

    @property
    def isCalculable(self):
        self.status = 0
        if sum(self.kwargs["caudalSolido"]) > 0:
            if self.kwargs["distribucion_fraccion"] and \
                    self.kwargs["distribucion_diametro"]:
                self.status = 2
            elif self.kwargs["diametroMedio"]:
                self.status = 1
        return self.status

    def calculo(self):
        Config = getMainWindowConfig()
        txt = Config.get("Components", "Solids")
        if isinstance(txt, str):
            self.ids = eval(txt)
        else:
            self.ids = txt
        self.componente = [Componente(int(i)) for i in self.ids]

        caudal = self.kwargs.get("caudalSolido", [])
        diametro_medio = self.kwargs.get("diametroMedio", 0.0)
        fraccion = self.kwargs.get("distribucion_fraccion", [])
        diametros = self.kwargs.get("distribucion_diametro", [])

        if self.status == 0:
            self._bool = False
            return
        else:
            self._bool = True

        self.caudalUnitario = [MassFlow(i) for i in caudal]
        self.caudal = MassFlow(sum(self.caudalUnitario))
        self.diametros = diametros
        self.fracciones = fraccion
        if self.status == 2:
            self.diametros = [Length(i, magnitud="ParticleDiameter")
                              for i in diametros]
            self.fracciones = fraccion
            diametro_medio = 0
            self.fracciones_acumuladas = [0]
            for di, xi in zip(diametros, fraccion):
                diametro_medio += di*xi
                self.fracciones_acumuladas.append(
                    xi+self.fracciones_acumuladas[-1])
            del self.fracciones_acumuladas[0]
        self.diametro_medio = Length(diametro_medio,
                                     magnitud="ParticleDiameter")

    def RhoS(self, T):
        densidad = 0
        for i in range(len(self.ids)):
            densidad += self.caudalUnitario[i]/self.caudal * \
                self.componente[i].RhoS(T)
        self.rho = Density(densidad)

    def __repr__(self):
        if self.status:
            return "Solid with %s and dm %s" % (self.caudal.str,
                                                self.diametro_medio.str)
        else:
            return "%s empty" % (self.__class__)

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

        elif eq == 5:
            model = "Harris"
            inicio = r_[d[-1]*2, 1, 1]

            def function(p, d):
                return 1-(1-d/p[0]**p[1])**p[2]

        def residuo(p, d, y):
            return function(p, d) - y
        ajuste, exito = leastsq(residuo, inicio, args=(d, y))
        return d, function(ajuste, d), model

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
        elif rendimiento_global == 0:
            return self, None
        else:
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

    def writeStatetoJSON(self, solid):
        if self.status:
            solid["ids"] = self.ids
            solid["unitFlow"] = self.caudalUnitario
            solid["caudal"] = self.caudal
            solid["diametros"] = self.diametros
            solid["fracciones"] = self.fracciones
            solid["fracciones_acumuladas"] = self.fracciones_acumuladas
            solid["diametro_medio"] = self.diametro_medio

    def readStatefromJSON(self, solid):
        if solid:
            self._bool = True
            self.status = solid["status"]
            self.ids = solid["ids"]
            self.componente = [Componente(int(i)) for i in self.ids]
            self.caudalUnitario = [MassFlow(q) for q in solid["unitFlow"]]
            self.caudal = MassFlow(solid["caudal"])
            self.diametros = [Length(d) for d in solid["diametros"]]
            self.fracciones = solid["fracciones"]
            self.fracciones_acumuladas = solid["fracciones_acumuladas"]
            self.diametro_medio = Length(solid["dm"])
        else:
            self._bool = False
            self.status = False


if __name__ == '__main__':
    distribucion = [[17.5, 0.02],
                    [22.4, 0.03],
                    [26.2,  0.05],
                    [31.8,  0.1],
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
    for diametro, fraccion in distribucion:
        diametros.append(diametro)
        fracciones.append(fraccion)

    solido = Solid(caudalSolido=[5], distribucion_diametro=diametros,
                   distribucion_fraccion=fracciones)
    print solido._def
