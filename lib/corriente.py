#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.



Module with stream definition:

  * :class:`Corriente`: Stream general class model
  * :class:`PsyStream`: Stream specified as psychrometric state

"""

import logging
import os

from scipy.constants import R
from tools.qt import QtWidgets, translate

from lib.physics import R_atml
from lib import unidades, config
from lib import EoS, mEoS, gerg, iapws97, freeSteam, refProp, coolProp
from lib.solids import Solid
from lib.mezcla import Mezcla, mix_molarflow_molarfraction
from lib.psycrometry import PsychroState
from lib.thermo import ThermoWater, ThermoAdvanced, ThermoRefProp


class Corriente(config.Entity):
    """Class to model a stream object
    Parameters:
        -T: temperature, Kelvin
        -P: Pressure, Pa
        -x: quality

        -caudalMasico: mass flow in kg/s (solid component excluded)
        -caudalMolar: molar flow in kmol/s (solid component excluded)
        -caudalVolumetrico: volumetric flow in m3/s
        -caudalUnitarioMolar: Array with molar flows for each component
        -caudalUnitarioMasico: Array with mass flows for each component
        -fraccionMolar: Array with molar fractions for each component
        -fraccionMasica: Array with mass fractions for each component
        -mezcla: instance of class Mezcla to define all mixture variables

        -solido: Instance of class Solid to define all solid variables
        -caudalSolido: mass flow, in kg/h for solids component, array for each
            solid component
        -diametroMedio: mean diameter for particules, in micrometer
        -distribucion_fraccion: Array with mass fraccion of solid particle
            distribution
        -distribucion_diametro: Array with particle diameter of solid particle
            distribution, in micrometer

        -notas: Description text for stream

    Aditional parameter for define a corriente with different properties than
    the  project config, components for debug, thermodynamic method for add per
    stream configuration:
        -ids: Array with index of components
        -solids: Array with index of solids components

        -K: Name of method for K values calculation, ej: "SRK", "Lee-Kesler"
        -alfa: Name of method for alpha calculation, available only for some K
            methods
        -mix: Mixing rule for calculate mix properties
        -H: Name of method for enthalpy calculation, ej: "SRK", "Lee-Kesler"
        -Cp_ideal: Bool to choose method for ideal cp:
            0 - Ideal parameters
            1 - DIPPR parameters
        -MEoS: Use meos equation if is available
        -iapws: Use iapws97 standard for water
        -GERG: Use GERG-2008 equation if is available
        -freesteam: Use freesteam external library for water
        -coolProp: Use coolProp external library if is available
        -refprop: Use refProp external library if is available
    """
    kwargs = {"T": 0.0,
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
              "distribucion_diametro": [],

              "ids": [],
              "solids": None,
              "K": "",
              "alfa": "",
              "mix": "",
              "H": "",
              "Cp_ideal": None,
              "MEoS": None,
              "iapws": None,
              "GERG": None,
              "freesteam": None,
              "coolProp": None,
              "refprop": None}

    status = 0
    msg = translate("Corriente", "Unknown variables")
    kwargs_forbidden = ["entrada", "mezcla", "solido"]
    solido = None

    def __init__(self, **kwargs):
        self.kwargs = Corriente.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        # Clean input parameters to get only the defined in the last run
        if kwargs.get("mezcla", None):
            kwargs.update(kwargs["mezcla"].kwargs)
        if kwargs.get("solido", None):
            kwargs.update(kwargs["solido"].kwargs)

        if kwargs.get("caudalUnitarioMasico", []):
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMolar"] = []
            self.kwargs["fraccionMasica"] = []
            self.kwargs["caudalMasico"] = None
            self.kwargs["caudalVolumetrico"] = None
            self.kwargs["caudalMolar"] = None

        elif kwargs.get("caudalUnitarioMolar", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["fraccionMolar"] = []
            self.kwargs["fraccionMasica"] = []
            self.kwargs["caudalMasico"] = None
            self.kwargs["caudalVolumetrico"] = None
            self.kwargs["caudalMolar"] = None

        elif kwargs.get("caudalMasico", None) and \
                kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMasica"] = []
            self.kwargs["caudalVolumetrico"] = None
            self.kwargs["caudalMolar"] = None

        elif kwargs.get("caudalMasico", None) and \
                kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMolar"] = []
            self.kwargs["caudalVolumetrico"] = None
            self.kwargs["caudalMolar"] = None

        elif kwargs.get("caudalMasico", None):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["caudalVolumetrico"] = None
            self.kwargs["caudalMolar"] = None

        elif kwargs.get("caudalMolar", None) and \
                kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMasica"] = []
            self.kwargs["caudalMasico"] = None
            self.kwargs["caudalVolumetrico"] = None

        elif kwargs.get("caudalMolar", None) and \
                kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMolar"] = []
            self.kwargs["caudalMasico"] = None
            self.kwargs["caudalVolumetrico"] = None

        elif kwargs.get("caudalMolar", None):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["caudalMasico"] = None
            self.kwargs["caudalVolumetrico"] = None

        elif kwargs.get("caudalVolumetrico", None) and \
                kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMasica"] = []
            self.kwargs["caudalMasico"] = None

        elif kwargs.get("caudalVolumetrico") and \
                kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMolar"] = []
            self.kwargs["caudalMasico"] = None

        elif kwargs.get("caudalVolumetrico", None):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["caudalMolar"] = []
            self.kwargs["caudalMasico"] = None

        elif kwargs.get("fraccionMasica", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMolar"] = []

        elif kwargs.get("fraccionMolar", []):
            self.kwargs["caudalUnitarioMasico"] = []
            self.kwargs["caudalUnitarioMolar"] = []
            self.kwargs["fraccionMasica"] = []

        elif kwargs.get("x", None) and self.kwargs["T"] and self.kwargs["P"]:
            self.kwargs["T"] = 0.0
        elif kwargs.get("T", 0.0) and self.kwargs["x"] and self.kwargs["P"]:
            self.kwargs["P"] = 0.0
        elif kwargs.get("P", 0.0) and self.kwargs["T"] and self.kwargs["x"]:
            self.kwargs["x"] = None

        self.kwargs.update(kwargs)

        for key, value in list(self.kwargs.items()):
            if value:
                self._bool = True
                break

        logging.info('Calculate STREAM')
        kw_new = {}
        for key, value in list(kwargs.items()):
            if self.__class__.kwargs[key] != value:
                kw_new[key] = value
        logging.debug('kwarg; %s' % kw_new)
        if self.calculable:
            statusmsg = (
                translate("Corriente", "Underspecified"),
                translate("Corriente", "Solved"),
                translate("Corriente", "Ignored"),
                translate("Corriente", "Warning"),
                translate("Corriente", "Calculating..."),
                translate("Corriente", "Error"))
            status = statusmsg[self.status]
            logging.debug('%s %s' % (status, self.msg))
            QtWidgets.QApplication.processEvents()

            self.status = 1
            self.calculo()
            self.msg = ""

        elif self.tipoFlujo:
            if self.kwargs["mezcla"]:
                self.mezcla = self.kwargs["mezcla"]
            else:
                self.mezcla = Mezcla(self.tipoFlujo, **self.kwargs)

        elif self.tipoSolido:
            if self.kwargs["solido"]:
                self.solido = self.kwargs["solido"]
            else:
                self.solido = Solid(**self.kwargs)
            if self.solido:
                self.solido.RhoS(self.kwargs["T"])

    @property
    def calculable(self):
        # Thermo definition
        self.tipoTermodinamica = ""
        if self.kwargs["T"] and self.kwargs["P"]:
            self.tipoTermodinamica = "TP"
        elif self.kwargs["T"] and self.kwargs["x"]:
            self.tipoTermodinamica = "Tx"
        elif self.kwargs["P"] and self.kwargs["x"]:
            self.tipoTermodinamica = "Px"

        # Mix definition
        self.tipoFlujo = 0
        if self.kwargs["caudalUnitarioMasico"]:
            self.tipoFlujo = 1

        elif self.kwargs["caudalUnitarioMolar"]:
            self.tipoFlujo = 2

        elif self.kwargs["caudalMasico"] and self.kwargs["fraccionMolar"]:
            self.tipoFlujo = 3

        elif self.kwargs["caudalMasico"] and self.kwargs["fraccionMasica"]:
            self.tipoFlujo = 4

        elif self.kwargs["caudalMolar"] and self.kwargs["fraccionMolar"]:
            self.tipoFlujo = 5

        elif self.kwargs["caudalMolar"] and self.kwargs["fraccionMasica"]:
            self.tipoFlujo = 6

        elif self.kwargs["caudalVolumetrico"] and self.kwargs["fraccionMolar"]:
            self.kwargs["caudalMolar"] = 1
            self.tipoFlujo = 5

        elif self.kwargs["caudalVolumetrico"] and \
                self.kwargs["fraccionMasica"]:
            self.kwargs["caudalMolar"] = 1
            self.tipoFlujo = 6

        elif self.kwargs["mezcla"]:
            self.tipoFlujo = 7

        # Solid definition
        self.tipoSolido = 0
        if sum(self.kwargs["caudalSolido"]) > 0:
            if self.kwargs["distribucion_fraccion"] and \
                    self.kwargs["distribucion_diametro"]:
                self.tipoSolido = 2
            elif self.kwargs["diametroMedio"]:
                self.tipoSolido = 1
        if self.kwargs["solido"]:
            self.tipoSolido = self.kwargs["solido"].status
        return self.tipoTermodinamica and self.tipoFlujo

    def calculo(self):
        Config = config.getMainWindowConfig()
        if self.kwargs["mezcla"]:
            self.mezcla = self.kwargs["mezcla"]
        else:
            self.mezcla = Mezcla(self.tipoFlujo, **self.kwargs)

        self.ids = self.mezcla.ids
        self.componente = self.mezcla.componente
        self.fraccion = self.mezcla.fraccion
        self.caudalmasico = self.mezcla.caudalmasico
        self.caudalmolar = self.mezcla.caudalmolar
        self.fraccion_masica = self.mezcla.fraccion_masica
        self.caudalunitariomasico = self.mezcla.caudalunitariomasico
        self.caudalunitariomolar = self.mezcla.caudalunitariomolar

        T = unidades.Temperature(self.kwargs.get("T", None))
        P = unidades.Pressure(self.kwargs.get("P", None))
        x = self.kwargs.get("x", None)

        self._method()
        setData = True

        if self._thermo == "freesteam":
            compuesto = freeSteam.Freesteam(**self.kwargs)
        elif self._thermo == "iapws":
            compuesto = iapws97.IAPWS97(**self.kwargs)
        elif self._thermo == "refprop":
            if not self.kwargs["ids"]:
                self.kwargs["ids"] = self.ids

            # Avoid overwrite refprop H parameter
            kwargs = self.kwargs.copy()
            del kwargs["H"]

            compuesto = refProp.RefProp(**kwargs)
        elif self._thermo == "gerg":
            ids = []
            for id in self.ids:
                ids.append(gerg.id_GERG.index(id))
            kwargs = self.kwargs
            kwargs["mezcla"] = self.mezcla
            compuesto = gerg.GERG(
                    componente=ids, fraccion=self.fraccion, **kwargs)
        elif self._thermo == "coolprop":
            if not self.kwargs["ids"]:
                self.kwargs["ids"] = self.ids
            compuesto = coolProp.CoolProp(**self.kwargs)
        elif self._thermo == "meos":
            class_ = mEoS.__all__[mEoS.id_mEoS.index(self.ids[0])]
            if self.tipoTermodinamica == "TP":
                compuesto = class_(T=T, P=P)
            elif self.tipoTermodinamica == "Tx":
                compuesto = class_(T=T, x=x)
            elif self.tipoTermodinamica == "Px":
                compuesto = class_(P=P, x=x)
        elif self._thermo == "eos":
            if self.kwargs["K"]:
                index = EoS.K_status.index(self.kwargs["K"])
                K = EoS.K[index]
            else:
                K = EoS.K[Config.getint("Thermo", "K")]
            if self.kwargs["H"]:
                index = EoS.H_status.index(self.kwargs["H"])
                H = EoS.H[index]
            else:
                H = EoS.H[Config.getint("Thermo", "H")]

            setData = False
            self.M = unidades.Dimensionless(self.mezcla.M)
            self.Tc = self.mezcla.Tc
            self.Pc = self.mezcla.Pc
            self.SG = unidades.Dimensionless(self.mezcla.SG)

            if self.tipoTermodinamica == "TP":
                self.T = unidades.Temperature(T)
                self.P = unidades.Pressure(P)
                eos = K(self.T, self.P, self.mezcla)
                self.eos = eos
                self.x = unidades.Dimensionless(eos.x)
            else:
                self.x = unidades.Dimensionless(x)

#             self.mezcla.recallZeros(eos.xi)
#             self.mezcla.recallZeros(eos.yi)
#             self.mezcla.recallZeros(eos.Ki, 1.)

            if 0. < self.x < 1.:
                self.Liquido = Mezcla(
                    tipo=5, ids=self.mezcla.ids, fraccionMolar=eos.xi,
                    caudalMolar=self.caudalmolar*(1-self.x))
                self.Gas = Mezcla(
                    tipo=5, ids=self.mezcla.ids,  fraccionMolar=eos.yi,
                    caudalMolar=self.caudalmolar*self.x)
            elif self.x <= 0:
                self.Liquido = self.mezcla
                self.Gas = Mezcla()
            else:
                self.Liquido = Mezcla()
                self.Gas = self.mezcla
#             self.Gas.Z = unidades.Dimensionless(eos.Zg)
#             self.Liquido.Z = unidades.Dimensionless(eos.Zl)

            # if H == K:
                # eosH = eos
            # else:
                # eosH = H(self.T, self.P.atm, self.mezcla)
            # self.H_exc = eosH.H_exc
            self.H_exc = (0, 0)

            self.Liquido.Q = unidades.VolFlow(0)
            self.Gas.Q = unidades.VolFlow(0)
            self.Liquido.h = unidades.Power(0)
            self.Gas.h = unidades.Power(0)

            if self.x < 1:
                # There is liquid phase

                if Config.getboolean("Transport", "RhoLEoS") and eos.rhoL:
                    self.Liquido.rho = eos.rhoL
                    self.Liquido.Z = unidades.Dimensionless(float(eos.Zl))
                else:
                    self.Liquido.rho = self.Liquido.RhoL(T, P)
                    Z = self.Liquido.rho*R*T/P
                    self.Liquido.Z = unidades.Dimensionless(Z)

                # Hl = (self.Liquido._Ho(self.T).Jg-self.Liquido.Hv_DIPPR(self.T).Jg)*self.Liquido.caudalmasico.gh
                Hl = 0
                self.Liquido.h = unidades.Power(Hl-R*self.T/self.M*self.H_exc[1]*(1-self.x)*self.Liquido.caudalmasico.gh, "Jh")
                self.Liquido.cp = self.Liquido.Cp_Liquido(T)


                self.Liquido.Q = unidades.VolFlow(self.Liquido.caudalmasico/self.Liquido.rho)

                # Transport properties
                self.Liquido.mu = self.Liquido.Mu_Liquido(T, P)
                k = self.Liquido.ThCond_Liquido(T, P, self.Liquido.rho)
                self.Liquido.k = k
                self.Liquido.sigma = self.Liquido.Tension(T)
                self.Liquido.Prandt = self.Liquido.cp*self.Liquido.mu/self.Liquido.k

            if self.x > 0:
                # There is gas phase

                if eos.rhoG:
                    self.Gas.Z = unidades.Dimensionless(float(eos.Zg))
                    self.Gas.rho = eos.rhoG
                else:
                    self.Gas.rho = self.Liquido.RhoL(T, P)
                    Z = self.Gas.rho*R*T/P
                    self.Gas.Z = unidades.Dimensionless(Z)
                self.Gas.rho = unidades.Density(self.P/self.Gas.Z/R/self.T*self.M, "gl")
                self.Gas.rhoSd = unidades.Density(1./self.Gas.Z/R_atml/298.15*self.M, "gl")

                Hg = self.Gas._Ho(self.T).Jg*self.Gas.caudalmasico.gh
                self.Gas.h = unidades.Power(Hg-R*self.T/self.M*self.H_exc[0]*self.x*self.Gas.caudalmasico.gh, "Jh")
                self.Gas.cp = self.Gas.Cp_Gas(T, self.P)
                self.Gas.mu = self.Gas.Mu_Gas(T, self.P, self.Gas.rho)
                self.Gas.k = self.Gas.ThCond_Gas(T, self.P, self.Gas.rho)
                self.Gas.Q = unidades.VolFlow(self.Gas.caudalmasico/self.Gas.rho)
                self.Gas.Prandt = self.Gas.cp*self.Gas.mu/self.Gas.k

            self.Q = unidades.VolFlow(self.Liquido.Q+self.Gas.Q)
            self.h = unidades.Power(self.Liquido.h+self.Gas.h)
            self.Molaridad = [caudal/self.Q.m3h for caudal in self.caudalunitariomolar]

            # print(self.Liquido.mu, self.Liquido.k)
#             self.Liquido.fraccion = self.Liquido.recallZeros(self.Liquido.fraccion)
#             print(self.Liquido.fraccion)
            # self.Liquido.recallZeros(eos.yi)
            # self.Liquido.recallZeros(eos.Ki, 1.)

            # TODO:
            self.cp_cv = 0.5
            self.cp_cv_ideal = 0.5
            self.s = 0
            self.rho = 0

        if setData:
            # Asignación de valores comun
            self.cmp = compuesto
            self.T = compuesto.T
            self.P = compuesto.P
            self.x = compuesto.x
            self.M = unidades.Dimensionless(compuesto.M)
            self.Tc = compuesto.Tc
            self.Pc = compuesto.Pc
            self.h = unidades.Power(compuesto.h*self.caudalmasico)
            self.s = unidades.Entropy(compuesto.s*self.caudalmasico)
            self.rho = compuesto.rho
            self.Q = unidades.VolFlow(compuesto.v*self.caudalmasico)

            # harcoded thermo derived global properties in corriente for avoid
            # losing data
            if self._thermo != "eos":
                compuesto._fillCorriente(self)

#        if self.__class__!=mEoS.H2O:
#            agua=mEoS.H2O(T=self.T, P=self.P.MPa)
#            self.SG=unidades.Dimensionless(self.rho/agua.rho)
#        else:
            self.SG = unidades.Dimensionless(1.)

            self.Liquido = compuesto.Liquido
            self.Gas = compuesto.Gas

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

            if self.x > 0:
                self.Gas.Q = unidades.VolFlow(self.Q*(1-self.x))
                self.Gas.caudalmasico = unidades.MassFlow(self.caudalmasico*self.x)
                self.Gas.caudalmolar = unidades.MolarFlow(self.caudalmolar*self.x)
                kw = mix_molarflow_molarfraction(
                    self.Gas.caudalmolar, self.Gas.fraccion, self.componente)
                self.Gas.caudalunitariomasico = [
                    unidades.MassFlow(f) for f in kw["unitMassFlow"]]
                self.Gas.caudalunitariomolar = [
                    unidades.MolarFlow(f) for f in kw["unitMolarFlow"]]
                self.Gas.ids = self.ids

            if self.x < 1:
                self.Liquido.Q = unidades.VolFlow(self.Q*(1-self.x))
                self.Liquido.caudalmasico = unidades.MassFlow(self.caudalmasico*(1-self.x))
                self.Liquido.caudalmolar = unidades.MolarFlow(self.caudalmolar*(1-self.x))
                kw = mix_molarflow_molarfraction(
                    self.Liquido.caudalmolar, self.Liquido.fraccion, self.componente)
                self.Liquido.caudalunitariomasico = [
                    unidades.MassFlow(f) for f in kw["unitMassFlow"]]
                self.Liquido.caudalunitariomolar = [
                    unidades.MolarFlow(f) for f in kw["unitMolarFlow"]]
                self.Liquido.sigma = compuesto.sigma
                self.Liquido.ids = self.ids

        if Config.get("Components", "Solids"):
            if self.kwargs["solido"]:
                self.solido = self.kwargs["solido"]
            else:
                self.solido = Solid(**self.kwargs)
            if self.solido.status:
                self.solido.RhoS(T)
        else:
            self.solido = None

        if self.kwargs["caudalVolumetrico"]:
            self.kwargs["caudalMolar"] *= self.kwargs["caudalVolumetrico"]/self.Q
            Q = self.kwargs["caudalVolumetrico"]
            self.kwargs["caudalVolumetrico"] = None
            self.calculo()
            self.kwargs["caudalVolumetrico"] = Q
            self.kwargs["caudalMolar"] = None

    def _method(self):
        """Find the thermodynamic method to use
        Define internal variables to know the definition:
            _thermo: Method of definition, one of this values:
                eos, iapws, freesteam, meos, coolprop, refprop, gerg
            _dependence: In advanced method with external dependencias define
                the neccesary dependences to calculate
        """
        Config = config.getMainWindowConfig()

        # MEoS availability,
        if self.kwargs["MEoS"] is not None:
            _meos = self.kwargs["MEoS"]
        else:
            _meos = Config.getboolean("Thermo", "MEoS")
        mEoS_available = self.ids[0] in mEoS.id_mEoS
        MEoS = _meos and len(self.ids) == 1 and mEoS_available

        # iapws availability
        if self.kwargs["iapws"] is not None:
            _iapws = self.kwargs["iapws"]
        else:
            _iapws = Config.getboolean("Thermo", "iapws")
        IAPWS = _iapws and len(self.ids) == 1 and self.ids[0] == 62

        # freesteam availability
        if self.kwargs["freesteam"] is not None:
            _freesteam = self.kwargs["freesteam"]
        else:
            _freesteam = Config.getboolean("Thermo", "freesteam")
        FREESTEAM = _freesteam and len(self.ids) == 1 and \
            self.ids[0] == 62 and os.environ["freesteam"]

        COOLPROP_available = True
        GERG_available = True
        REFPROP_available = True
        for id in self.ids:
            if id not in coolProp.__all__:
                COOLPROP_available = False
            if id not in refProp.__all__:
                REFPROP_available = False
            if id not in gerg.id_GERG:
                GERG_available = False

        # coolprop availability
        if self.kwargs["coolProp"] is not None:
            _coolprop = self.kwargs["coolProp"]
        else:
            _coolprop = Config.getboolean("Thermo", "coolprop")
        COOLPROP = _coolprop and os.environ["CoolProp"] and COOLPROP_available

        # refprop availability
        if self.kwargs["refprop"] is not None:
            _refprop = self.kwargs["refprop"]
        else:
            _refprop = Config.getboolean("Thermo", "refprop")
        REFPROP = _refprop and os.environ["refprop"] and REFPROP_available

        # GERG availability
        if self.kwargs["GERG"] is not None:
            _gerg = self.kwargs["GERG"]
        else:
            _gerg = Config.getboolean("Thermo", "GERG")
        GERG = _gerg and GERG_available

        # Final selection
        if IAPWS and FREESTEAM:
            self._thermo = "freesteam"
            self._dependence = "freesteam"
        elif IAPWS:
            self._thermo = "iapws"
        elif _meos and REFPROP:
            self._thermo = "refprop"
            self._dependence = "refprop"
        elif _meos and COOLPROP:
            self._thermo = "coolprop"
            self._dependence = "CoolProp"
        elif MEoS and GERG:
            self._thermo = "gerg"
        elif MEoS:
            self._thermo = "meos"
        else:
            self._thermo = "eos"

    def setSolid(self, solid):
        self.solido = solid

    @property
    def psystream(self):
        xw = self.fraccion_masica[self.ids.index(62)]
        xa = self.fraccion_masica[self.ids.index(475)]
        psystream = PsyStream(caudalMasico=self.caudalMasico, P=self.P,
                              tdb=self.T, w=xw/xa)
        return psystream

    def clone(self, **kwargs):
        """Create a new stream instance with change only kwags new values"""
        old_kwargs = self.kwargs.copy()
        if "split" in kwargs:
            split = kwargs["split"]
            del kwargs["split"]
            if self.kwargs["caudalUnitarioMasico"]:
                kwargs["caudalUnitarioMasico"] = []
                for caudal in self.kwargs["caudalUnitarioMasico"]:
                    kwargs["caudalUnitarioMasico"].append(split*caudal)
            if self.kwargs["caudalUnitarioMolar"]:
                kwargs["caudalUnitarioMolar"] = []
                for caudal in self.kwargs["caudalUnitarioMolar"]:
                    kwargs["caudalUnitarioMolar"].append(split*caudal)
            if self.kwargs["caudalMasico"]:
                kwargs["caudalMasico"] = split*self.kwargs["caudalMasico"]
            if self.kwargs["caudalMolar"]:
                kwargs["caudalMolar"] = split*self.kwargs["caudalMolar"]
        if "x" in kwargs:
            del old_kwargs["T"]
        if "mezcla" in kwargs:
            old_kwargs.update(kwargs["mezcla"].kwargs)
            del kwargs["mezcla"]
        old_kwargs.update(kwargs)
        return Corriente(**old_kwargs)

    def __repr__(self):
        """String representation of any instance"""
        if self.status:
            return "Corriente at %0.2fK and %0.2fatm" % (self.T, self.P.atm)
        else:
            return "%s empty" % (self.__class__)

    def txt(self):
        txt = str(self.notasPlain)+os.linesep+os.linesep
        txt += "#---------------"
        txt += translate("Corriente", "Input properties")
        txt += "-----------------#"+os.linesep
        for key, value in list(self.kwargs.items()):
            if value:
                txt += key+": "+str(value)+os.linesep

        if self.calculable:
            txt += os.linesep + "#---------------"
            txt += translate("Corriente", "Global stream")
            txt += "-------------------#"+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "Temperature"),
                self.T.str)+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "Pressure"),
                self.P.str)+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "Vapor Fraction"),
                self.x.str)+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "Molar Flow"),
                self.caudalmasico.str)+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "Mass Flow"),
                self.caudalmolar.str)+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "Volumetric Flow"),
                self.Q.str)+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "Enthalpy"),
                self.h.str)+os.linesep
            txt += "%-25s\t%s" % ("Tc", self.Tc.str)+os.linesep
            txt += "%-25s\t%s" % ("Pc", self.Pc.str)+os.linesep
            txt += "%-25s\t%s" % (
                translate("Corriente", "SG, water=1"),
                self.SG.str)+os.linesep
            txt += os.linesep+"%-25s\t%s" % (
                translate("Corriente", "Molecular weight"),
                self.M.str)+os.linesep
            txt += "#"+translate("Corriente", "Molar Composition")
            txt += os.linesep
            for cmp, xi in zip(self.componente, self.fraccion):
                txt += "%-25s\t %0.4f" % (cmp.name, xi)+os.linesep

            if self.x > 0:
                txt += os.linesep+"#---------------"
                txt += translate("Corriente", "Vapor Only")
                txt += "--------------------#"+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Molar Flow"),
                    self.Gas.caudalmasico.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Mass Flow"),
                    self.Gas.caudalmolar.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Volumetric Flow"),
                    self.Gas.Q.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Molecular weight"),
                    self.Gas.M.str)+os.linesep
                txt += os.linesep+"#"
                txt += translate("Corriente", "Molar Composition")
                txt += os.linesep
                for cmp, xi in zip(self.componente, self.Gas.fraccion):
                    txt += "%-25s\t %0.4f" % (cmp.name, xi)+os.linesep

                txt += os.linesep+"%-25s\t%s" % (
                    translate("Corriente", "Density"),
                    self.Gas.rho.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Compresibility"),
                    self.Gas.Z.str)+os.linesep

                txt += "%-25s\t%s" % (
                    translate("Corriente", "Enthalpy"),
                    self.Gas.h.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Heat Capacity"),
                    self.Gas.cp.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Viscosity"),
                    self.Gas.mu.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Thermal conductivity"),
                    self.Gas.k.str)+os.linesep

            if self.x < 1:
                txt += os.linesep+"#---------------"
                txt += translate("Corriente", "Liquid Only")
                txt += "-------------------#"+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Molar Flow"),
                    self.Liquido.caudalmasico.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Mass Flow"),
                    self.Liquido.caudalmolar.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Volumetric Flow"),
                    self.Liquido.Q.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Molecular weight"),
                    self.Liquido.M.str)+os.linesep
                txt += os.linesep+"#"
                txt += translate("Corriente", "Molar Composition")
                txt += os.linesep
                for cmp, xi in zip(self.componente, self.Liquido.fraccion):
                    txt += "%-25s\t %0.4f" % (cmp.name, xi)+os.linesep

                txt += os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Density"),
                    self.Liquido.rho.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Compresibility"),
                    self.Liquido.Z.str)+os.linesep

                txt += "%-25s\t%s" % (
                    translate("Corriente", "Enthalpy"),
                    self.Liquido.h.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Heat Capacity"),
                    self.Liquido.cp.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Viscosity"),
                    self.Liquido.mu.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Thermal Conductivity"),
                    self.Liquido.k.str)+os.linesep
                txt += "%-25s\t%s" % (
                    translate("Corriente", "Surface Tension"),
                    self.Liquido.sigma.str)+os.linesep

        else:
            txt += os.linesep+"#---------------"
            txt += translate("Corriente", "No Fluid Stream")
            txt += "-------------------#"+os.linesep

        if self.solido.status:
            txt += os.linesep+"#---------------"
            txt += translate("Corriente", "Solid")
            txt += "-------------------#"+os.linesep
            for cmp, G in zip(self.solido.componente, self.solido.caudalUnitario):
                txt += "%-25s\t%s" % (cmp.name, G.str)+os.linesep
            txt += os.linesep
            txt += "%-25s\t%s" % (translate("Corriente", "Density"),
                                  self.solido.rho.str)+os.linesep
            txt += "%-25s\t%s" % (translate("Corriente", "Mean Diameter"),
                                  self.solido.diametro_medio.str)+os.linesep
            if self.solido.diametros:
                txt += os.linesep + "#"
                txt += translate("Corriente", "Particle Size Distribution")
                txt += os.linesep
                txt += "%s, %s \t%s" % (translate("Corriente", "Diameter"), unidades.Length.text("ParticleDiameter"), translate("Corriente", "Fraction"))+os.linesep
                for di, xi in zip(self.solido.diametros, self.solido.fracciones):
                    txt += "%10.4f\t%0.4f\t" % (di.config("ParticleDiameter"),
                                                xi)+os.linesep

        # Add advanced properties to report
        if self.calculable and self._thermo != "eos":
            doc = self._doc()

            if 0 < self.x < 1:
                param = "%-40s\t%-20s\t%-20s"
            else:
                param = "%-40s\t%s"
            if self.x == 0:
                txtphases = "%60s" % translate("Corriente", "Liquid")+os.linesep
                phases = [self.Liquido]
            elif self.x == 1:
                txtphases = "%60s" % translate("Corriente", "Gas")+os.linesep
                phases = [self.Gas]
            else:
                txtphases = "%60s\t%20s" % (translate("Corriente", "Liquid"), translate("Corriente", "Gas"))+os.linesep
                phases = [self.Liquido, self.Gas]

            complejos = ""
            data = self.cmp.properties()
            for propiedad, key, unit in data:
                if key == "sigma":
                    if self.x < 1:
                        complejos += "%-40s\t%s" % (propiedad, self.Liquido.sigma.str)
                        complejos += os.linesep
                elif key in ["f", "fi"]:
                    complejos += propiedad + os.linesep
                    for i, cmp in enumerate(self.componente):
                        values = ["  " + cmp.name]
                        for phase in phases:
                            values.append(phase.__getattribute__(key)[i].str)
                        complejos += param % tuple(values) + os.linesep
                elif key in ["K", "csat", "dpdt_sat", "cv2p", "chempot"]:
                    complejos += propiedad + os.linesep
                    for i, cmp in enumerate(self.componente):
                        values = ["  " + cmp.name]
                        values.append(self.__getattribute__(key)[i].str)
                        complejos += "%-40s\t%s" % tuple(values)
                        complejos += os.linesep
                elif key in self.Gas.__dict__ or key in self.Liquido.__dict__:
                    values = [propiedad]
                    for phase in phases:
                        values.append(phase.__getattribute__(key).str)
                    complejos += param % tuple(values) + os.linesep
                else:
                    complejos += "%-40s\t%s" % (propiedad, self.__getattribute__(key).str)
                    complejos += os.linesep

            txt += doc + os.linesep + txtphases + complejos

        return txt

    def _doc(self):
        """Return a text repr of class with all properties"""
        if self._thermo == "meos":
            title = translate("Corriente", "Advanced MEoS properties")
            doc_param = [self.cmp._constants["__doi__"]]
        else:
            title = translate("Corriente", "Advanced thermo properties")
            doc_param = self.cmp.__doi__
        doc = ""
        for doi in doc_param:
            doc += doi["autor"] + "; " + doi["title"] + "; " + doi["ref"]
            doc += os.linesep

        txt = os.linesep + os.linesep + "#---------------"
        txt += title + "-------------------#" + os.linesep
        txt += doc
        return txt

    @classmethod
    def propertiesNames(cls):
        list = [
            (translate("Corriente", "Temperature"), "T", unidades.Temperature),
            (translate("Corriente", "Pressure"), "P", unidades.Pressure),
            (translate("Corriente", "Vapor Fraction"), "x", unidades.Dimensionless),
            (translate("Corriente", "Molar Flow"), "caudalmolar", unidades.MolarFlow),
            (translate("Corriente", "Mass Flow"), "caudalmasico", unidades.MassFlow),
            (translate("Corriente", "Volumetric Flow"), "Q", unidades.VolFlow),
            (translate("Corriente", "Enthalpy"), "h", unidades.Enthalpy),
            (translate("Corriente", "Critic Temperature"), "Tc", unidades.Temperature),
            (translate("Corriente", "Critic Pressure"), "Pc", unidades.Pressure),
            (translate("Corriente", "SG, water=1"), "SG", unidades.Dimensionless),
            (translate("Corriente", "Molecular weight"), "M", unidades.Dimensionless),
            (translate("Corriente", "Molar Composition"), "fraccion", unidades.Dimensionless),
            (translate("Corriente", "Mass Composition"), "fraccion_masica", unidades.Dimensionless),
            (translate("Corriente", "Molar Component Flow"), "caudalunitariomolar", unidades.MolarFlow),
            (translate("Corriente", "Mass Component Flow"),  "caudalunitariomasico", unidades.MassFlow),
            (translate("Corriente", "Notes"), "notasPlain", str)]
        return list

    def propertiesListTitle(self, index):
        """Define los titulos para los popup de listas"""
        lista = [comp.name for comp in self.componente]
        return lista

    def writeToJSON(self, data):
        """Read entity from file"""
        config.Entity.writeToJSON(self, data)

        # Solid state, can be defined without a thermo status
        solid = {}
        if self.solido is not None:
            self.solido.writeStatetoJSON(solid)
        data["solid"] = solid

    def writeStatetoJSON(self, state):
        state["thermo"] = self._thermo
        state["bool"] = self._bool
        state["thermoType"] = self.tipoTermodinamica
        state["T"] = self.T
        state["P"] = self.P
        state["x"] = self.x
        state["M"] = self.M
        state["Tc"] = self.Tc
        state["Pc"] = self.Pc
        state["h"] = self.h
        state["s"] = self.s
        state["rho"] = self.rho
        state["Q"] = self.Q
        state["SG"] = self.SG

        if self._thermo != "eos":
            self.cmp._writeGlobalState(self, state)

        state["fluxType"] = self.tipoFlujo
        self.mezcla.writeStatetoJSON(state)

        self.Liquido.writeStatetoJSON(state, "liquid")
        self.Gas.writeStatetoJSON(state, "gas")
        if state["liquid"]:
            state["liquid"]["sigma"] = self.Liquido.sigma

        if self._thermo == "meos":
            state["meos_eq"] = self.cmp.kwargs["eq"]

    def readFromJSON(self, data):
        """Read entity from file"""
        config.Entity.readFromJSON(self, data)

        # Read solid
        self.solido = Solid()
        self.solido.readStatefromJSON(data["solid"])

    def readStatefromJSON(self, state):
        self._thermo = state["thermo"]
        self._bool = state["bool"]
        self.tipoTermodinamica = state["thermoType"]
        self.T = unidades.Temperature(state["T"])
        self.P = unidades.Pressure(state["P"])
        self.x = unidades.Dimensionless(state["x"])
        self.M = unidades.Dimensionless(state["M"])
        self.Tc = unidades.Temperature(state["Tc"])
        self.Pc = unidades.Pressure(state["Pc"])
        self.h = unidades.Power(state["h"])
        self.s = unidades.Entropy(state["s"])
        self.rho = unidades.Density(state["rho"])
        self.Q = unidades.VolFlow(state["Q"])
        self.SG = unidades.Dimensionless(state["SG"])

        self.tipoFlujo = state["fluxType"]
        self.mezcla = Mezcla()
        self.mezcla.readStatefromJSON(state["mezcla"])
        self.ids = self.mezcla.ids
        self.componente = self.mezcla.componente
        self.fraccion = self.mezcla.fraccion
        self.caudalmasico = self.mezcla.caudalmasico
        self.caudalmolar = self.mezcla.caudalmolar
        self.fraccion_masica = self.mezcla.fraccion_masica
        self.caudalunitariomasico = self.mezcla.caudalunitariomasico
        self.caudalunitariomolar = self.mezcla.caudalunitariomolar

        # Define the advanced thermo component if is necessary
        if self._thermo == "freesteam":
            self.cmp = freeSteam.Freesteam()
        elif self._thermo == "iapws":
            self.cmp = iapws97.IAPWS97()
        elif self._thermo == "refprop":
            self.cmp = refProp.RefProp(ids=self.ids)
        elif self._thermo == "coolprop":
            self.cmp = coolProp.CoolProp(ids=self.ids)
        elif self._thermo == "meos":
            eq = state["meos_eq"]
            self.cmp = mEoS.__all__[mEoS.id_mEoS.index(self.ids[0])](eq=eq)

        # Load advanced properties from global phase
        if self._thermo != "eos":
            self.cmp._readGlobalState(self, state)

        # Load state for phases
        if self._thermo in ["iapws", "freesteam"]:
            self.Liquido = ThermoWater()
            self.Gas = ThermoWater()
        elif self._thermo in ["coolprop", "meos"]:
            self.Liquido = ThermoAdvanced()
            self.Gas = ThermoAdvanced()
        elif self._thermo == "refprop":
            self.Liquido = ThermoRefProp()
            self.Gas = ThermoRefProp()
        else:
            self.Liquido = Mezcla()
            self.Gas = Mezcla()
        self.Liquido.readStatefromJSON(state["liquid"])
        self.Gas.readStatefromJSON(state["gas"])

        if state["liquid"]:
            self.Liquido.sigma = unidades.Tension(state["liquid"]["sigma"])


class PsyStream(config.Entity):
    """
    Class to model a stream as psychrometric state
    kwargs definition parameters:
        P: Pressure, Pa
        z: altitude, m

        tdp: dew-point temperature
        tdb: dry-bulb temperature
        twb: web-bulb temperature
        w: Humidity Ratio [kg water/kg dry air]
        HR: Relative humidity
        h: Mixture enthalpy
        v: Mixture specified volume
        state: opcional for predefined psychrometric state

    P: input for barometric pressure, z is an alternate pressure input

    For flow definition, one of:
        caudalMasico
        caudalVolumetrico
        caudalMolar
    """
    kwargs = {"z": 0.0,
              "P": 0.0,

              "tdb": 0.0,
              "tdb": 0.0,
              "twb": 0.0,
              "w": None,
              "HR": None,
              "h": None,
              "v": 0.0,

              "caudalMasico": 0.0,
              "caudalVolumetrico": 0.0,
              "caudalMolar": 0.0,
              "state": None}
    status = 0
    msg = "undefined"

    def __init__(self, **kwargs):
        self.kwargs = PsyStream.kwargs.copy()
        self.__call__(**kwargs)

    def __call__(self, **kwargs):
        if kwargs.get("state", None):
            kwargs.update(kwargs["state"].kwargs)

        self.kwargs.update(kwargs)

        for key, value in list(self.kwargs.items()):
            if value:
                self._bool = True
                break

        if self.calculable:
            self.status = 1
            self.calculo()
            self.msg = ""

    @property
    def calculable(self):
        # State definition
        tdp = self.kwargs.get("tdp", 0)
        tdb = self.kwargs.get("tdb", 0)
        twb = self.kwargs.get("twb", 0)
        w = self.kwargs.get("w", None)
        HR = self.kwargs.get("HR", None)
        h = self.kwargs.get("h", None)
        v = self.kwargs.get("v", 0)

        self.mode = -1
        if tdb and w is not None:
            self.mode = 0
        elif tdb and HR is not None:
            self.mode = 1
        elif tdb and twb:
            self.mode = 2
        elif tdb and tdp:
            self.mode = 3
        elif tdp and HR is not None:
            self.mode = 4
        elif self.kwargs["state"]:
            self.mode = 5

        # Flow definition
        caudal = self.kwargs["caudalVolumetrico"] \
            or self.kwargs["caudalMolar"] or self.kwargs["caudalMasico"]

        return bool(self.mode+1) and caudal

    def calculo(self):
        if self.kwargs["state"]:
            self.state = self.kwargs["state"]
        else:
            self.state = PsychroState(**self.kwargs)

        kwargs = self.kwargs
        self.__dict__.update(self.state.__dict__)
        self.kwargs = kwargs

        self.updateFlow()

    def updatekwargsFlow(self, key, value):
        self.kwargs[key] = value
        if key == "caudalMasico" and value:
            self.kwargs["caudalVolumetrico"] = 0.0
            self.kwargs["caudalMolar"] = 0.0
        elif key == "caudalMolar" and value:
            self.kwargs["caudalVolumetrico"] = 0.0
            self.kwargs["caudalMasico"] = 0.0
        elif key == "caudalVolumetrico" and value:
            self.kwargs["caudalMasico"] = 0.0
            self.kwargs["caudalMolar"] = 0.0

        if self.status:
            self.updateFlow()
        elif self.calculable:
            self.status = 1
            self.calculo()
            self.msg = ""

    def updateFlow(self):
        if self.kwargs["caudalMasico"]:
            G = self.kwargs["caudalMasico"]
            M = G/(self.Xw*18.01528+self.Xa*28.9645)
            Q = G*self.v
        elif self.kwargs["caudalMolar"]:
            M = self.kwargs["caudalMolar"]
            G = M*self.Xw*18.01528+M*self.Xa*28.9645
            Q = G*self.v
        elif self.kwargs["caudalVolumetrico"]:
            Q = self.kwargs["caudalVolumetrico"]
            G = Q*self.rho
            M = G/(self.Xw*18.01528+self.Xa*28.9645)
        self.caudalMasico = unidades.MassFlow(G)
        self.caudalVolumetrico = unidades.VolFlow(Q)
        self.caudalMolar = unidades.MolarFlow(M)

    @property
    def corriente(self):
        corriente = Corriente(T=self.state.twb, P=self.state.P,
                              caudalMasico=self.caudalMasico, ids=[62, 475],
                              fraccionMolar=[self.state.Xw, self.state.Xa])
        return corriente


if __name__ == '__main__':
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

#    corr=Corriente(T=300, P=101325.)
#    if corr:
#        print bool(corr)
#    aire=Corriente(T=350, P=101325, caudalMasico=0.01, ids=[475, 62], fraccionMolar=[0.99, 0.01])
#    agua=Corriente(T=300, P=101325, caudalMasico=0.1, ids=[62], fraccionMolar=[1.])


#    aire=PsyStream(caudal=5, tdb=300, HR=50)

#     agua=Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.], MEoS=True)
#     agua2=Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.], iapws=True)
#     from pprint import pprint
#     pprint(agua.__dict__)
#     pprint(agua2.__dict__)
#     print(agua.rho, agua2.rho)

    T = unidades.Temperature(120, "F")
    P = unidades.Pressure(485, "psi")
    stream = Corriente(T=T, P=P, caudalMolar=5.597, ids=[1, 2, 40, 41],
                       fraccionMolar=[0.3177, 0.5894, 0.0715, 0.0214],
                       K="Peng-Robinson", H="Benedict-Webb-Rubin-Starling")
    # from lib.mEoS import H2O
    # agua = H2O(T=300, P=101325)
    # agua1 = Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.], MEoS=True)
    # # agua2=Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.], iapws=True)
    # from pprint import pprint
    # pprint(agua.__dict__)
    # pprint(agua1.__dict__)

    # agua2=Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.], MEoS=True, refprop=True)
    # print(agua.rho, agua1.rho, agua2.rho)

    # agua3 = Corriente(T=300, P=101325, caudalMasico=1., ids=[62], fraccionMolar=[1.], MEoS=False, K="BWRS", H="BWRS")
    # print(agua.rho, agua1.rho, agua3.Gas.rho)

    # P = unidades.Pressure(410.3, "psi")
    # T = unidades.Temperature(400, "F")
    # mix = Corriente(T=T, P=P, caudalMasico=1., ids=[4, 40],
                    # caudalUnitarioMasico=[26.92, 73.08],  K="GS", H="BWRS")
    # print(mix.x)

    zi = [0.31767, 0.58942, 0.07147, 0.02144]
    P = unidades.Pressure(485, "psi")
    T = unidades.Temperature(100, "F")
    mix = Corriente(T=T, P=P, ids=[1, 2, 40, 41], caudalUnitarioMolar=zi,
                    K="GS", H="BWRS")
    print(mix.x)
    print(mix.Liquido.fraccion)
    print(mix.Gas.fraccion)

    zi = [0.31767, 0.58942, 0.07147, 0.02144, 0]
    mix = Corriente(T=T, P=P, ids=[1, 2, 40, 41, 34], caudalUnitarioMolar=zi,
                    K="GS", H="BWRS")
    print(mix.x)
    print(mix.Liquido.fraccion)
    print(mix.Gas.fraccion)

    # P = unidades.Pressure(485, "psi")
    # T = unidades.Temperature(100, "F")
    # mix = Corriente(T=T, P=P, caudalMasico=1., ids=[1, 2, 40, 41],
                    # caudalUnitarioMolar=[0.31767, 0.58942, 0.07147, 0.02144],
                    # K="PR76", H="PR76")

    # mix = Corriente(T=293.15, P=5e6, caudalMasico=1.,
                    # ids=[2, 3, 4, 6, 5, 8, 46, 49, 50, 22],
                    # caudalUnitarioMolar=[1]*10, K="PR76", H="PR76")
    # print(mix.x)
