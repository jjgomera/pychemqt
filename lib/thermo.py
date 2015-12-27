#!/usr/bin/python
# -*- coding: utf-8 -*-

#################################################################################
# Module with common thermal utilities:
#   - Fluid: Common class for thermodynamic method
#      Implement properites save/load to stream. Null properties state
#################################################################################


from . import unidades


class Fluid(object):
    """Custom object with null parameter to model a fluid with properties"""
    _bool = False
    h = 0
    s = 0
    cp = 0
    cv = 0
    cp_cv = 0
    cp0 = 0

    def writeStatetoStream(self, stream):
        stream.writeBool(self._bool)
        if self._bool:
            stream.writeFloat(self.M)
            stream.writeFloat(self.v)

            stream.writeFloat(self.h)
            stream.writeFloat(self.s)
            stream.writeFloat(self.u)
            stream.writeFloat(self.a)
            stream.writeFloat(self.g)

            stream.writeFloat(self.cv)
            stream.writeFloat(self.cp)
            stream.writeFloat(self.cp_cv)
            stream.writeFloat(self.w)

            stream.writeFloat(self.Z)
            stream.writeFloat(self.alfav)
            stream.writeFloat(self.xkappa)

            stream.writeFloat(self.mu)
            stream.writeFloat(self.k)
            stream.writeFloat(self.nu)
            stream.writeFloat(self.epsilon)
            stream.writeFloat(self.Prandt)
            stream.writeFloat(self.n)

            stream.writeFloat(self.alfa)
            stream.writeFloat(self.joule)
            stream.writeFloat(self.deltat)
            stream.writeFloat(self.gamma)

            stream.writeFloat(self.alfap)
            stream.writeFloat(self.betap)

            stream.writeFloat(self.v0)
            stream.writeFloat(self.h0)
            stream.writeFloat(self.u0)
            stream.writeFloat(self.s0)
            stream.writeFloat(self.a0)
            stream.writeFloat(self.g0)

            stream.writeFloat(self.cp0)
            stream.writeFloat(self.cv0)
            stream.writeFloat(self.cp0_cv)
            stream.writeFloat(self.w0)
            stream.writeFloat(self.gamma0)
            stream.writeFloat(self.f)

            stream.writeFloat(self.Q)
            stream.writeFloat(self.caudalmasico)
            stream.writeFloat(self.caudalmolar)
            stream.writeInt32(len(self.fraccion))
            for x in self.fraccion:
                stream.writeFloat(x)
            for x in self.fraccion_masica:
                stream.writeFloat(x)
            for x in self.caudalunitariomasico:
                stream.writeFloat(x)
            for x in self.caudalunitariomolar:
                stream.writeFloat(x)

    def readStatefromStream(self, stream):
        self._bool = stream.readBool()
        if self._bool:
            self.M = unidades.Dimensionless(stream.readFloat())
            self.v = unidades.SpecificVolume(stream.readFloat())
            self.rho = unidades.Density(1/self.v)

            self.h = unidades.Enthalpy(stream.readFloat())
            self.s = unidades.SpecificHeat(stream.readFloat())
            self.u = unidades.Enthalpy(stream.readFloat())
            self.a = unidades.Enthalpy(stream.readFloat())
            self.g = unidades.Enthalpy(stream.readFloat())

            self.cv = unidades.SpecificHeat(stream.readFloat())
            self.cp = unidades.SpecificHeat(stream.readFloat())
            self.cp_cv = unidades.Dimensionless(stream.readFloat())
            self.w = unidades.Speed(stream.readFloat())

            self.Z = unidades.Dimensionless(stream.readFloat())
            self.alfav = unidades.InvTemperature(stream.readFloat())
            self.xkappa = unidades.InvPressure(stream.readFloat())

            self.mu = unidades.Viscosity(stream.readFloat())
            self.k = unidades.ThermalConductivity(stream.readFloat())
            self.nu = unidades.Diffusivity(stream.readFloat())
            self.epsilon = unidades.Dimensionless(stream.readFloat())
            self.Prandt = unidades.Dimensionless(stream.readFloat())
            self.n = unidades.Dimensionless(stream.readFloat())

            self.alfa = unidades.Diffusivity(stream.readFloat())
            self.joule = unidades.TemperaturePressure(stream.readFloat())
            self.deltat = unidades.EnthalpyPressure(stream.readFloat())
            self.gamma = unidades.Dimensionless(stream.readFloat())

            self.alfap=stream.readFloat()
            self.betap=stream.readFloat()

            self.v0 = unidades.SpecificVolume(stream.readFloat())
            self.h0 = unidades.Enthalpy(stream.readFloat())
            self.u0 = unidades.Enthalpy(stream.readFloat())
            self.s0 = unidades.SpecificHeat(stream.readFloat())
            self.a0 = unidades.Enthalpy(stream.readFloat())
            self.g0 = unidades.Enthalpy(stream.readFloat())

            self.cp0 = unidades.SpecificHeat(stream.readFloat())
            self.cv0 = unidades.SpecificHeat(stream.readFloat())
            self.cp0_cv = unidades.Dimensionless(stream.readFloat())
            self.w0 = unidades.Speed(stream.readFloat())
            self.gamma0 = unidades.Dimensionless(stream.readFloat())
            self.f = unidades.Pressure(stream.readFloat())

            self.Q = unidades.VolFlow(stream.readFloat())
            self.caudalmasico = unidades.MassFlow(stream.readFloat())
            self.caudalmolar = unidades.MolarFlow(stream.readFloat())
            num = stream.readInt32()
            self.fraccion = []
            for i in range(num):
                self.fraccion.append(stream.readFloat())
            self.fraccion_masica = []
            for i in range(num):
                self.fraccion_masica.append(stream.readFloat())
            self.caudalunitariomasico = []
            for i in range(num):
                self.caudalunitariomasico.append(stream.readFloat())
            self.caudalunitariomolar = []
            for i in range(num):
                self.caudalunitariomolar.append(stream.readFloat())
