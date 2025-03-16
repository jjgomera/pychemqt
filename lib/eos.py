#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


from numpy import exp
from numpy.lib.scimath import log
from scipy.optimize import fsolve

from lib import unidades
from lib.utilities import refDoc


__doi__ = {
    1:
        {"autor": "Rachford, H.H., Rice, J.D.",
         "title": "Procedure for Use of Electronic Digital Computers in "
                  "Calculating Flash Vaporization Hydrocarbon Equilibrium",
         "ref": "Petroleum Transactions, AIME 195 (1952) 327-328",
         "doi": "10.2118/952327-G"},
    2:
        {"autor": "Peng, D.-Y.",
         "title": "Accelerated Successive Substitution Schemes for "
                  "Bubble-Point and Dew-Point Calculations",
         "ref": "Can. J. Chem. Eng. 69(4) (1991) 978-985",
         "doi": "10.1002/cjce.5450690421"},
    3:
        {"autor": "",
         "title": "",
         "ref": "",
         "doi": ""},
        }


# TODO: Add UI configuration support for specific parameters
# Grayson-Streed: Flory option
# SRK: alpha option
# PR: alpha option
# PRMathiasCopeman: alpha option
# BWRS: extended option
# virial: B, C, method configuration


@refDoc(__doi__, [1])
class EoS(object):
    """Base class for equation of state modeling, define common functionality
    as LV flash algorithm"""

    def __init__(self, T, P, mezcla, **kwargs):
        self.T = unidades.Temperature(T)
        self.P = unidades.Pressure(P)
        self.mezcla = mezcla
        self.componente = mezcla.componente
        self.zi = mezcla.fraccion
        self.kwargs = kwargs

    def _fug(self, xi, yi, T, P):
        """Each child class with liquid-vapor equilibrium support must define
        this procedure to do the flash iteration, it return fugacity
        coefficient for both phases"""
        # If this method if executed the child class don't define it
        # So raise NotImplementedError to let user to know that
        msg = "Derived class %s don't define Liquid-Vapor fugacities" % (
            self.__class__.__name__)
        raise NotImplementedError(msg)

    def _Flash(self):
        """Calculation K values for liquid-vapour phase equilibrium
        Naji, H.S.
        Conventional and Rapid Flash Calculations for the Soave-Redlich-Kwong
        and Peng-Robinson Equations of State
        Emirates J. Eng. Res. 13(3) (2008) 81-91
        """
        # Initial estimation using Wilson correlation, Eq 19
        Ki = []
        for c in self.componente:
            Ki.append(c.Pc/self.P*exp(5.37*(1.+c.f_acent)*(1.-c.Tc/self.T)))

        # Rachford-Rice equation, [1]_
        def RR(q):
            f = 0
            for zi, ki in zip(self.zi, Ki):
                f += zi*(1-ki)/(1+q*(ki-1))
            return f

        if RR(0) > 0 and RR(1) > 0:
            # q>1, superheated gas
            xi = self.zi
            yi = self.zi
            q = 1
            Zv = self._Z(self.zi, self.T, self.P)[-1]
            Zl = None
        elif RR(0) < 0 and RR(1) < 0:
            # q<0, subcooled liquid
            xi = self.zi
            yi = self.zi
            q = 0
            Zl = self._Z(self.zi, self.T, self.P)[0]
            Zv = None
        else:
            q = 0.5
            while True:
                qo = q
                solucion = fsolve(RR, q, full_output=True)
                if solucion[2] != 1:
                    break
                else:
                    q = solucion[0][0]
                    xi = []
                    yi = []
                    for zi, ki in zip(self.zi, Ki):
                        xi.append(zi/(1+q*(ki-1)))
                        yi.append(zi*ki/(1+q*(ki-1)))

                    tital, titav = self._fug(xi, yi, self.T, self.P)

                    fiv = [z*t*self.P for z, t in zip(yi, titav)]
                    fil = [z*t*self.P for z, t in zip(xi, tital)]

                    # criterio de convergencia Eq 21
                    err = sum([abs(l/v-1) for l, v in zip(fil, fiv)])
                    if err < 1e-12 and (q-qo)**2 < 1e-15:
                        break
                    else:
                        Ki = [l/v for l, v in zip(tital, titav)]

            Zl = self._Z(xi, self.T, self.P)[0]
            Zv = self._Z(yi, self.T, self.P)[-1]

        return q, Zl, Zv, xi, yi, Ki

    def _Bubble_T(self, P):
        """Calculation Bubble Point Temperature"""
        # Initial estimation of temperature
        T = 0
        for xi, cmp in zip(self.zi, self.componente):
            Ti = cmp.Tc/(1-3*log(P/cmp.Pc)/(log(10)*(7+7*cmp.f_acent)))
            T += Ti*xi

        # Initial estimation of K using inverted Wilson correlation
        Ki = []
        for c in self.componente:
            Ki.append(c.Pc/P*exp(5.37*(1.+c.f_acent)*(1.-c.Tc/T)))

        yi = [k*x for k, x in zip(Ki, self.zi)]
        ym = sum(yi)
        yi = [y/ym for y in yi]

        def f():
            val = 0
            for zi, ki in zip(self.zi, Ki):
                val += zi*ki
            return val - 1

        c = 0
        while True:
            c += 1
            tital, titav = self._fug(self.zi, yi, T, P)
            Ki = [l/v for l, v in zip(tital, titav)]
            if abs(f()) <= 1e-9:
                find = 1
                break
            else:
                yi = [k*x for k, x in zip(Ki, self.zi)]
                ym = sum(yi)
                yi = [y/ym for y in yi]

                # Estimating next temperature value from [2]_, eq 11
                s = 0
                for cmp, x, k in zip(self.componente, self.zi, Ki):
                    s += (1+cmp.f_acent)*cmp.Tc*x*k
                ed = T/5.373/s
                T *= (1-ed*f())

            if c > 500:
                print("reach limit iteration count")
                find = 0
                break

        return T, find

    def _Dew_T(self, P):
        """Calculation Dew Point Temperature"""
        # Initial estimation of temperature
        T = 0
        for xi, cmp in zip(self.zi, self.componente):
            Ti = cmp.Tc/(1-3*log(P/cmp.Pc)/(log(10)*(7+7*cmp.f_acent)))
            T += Ti*xi

        # Initial estimation of K using inverted Wilson correlation
        Ki = []
        for c in self.componente:
            Ki.append(c.Pc/P*exp(5.37*(1.+c.f_acent)*(1.-c.Tc/T)))

        xi = [k/y for k, y in zip(Ki, self.zi)]
        xm = sum(xi)
        xi = [x/xm for x in xi]

        def f():
            val = 1
            for zi, ki in zip(self.zi, Ki):
                val -= zi/ki
            return val

        c = 0
        while True:
            c += 1
            tital, titav = self._fug(xi, self.zi, T, P)
            # print(xi, tital, titav)
            Ki = [l/v for l, v in zip(tital, titav)]
            if abs(f()) <= 1e-9:
                find = 1
                break
            else:
                xi = [k/y for k, y in zip(Ki, self.zi)]
                xm = sum(xi)
                xi = [x/xm for x in xi]

                # Estimating next temperature value from [2]_, eq 13
                s = 0
                for cmp, x, k in zip(self.componente, self.zi, Ki):
                    s += (1+cmp.f_acent)*cmp.Tc*x/k
                ed = T/5.373/s
                T *= (1-ed*f())

            if c > 500:
                print("reach limit iteration count")
                find = 0
                break

        return T, find

    def _Bubble_P(self, T):
        """Calculation Bubble Point Pressure"""
        # Initial estimation of pressure
        P = 0
        for xi, cmp in zip(self.zi, self.componente):
            Pi = cmp.Pc*exp(5.373*(1+cmp.f_acent)*(1.-cmp.Tc/T))
            if T <= cmp.Tc:
                P += Pi*xi
            else:
                P += (Pi*cmp.Pc)**0.5*xi

        # Initial estimation of K using inverted Wilson correlation
        Ki = []
        for c in self.componente:
            Ki.append(c.Pc/P*exp(5.37*(1.+c.f_acent)*(1.-c.Tc/T)))

        yi = [k*x for k, x in zip(Ki, self.zi)]

        def f():
            val = 0
            for zi, ki in zip(self.zi, Ki):
                val += zi*ki
            return val - 1

        c = 0
        while True:
            c += 1
            tital, titav = self._fug(self.zi, yi, T, P)
            Ki = [l/v for l, v in zip(tital, titav)]
            if abs(f()) <= 1e-9:
                find = 1
                break
            else:
                yi = [k*x for k, x in zip(Ki, self.zi)]
                ym = sum(yi)
                yi = [y/ym for y in yi]

                # Estimating next temperature value from [2]_, eq 8
                s = sum([zi*ki for zi, ki in zip(self.zi, Ki)])
                P *= 2-1*s

            if c > 100:
                print("reach limit iteration count")
                find = 0
                break

        return P, find

    def _Dew_P(self, T):
        """Calculation Dew Point Pressure"""
        # Initial estimation of pressure
        P = 0
        for xi, cmp in zip(self.zi, self.componente):
            Pi = cmp.Pc*exp(5.373*(1+cmp.f_acent)*(1-cmp.Tc/T))
            if T <= cmp.Tc:
                P += xi/Pi
            else:
                P += xi/(Pi*cmp.Pc)**0.5
        P = 1/P

        # Initial estimation of K using inverted Wilson correlation
        Ki = []
        for c in self.componente:
            Ki.append(c.Pc/P*exp(5.373*(1+c.f_acent)*(1-c.Tc/T)))

        xi = [k/y for k, y in zip(Ki, self.zi)]
        xm = sum(xi)
        xi = [x/xm for x in xi]

        def f():
            val = 1
            for zi, ki in zip(self.zi, Ki):
                val -= zi/ki
            return val

        c = 0
        while True:
            c += 1
            tital, titav = self._fug(xi, self.zi, T, P)
            # print(xi, tital, titav)
            Ki = [l/v for l, v in zip(tital, titav)]
            if abs(f()) <= 1e-9:
                find = 1
                break
            else:
                xi = [k/y for k, y in zip(Ki, self.zi)]
                xm = sum(xi)
                xi = [x/xm for x in xi]

                # Estimating next temperature value from [2]_, eq 9
                s = sum([zi/ki for zi, ki in zip(self.zi, Ki)])
                P /= s

            if c > 100:
                print("reach limit iteration count")
                find = 0
                break

        return P, find

    def envelope(self):
        # TODO: Implement algorithm
        from numpy import logspace
        P = logspace(6, 7, 50)
        Pd = []
        dew = []
        Pb = []
        bubble = []
        for p in P:
            try:
                print(p)
                Td, success = self._Dew_T(p)
                if success:
                    Pd.append(p)
                    dew.append(Td)
                Tb, success = self._Bubble_T(p)
                if success:
                    Pb.append(p)
                    bubble.append(Tb)
            except:
                pass

        from pylab import plot, show
        plot(dew, Pd, "o")
        plot(bubble, Pb, "x")
        # yscale("log")
        show()
