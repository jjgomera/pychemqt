#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2017, Juan José Gómez Romera <jjgomera@gmail.com>

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


from scipy import exp
from scipy.optimize import fsolve

from lib import unidades


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

    def _fug(self, xi, yi):
        """Each child class with liquid-vapor equilibrium support must define
        this procedure to do the flash iteration, it return fugacity
        coefficient for both phases"""
        # If this method if executed the child class don't define it
        print("ERROR: Liquid-Vapor fugacities unimplemented")
        return

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

        def RR(q):
            """Rachford-Rice equation"""
            f = 0
            for zi, ki in zip(self.zi, Ki):
                f += zi*(1-ki)/(1+q*(ki-1))
            return f

        if RR(0) > 0 and RR(1) > 0:
            # q>1, superheated gas
            xi = self.zi
            yi = self.zi
            q = 1
            Zv = self._Z(self.zi)[-1]
            Zl = None
        elif RR(0) < 0 and RR(1) < 0:
            # q<0, subcooled liquid
            xi = self.zi
            yi = self.zi
            q = 0
            Zl = self._Z(self.zi)[0]
            Zv = None
        else:
            q = 0.5
            while True:
                qo = q
                solucion = fsolve(RR, q, full_output=True)
                if solucion[2] != 1:
                    print(solucion)
                    break
                else:
                    q = solucion[0][0]
                    xi = []
                    yi = []
                    for zi, ki in zip(self.zi, Ki):
                        xi.append(zi/(1+q*(ki-1)))
                        yi.append(zi*ki/(1+q*(ki-1)))

                    tital, titav = self._fug(xi, yi)

                    # print("tital", tital)
                    # print("titav", titav)

                    # print("tital2", self._fugl2(Zl, xi, Al, Bl))
                    # print("titav2", self._fugl2(Zv, yi, Av, Bv))

                    # tital2 = self._fugl2(Zl, xi)
                    # titav2 = self._fugl2(Zg, yi)

                    fiv = [z*t*self.P for z, t in zip(yi, titav)]
                    fil = [z*t*self.P for z, t in zip(xi, tital)]

                    # criterio de convergencia Eq 21
                    err = sum([abs(l/v-1) for l, v in zip(fil, fiv)])
                    if err < 1e-12 and (q-qo)**2 < 1e-15:
                        break
                    else:
                        Ki = [l/v for l, v in zip(tital, titav)]

            Zl = self._Z(xi)[0]
            Zv = self._Z(yi)[-1]

        return q, Zl, Zv, xi, yi, Ki

    def _Bubble_T(self):
        def f(T):
            eq=self.__class__(T, self.P.atm, self.mezcla)
            return sum([k*x for k, x in zip(eq.Ki, self.zi)])-1.

        T=fsolve(f, self.T)
        return unidades.Temperature(T)

    def _Bubble_P(self):
        def f(P):
            eq=self.__class__(self.T, P, self.mezcla)
            return sum([k*x for k, x in zip(eq.Ki, self.zi)])-1.

        P=fsolve(f, self.P.atm)
        return unidades.Pressure(P, "atm")

    def _Dew_T(self):
        def f(T):
            eq=self.__class__(T, self.P.atm, self.mezcla)
            return 1./sum([x/k for k, x in zip(eq.Ki, self.zi)])-1.

        T=fsolve(f, self.T)
        return unidades.Temperature(T)

    def _Dew_P(self):
        def f(P):
            eq=self.__class__(self.T, P, self.mezcla)
            return sum([x/k for k, x in zip(eq.Ki, self.zi)])-1.

        P=fsolve(f, self.P.atm)
        return unidades.Pressure(P, "atm")

    # @property
    # def H_exc(self):
