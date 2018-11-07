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
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Library to implement solid-liquid interaction equipment:

    * :class:`Centrifuge`: Centrifuge separator equipment
    * :class:`Crystallizer`: Crystallizer separator equipment
    * :class:`Filter`: Filter separator equipment
'''


from math import log, exp

from lib.unidades import Currency
from lib.corriente import Corriente
from .parents import equipment


class Crystallizer(equipment):
    """Clase que define un cristalizador"""
    title = "Cristalizador"
    help = ""

    def coste(self, *args, **kwargs):
        """
        tipo:
            0   -   External forced circulation
            1   -   Internal draft tube
            2   -   Batch vacuum
        material:
            Para tipos 0 y 1:
                0   -   Acero inoxidable
                1   -   Acero dulce
            Para tipo 2:
                0   -   Acero dulce
                1   -   Acero inoxidable
                2   -   Cubierto de goma
        """
        self._indicesCoste(*args)

        self.tipo = kwargs.get("tipo", 0)
        self.material = kwargs.get("material", 0)
        V = self.V.ft3

        W = entrada.caudal_solido.lbh/1000.

        if self.tipo == 0:
            # External forced circulation
            if self.material == 0:
                Fm = 1
            else:
                Fm = 2.5
            C = Fm*exp(4.868+0.3092*log(W)+0.0548*log(W)**2)*1000
        elif self.tipo == 1:
            # Internal draft tube
            if self.material == 0:
                Fm = 1
            else:
                Fm = 2.5
            C = 178.*Fm*W**0.58*1000
        else:
            # Batch vacuum
            if self.material == 0:
                Fm = 1.
            elif self.material == 1:
                Fm = 2.
            else:
                Fm = 1.3
            C = 8.16*Fm*V**0.47*1000

        self.C_adq = Currency(C*self.Current_index/self.Base_index)
        self.C_inst = Currency(self.C_adq*self.f_install)


class Filter(equipment):
    """Clase que define un filtro, ya sea de vacío o de presión"""
    title = "Filtro"
    help = ""

    def coste(self, *args, **kwargs):
        """
        tipo:
            0   -   Rotary vacuum belt discharge
            1   -   Rotary vacuum drum scraper discharge
            2   -   Rotary vacuum disk
            3   -   Horizontal vacuum belt
            4   -   Pressure leaf
            5   -   Plate and frame
        """
        self._indicesCoste(*args)

        self.tipo = kwargs.get("tipo", 0)

        A = self.area.ft2

        if self.tipo == 0:
            # Rotary vacuum belt discharge
            C = exp(11.20 - 1.2252 * log(A) + 0.0587 * log(A)**2)*A
        elif self.tipo == 1:
            # Rotary vacuum drum scraper discharge
            C = exp(11.27 - 1.3408 * log(A) + 0.0709 * log(A)**2)*A
        elif self.tipo == 2:
            # Rotary vacuum disk
            C = exp(10.50 - 1.008 * log(A) + 0.0344 * log(A)**2)*A
        elif self.tipo == 3:
            # Horizontal vacuum belt
            C = 28300./A**0.5*A
        elif self.tipo == 4:
            # Pressure leaf
            C = 695./A**0.29*A
        else:
            # Plate and frame
            C = 460./A**0.45*A

        self.C_adq = Currency(C*self.Current_index/self.Base_index)
        self.C_inst = Currency(self.C_adq*self.f_install)


class Centrifuge(equipment):
    """ Clase que define un filtro centrifgo"""
    title = "Centrifuga"
    help = ""

    def coste(self, *args, **kwargs):
        """
        modelo:
            0   -   Proceso inorgánico
            1   -   Proceso orgánico
        material:
            0   -   Acero al carbón     Procesos organicos, Hastelloy, Ni-Fe-Mo
            1   -   Acero Inoxidable
            2   -   Monel (Ni-Cu)
            3   -   Niquel
        """
        self._indicesCoste(*args)

        self.modelo = kwargs.get("modelo", 0)
        self.material = kwargs.get("material", 0)

        if self.modelo == 0:
            # inorganic process
            if self.material == 0:
                # Carbon steel
                a, b = 42., 1.63
            elif self.material == 1:
                # Acero Inoxidable
                a, b = 65., 3.50
            elif self.material == 2:
                # Monel
                a, b = 70., 5.50
            else:
                # Nickel
                a, b = 84.4, 6.56
        else:
            # Inorganic process
            if self.material == 0:
                # Hastelloy
                a, b = 300., 10.
            elif self.material == 1:
                # Acero Inoxidable
                a, b = 98., 5.06
            elif self.material == 2:
                # Monel
                a, b = 114., 7.14
            else:
                # Nickel
                a, b = 143., 9.43

        W = self.entrada.caudal.Tmh
        C = (a+b*W)*1000
        self.C_adq = Currency(C*self.Current_index/self.Base_index)
        self.C_inst = Currency(self.C_adq*self.f_install)


if __name__ == '__main__':

    entrada = Corriente(423.15, 3, 9000,  [[46, 47, 49, 64],
                        [0.78, 0.21, 0.0001, 0.0009]], [64], [[1e-6, 1]])
    cristalizador = Crystallizer(entrada, 0.1)

#    agua=Corriente(300, 1, 4.3756, [[62], [1]])
#    scrubber=Scrubber(entrada, agua, 0.1)
