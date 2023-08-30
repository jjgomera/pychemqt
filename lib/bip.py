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


from numpy import zeros

from lib import config
from lib.sql import databank


EoSBIP = ["SRK", "PR", "APISRK", "BWRS", "NRTL", "UNIQUAC", "WILSON"]


def Kij(ids, EOS=None, dat=None):
    """Calculate binary interaction matrix for component of mixture,
    use bip data from database

    Parameters
    ----------
    ids : list
        Index of components in database, [-]
    EOS : string
        Code of equation of state: SRK, APISRK, PR, BWRS, NRTL, UNIQUAC, WILSON
    dat : dict
        For EoS with interaction parameter don't defined in database
    """
    # Return null bip if EOS is not specified
    if EOS is None or (EOS not in EoSBIP and dat is None):
        kij = zeros((len(ids), len(ids)))
        return kij

    # Continue with real procedure
    # Simple case with only a symetric parameter
    if EOS in ["SRK", "APISRK", "PR", "BWRS"]:
        kij = []
        query = "SELECT kij FROM %sbip WHERE i=? AND j=?" % EOS
        for i in ids:
            kiji = []
            for j in ids:
                i_ = min(i, j)
                j_ = max(i, j)
                databank.execute(query, (i_, j_))
                k = databank.fetchone()
                if k:
                    kiji.append(float(k[0]))
                else:
                    kiji.append(0)

            kij.append(kiji)
        return kij

    elif EOS == "LK":
        kij = []
        for i in ids:
            kiji = []
            for j in ids:
                key = f"{i}-{j}"
                key2 = f"{j}-{i}"
                if key in dat:
                    kiji.append(dat[key])
                elif key2 in dat:
                    kiji.append(dat[key2])
                else:
                    kiji.append(1)
            kij.append(kiji)
        return kij

    # Asymetric BIP
    else:
        kij = []
        query = "SELECT * FROM %sbip WHERE i=? AND j=?" % EOS
        for i in ids:
            kiji = []
            for j in ids:
                i_ = min(i, j)
                j_ = max(i, j)
                databank.execute(query, (i_, j_))
                k = databank.fetchone()
                if i == i_ and k:
                    kiji.append(float(k[3]))
                elif k:
                    kiji.append(float(k[4]))
                else:
                    kiji.append(0)
            kij.append(kiji)

        # Get second parameter for NRTL
        if EOS == "NRTL":
            alpha = []
            query = "SELECT alpha FROM NRTLbip WHERE i=? AND j=?"
            for i in ids:
                alphai = []
                for j in ids:
                    i_ = min(i, j)
                    j_ = max(i, j)
                    databank.execute(query, (i_, j_))
                    a = databank.fetchone()
                    if a:
                        alphai.append(float(a[0]))
                    else:
                        alphai.append(0)
                alpha.append(alphai)

            return kij, alpha
        else:
            return kij


# Mixing Rules
def Mix_vdW1f(xi, parameters, kij):
    """Mixing rules of van der Waals"""
    ai = parameters[0]
    bi = parameters[1:]

    # Geometric mean rule for croos-energy parameter
    a = 0
    for x_i, a_i, kiji in zip(xi, ai, kij):
        for x_j, a_j, k in zip(xi, ai, kiji):
            a += x_i * x_j * (a_i*a_j)**0.5 * (1-k)

    # Arithmetic mean rule for the aditional parameters
    b = []
    for b_i in bi:
        suma = 0
        for x_i, bii in zip(xi, b_i):
            suma += x_i*bii
        b.append(suma)

    return tuple([a]+b)


def Mix_Stryjek_Vera(self, parameters, kij):
    """Mixing rules of Stryjek and Vera (1986)"""
    ai = parameters[0]
    bi = parameters[1:]
    b = [0]*len(bi)
    a = 0
    for i in range(len(self.componente)):
        for j in range(len(bi)):
            b[j] += self.fraccion[i]*bi[j][i]
        for j in range(len(self.componente)):
            if kij[i][j] == 0 and kij[j][i] == 0:
                k = 0.
            else:
                k = kij[i][j]*kij[j][i]/(self.fraccion[i]*kij[i][j]+self.fraccion[j]*kij[j][i])
            a += self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-k)
    return tuple([a]+b)

def Mix_Panagiotopoulos(parameters, kij):
    """Mixing Rules of Panagiotopoulos (1985)"""
    ai = parameters[0]
    bi = parameters[1:]
    b = [0]*len(bi)
    a = 0
    for i in range(len(self.componente)):
        for j in range(len(bi)):
            b[j] += self.fraccion[i]*bi[j][i]
        for j in range(len(self.componente)):
            a += self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j]+(kij[i][j]-kij[j][i])*self.fraccion[i])
    return tuple([a]+b)

def Mix_Melhem(parameters, kij):
    """Mixing Rules of Melhem (1991)"""
    ai = parameters[0]
    bi = parameters[1:]
    b = [0]*len(bi)
    a = 0
    for i in range(len(self.componente)):
        for j in range(len(bi)):
            b[j] += self.fraccion[i]*bi[j][i]
        for j in range(len(self.componente)):
            a += self.fraccion[i]*self.fraccion[j]*(ai[i]*ai[j])**0.5*(1-kij[i][j]+(kij[i][j]-kij[j][i])*self.fraccion[i]/(self.fraccion[i]+self.fraccion[j]))
    return tuple([a]+b)


mixing = [Mix_vdW1f, Mix_Stryjek_Vera, Mix_Panagiotopoulos, Mix_Melhem]
conf = config.getMainWindowConfig().getint("Thermo", "Mixing")
Mixing_Rule = mixing[conf]
