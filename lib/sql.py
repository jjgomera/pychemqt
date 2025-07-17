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
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Module with actions to manipulate compound database

    * :func:`transformElement`: Transform values from dialog in valid list to \
save in databank
    * :func:`inserElementsFromArray`: Insert element to a database
    * :func:`updateElement`: Update element with indice in database
    * :func:`deleteElement`: Delete Element with indice from custom Database
    * :func:`getElement`: Get element from database
    * :func:`copyElement`: Create a copy of element of indice in custom Database

API reference
-------------

'''


import os
import sqlite3


databank_name = os.path.join(os.environ["pychemqt"], 'dat', 'databank.db')
databank = sqlite3.connect(databank_name).cursor()
databank.execute("SELECT COUNT(*) AS Total FROM compuestos")
N_comp = databank.fetchone()[0]

conf_dir = os.path.join(os.path.expanduser('~'), ".pychemqt")
databank_Custom_name = conf_dir + os.sep + 'databank.db'
if os.path.isfile(databank_Custom_name):
    databank_Custom = sqlite3.connect(databank_Custom_name).cursor()
    databank_Custom.execute("SELECT COUNT(*) AS Total FROM compuestos")
    N_comp_Custom = databank_Custom.fetchone()[0]
else:
    N_comp_Custom = 0


def transformElement(elemento):
    """Transform list generated in dialog in valid list to save in database"""
    vals = []
    vals.append(str(elemento[0]))   # formula
    vals.append(str(elemento[1]))   # name
    vals.append(elemento[2])        # M
    vals.append(elemento[3])        # tc
    vals.append(elemento[4])        # pc
    vals.append(elemento[5])        # vc
    vals.append(elemento[6])        # API

    if elemento[7]:                 # Cp_ideal
        vals += elemento[7]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[8]:                 # Antoine
        vals += elemento[8]
    else:
        vals += [0.0, 0.0, 0.0]

    if elemento[9]:                 # Henry
        vals += elemento[9]
    else:
        vals += [0.0, 0.0, 0.0, 0.0]

    if elemento[10]:                # Visco
        vals += elemento[10]
    else:
        vals += [0.0, 0.0]

    if elemento[11]:                # tension
        vals += elemento[11]
    else:
        vals += [0.0, 0.0]

    if elemento[12]:                # rhoS_DIPPR
        vals += elemento[12][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[13]:                # rhoL_DIPPR
        vals += elemento[13][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[14]:                # Pv_DIPPR
        vals += elemento[14][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[15]:                # Hv_DIPPR
        vals += elemento[15][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[16]:                # Cps_DIPPR
        vals += elemento[16][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[17]:                # CpL_DIPPR
        vals += elemento[17][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[18]:                # CpG_DIPPR
        vals += elemento[18][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[19]:                # muL_DIPPR
        vals += elemento[19][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[20]:                # muG_DIPPR
        vals += elemento[20][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[21]:                # ThCondL_DIPPR
        vals += elemento[21][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[22]:                # ThCondG_DIPPR
        vals += elemento[22][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    if elemento[23]:                # tension_DIPPR
        vals += elemento[23][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    vals.append(elemento[24])       # Dipole
    vals.append(elemento[25])       # V_liq
    vals.append(elemento[26])       # Rackett
    vals.append(elemento[27])       # SG
    vals.append(elemento[28])       # f_acent
    vals.append(elemento[29])       # SolubilityParameter
    vals.append(elemento[30])       # watson

    # if elemento[31]:                # MSRK
        # vals += elemento[31][0:2]
    # else:
        # vals += [0.0, 0.0]

    vals.append(elemento[31])       # Stiehl
    vals.append(elemento[32])       # Tb
    vals.append(elemento[33])       # Tf
    vals.append(str(elemento[34]))  # CAS
    vals.append(str(elemento[35]))  # alternateFormula
    vals.append(elemento[36])       # UNIFAC
    vals.append(elemento[37])       # Dm
    vals.append(elemento[38])       # Eps_k
    vals.append(elemento[39])       # UNIQUAC_area
    vals.append(elemento[40])       # UNIQUAC_volumen
    vals.append(elemento[41])       # f_acent
    vals.append(elemento[42])       # Hf
    vals.append(elemento[43])       # Gf
    vals.append(elemento[44])       # volumen_wilson
    vals.append(elemento[45])       # NetHeating
    vals.append(elemento[46])       # GrossHeating
    vals.append(str(elemento[47]))  # Synonyms
    vals.append(elemento[48])       # volumen_caracteristico
    vals.append(elemento[49])       # calor_formacion_solido
    vals.append(elemento[50])       # energia_libre_solido
    vals.append(elemento[51])       # PolarParameter
    vals.append(str(elemento[52]))  # smile

    if elemento[53]:                # antoine_extended
        vals += elemento[53][0:4]
    else:
        vals += [0.0, 0.0, 0.0, 0.0]

    if elemento[54]:                # wagner
        vals += elemento[54][0:4]
    else:
        vals += [0.0, 0.0, 0.0, 0.0]

    return vals


def inserElementsFromArray(name, lista):
    """Insert element to a database
    lista: array with component data in text format"""
    conn = sqlite3.connect(name)
    curs = conn.cursor()
    curs.execute("SELECT COUNT(*) AS Total FROM compuestos")
    numero = curs.fetchone()[0]
    if name == databank_Custom_name:
        numero += 10000
    query = "INSERT INTO compuestos VALUES "
    for indice, elemento in enumerate(lista):
        vals = transformElement(elemento)
        vals.insert(0, numero+indice+1)
        curs.execute(query+str(tuple(vals)))
    conn.commit()
    conn.close()


def updateElement(elemento, indice):
    """Update element with indice in custom Database"""
    variables = [
        "formula", "name", "M", "tc", "pc", "vc", "API",
        "Cp_ideal_A", "Cp_ideal_B", "Cp_ideal_C", "Cp_ideal_D", "Cp_ideal_E",
        "Cp_ideal_F", "antoine_A", "antoine_B", "antoine_C", "henry_A",
        "henry_B", "henry_C", "henry_D", "visco_A", "visco_B", "tension_A",
        "tension_B", "rhoS_DIPPR_EQ", "rhoS_DIPPR_A", "rhoS_DIPPR_B",
        "rhoS_DIPPR_C", "rhoS_DIPPR_D", "rhoS_DIPPR_E", "rhoS_DIPPR_tmin",
        "rhoS_DIPPR_tmax", "rhoL_DIPPR_EQ", "rhoL_DIPPR_A", "rhoL_DIPPR_B",
        "rhoL_DIPPR_C", "rhoL_DIPPR_D", "rhoL_DIPPR_E", "rhoL_DIPPR_tmin",
        "rhoL_DIPPR_tmax", "Pv_DIPPR_EQ", "Pv_DIPPR_A", "Pv_DIPPR_B",
        "Pv_DIPPR_C", "Pv_DIPPR_D", "Pv_DIPPR_E", "Pv_DIPPR_tmin",
        "Pv_DIPPR_tmax", "Hv_DIPPR_EQ", "Hv_DIPPR_A", "Hv_DIPPR_B",
        "Hv_DIPPR_C", "Hv_DIPPR_D", "Hv_DIPPR_E", "Hv_DIPPR_tmin",
        "Hv_DIPPR_tmax", "CpS_DIPPR_EQ", "CpS_DIPPR_A", "CpS_DIPPR_B",
        "CpS_DIPPR_C", "CpS_DIPPR_D", "CpS_DIPPR_E", "CpS_DIPPR_tmin",
        "CpS_DIPPR_tmax", "CpL_DIPPR_EQ", "CpL_DIPPR_A", "CpL_DIPPR_B",
        "CpL_DIPPR_C", "CpL_DIPPR_D", "CpL_DIPPR_E", "CpL_DIPPR_tmin",
        "CpL_DIPPR_tmax", "CpG_DIPPR_EQ", "CpG_DIPPR_A", "CpG_DIPPR_B",
        "CpG_DIPPR_C", "CpG_DIPPR_D", "CpG_DIPPR_E", "CpG_DIPPR_tmin",
        "CpG_DIPPR_tmax", "muL_DIPPR_EQ", "muL_DIPPR_A", "muL_DIPPR_B",
        "muL_DIPPR_C", "muL_DIPPR_D", "muL_DIPPR_E", "muL_DIPPR_tmin",
        "muL_DIPPR_tmax", "muG_DIPPR_EQ", "muG_DIPPR_A", "muG_DIPPR_B",
        "muG_DIPPR_C", "muG_DIPPR_D", "muG_DIPPR_E", "muG_DIPPR_tmin",
        "muG_DIPPR_tmax", "ThcondL_DIPPR_EQ", "ThcondL_DIPPR_A",
        "ThcondL_DIPPR_B", "ThcondL_DIPPR_C", "ThcondL_DIPPR_D",
        "ThcondL_DIPPR_E", "ThcondL_DIPPR_tmin", "ThcondL_DIPPR_tmax",
        "ThcondG_DIPPR_EQ", "ThcondG_DIPPR_A", "ThcondG_DIPPR_B",
        "ThcondG_DIPPR_C", "ThcondG_DIPPR_D", "ThcondG_DIPPR_E",
        "ThcondG_DIPPR_tmin", "ThcondG_DIPPR_tmax", "tension_DIPPR_EQ",
        "tension_DIPPR_A", "tension_DIPPR_B", "tension_DIPPR_C",
        "tension_DIPPR_D", "tension_DIPPR_E", "tension_DIPPR_tmin",
        "tension_DIPPR_tmax", "dipole", "V_liq",
        "Rackett", "SG", "f_acent",
        "SolubilityParameter", "watson", "Stiehl",
        "Tb", "Tf", "CAS", "alternateFormula", "UNIFAC",
        "Dm", "Eps_k", "UNIQUAC_area", "UNIQUAC_volumen",
        "f_acent_MSRK", "Hf",
        "Gf", "volumen_wilson", "NetHeating",
        "GrossHeating", "Synonyms",
        "volumen_caracteristico", "calor_formacion_solido",
        "energia_libre_solido", "PolarParameter", "smile", "antoine_to",
        "antoine_n", "antoine_E", "antoine_F", "wagner_a", "wagner_b",
        "wagner_c", "wagner_d"]
    vals = transformElement(elemento)

    conn = sqlite3.connect(databank_Custom_name)
    curs = conn.cursor()
    for variable, valor in zip(variables, vals):
        if isinstance(valor, int):
            curs.execute(
                f'UPDATE compuestos SET {variable}={valor} WHERE id=={indice}')
        elif isinstance(valor, float):
            curs.execute(
                f'UPDATE compuestos SET {variable}={valor} WHERE id=={indice}')
        elif isinstance(valor, str):
            valor = '"'+valor+'"'
            curs.execute(
                f'UPDATE compuestos SET {variable}={valor} WHERE id=={indice}')
    conn.commit()
    conn.close()


def deleteElement(indice):
    """Delete Element with indice from custom Database"""
    conn = sqlite3.connect(databank_Custom_name)
    curs = conn.cursor()
    curs.execute(f"DELETE FROM compuestos WHERE id={indice}")
    conn.commit()
    conn.close()


def getElement(indice):
    """Get element from database
    indice: index in databank of element"""
    if indice > 10000:
        db_Custom = sqlite3.connect(databank_Custom_name).cursor()
        db_Custom.execute(f"select * from compuestos where id=={indice}")
        componente = db_Custom.fetchone()
    else:
        db = sqlite3.connect(databank_name).cursor()
        db.execute(f"select * from compuestos where id=={indice}")
        componente = db.fetchone()

    # Change none values for 0s
    componente = list(componente)
    while None in componente:
        componente[componente.index(None)] = 0

    return tuple(componente)


def copyElement(indice):
    """Create a copy of element of indice in custom Database"""
    elemento = getElement(indice)
    vals = elemento[1:]
    conn = sqlite3.connect(databank_Custom_name)
    curs = conn.cursor()
    curs.execute("INSERT INTO compuestos VALUES"
                 + str((10001+N_comp_Custom, ) + vals))
    conn.commit()
    conn.close()
