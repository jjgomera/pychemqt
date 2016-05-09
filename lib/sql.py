#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
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
along with this program.  If not, see <http://www.gnu.org/licenses/>.'''


###############################################################################
# Module for properties database function
#   -createDatabase: Create empty database
#   -transformElement
#   -inserElementsFromArray: Insert element to a database
#   -updateElement: Update element with indice in database
#   -deleteElement: Delete Element with indice from custom Database
#   -getElement: Get element from database
#   -copyElement: Create a copy of element of indice in custom Database
###############################################################################


import os
import sqlite3


databank_name = os.environ["pychemqt"] + 'dat'+os.sep+'databank.db'
databank = sqlite3.connect(databank_name).cursor()
databank.execute("SELECT COUNT(*) AS Total FROM compuestos")
N_comp = databank.fetchone()[0]

conf_dir = os.path.expanduser('~') + os.sep+".pychemqt"+os.sep
databank_Custom_name = conf_dir + 'databank.db'
databank_Custom = sqlite3.connect(databank_Custom_name).cursor()
databank_Custom.execute("SELECT COUNT(*) AS Total FROM compuestos")
N_comp_Custom = databank_Custom.fetchone()[0]


def createDatabase(name):
    """Create empty database"""
    conn = sqlite3.connect(name)
    curs = conn.cursor()
    curs.execute("""
                 CREATE TABLE compuestos (
                 id  INTEGER PRIMARY KEY,
                 formula TEXT,
                 nombre  TEXT,
                 peso_molecular FLOAT,
                 tc          FLOAT,
                 pc          FLOAT,
                 vc          FLOAT,
                 API         FLOAT,
                 Cp_ideal_A    FLOAT,
                 Cp_ideal_B    FLOAT,
                 Cp_ideal_C    FLOAT,
                 Cp_ideal_D    FLOAT,
                 Cp_ideal_E    FLOAT,
                 Cp_ideal_F    FLOAT,
                 antoine_A   FLOAT,
                 antoine_B   FLOAT,
                 antoine_C   FLOAT,
                 henry_A     FLOAT,
                 henry_B     FLOAT,
                 henry_C     FLOAT,
                 henry_D     FLOAT,
                 visco_A     FLOAT,
                 visco_B     FLOAT,
                 tension_A       FLOAT,
                 tension_B       FLOAT,
                 rhoS_DIPPR_EQ   INTEGER,
                 rhoS_DIPPR_A   FLOAT,
                 rhoS_DIPPR_B   FLOAT,
                 rhoS_DIPPR_C   FLOAT,
                 rhoS_DIPPR_D   FLOAT,
                 rhoS_DIPPR_E  FLOAT,
                 rhoS_DIPPR_tmin   FLOAT,
                 rhoS_DIPPR_tmax   FLOAT,
                 rhoL_DIPPR_EQ   INTEGER,
                 rhoL_DIPPR_A   FLOAT,
                 rhoL_DIPPR_B   FLOAT,
                 rhoL_DIPPR_C   FLOAT,
                 rhoL_DIPPR_D   FLOAT,
                 rhoL_DIPPR_E  FLOAT,
                 rhoL_DIPPR_tmin   FLOAT,
                 rhoL_DIPPR_tmax   FLOAT,
                 Pv_DIPPR_EQ   INTEGER,
                 Pv_DIPPR_A   FLOAT,
                 Pv_DIPPR_B   FLOAT,
                 Pv_DIPPR_C   FLOAT,
                 Pv_DIPPR_D   FLOAT,
                 Pv_DIPPR_E  FLOAT,
                 Pv_DIPPR_tmin   FLOAT,
                 Pv_DIPPR_tmax   FLOAT,
                 Hv_DIPPR_EQ   INTEGER,
                 Hv_DIPPR_A   FLOAT,
                 Hv_DIPPR_B   FLOAT,
                 Hv_DIPPR_C   FLOAT,
                 Hv_DIPPR_D   FLOAT,
                 Hv_DIPPR_E  FLOAT,
                 Hv_DIPPR_tmin   FLOAT,
                 Hv_DIPPR_tmax   FLOAT,
                 CpS_DIPPR_EQ   INTEGER,
                 CpS_DIPPR_A   FLOAT,
                 CpS_DIPPR_B   FLOAT,
                 CpS_DIPPR_C   FLOAT,
                 CpS_DIPPR_D   FLOAT,
                 CpS_DIPPR_E  FLOAT,
                 CpS_DIPPR_tmin   FLOAT,
                 CpS_DIPPR_tmax   FLOAT,
                 CpL_DIPPR_EQ   INTEGER,
                 CpL_DIPPR_A   FLOAT,
                 CpL_DIPPR_B   FLOAT,
                 CpL_DIPPR_C   FLOAT,
                 CpL_DIPPR_D   FLOAT,
                 CpL_DIPPR_E  FLOAT,
                 CpL_DIPPR_tmin   FLOAT,
                 CpL_DIPPR_tmax   FLOAT,
                 CpG_DIPPR_EQ   INTEGER,
                 CpG_DIPPR_A   FLOAT,
                 CpG_DIPPR_B   FLOAT,
                 CpG_DIPPR_C   FLOAT,
                 CpG_DIPPR_D   FLOAT,
                 CpG_DIPPR_E  FLOAT,
                 CpG_DIPPR_tmin   FLOAT,
                 CpG_DIPPR_tmax   FLOAT,
                 muL_DIPPR_EQ   INTEGER,
                 muL_DIPPR_A   FLOAT,
                 muL_DIPPR_B   FLOAT,
                 muL_DIPPR_C   FLOAT,
                 muL_DIPPR_D   FLOAT,
                 muL_DIPPR_E  FLOAT,
                 muL_DIPPR_tmin   FLOAT,
                 muL_DIPPR_tmax   FLOAT,
                 muG_DIPPR_EQ   INTEGER,
                 muG_DIPPR_A   FLOAT,
                 muG_DIPPR_B   FLOAT,
                 muG_DIPPR_C   FLOAT,
                 muG_DIPPR_D   FLOAT,
                 muG_DIPPR_E  FLOAT,
                 muG_DIPPR_tmin   FLOAT,
                 muG_DIPPR_tmax   FLOAT,
                 ThcondL_DIPPR_EQ   INTEGER,
                 ThcondL_DIPPR_A   FLOAT,
                 ThcondL_DIPPR_B   FLOAT,
                 ThcondL_DIPPR_C   FLOAT,
                 ThcondL_DIPPR_D   FLOAT,
                 ThcondL_DIPPR_E  FLOAT,
                 ThcondL_DIPPR_tmin   FLOAT,
                 ThcondL_DIPPR_tmax   FLOAT,
                 ThcondG_DIPPR_EQ   INTEGER,
                 ThcondG_DIPPR_A   FLOAT,
                 ThcondG_DIPPR_B   FLOAT,
                 ThcondG_DIPPR_C   FLOAT,
                 ThcondG_DIPPR_D   FLOAT,
                 ThcondG_DIPPR_E  FLOAT,
                 ThcondG_DIPPR_tmin   FLOAT,
                 ThcondG_DIPPR_tmax   FLOAT,
                 tension_DIPPR_EQ   INTEGER,
                 tension_DIPPR_A   FLOAT,
                 tension_DIPPR_B   FLOAT,
                 tension_DIPPR_C   FLOAT,
                 tension_DIPPR_D   FLOAT,
                 tension_DIPPR_E  FLOAT,
                 tension_DIPPR_tmin   FLOAT,
                 tension_DIPPR_tmax   FLOAT,
                 momento_dipolar     FLOAT,
                 constante_volumen_liquido   FLOAT,
                 constante_rackett   FLOAT,
                 densidad_especifica     FLOAT,
                 factor_acentrico    FLOAT,
                 parametro_solubilidad   FLOAT,
                 watson      FLOAT,
                 MSRK_A    FLOAT,
                 MSRK_B    FLOAT,
                 Stiehl  FLOAT,
                 t_ebullicion  FLOAT,
                 t_fusion  FLOAT,
                 CAS_id  TEXT,
                 formula_alternativa TEXT,
                 UNIFAC  TEXT,
                 diametro_molecular  FLOAT,
                 Eps_k   FLOAT,
                 UNIQUAC_area    FLOAT,
                 UNIQUAC_volumen FLOAT,
                 factor_acentrico_modificado FLOAT,
                 calor_formacion_gas FLOAT,
                 energia_libre_gas   FLOAT,
                 volumen_wilson  FLOAT,
                 calor_combustion_neto FLOAT,
                 calor_combustion_bruto   FLOAT,
                 nombre_alternativo TEXT,
                 volumen_caracteristico FLOAT,
                 calor_formacion_solido  FLOAT,
                 energia_libre_solido    FLOAT,
                 parametro_polar FLOAT,
                 smile   TEXT)
                 """)
    conn.commit()
    conn.close()


def transformElement(elemento):
    vals = []
    vals.append(str(elemento[0]))
    vals.append(str(elemento[1]))
    vals.append(elemento[2])
    vals.append(elemento[3])
    vals.append(elemento[4])
    vals.append(elemento[5])
    vals.append(elemento[6])

    if elemento[7]:
        vals += elemento[7]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[8]:
        vals += elemento[8]
    else:
        vals += [0.0, 0.0, 0.0]

    if elemento[9]:
        vals += elemento[9]
    else:
        vals += [0.0, 0.0, 0.0, 0.0]
    if elemento[10]:
        vals += elemento[10]
    else:
        vals += [0.0, 0.0]
    if elemento[11]:
        vals += elemento[11]
    else:
        vals += [0.0, 0.0]
    if elemento[12]:
        vals += elemento[12][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[13]:
        vals += elemento[13][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[14]:
        vals += elemento[14][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[15]:
        vals += elemento[15][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[16]:
        vals += elemento[16][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[17]:
        vals += elemento[17][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[18]:
        vals += elemento[18][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[19]:
        vals += elemento[19][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[20]:
        vals += elemento[20][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[21]:
        vals += elemento[21][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[22]:
        vals += elemento[22][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if elemento[23]:
        vals += elemento[23][0:8]
    else:
        vals += [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    vals.append(elemento[24])
    vals.append(elemento[25])
    vals.append(elemento[26])
    vals.append(elemento[27])
    vals.append(elemento[28])
    vals.append(elemento[29])
    vals.append(elemento[30])
    if elemento[31]:
        vals += elemento[31]
    else:
        vals += [0.0, 0.0]
    vals.append(elemento[32])
    vals.append(elemento[33])
    vals.append(elemento[34])
    vals.append(str(elemento[35]))
    vals.append(str(elemento[36]))
    vals.append(str(elemento[37]))
    vals.append(elemento[38])
    vals.append(elemento[39])
    vals.append(elemento[40])
    vals.append(elemento[41])
    vals.append(elemento[42])
    vals.append(elemento[43])
    vals.append(elemento[44])
    vals.append(elemento[45])
    vals.append(elemento[46])
    vals.append(elemento[47])
    vals.append(str(elemento[48]))
    vals.append(elemento[49])
    vals.append(elemento[50])
    vals.append(elemento[51])
    vals.append(elemento[53])
    try:
        vals.append(elemento[59])
    except:
        vals.append("")
    return vals


def inserElementsFromArray(name, lista):
    """Insert element to a database
    lista: array with component data in text format"""
    conn = sqlite3.connect(name)
    curs = conn.cursor()
    curs.execute("SELECT COUNT(*) AS Total FROM compuestos")
    numero = curs.fetchone()[0]
    if name == databank_Custom_name:
        numero += 1000
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
        "formula", "nombre", "peso_molecular", "tc", "pc", "vc", "API",
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
        "tension_DIPPR_tmax", "momento_dipolar", "constante_volumen_liquido",
        "constante_rackett", "densidad_especifica", "factor_acentrico",
        "parametro_solubilidad", "watson", "MSRK_A", "MSRK_B", "Stiehl",
        "t_ebullicion", "t_fusion", "CAS_id", "formula_alternativa", "UNIFAC",
        "diametro_molecular", "Eps_k", "UNIQUAC_area", "UNIQUAC_volumen",
        "factor_acentrico_modificado", "calor_formacion_gas",
        "energia_libre_gas", "volumen_wilson", "calor_combustion_neto",
        "calor_combustion_bruto", "nombre_alternativo",
        "volumen_caracteristico", "calor_formacion_solido",
        "energia_libre_solido", "parametro_polar", "smile"]
    vals = transformElement(elemento)

    conn = sqlite3.connect(databank_Custom_name)
    curs = conn.cursor()
    for variable, valor in zip(variables, vals):
        if isinstance(valor, int):
            curs.execute('UPDATE compuestos SET %s=%i WHERE id==%i'
                         % (variable, valor, indice))
        elif isinstance(valor, float):
            curs.execute('UPDATE compuestos SET %s=%f WHERE id==%i'
                         % (variable, valor, indice))
        elif isinstance(valor, str):
            valor = '"'+valor+'"'
            curs.execute('UPDATE compuestos SET %s=%s WHERE id==%i'
                         % (variable, valor, indice))
    conn.commit()
    conn.close()


def deleteElement(indice):
    """Delete Element with indice from custom Database"""
    conn = sqlite3.connect(databank_Custom_name)
    curs = conn.cursor()
    curs.execute("DELETE FROM compuestos WHERE id=%i" % indice)
    conn.commit()
    conn.close()


def getElement(indice):
    """Get element from database
    indice: index in databank of element"""
    if indice > 1000:
        db_Custom = sqlite3.connect(databank_Custom_name).cursor()
        db_Custom.execute("select * from compuestos where id==%i" % indice)
        componente = db_Custom.fetchone()
    else:
        db = sqlite3.connect(databank_name).cursor()
        db.execute("select * from compuestos where id==%i" % indice)
        componente = db.fetchone()
    return componente


def copyElement(indice):
    """Create a copy of element of indice in custom Database"""
    elemento = getElement(indice)
    vals = elemento[1:]
    conn = sqlite3.connect(databank_Custom_name)
    curs = conn.cursor()
    curs.execute("INSERT INTO compuestos VALUES" +
                 str((1001+N_comp_Custom, ) + vals))
    conn.commit()
    conn.close()


if __name__ == "__main__":
    copyElement(5)
