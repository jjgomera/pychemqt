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



##############################
###   librería de definición de spreadsheet   ###
##############################

import os
import string

try:
    import ezodf
    import openpyxl
except:
    pass

from PyQt5.QtWidgets import  QApplication

from .parents import equipment


class Spreadsheet(equipment):
    """Clase que define un la interaccion con un hoja de calculo de libreoffice/openoffice

    Parámetros:
        project: instancia project
        input: entity de entrada
        output: entity de salida
        filename: Path del archivo ods
        datamap: Array con la estructura de datos a traspasar,
            cada elemento es un dicciontario con los datos de conversion de valores
                entity: corriente o equipo del que exportar el valor
                variable: nombre del valor de la variable de la corriente
                unidad: unidad del valor a pasar en el caso de magitudes
                hoja: Nombre de la hoja
                celda: Celda en la que colocar el dato

    """
    title = QApplication.translate("pychemqt", "Spreadsheet")
    help = ""
    kwargs = {
        "project": None,
        "input": "",
        "output": "",
        "filename": "",
        "datamap": []}
    kwargs_forbidden = ["project", ]

    @property
    def isCalculable(self):
        self.msg = ""
        self.status = 1
        if not self.kwargs["filename"] or \
                not os.path.isfile(self.kwargs["filename"]):
            self.msg = QApplication.translate(
                "pychemqt", "undefined spreadsheet filename")
            self.status = 0
            return
        if not self.kwargs["datamap"]:
            self.msg = QApplication.translate(
                "pychemqt", "undefined spreadsheet data map")
            self.status = 3
        return True

    def cleanOldValues(self, **kwargs):
        """Si se cambia la ruta de la hoja de cálculo se reinicia el datamap"""
        if kwargs.get("filename", "") and \
                kwargs.get("filename", "") != self.kwargs["filename"]:
            self.kwargs["datamap"] = []
        self.kwargs.update(kwargs)

    def calculo(self):
        ext = self.kwargs["filename"].split(".")[-1]
        if ext == "ods":
            spreadsheet = ezodf.opendoc(self.kwargs["filename"])
            self.sheets = [name for name in spreadsheet.sheets.names()]
            if self.kwargs["datamap"]:
                for data in self.kwargs["datamap"]:
                    entity = self.kwargs["project"].getObject(data["entity"])
                    sheet = spreadsheet.sheets[data["sheet"]]
                    indProp = entity.propertiesTitle().index(data["property"])
                    if entity.propertiesUnit()[indProp] == str:
                        value = entity.__getattribute__(
                            entity.propertiesAttribute()[indProp])
                    else:
                        indUnit = entity.propertiesUnit()[indProp].__text__.index(data["unit"])
                        units = entity.propertiesUnit()[indProp].__units__
                        value = entity.__getattribute__(entity.propertiesAttribute()[indProp]).__getattribute__(units[indUnit])

                    # Chequear
                    celda = list(data["cell"])
                    column = []
                    while celda[0] in string.ascii_uppercase:
                        column.append(celda.pop(0))
                    base = len(string.ascii_uppercase)
                    exponente = 0
                    columna = 0
                    while column:
                        ordinal = ord(column.pop())-64
                        columna += ordinal*base**exponente
                        exponente += 1
                    fila = int("".join(celda))
                    if fila > sheet.nrows():
                        sheet.append_rows(fila-sheet.nrows())
                    if columna > sheet.ncols():
                        sheet.append_columns(columna-sheet.ncols())

                    sheet[data["cell"]].set_value(value)
                spreadsheet.save()

        elif ext == "xlsx":
            spreadsheet = openpyxl.load_workbook(self.kwargs["filename"])
            self.sheets = spreadsheet.get_sheet_names()
            if self.kwargs["datamap"]:
                for data in self.kwargs["datamap"]:
                    entity = self.kwargs["project"].getObject(data["entity"])
                    sheet = spreadsheet[data["sheet"]]
                    indProp = entity.propertiesTitle().index(data["property"])
                    if entity.propertiesUnit()[indProp] == str:
                        value = entity.__getattribute__(entity.propertiesAttribute()[indProp])
                    else:
                        indUnit = entity.propertiesUnit()[indProp].__text__.index(data["unit"])
                        units = entity.propertiesUnit()[indProp].__units__
                        value = entity.__getattribute__(entity.propertiesAttribute()[indProp]).__getattribute__(units[indUnit])
                    sheet[data["cell"]] = value
                    comentario = openpyxl.comments.Comment("{0[entity]}.{0[property]}.{0[unit]} ---> {0[sheet]}.{0[cell]}".format(data), 'pychemqt')
                    sheet[data["cell"]].comment = comentario
                spreadsheet.save(".".join(self.kwargs["filename"].split(".")[:-1])+"-bak"+".xlsx")

            elif ext == "xls":
                # TODO: Implement old office support
                pass

        self.salida = []

    def writeListtoJSON(self, kwarg, key, value):
        """Personalizar en el caso de equipos con listas complejas"""
        kwarg_list = {}
        if key == "datamap":
            for i, data in enumerate(value):
                kwarg_list[i] = data
        kwarg[key] = kwarg_list

    def readListFromJSON(self, data, key):
        """Read list from file, customize in entities with complex list"""
        kwarg = []
        if key == "datamap":
            for i, data in data[key].items():
                kwarg.append(data)
        return kwarg

    def propTxt(self):
        txt = "#---------------"
        txt += QApplication.translate("pychemqt", "Data map")
        txt += "-----------------#" + os.linesep
        txt += self.propertiesToText(0)
        if self.kwargs["datamap"]:
            for data in self.kwargs["datamap"]:
                txt += "{0[entity]}.{0[property]}.{0[unit]} ---> {0[sheet]}.{0[cell]}".format(data)+os.linesep
        else:
            txt += QApplication.translate("pychemqt", "Undefined")+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        l = [(QApplication.translate("pychemqt", "Spreadsheet path"),
              "filename", str),
             (QApplication.translate("pychemqt", "Data map"), "datamap", None)]
        return l

    def propertiesListTitle(self, index):
        lista = []
        for data in self.kwargs["datamap"]:
            lista.append("{0[entity]}.{0[property]}.{0[unit]} ---> {0[sheet]}.{0[cell]}".format(data))
        return lista

    def writeStatetoJSON(self, state):
        """Write instance parameter to file"""
        state["sheets"] = self.sheets

    def readStatefromJSON(self, state):
        """Load instance parameter from saved file"""
        self.sheets = state["sheets"]
        self.salida = [None]


if __name__ == '__main__':
    spreadsheet=Spreadsheet(filename="/media/datos/ejemplo.ods")
#    spreadsheet = ezodf.opendoc("/media/datos/ejemplo.ods")
#    hoja=spreadsheet.sheets["Prueba"]
#    hoja["B6"].set_value(5.)
#    hoja["B7"].set_value(4.)
#    spreadsheet.save()
#    print hoja["D6"].value

#    spreadsheet=Spreadsheet(filename="/media/datos/ejemplo.xlsx")
#    ws=wb["Calculos"]
#    ws['I4']=32
#    wb.save('/media/datos/ejemplo2.xlsx')

#    import xlrd
#
#    book = xlrd.open_workbook("/media/datos/ejemplo.xls")
#
#    for sheet_name in book.sheet_names():
#        sheet = book.sheet_by_name(sheet_name)
##       print sheet.row_values(0)[0]
#        print sheet_name
#    print dir(book)
#    book.save()
