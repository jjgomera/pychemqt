#!/usr/bin/python
# -*- coding: utf-8 -*-

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

from PyQt4.QtGui import QApplication

from parents import equipment


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
    title=QApplication.translate("pychemqt", "Spreadsheet")
    help=""
    kwargs={"project": None,
                    "input": "",
                    "output": "",
                    "filename": "",
                    "datamap": []}
    kwargs_forbidden=["project", ]

    @property
    def isCalculable(self):
        self.msg=""
        self.status=1
        if not self.kwargs["filename"] or not os.path.isfile(self.kwargs["filename"]):
            self.msg=QApplication.translate("pychemqt", "undefined spreadsheet filename")
            self.status=0
            return
        if not self.kwargs["datamap"]:
            self.msg=QApplication.translate("pychemqt", "undefined spreadsheet data map")
            self.status=3
        return True

    def cleanOldValues(self, **kwargs):
        """Si se cambia la ruta de la hoja de cálculo se reinicia el datamap"""
        if kwargs.get("filename", "") and kwargs.get("filename", "")!=self.kwargs["filename"]:
            self.kwargs["datamap"]=[]
        self.kwargs.update(kwargs)

    def calculo(self):
        ext=self.kwargs["filename"].split(".")[-1]
        if ext=="ods":
            spreadsheet = ezodf.opendoc(self.kwargs["filename"])
            self.sheets=[name for name in spreadsheet.sheets.names()]
            if self.kwargs["datamap"]:
                for data in self.kwargs["datamap"]:
                    entity=self.kwargs["project"].getObject(data["entity"])
                    sheet=spreadsheet.sheets[data["sheet"]]
                    indProp=entity.propertiesTitle().index(data["property"])
                    if entity.propertiesUnit()[indProp]==str:
                        value=entity.__getattribute__(entity.propertiesAttribute()[indProp])
                    else:
                        indUnit=entity.propertiesUnit()[indProp].__text__.index(data["unit"])
                        units=entity.propertiesUnit()[indProp].__units__
                        value=entity.__getattribute__(entity.propertiesAttribute()[indProp]).__getattribute__(units[indUnit])

                    #Chequear
                    celda=list(data["cell"])
                    column=[]
                    while celda[0] in string.uppercase:
                        column.append(celda.pop(0))
                    base=len(string.uppercase)
                    exponente=0
                    columna=0
                    while column:
                        ordinal=ord(column.pop())-64
                        columna+=ordinal*base**exponente
                        exponente+=1
                    fila=int("".join(celda))
                    if fila>sheet.nrows():
                        sheet.append_rows(fila-sheet.nrows())
                    if columna>sheet.ncols():
                        sheet.append_columns(columna-sheet.ncols())

                    sheet[data["cell"]].set_value(value)
                spreadsheet.save()

        elif ext=="xlsx":
            spreadsheet = openpyxl.load_workbook(self.kwargs["filename"])
            self.sheets=spreadsheet.get_sheet_names()
            if self.kwargs["datamap"]:
                for data in self.kwargs["datamap"]:
                    entity=self.kwargs["project"].getObject(data["entity"])
                    sheet=spreadsheet[data["sheet"]]
                    indProp=entity.propertiesTitle().index(data["property"])
                    if entity.propertiesUnit()[indProp]==str:
                        value=entity.__getattribute__(entity.propertiesAttribute()[indProp])
                    else:
                        indUnit=entity.propertiesUnit()[indProp].__text__.index(data["unit"])
                        units=entity.propertiesUnit()[indProp].__units__
                        value=entity.__getattribute__(entity.propertiesAttribute()[indProp]).__getattribute__(units[indUnit])
                    sheet[data["cell"]] = value
                    comentario = openpyxl.comments.Comment("{0[entity]}.{0[property]}.{0[unit]} ---> {0[sheet]}.{0[cell]}".format(data), 'pychemqt')
                    sheet[data["cell"]].comment=comentario
                spreadsheet.save(".".join(self.kwargs["filename"].split(".")[:-1])+"-bak"+".xlsx")

        self.salida=[]


    def writeListtoStream(self, stream, key, value):
        """Personalizar en el caso de equipos con listas complejas"""
        if key=="datamap":
            stream.writeInt32(len(value))
            for data in value:
                stream.writeString(data["entity"])
                stream.writeString(data["property"])
                stream.writeString(data["unit"])
                stream.writeString(data["sheet"])
                stream.writeString(data["cell"])

    def readListFromStream(self, stream, key):
        """Personalizar en el caso de equipos con listas complejas"""
        valor=[]
        if key=="datamap":
            for i in range(stream.readInt32()):
                data={}
                data["entity"]=stream.readString()
                data["property"]=stream.readString()
                data["unit"]=stream.readString()
                data["sheet"]=stream.readString()
                data["cell"]=stream.readString()
                valor.append(data)
        return valor

    def propTxt(self):
        txt="#---------------"+QApplication.translate("pychemqt", "Data map")+"-----------------#"+os.linesep
        txt+=QApplication.translate("pychemqt", "Spreadsheet path")+": "+self.kwargs["filename"]+os.linesep
        if self.kwargs["datamap"]:
            for data in self.kwargs["datamap"]:
                txt+="{0[entity]}.{0[property]}.{0[unit]} ---> {0[sheet]}.{0[cell]}".format(data)+os.linesep
        else:
            txt+=QApplication.translate("pychemqt", "Undefined")+os.linesep

        return txt

    @classmethod
    def propertiesEquipment(cls):
        list=[(QApplication.translate("pychemqt", "Filename"), "filename", str),
                (QApplication.translate("pychemqt", "Data map"), "datamap", None)]
        return list

    def propertiesListTitle(self, index):
        lista=[]
        for data in self.kwargs["datamap"]:
            lista.append("{0[entity]}.{0[property]}.{0[unit]} ---> {0[sheet]}.{0[cell]}".format(data))
        return lista


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
