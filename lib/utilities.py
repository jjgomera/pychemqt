#!/usr/bin/python
# -*- coding: utf-8 -*-

#################################################################################
# Module with utilities:
#   - format2txt: Function to convert dict format config in a string value
#   - representacion: Function for string representation of float values
#   - colors: Function to generate colors
#   - exportTable; Save data to a file
#################################################################################


import random
import os
from string import maketrans

from PyQt4.QtGui import QApplication


def format2txt(formato):
    """Function to convert dict format config in a string equivalent"""
    if formato["signo"]:
        txt = "+"
    else:
        txt = ""
    if formato["format"] == 0:
        txt += "{total}.{decimales} fixed".format(**formato)
    elif formato["format"] == 1:
        txt += "{decimales} sign".format(**formato)
    elif formato["format"] == 2:
        txt += "{decimales} exp".format(**formato)
    if formato.get("exp", False):
        txt += " ({tol} exp)".format(**formato)
    return txt


def representacion(float, format=0, total=0, decimales=4, exp=False, tol=4, signo=False, thousand=False):
    """Function for string representation of float values
    float: number to transform
    format: mode
        0   -   fixed point
        1   -   Significant figures
        2   -   Engineering format
    total: total number of digits
    decimales: decimal number
    exp: boolean to use engineering repr for big or small number
    tol: exponent limit tu use engineering repr over float normal repr
    signo: show sign
    thousand: use thousan separator point
    """
    if type(float) is str:
        return float

    if signo:
        start="{:+"
    else:
        start="{: "

    if thousand:
        coma=",."
    else:
        coma="."

    if exp:
        if -10**tol > float or -10**-tol < float < 10**-tol or float > 10**tol:
            format=2

    if format==1:
        string=start+"{}{:d}g".format(coma, decimales)+"}"
    elif format==2:
        string=start+"{:d}{}{:d}e".format(total, coma, decimales)+"}"
    else:
        string=start+"{:d}{}{:d}f".format(total, coma, decimales)+"}"

    return string.format(float)


def colors(number, mix=""):
    """Function to generate colors
    Input:
        number: number of required colors
        mix: string name color for mix
    Output:
        Array with color hex repr
    """
    colors = []
    for i in number:
        red = random.randint(0, 255)
        green = random.randint(0, 255)
        blue = random.randint(0, 255)
        if mix:
            red_mix = int(mix[1:3], base=16)
            red = (red + red_mix) / 2
            green_mix = int(mix[3:5], base=16)
            green = (green + green_mix) / 2
            blue_mix = int(mix[5:], base=16)
            blue = (blue + blue_mix) / 2
        
        colors.append(('#%02X%02X%02X' % (r(),r(),r())))
    return
    
    
def exportTable(matrix, fname, format, title=None):
    """Save data to a file
    Inputs
        matrix: array with data to save
        fname: name of file to save
        format: name of format to save
            csv | ods | xls | xlsx
        title: column title array, optional
    """
    sheetTitle = unicode(QApplication.translate("pychemqt", "Table"))
    if fname.split(".")[-1] != format:
        fname+=".%s" % format

    # Format title
    if title:
        header = []
        for ttl in title:
            line = unicode(ttl).split(os.linesep)
            if line[-1] != "[-]":
                line[-1] = "["+line[-1]+"]"
            header.append(" ".join(line))
        c_newline=maketrans(os.linesep, " ")

    if format == "csv":
        import csv
        with open(fname, "w") as archivo:
            writer = csv.writer(archivo, delimiter='\t', quotechar='"', quoting=csv.QUOTE_NONE)
            
            # Add Data
            if title:
                writer.writerow([ttl.translate(c_newline) for ttl in header])
            c_float=maketrans(".", ",")
            for row in matrix:
                writer.writerow([str(data).translate(c_float) for data in row])

    elif format == "ods":
        import ezodf
        spreadsheet = ezodf.newdoc("ods", fname)
        sheets = spreadsheet.sheets
        sheet=ezodf.Table(sheetTitle)
        sheets+=sheet
        sheet.reset(size=(len(matrix)+1, len(matrix[0])))
        
        # Add Data
        if title:
            for i, ttl in enumerate(header):
                sheet["%s%i"%(spreadsheetColumn(i), 1)].set_value(ttl)
        for j, row in enumerate(matrix):
            for i, data in enumerate(row):
                sheet["%s%i"%(spreadsheetColumn(i), j+2)].set_value(data)
        spreadsheet.save()

    elif format == "xls":
        import xlwt
        spreadsheet = xlwt.Workbook()
        sheet = spreadsheet.add_sheet(sheetTitle)
        
        font = xlwt.Font()
        font.bold = True
        style = xlwt.XFStyle()
        style.font = font

        # Add Data
        if title:
            for i, ttl in enumerate(header):
                sheet.write(0, i, ttl, style)
        for j, row in enumerate(matrix):
            for i, data in enumerate(row):
                sheet.write(j+1, i, data)
        spreadsheet.save(fname)

    elif format == "xlsx":
        import openpyxl
        from openpyxl.styles import Style, Font
        spreadsheet = openpyxl.Workbook()
        sheet = spreadsheet.active
        sheet.title = sheetTitle
        
        font1 = Font()
        font1.size = 9
        font1.bold = True
        font2 = Font()
        font2.size = 9
        
        # Add Data
        if title:
            for i, ttl in enumerate(header):
                sheet["%s%i"%(spreadsheetColumn(i), 1)] = ttl
                sheet["%s%i"%(spreadsheetColumn(i), 1)].style.font= font1
        for j, row in enumerate(matrix):
            for i, data in enumerate(row):
                sheet["%s%i"%(spreadsheetColumn(i), j+2)] = data
                sheet["%s%i"%(spreadsheetColumn(i), j+2)].style.font = font2
        spreadsheet.save(filename=fname)
    
    else:
        raise ValueError(QApplication.translate(
            "pychemqt", "Unsopported format") + " " + format)


def spreadsheetColumn(index):
    """Procedure to convert index column in AAA spreadsheet column namestyle
    Input:
        index: index of column start with 0
    """
    index += 1
    letters = ""
    while index:
        mod = index % 26
        index = index // 26
        letters += chr(mod + 64)
    return "".join(reversed(letters))


if __name__ == "__main__":
#    import math
#    print representacion(math.pi, decimales=6, tol=1)
#    print repr(Configuracion("Density", "DenGas").text())
#    print representacion("3232326262")

    print spreadsheetColumn(55)
