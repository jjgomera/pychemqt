#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=import-outside-toplevel

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


###############################################################################
# Module with utilities:
#   - representacion: Function for string representation of float values
#   - colors: Function to generate colors
#   - exportTable; Save data to a file
#   - formatLine: Return a matplotlib line formatting kw
###############################################################################


import os
import random
from math import exp

from tools.qt import translate


def representacion(number, fmt=0, total=0, decimales=4, eng=False, tol=5,
                   signo=False, thousand=False):
    """Function for string representation of float values

    Parameters
    ----------
    number : float
        number to transform
    fmt : mode
        0   -   fixed point
        1   -   Significant figures
        2   -   Engineering format
    total : integer
        total number of digits
    decimales : integer
        decimal number
    eng : boolean
        Use engineering repr for big or small number
    tol : integer
        Exponent limit to use engineering repr over float normal repr
    signo : boolean
        Show sign, use to force positive sign
    thousand : boolean
        Use thousand separator point

    Returns
    -------
    repr : str
        String representation

    Examples
    --------
    >>> import math
    >>> representacion(math.pi, signo=True)
    '+3.1416'
    """
    if isinstance(number, str):
        return number

    if signo:
        start = "{:+"
    else:
        start = "{: "

    if thousand:
        coma = ",."
    else:
        coma = "."

    if eng and (-10**tol > number
                or (-10**-tol < number < 10**(-tol+1) and number != 0)
                or number > 10**tol):
        fmt = 2
    if number == 0:
        decimales = 1

    if fmt == 1:
        string = start+"{}{:d}g".format(coma, decimales)+"}"
    elif fmt == 2:
        string = start+"{:d}{}{:d}e".format(total, coma, decimales)+"}"
    else:
        string = start+"{:d}{}{:d}f".format(total, coma, decimales)+"}"

    return string.format(number)


def colors(number, mix="", scale=False):
    """Function to generate colors

    Parameters
    ----------
    number: integer
        Number of required colors
    mix : str
        Name color for mix
    scale : boolean
        Define  output in a scale red-blue

    Returns
    -------
    rgb : list
        List with color hexadecimal code
    """
    rgb = []
    for i in range(number):
        if scale:
            red = int(255*(i/number))
            green = 0
            blue = int(255*((number-i)/number))
        else:
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

        rgb.append(('#%02X%02X%02X' % (red, green, blue)))
    return rgb


def exportTable(matrix, fname, ext, title=None):
    """Save data to a spreadsheet file

    Parameters
    ----------
    matrix : list
        array with data to save
    ext : str
        name of format to save, with supported csv | ods | xls | xlsx
    title : list, optional
        list with column title
    """
    sheetTitle = translate("Table")
    if fname.split(".")[-1] != ext:
        fname += ".%s" % ext

    # Format title
    header = []
    if title:
        for ttl in title:
            line = str(ttl).split(os.linesep)
            if line[-1] != "[-]":
                line[-1] = "["+line[-1]+"]"
            header.append(" ".join(line))

    if ext == "csv":
        import csv
        with open(fname, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=header)

            writer.writeheader()
            for row in matrix:
                kw = {}
                for ttl, value in zip(header, row):
                    kw[ttl] = value
                writer.writerow(kw)

    elif ext == "ods":
        import ezodf
        spreadsheet = ezodf.newdoc("ods", fname)
        sheets = spreadsheet.sheets
        sheet = ezodf.Table(sheetTitle)
        sheets += sheet
        sheet.reset(size=(len(matrix)+1, len(matrix[0])))

        # Add Data
        if title:
            for i, ttl in enumerate(header):
                sheet["%s%i" % (spreadsheetColumn(i), 1)].set_value(ttl)
        for j, row in enumerate(matrix):
            for i, data in enumerate(row):
                sheet["%s%i" % (spreadsheetColumn(i), j+2)].set_value(data)
        spreadsheet.save()

    elif ext == "xls":
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

    elif ext == "xlsx":
        import openpyxl
        from openpyxl.styles import Font
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
                sheet["%s%i" % (spreadsheetColumn(i), 1)] = ttl
                sheet["%s%i" % (spreadsheetColumn(i), 1)].style.font = font1
        for j, row in enumerate(matrix):
            for i, data in enumerate(row):
                sheet["%s%i" % (spreadsheetColumn(i), j+2)] = data
                sheet["%s%i" % (spreadsheetColumn(i), j+2)].style.font = font2
        spreadsheet.save(filename=fname)

    else:
        raise ValueError(translate("Unsopported format") + " " + ext)


def spreadsheetColumn(index):
    """Procedure to convert index column in AAA spreadsheet column namestyle

    Parameters
    ----------
    index : integer
        index of column start with 0

    Returns
    -------
    key : str
        column letter code, ej: A, C

    Examples
    --------
    >>> spreadsheetColumn(5)
    'F'
    >>> spreadsheetColumn(75)
    'BX'
    """
    index += 1
    letters = ""
    while index:
        mod = index % 26
        index = index // 26
        letters += chr(mod + 64)
    return "".join(reversed(letters))


def formatLine(config, section, name):
    """Return a matplotlib line formatting kw

    Parameters
    ----------
    config : configparser.Configparser instance
        configuration
    section : str
        Section name in Configparser
    name : str
        Line name prefix

    Returns
    -------
    kw : dict
        matplotlib kwargs formatting dictionary
    """
    kw = {}
    kw["ls"] = config.get(section, name+"lineStyle")
    kw["lw"] = config.getfloat(section, name+"lineWidth")
    kw["color"] = config.get(section, name+"Color")
    kw["alpha"] = config.getfloat(section, name+"alpha")/255

    kw["marker"] = config.get(section, name+"marker")
    kw["ms"] = config.getfloat(section, name+"markersize")
    kw["mfc"] = config.get(section, name+"markerfacecolor")
    kw["mew"] = config.getfloat(section, name+"markeredgewidth")
    kw["mec"] = config.get(section, name+"markeredgecolor")

    # Additional options for grid
    if config.has_option(section, name+"which"):
        kw["which"] = config.get(section, name+"which")
    if config.has_option(section, name+"which"):
        kw["axis"] = config.get(section, name+"axis")

    return kw


def SimpleEq(Tc, T, coef):
    r"""Common procedure for calculation of simple properties like the
    ancillary equation for vapor pressure, saturated densities of liquid
    and vapor

    The procedure define several equations:

    .. math::
        \frac{x}{x_r} = 1 + \sum_i N_i\theta^{t_i}

    .. math::
        \ln\left(\frac{x}{x_r}\right) = \sum_i N_i\theta^{t_i}

    .. math::
        \ln\left(\frac{x}{x_r}\right) = \frac{T_c}{T}\sum_i N_i\theta^{t_i}

    .. math::
        \theta = 1 - \frac{T}{Tc}

    The coef must define the parameters:

        * eq: Index of equation to use
        * n: Polinomial parameter
        * t: Exponent of tita


    Parameters
    ----------
    Tc : float
        Critical Temperature, [K]
    T : float
        Temperature, [K]
    coef: dict
        Coefficient for correlation
    """

    Tita = 1-T/Tc

    pr = 0
    for n, t in zip(coef["n"], coef["t"]):
        pr += n*Tita**t

    if coef["eq"] == 1:
        pr += 1

    if coef["eq"] == 3:
        pr *= Tc/T

    if coef["eq"] in [2, 3]:
        pr = exp(pr)

    return pr


def refDoc(doi, refs, tab=4):
    """Function decorator used to automatic addiction of References section
    to documentation of procedures

    Parameters
    ----------
    doi : dict
        Dictionary with library references
    refs : list
        List of number to report in References section
    """
    def decorator(f):
        f.__doc__ += os.linesep + os.linesep
        f.__doc__ += " "*tab + "References" + os.linesep
        f.__doc__ += " "*tab + "----------" + os.linesep
        for ref in refs:
            rf = doi[ref]
            f.__doc__ += " "*tab + "[%i]_ %s; %s. %s" % (
                ref, rf["autor"], rf["title"], rf["ref"]) + os.linesep
            f.__doc__ += os.linesep

        return f
    return decorator


if __name__ == "__main__":
    import math
    print(representacion(math.pi*1000, decimales=6, tol=5))
    print(representacion(0, decimales=6, tol=1))
    # print repr(Configuracion("Density", "DenGas").text())
    # print representacion("3232326262")

    # print(spreadsheetColumn(55))
    # print(colors(5))
