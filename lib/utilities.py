#!/usr/bin/python
# -*- coding: utf-8 -*-

#################################################################################
# Module with utilities:
#   - format2txt: function to convert dict format config in a string value
#   - representacion
#   - colors
#################################################################################


import random


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
    """Función que expresa un valor de tipo float en forma de string
    float: numero a representar
    format: tipo de modo
        0   -   fixed point
        1   -   Significant figures
        2   -   Engineering format
    total: numero total de digitos
    decimales: numero de decimales
    exp: boolean que indica si se usa notacion exponencial para numeros grandes y pequeños
    tol: potencia por encima de la cual se usa notacion exponencial
    signo: mostrar signo positivo
    thousand: usa separador para miles
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

#public Color generateRandomColor(Color mix) {
#    Random random = new Random();
#    int red = random.nextInt(256);
#    int green = random.nextInt(256);
#    int blue = random.nextInt(256);
#
#    // mix the color
#    if (mix != null) {
#        red = (red + mix.getRed()) / 2;
#        green = (green + mix.getGreen()) / 2;
#        blue = (blue + mix.getBlue()) / 2;
#    }
#
#    Color color = new Color(red, green, blue);
#    return color;
#}

def colors(number, mix=""):
    """Function to generate colors
    number: number of required colors
    mix: string name color for mix"""
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
    

if __name__ == "__main__":
    import math
    print representacion(math.pi, decimales=6, tol=1)
#    print repr(Configuracion("Density", "DenGas").text())
    print representacion("3232326262")
