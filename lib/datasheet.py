#!/usr/bin/python3
# -*- coding: utf-8 -*-

'''Pychemqt, Chemical Engineering Process simulator
Copyright (C) 2009-2022, Juan José Gómez Romera <jjgomera@gmail.com>

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
# Module for pdf integration
# TODO: to improve studying reportlab library
###############################################################################

import os

from lib import config

if os.environ["reportlab"] == "True":
    from reportlab.pdfgen import canvas
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import cm


class pdf():
    def __init__(self, titulo):
        self.c = canvas.Canvas("datasheet.pdf", pagesize=A4)
        self.c.setTitle(titulo)

    def dibujar(self):
        self.c.showPage()
        self.c.save()

    def marco(self, titulo):
        # marco
        self.c.setLineWidth(2)
        self.c.line(3.5*cm, 2.5*cm, 19.5*cm,2.5*cm)
        self.c.line(3.5*cm, 2.5*cm, 3.5*cm, 27.2*cm)
        self.c.line(3.5*cm, 27.2*cm,19.5*cm,27.2*cm)
        self.c.line(19.5*cm,2.5*cm,19.5*cm,27.2*cm)

        # encabezado
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 660, 19.5*cm,660)
        self.c.line(10*cm, 27.2*cm, 10*cm,660)
        self.c.line(3.5*cm, 741, 19.5*cm,741)
        self.c.setFont("Helvetica-Bold", 12)
        self.c.drawCentredString(14.75*cm, 753, titulo)

        self.c.setLineWidth(0.5)
        self.c.line(3.5*cm, 676, 10*cm,676)         #revisiones
        self.c.line(3.5*cm, 692, 10*cm,692)
        self.c.line(3.5*cm, 708, 10*cm,708)
        self.c.line(3.5*cm, 724, 10*cm,724)
        self.c.line(4.6*cm, 741, 4.6*cm, 660)
        self.c.line(6.4*cm, 741, 6.4*cm, 660)
        self.c.line(8.2*cm, 741, 8.2*cm, 660)
        self.c.setFont("Helvetica", 10)
        self.c.drawString(3.6*cm, 728.5, "Rev.")
        self.c.drawRightString(4.3*cm, 712.5, "1")
        self.c.drawRightString(4.3*cm, 696.5, "2")
        self.c.drawRightString(4.3*cm, 680.5, "3")
        self.c.drawRightString(4.3*cm, 664.5, "4")
        self.c.drawCentredString(5.5*cm, 728.5, "Sustit por")
        self.c.drawCentredString(7.3*cm, 728.5, "Fecha")
        self.c.drawCentredString(9.1*cm, 728.5, "Sustit a")

        self.c.line(10*cm, 714, 19.5*cm,714)        #propietario
        self.c.line(10*cm, 687, 19.5*cm,687)
        self.c.line(13.17*cm, 741, 13.17*cm, 660)
        self.c.line(16.33*cm, 741, 16.33*cm, 660)
        self.c.drawString(10.1*cm, 731, "Cliente")
        self.c.drawString(13.2*cm, 731, "Nº Equipo")
        self.c.drawString(16.4*cm, 731, "Página")
        self.c.drawString(10.1*cm, 705, "W.O.")
        self.c.drawString(13.2*cm, 705, "Nº Requerimiento")
        self.c.drawString(16.4*cm, 705, "Nº Especificación")
        self.c.drawString(10.1*cm, 678, "Area")
        self.c.drawString(13.2*cm, 678, "Suministrado por")
        self.c.drawString(16.4*cm, 678, "Instalado por")


    def bomba(self, bomba):
        titulo = "BOMBA CENTRÍFUGA"
        nombre = "Bomba impulsión tubería 1A"
        fluido = "Líquido madre"
        necesarios = 1
        tipo = "ANSI AABack PullOut"
        situacion = "Interior"
        fabricante = "Goulds"
        modelo = "3196STD"
        tamano = "1 x1-1/2 x8"
        iscorrosivo = "No corrosivo"
        corrosivos = "N/A"
        solidos = "No"
        peligros = ""
        caudal_diseno = 50
        caudal_normal = 50
        temperatura = 95
        presion = bomba.entrada.P
        viscosidad = bomba.entrada.Liquido.mu
        gravedad = 1.1

        self.marco(titulo)

        # General
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 640, 19.5*cm,640)
        self.c.line(3.5*cm, 625, 19.5*cm,625)
        self.c.setLineWidth(0.5)
        self.c.line(3.5*cm, 610, 19.5*cm,610)
        self.c.line(3.5*cm, 595, 19.5*cm,595)
        self.c.line(3.5*cm, 580, 19.5*cm,580)
        self.c.line(11.5*cm, 565, 11.5*cm,625)
        self.c.drawCentredString(11.5*cm, 647, nombre)
        self.c.drawCentredString(11.5*cm, 629, "General")
        self.c.setFont("Helvetica", 8)
        self.c.drawString(3.6*cm, 615, "Fluido:  "+fluido)
        self.c.drawString(3.6*cm, 600, "Bombas necesarias:  "+ str(necesarios))
        self.c.drawString(3.6*cm, 585, "Tipo:  "+ tipo)
        self.c.drawString(3.6*cm, 570, "Situación:  "+ situacion)
        self.c.drawString(11.6*cm, 615, "Fabricante:  "+ fabricante)
        self.c.drawString(11.6*cm, 600, "Modelo:  "+ modelo)
        self.c.drawString(11.6*cm, 585, "Tamaño:  "+ tamano)

        # Datos de proceso
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 565, 19.5*cm,565)
        self.c.line(3.5*cm, 550, 19.5*cm,550)
        self.c.setLineWidth(0.5)
        self.c.line(3.5*cm, 535, 19.5*cm,535)
        self.c.line(3.5*cm, 520, 19.5*cm,520)
        self.c.line(3.5*cm, 505, 19.5*cm,505)
        self.c.line(3.5*cm, 490, 19.5*cm,490)
        self.c.line(11.5*cm, 475, 19.5*cm,475)
        self.c.line(11.5*cm, 550, 11.5*cm,460)
        self.c.setFont("Helvetica", 10)
        self.c.drawCentredString(11.5*cm, 555, "Datos de proceso")
        self.c.setFont("Helvetica", 8)
        self.c.drawString(3.6*cm, 540, "Fluido:  "+ fluido)
        self.c.drawString(3.6*cm, 525, "Características:  "+ iscorrosivo)
        self.c.drawString(3.6*cm, 510, "Compuestos corrosivos:  "+ corrosivos)
        self.c.drawString(3.6*cm, 495, "Sólidos:  "+ solidos)
        self.c.drawString(3.6*cm, 480, "Peligros:  "+ peligros)
        self.c.drawString(11.6*cm, 540, "Caudal normal:  "+ str(caudal_normal) + "gdm")
        self.c.drawString(11.6*cm, 525, "Caudal de diseño:  "+ str(caudal_diseno) + "gdm")
        self.c.drawString(11.6*cm, 510, "Temperatura de bombeo:  "+ str(temperatura) + "ºF")
        self.c.drawString(11.6*cm, 495, "Presión de vapor @ T.B.:  "+ representacion(presion) + config.Configuracion("Pressure").text())
        self.c.drawString(11.6*cm, 480, "Viscosidad @ T.B.:  "+ representacion(viscosidad) + config.Configuracion("Viscosity").text())
        self.c.drawString(11.6*cm, 465, "Gravedad específica @ T.B.:  "+ str(gravedad))

        # Condiciones de bombeo
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 445, 19.5*cm,445)
        self.c.line(3.5*cm, 460, 19.5*cm,460)
        self.c.setLneWidth(0.5)
        self.c.line(3.5*cm, 445, 19.5*cm,445)
        self.c.line(3.5*cm, 430, 19.5*cm,430)
        self.c.line(3.5*cm, 415, 19.5*cm,415)
        self.c.line(3.5*cm, 400, 19.5*cm,400)
        self.c.line(3.5*cm, 385, 19.5*cm,385)
        self.c.line(3.5*cm, 370, 19.5*cm,370)
        self.c.line(3.5*cm, 355, 19.5*cm,355)
        self.c.line(3.5*cm, 340, 19.5*cm,340)
        self.c.line(3.5*cm, 325, 19.5*cm,325)
        self.c.line(3.5*cm, 310, 19.5*cm,310)
        self.c.line(3.5*cm, 295, 19.5*cm,295)
        self.c.line(3.5*cm, 280, 19.5*cm,280)
        self.c.line(3.5*cm, 265, 19.5*cm,265)
        self.c.line(8.5*cm, 445, 8.5*cm,250)
        self.c.line(14*cm, 445, 14*cm,250)
        self.c.setFont("Helvetica", 10)
        self.c.drawCentredString(11.5*cm, 450, "Condiciones de bombeo")
        self.c.setFont("Helvetica", 8)
        self.c.drawCentredString(11.25*cm, 433, "Succión")
        self.c.drawCentredString(16.75*cm, 433, "Impulsión")
        self.c.drawString(3.6*cm, 420, "Presión, mmHg")


        # Datos mecánicos
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 250, 19.5*cm,250)
        self.c.line(3.5*cm, 235, 19.5*cm,235)
        self.c.setLineWidth(0.5)
        self.c.line(3.5*cm, 220, 19.5*cm,220)
        self.c.line(3.5*cm, 205, 19.5*cm,205)
        self.c.line(3.5*cm, 190, 19.5*cm,190)
        self.c.line(3.5*cm, 175, 19.5*cm,175)
        self.c.line(3.5*cm, 160, 19.5*cm,160)
        self.c.line(3.5*cm, 145, 19.5*cm,145)

        # Datos eléctricos
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 130, 19.5*cm,130)
        self.c.line(3.5*cm, 115, 19.5*cm,115)
        self.c.setLineWidth(0.5)
        self.c.line(3.5*cm, 100, 19.5*cm,100)
        self.c.line(3.5*cm, 85, 19.5*cm,85)
        self.c.line(3.5*cm, 70, 19.5*cm,70)

    def ciclon(self, nombre, ciclon):
#               numero, densidad_gas, densidad_solido, viscosidad, caudal, solidos_entrada, solidos_salida, rendimiento, delta_P, Dc, Bc, Hc, Jc, Lc, Zc, De, Sc, caudal_unitario, modelo, f_instal, indice, c_adqui, c_instal):
        titulo = "CICLÓN"
        self.marco(titulo)

        # variables
        temperatura = 45
        fluido = "Aire"
        solidos = "Polvo"

        # General
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 640, 19.5*cm,640)
        self.c.line(3.5*cm, 625, 19.5*cm,625)
        self.c.setLineWidth(0.5)
        self.c.line(3.5*cm, 610, 19.5*cm,610)
        self.c.line(3.5*cm, 595, 19.5*cm,595)
        self.c.line(3.5*cm, 580, 19.5*cm,580)
        self.c.line(3.5*cm, 565, 19.5*cm,565)
        self.c.line(3.5*cm, 550, 19.5*cm,550)
        self.c.line(9*cm, 580, 9*cm,610)
        self.c.line(14*cm, 580, 14*cm,610)
        self.c.line(11.5*cm, 550, 11.5*cm,580)
        self.c.drawCentredString(11.5*cm, 647, nombre)
        self.c.drawCentredString(11.5*cm, 629, "General")
        self.c.setFont("Helvetica", 8)
        self.c.drawString(3.6*cm, 615, "Temperatura:  " + representacion(ciclon.entrada.T.config())+config.Configuracion("Temperature").text())
        self.c.drawString(3.6*cm, 600, "Fluido:  " + fluido)
        self.c.drawString(3.6*cm, 585, "Sólidos:  " + solidos)
        self.c.drawString(9.2*cm, 600, "Densidad @ T.P.:  " + representacion(ciclon.entrada.Gas.rho.config("DenGas"))+config.Configuracion("Density", "DenGas").text())
        self.c.drawString(14.2*cm, 600, "Viscosidad @ T.P.:  " + representacion(ciclon.entrada.Gas.mu.config())+config.Configuracion("Viscosity").text())
        self.c.drawString(9.2*cm, 585, "Densidad @ T.P.:  " + representacion(ciclon.entrada.solido.rho.config("DenLiq"))+config.Configuracion("Density", "DenLiq").text())
        self.c.drawString(3.6*cm, 570, "Caudal:  " + representacion(ciclon.entrada.Q.config("QGas"))+config.Configuracion("VolFlow", "QGas").text())
        self.c.drawString(11.6*cm, 570, "Eficiencia:  " + representacion(ciclon.rendimiento))
        self.c.drawString(3.6*cm, 555, "Concentración sólidos entrada:  " + representacion(ciclon.entrada.solido.caudal.config())+config.Configuracion("MassFlow").text())
        self.c.drawString(11.6*cm, 555, "Concentración sólidos salida:  " + representacion(ciclon.SalidaAire.solido.caudal.config())+config.Configuracion("MassFlow").text())
        self.c.drawString(3.6*cm, 540, "Perdida de presión:  " + representacion(ciclon.DeltaP.config())+config.Configuracion("Pressure").text())

        # dimensiones
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 535, 19.5*cm,535)
        self.c.line(3.5*cm, 520, 19.5*cm,520)
        self.c.setLineWidth(0.5)
        self.c.setFont("Helvetica", 10)
        self.c.drawCentredString(11.5*cm, 525, "Dimensiones")
        self.c.setFont("Helvetica", 8)
        self.c.drawString(3.7*cm, 510, "Ciclones necesarios:  " + str(ciclon.num_ciclones))
        self.c.line(3.5*cm, 505, 10*cm,505)
        self.c.line(3.5*cm, 490, 10*cm,490)
        self.c.line(10*cm, 520, 10*cm,490)
        self.c.drawString(3.7*cm, 495, "Caudal Unitario:  " + representacion(ciclon.entrada.Q.config("QGas"))+config.Configuracion("VolFlow", "QGas").text())
        self.c.drawString(6.8*cm, 465, "Dc:  " + representacion(ciclon.Dc.config())+config.Configuracion("Length").text())
        self.c.drawString(6.8*cm, 450, "Bc:  " + representacion(ciclon.Bc.config())+config.Configuracion("Length").text())
        self.c.drawString(6.8*cm, 435, "Hc:  " + representacion(ciclon.Hc.config())+config.Configuracion("Length").text())
        self.c.drawString(6.8*cm, 420, "Jc:  " + representacion(ciclon.Jc.config())+config.Configuracion("Length").text())
        self.c.drawString(6.8*cm, 405, "Lc:  " + representacion(ciclon.Lc.config())+config.Configuracion("Length").text())
        self.c.drawString(6.8*cm, 390, "Zc:  " + representacion(ciclon.Zc.config())+config.Configuracion("Length").text())
        self.c.drawString(6.8*cm, 375, "De:  " + representacion(ciclon.De.config())+config.Configuracion("Length").text())
        self.c.drawString(6.8*cm, 360, "Sc:  " + representacion(ciclon.Sc.config())+config.Configuracion("Length").text())
        self.c.drawImage(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equip/ciclon_datasheet.gif"), 11*cm, 260, width=180,preserveAspectRatio=True)

        # Costes
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 300, 19.5*cm,300)
        self.c.line(3.5*cm, 285, 19.5*cm,285)
        self.c.setLineWidth(0.5)
        self.c.setFont("Helvetica", 10)
        self.c.drawCentredString(11.5*cm, 290, "Costes")
        self.c.setFont("Helvetica", 8)
        self.c.line(3.5*cm, 270, 19.5*cm,270)
        self.c.line(3.5*cm, 255, 19.5*cm,255)
        self.c.line(11.5*cm, 240, 11.5*cm,270)
        self.c.drawString(3.7*cm, 275, "Modelo: " + str(ciclon.tipoCosto))
        self.c.drawString(3.7*cm, 260, "Factor Instalación: " + str(ciclon.f_install))
        self.c.drawString(11.6*cm, 260, "Indice de costes: " + str(ciclon.Current_index))
        self.c.drawString(3.7*cm, 245, "Coste Adquisición: $" + str(ciclon.C_adqTotal))
        self.c.drawString(11.6*cm, 245, "Coste Instalación: $" + str(ciclon.C_instTotal))

        # Notas
        self.c.setLineWidth(1)
        self.c.line(3.5*cm, 240, 19.5*cm,240)
        self.c.setLineWidth(0.5)
        self.c.setFont("Helvetica", 10)
        self.c.drawString(3.6*cm, 225, "Notas:")
        self.c.setFont("Helvetica", 8)


if __name__ == "__main__":
    pdf = pdf("CICLÓN")
    pdf.ciclon("Ciclon limpieza polvo", 12, 1.2, 1200, 1.84e-5, 126000, 93.6, 3.5, 99.56, 6500, 1, 0.2, 0.5, 0.375, 1.5, 2.5, 0.5, 0.5, 10.2, "Industrial", 1.4, 549.2, 157088, 219924)
    pdf.dibujar()

    #Si queremos abrirlo con algún visor de pdf
    os.system("evince datasheet.pdf")
