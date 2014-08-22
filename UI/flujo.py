#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Module for graphics elements of PFD
###############################################################################


from functools import partial
from datetime import datetime
import tempfile
import os
import subprocess

from PyQt4 import QtCore, QtGui, QtXml, QtSvg

from lib import unidades
from lib.project import Project
from lib.thread import WaitforClick
from lib.config import Preferences
from lib.corriente import Corriente
from UI import texteditor, UI_corriente
from UI.plots import Plot_Distribucion
from UI.widgets import createAction, Table_Graphics, PathConfig
from tools.UI_Preferences import ConfLine
from equipment import *


# Value for mouse wheel zoom
factor = 5


def loadFromFile(path):
    fh = QtCore.QFile(path)
    if not fh.open(QtCore.QIODevice.ReadOnly):
        raise IOError, unicode(fh.errorString())
    stream = QtCore.QDataStream(fh)

    magic = stream.readInt32()
    if magic != Project.MAGIC_NUMBER:
        raise IOError, "unrecognized file type"
    version = stream.readInt32()
    if version < Project.FILE_VERSION:
        raise IOError, "old and unreadable file format"
    elif version > Project.FILE_VERSION:
        raise IOError, "new and unreadable file format"
    stream.setVersion(QtCore.QDataStream.Qt_4_2)

    project=Project()
    project.loadFromStream(stream, False)
    return project

class SelectStreamProject(QtGui.QDialog):
    project=None
    def __init__(self, parent=None):
        super(SelectStreamProject, self).__init__(parent)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Select stream from file"))

        layout = QtGui.QGridLayout(self)
        label=QtGui.QApplication.translate("pychemqt", "Project path")+":"
        msg=QtGui.QApplication.translate("pychemqt", "Select pychemqt project file")
        patrones=QtCore.QStringList()
        patrones.append(QtGui.QApplication.translate("pychemqt", "pychemqt project file")+" (*.pcq)")
        patron=patrones.join(";;")
        self.filename=PathConfig(label, msg=msg, patron=patron)
        self.filename.valueChanged.connect(self.changeproject)
        layout.addWidget(self.filename,1,1,1,3)

        layout.addWidget(QtGui.QLabel(QtGui.QApplication.translate("pychemqt", "Streams")),2,1)
        self.stream=QtGui.QComboBox()
        layout.addWidget(self.stream,2,2)
        layout.addItem(QtGui.QSpacerItem(10,10,QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding),3,3)

        self.status=QtGui.QLabel()
        layout.addWidget(self.status,10,1)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.button(QtGui.QDialogButtonBox.Ok).setEnabled(False)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox,10,2,1,2)

    def changeproject(self, path):
        self.status.setText(QtGui.QApplication.translate("pychemqt", "Loading project..."))
        QtGui.QApplication.processEvents()
        try:
            self.project=loadFromFile(path)
            self.buttonBox.button(QtGui.QDialogButtonBox.Ok).setEnabled(True)
            self.status.setText(QtGui.QApplication.translate("pychemqt", "Project loaded succesfully"))
        except Exception as e:
            print e
            self.status.setText(QtGui.QApplication.translate("pychemqt", "Failed to loading project..."))
        self.stream.clear()
        for stream in self.project.streams.keys():
            self.stream.addItem("%i" %stream)



class ConfLineDialog(QtGui.QDialog, ConfLine):
    """Dialogo de definición de formatos de líneas"""
    def __init__(self, pen=None, parent=None):
        super(ConfLineDialog, self).__init__(pen)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Edit format line"))
        buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)
        self.layout().addWidget(buttonBox,10,1,1,5)


class TextItemDlg(QtGui.QDialog):
    """Dialogo de edición de textos del diagrama de flujo"""
    def __init__(self, text=None, parent=None):
        super(TextItemDlg, self).__init__(parent)
        layout = QtGui.QGridLayout(self)
        self.editor= texteditor.TextEditor()
        self.editor.notas.textChanged.connect(self.updateUi)
        layout.addWidget(self.editor, 1, 1, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(QtGui.QDialogButtonBox.Ok|QtGui.QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        layout.addWidget(self.buttonBox, 2, 1, 1, 1)
        self.editor.notas.setFocus()
        if text:
            self.editor.setText(text)
        self.setWindowTitle(QtGui.QApplication.translate("pychemqt", "Edit text"))
        self.updateUi()

    def updateUi(self):
        self.buttonBox.button(QtGui.QDialogButtonBox.Ok).setDisabled(self.editor.notas.toPlainText().isEmpty())



class GeometricItem(object):
    """Clase genérica con la funcionalidad comun de los elementos geométricos"""
    def __init__(self, parent=None):
        super(GeometricItem, self).__init__(parent)
        self.setPen(self._pen())
        self.setFlags(QtGui.QGraphicsItem.ItemIsSelectable|QtGui.QGraphicsItem.ItemIsMovable|QtGui.QGraphicsItem.ItemSendsGeometryChanges|QtGui.QGraphicsItem.ItemIsFocusable)
        self.setZValue(-1)

    def _pen(self):
        pen=QtGui.QPen(QtGui.QColor(Preferences.get("PFD", 'Color_Stream')))
        pen.setWidthF(Preferences.getfloat("PFD", 'Width'))
        pen.setJoinStyle([QtCore.Qt.MiterJoin, QtCore.Qt.BevelJoin, QtCore.Qt.RoundJoin][Preferences.getint("PFD", 'Union')])
        pen.setMiterLimit(Preferences.getfloat("PFD", 'Miter_limit'))
        pen.setCapStyle([QtCore.Qt.FlatCap, QtCore.Qt.RoundCap, QtCore.Qt.SquareCap][Preferences.getint("PFD", 'Punta')])
        pen.setStyle([QtCore.Qt.SolidLine, QtCore.Qt.DashLine, QtCore.Qt.DotLine, QtCore.Qt.DashDotLine, QtCore.Qt.DashDotDotLine][Preferences.getint("PFD", 'Guion')])
        pen.setDashOffset(Preferences.getfloat("PFD", 'Dash_offset'))
        return pen

    def delete(self):
        self.scene().delete(self)

    def format(self):
        dialog=ConfLineDialog(self.pen())
        if dialog.exec_():
            pen=dialog.pen()
            self.setPen(pen)
            self.itemChange(QtGui.QGraphicsItem.ItemPositionChange, 0)

    def contextMenu(self):
        contextMenu= QtGui.QMenu("%s Item" % self.type, self.scene().parent())
        contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(self.icon)))
        contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/editDelete.png"), QtGui.QApplication.translate("pychemqt", "Delete"), self.delete)
        contextMenu.addSeparator()
        contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Appearance"), self.format)
        return contextMenu

    def itemChange(self, change, variant):
        if self.scene():
            if change == QtGui.QGraphicsItem.ItemPositionChange:
                if self.scene().parent().dirty[self.scene().parent().idTab]==False:
                    self.scene().parent().dirty[self.scene().parent().idTab]=True
                    self.scene().parent().activeControl(True)
                    self.scene().parent().tabModified(self.scene().parent().idTab)
        return QtGui.QGraphicsItem.itemChange(self, change, variant)

    def keyPressEvent(self, event):
        if event.modifiers() & QtCore.Qt.ShiftModifier:
            if event.key()==QtCore.Qt.Key_Up:
                rect=self.rect()
                rect.setBottom(self.rect().bottom()-factor)
                self.setRect(rect)
            elif event.key()==QtCore.Qt.Key_Down:
                rect=self.rect()
                rect.setBottom(self.rect().bottom()+factor)
                self.setRect(rect)
            elif event.key()==QtCore.Qt.Key_Left:
                rect=self.rect()
                rect.setRight(self.rect().right()-factor)
                self.setRect(rect)
            elif event.key()==QtCore.Qt.Key_Right:
                rect=self.rect()
                rect.setRight(self.rect().right()+factor)
                self.setRect(rect)
        else:
            if event.key()==QtCore.Qt.Key_Delete or event.key()==QtCore.Qt.Key_Backspace:
                self.delete()
            elif event.key()==QtCore.Qt.Key_Escape:
                self.setSelected(False)
            elif event.key()==QtCore.Qt.Key_Return or event.key()==QtCore.Qt.Key_Enter:
                self.setCurrentCell(self.currentRow()-1, self.currentColumn())
            elif event.key()==QtCore.Qt.Key_Up:
                rect=self.rect()
                rect.moveTop(self.rect().y()-factor)
                self.setRect(rect)
            elif event.key()==QtCore.Qt.Key_Down:
                rect=self.rect()
                rect.moveTop(self.rect().y()+factor)
                self.setRect(rect)
            elif event.key()==QtCore.Qt.Key_Left:
                rect=self.rect()
                rect.moveLeft(self.rect().x()-factor)
                self.setRect(rect)
            elif event.key()==QtCore.Qt.Key_Right:
                rect=self.rect()
                rect.moveLeft(self.rect().x()+factor)
                self.setRect(rect)


class RectItem(GeometricItem, QtGui.QGraphicsRectItem):
    """Clase que define un rectangulo"""
    type="square"
    icon=os.environ["pychemqt"]+"/images/equipment/square.png"


class EllipseItem(GeometricItem, QtGui.QGraphicsEllipseItem):
    """Clase que define una circunferencia"""
    type="ellipse"
    icon=os.environ["pychemqt"]+"/images/equipment/circle.png"

class TextItem(QtGui.QGraphicsTextItem):
    """Clase que define los textos en el diagrama de flujo"""
    type="txt"

    def __init__(self, text, parent=None, position=QtCore.QPointF(0, 0), transform=QtGui.QTransform(), selectable=True):
        super(TextItem, self).__init__(parent=parent)
        if selectable:
            self.setFlags(QtGui.QGraphicsItem.ItemIsSelectable|QtGui.QGraphicsItem.ItemIsMovable|QtGui.QGraphicsItem.ItemSendsGeometryChanges|QtGui.QGraphicsItem.ItemIsFocusable)
        else:
            self.setFlags(QtGui.QGraphicsItem.ItemIsMovable)
        self.setHtml(text)
        self.setPos(position)
        self.setTransform(transform)
        self.selectable=selectable

    def delete(self):
        self.scene().delete(self)

    def mouseDoubleClickEvent(self, event=None):
        dialog = TextItemDlg(self.toHtml())
        if dialog.exec_():
            self.setHtml(dialog.editor.texto)
            self.itemChange(QtGui.QGraphicsItem.ItemPositionChange, 0)

    def contextMenu(self):
        if self.selectable:
            contextMenu= QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Text Item: %s" %self.toPlainText()), self.scene().parent())
            contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equipment/text.png")))
            contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/editDelete.png"), QtGui.QApplication.translate("pychemqt", "Delete"), self.delete)
            contextMenu.addSeparator()
            contextMenu.addAction("Edit", self.mouseDoubleClickEvent)
            return contextMenu

    def itemChange(self, change, variant):
        if self.scene():
            if change == QtGui.QGraphicsItem.ItemPositionChange:
                if self.scene().parent().dirty[self.scene().parent().idTab]==False:
                    self.scene().parent().dirty[self.scene().parent().idTab]=True
                    self.scene().parent().activeControl(True)
                    self.scene().parent().tabModified(self.scene().parent().idTab)
        return QtGui.QGraphicsItem.itemChange(self, change, variant)

    def keyPressEvent(self, event):
        if event.key()==QtCore.Qt.Key_Delete or event.key()==QtCore.Qt.Key_Backspace:
            self.delete()
        elif event.key()==QtCore.Qt.Key_Escape:
            self.setSelected(False)
        elif event.key()==QtCore.Qt.Key_Return or event.key()==QtCore.Qt.Key_Enter:
            self.mouseDoubleClickEvent()
        elif event.key()==QtCore.Qt.Key_Up:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()-factor))
        elif event.key()==QtCore.Qt.Key_Down:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()+factor))
        elif event.key()==QtCore.Qt.Key_Left:
            self.setPos(QtCore.QPointF(self.pos().x()-factor, self.pos().y()))
        elif event.key()==QtCore.Qt.Key_Right:
            self.setPos(QtCore.QPointF(self.pos().x()+factor, self.pos().y()))


class GraphicsEntity(object):
    """Clase que modela la funcionalidad comun a corrientes y equipos en el PFD"""

    def view(self):
        with tempfile.NamedTemporaryFile("w+r", delete=False, suffix=".txt") as temp:
            temp.write(QtGui.QApplication.translate("pychemqt","Project Name")+": "+self.scene().parent().currentFilename+os.linesep)
            if isinstance(self.entity, Corriente):
                temp.write(QtGui.QApplication.translate("pychemqt","Stream Id"))
            else:
                temp.write(QtGui.QApplication.translate("pychemqt","Equipment Id"))
            temp.write(": %i" %self.id+os.linesep)
            ahora=datetime.today()
            temp.write(QtGui.QApplication.translate("pychemqt","Report generated at")+ " %s - %s" %(ahora.strftime("%H:%M:%S"), ahora.strftime("%d/%m/%Y"))+os.linesep)
            temp.write(self.entity.txt())
            subprocess.Popen([Preferences.get("Applications", 'TextViewer'), temp.name])

    def exportExcel(self):
        msg=QtGui.QApplication.translate("pychemqt", "Select Spreadsheet")
        patrones=QtCore.QStringList()
        if os.environ["ezodf"]:
            patrones.append(QtGui.QApplication.translate("pychemqt", "Libreoffice spreadsheet files")+ " (*.ods)")
        if os.environ["openpyxl"]:
            patrones.append(QtGui.QApplication.translate("pychemqt", "Microsoft Excel 2007/2010 XML")+ " (*.xlsx)")
        if os.environ["xlwt"]:
            patrones.append(QtGui.QApplication.translate("pychemqt", "Microsoft Excel 97/2000/XP/2003 XMLL")+ " (*.xls)")
        patron=patrones.join(";;")
        dir=os.path.dirname(str(self.scene().parent().currentFilename))
        ruta = unicode(QtGui.QFileDialog.getSaveFileName(self.scene().parent(), msg, dir, patron))
        if ruta:
            name, ext=os.path.splitext(ruta)
            if not ext or ext not in (".ods", ".xlsx", ".xls"):
                ruta+="."+str(patrones[0]).split(".")[-1][:-1]

            if ruta[-3:]=="ods":
                import ezodf
                templatefile=os.environ["pychemqt"]+os.sep+"dat"+os.sep+"templates"+os.sep+self.entity.__class__.__name__.lower()+".ots"
                if os.path.isfile(templatefile):
                    spreadsheet = ezodf.newdoc("ods", ruta, templatefile)
                    sheet = spreadsheet.sheets[0]
                    for attr, type, cell in self.entity.datamap2xls():
                        prop=self.entity._prop(attr)
                        if type == "value":
                            value=prop.config()
                        else:
                            value=prop.text()
                        sheet[cell].set_value(value)
                else:
                    spreadsheet = ezodf.newdoc("ods", ruta)
                    sheets = spreadsheet.sheets
                    sheet=ezodf.Table('pychemqt - s%i'%self.id)
                    sheets+=sheet
                    propiedades=self.entity.properties()
                    sheet.reset(size=(len(propiedades)+1, 10))

                    for i, (name, attr, unit) in enumerate(self.entity.propertiesNames()):
                        value=self.entity._prop(attr)
                        txt=""
                        if unit in unidades._all and not isinstance(value, list):
                            txt=value.text()
                            value=value.config()
                        elif isinstance(value, list):
                            txt=value[0].text()

                        sheet["B%i"%(i+2)].set_value(name)
                        sheet["C%i"%(i+2)].set_value(value)
                        sheet["D%i"%(i+2)].set_value(txt)

                spreadsheet.save()
            elif ruta[-4:]=="xlsx":
                print ruta, "is xlsx"
            elif ruta[-3:]=="xls":
                print ruta, "is xls"

class StreamItem(GeometricItem, QtGui.QGraphicsPathItem, GraphicsEntity):
    """Clase que define una corriente gráficamente"""
    up=None
    down=None
    id=0
    free_id=[]
    type="stream"

    def __init__(self, parent=None):
        super(StreamItem, self).__init__()
        self.parent=parent
        self.setPen(self._pen())
        qp=QtGui.QPainterPath()
        self.setPath(qp)
        self.setFlags(QtGui.QGraphicsItem.ItemIsSelectable|QtGui.QGraphicsItem.ItemIsFocusable)
        if StreamItem.free_id:
            self.id=StreamItem.free_id.pop(0)
        else:
            self.id=StreamItem.id+1
            StreamItem.id+=1
        self.idLabel=TextItem("S%i" %self.id, self, selectable=False)
        self.idLabel.setZValue(2)
        self.setAcceptHoverEvents(True)

    @property
    def corriente(self):
        return self.scene().project.getStream(self.id)
    @property
    def entity(self):
        return self.corriente

    def setCorriente(self, corriente):
            self.scene().project.setStream(self.id, corriente)
#            kwargs={"entrada": corriente}
#            if isinstance(self.scene().project.getDownToStream(self.id), flux.Mixer):
#                kwargs["id_entrada"]=self.scene().project.streams[self.id][3]+1
#            self.scene().project.getDownToStream(self.id)(**kwargs)
            pen=self.pen()
            if corriente.status==1:
                pen.setColor(QtGui.QColor("blue"))
            else:
                pen.setColor(QtGui.QColor("red"))
            self.setPen(pen)
            self.itemChange(QtGui.QGraphicsItem.ItemPositionChange, 0)

    def mouseDoubleClickEvent(self, event=None):
        dialog = UI_corriente.Corriente_Dialog(self.corriente)
        if dialog.exec_():
            self.setCorriente(dialog.corriente)

    def copyFromProject(self):
        dialog=SelectStreamProject()
        if dialog.exec_():
            indice=int(dialog.stream.currentText())
            corriente=dialog.project.getStream(indice)
            self.setCorriente(corriente)


    def keyPressEvent(self, event):
        if event.key()==QtCore.Qt.Key_Delete or event.key()==QtCore.Qt.Key_Backspace:
            self.delete()
        elif event.key()==QtCore.Qt.Key_Escape:
            self.setSelected(False)
        elif event.key()==QtCore.Qt.Key_Return or event.key()==QtCore.Qt.Key_Enter:
            self.mouseDoubleClickEvent()

    def hoverEnterEvent(self, event):
        if not (self.scene().addObj and self.scene().addType=="stream"):
            self.tabla=Table_Graphics(self.corriente, self.id, self.scene().parent().Preferences)
            self.tabla.move(event.screenPos())
            self.tabla.show()

    def hoverLeaveEvent(self, event):
        if not (self.scene().addObj and self.scene().addType=="stream"):
            self.tabla.hide()
            self.tabla.deleteLater()

#    def hoverMoveEvent(self, event):
#        self.tabla.move(event.screenPos())
#        self.tabla.show()

    def contextMenu(self):
        ViewAction=createAction(QtGui.QApplication.translate("pychemqt", "View Properties"), slot=self.view, parent=self.scene())
        SolidDistributionAction=createAction(QtGui.QApplication.translate("pychemqt", "Solid Distribution Fit"), slot=self.solidFit, parent=self.scene())
        if self.corriente:
            if not self.corriente.solido:
                SolidDistributionAction.setEnabled(False)
        else:
            ViewAction.setEnabled(False)
            SolidDistributionAction.setEnabled(False)

        contextMenu= QtGui.QMenu("Stream %i" %self.id, self.scene().parent())
        contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equipment/stream.png")))
        contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Show/Hide Id Label"), self.idLabelVisibility)
        contextMenu.addSeparator()
        contextMenu.addAction(ViewAction)
        contextMenu.addAction(SolidDistributionAction)
        contextMenu.addSeparator()
        contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Copy from another project"), self.copyFromProject)
        contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Edit"), self.mouseDoubleClickEvent)
        contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Export to spreadsheet"), self.exportExcel)
        contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Appearance"), self.format)
        contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/editDelete.png"), QtGui.QApplication.translate("pychemqt", "Delete"), self.delete)
        return contextMenu


    def redraw(self, entrada=None, salida=None):
        if entrada:
            self.entrada=entrada
        if salida:
            self.salida=salida

        max_height=max(self.up.boundingRect().height(), self.down.boundingRect().height())
        max_width=max(self.up.boundingRect().width(), self.down.boundingRect().width())
        y_up=min(self.up.pos().y(), self.down.pos().y())
        y_down=max(self.up.pos().y(), self.down.pos().y())
        if self.up.pos().y()==y_up:
            height_sup=self.up.boundingRect().height()
        else:
            height_sup=self.down.boundingRect().height()
        x_mean=(self.entrada.x()+self.salida.x())/2.
        Xdist_entrada=0
        Ydist_entrada=0
        Xdist_salida=0
        Ydist_salida=0
        if self.Ang_entrada==0:
            Xdist_entrada=20
        elif self.Ang_entrada==180:
            Xdist_entrada=-20
        elif self.Ang_entrada==90:
            Ydist_entrada=20
        elif self.Ang_entrada==360:
            Ydist_entrada=-20

        if self.Ang_salida==0:
            Xdist_salida=-20
        elif self.Ang_salida==180:
            Xdist_salida=20

        qp=QtGui.QPainterPath()
        qp.moveTo(self.entrada)
        if self.Ang_entrada==self.Ang_salida:
            if self.salida.x() > self.entrada.x()+10:
                qp.lineTo(QtCore.QPointF(x_mean, self.entrada.y()))
                qp.lineTo(QtCore.QPointF(x_mean, self.salida.y()))
            else:
                if abs(self.entrada.y()-self.salida.y())>max_height: #cabe la linea entre ambos equipos
                    y_mean=(y_up+y_down+height_sup)/2.
                    qp.lineTo(QtCore.QPointF(self.entrada.x()+Xdist_entrada, self.entrada.y()+Ydist_entrada))
                    qp.lineTo(QtCore.QPointF(self.entrada.x()+Xdist_entrada, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x()+Xdist_salida, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x()+Xdist_salida, self.salida.y()+Ydist_entrada))
                else:   #sacamos la linea por encima de los equipos
                    y_mean=y_up-20
                    qp.lineTo(QtCore.QPointF(self.entrada.x()+Xdist_entrada, self.entrada.y()+Ydist_entrada))
                    qp.lineTo(QtCore.QPointF(self.entrada.x()+Xdist_entrada, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x()+Xdist_salida, y_mean))
                    qp.lineTo(QtCore.QPointF(self.salida.x()+Xdist_salida, self.salida.y()+Ydist_salida))
        else:
            x_mean=max(self.entrada.x(), self.salida.x())+Xdist_salida
            qp.lineTo(QtCore.QPointF(x_mean, self.entrada.y()))
            qp.lineTo(QtCore.QPointF(x_mean, self.salida.y()))

        qp.lineTo(self.salida)
        self.prepareGeometryChange()
        self.setPath(qp)
        if abs(self.entrada.y()-self.salida.y())<=30:
            self.idLabel.setPos(x_mean, max(self.entrada.y(), self.salida.y()))
        else:
            self.idLabel.setPos(x_mean, (self.entrada.y()+self.salida.y())/2.-10)

    def postDelete(self):
        StreamItem.free_id.append(self.id)
        self.up.down_used-=1
        self.down.up_used-=1
        self.up.down.remove(self)
        self.down.up.remove(self)

    def idLabelVisibility(self):
        self.idLabel.setVisible(not self.idLabel.isVisible())


    def solidFit(self):
        if self.corriente.solido:
            dialog=Plot_Distribucion(self.id, self.corriente.solido)
            self.scene().parent().currentMdi.addSubWindow(dialog)
            dialog.show()



class EquipmentItem(QtSvg.QGraphicsSvgItem, GraphicsEntity):
    """Clase que define los equipos en el diagrama de flujo"""
    up=[]
    down=[]
    up_used=0
    down_used=0
    id=0
    id_in=0
    id_out=0
    type="equip"
    def __init__(self, name, dialogoId, parent=None):
        self.name=name
        imagen=os.environ["pychemqt"]+"/images/equipment/%s.svg" % name
        super(EquipmentItem, self).__init__(imagen, parent=parent)
        self.dialogoId=dialogoId
        self.setFlags(QtGui.QGraphicsItem.ItemIsSelectable|QtGui.QGraphicsItem.ItemIsMovable|QtGui.QGraphicsItem.ItemSendsGeometryChanges|QtGui.QGraphicsItem.ItemIsFocusable)
        self.imagen=imagen
        self.angle=0
        self.setAcceptHoverEvents(True)

        if dialogoId!=None:
            self.dialogo=UI_equipments[dialogoId].UI_equipment
            EquipmentItem.id+=1
            self.id=EquipmentItem.id
            self.tipo="e"
            self.idLabel=TextItem("E%i" %self.id, self, selectable=False)
            self.idLabel.setPos(self.boundingRect().width()/3.,-20)

        else:
            self.dialogo=UI_corriente.Corriente_Dialog
            if name=="in":
                EquipmentItem.id_in+=1
                self.id=EquipmentItem.id_in
                self.tipo="i"
            else:
                EquipmentItem.id_out+=1
                self.id=EquipmentItem.id_out
                self.tipo="o"


        output=[]
        input=[]
        dom = QtXml.QDomDocument()
        fh = QtCore.QFile(imagen)
        if not fh.open(QtCore.QIODevice.ReadOnly):
            raise IOError, unicode(fh.errorString())
        if not dom.setContent(fh):
            raise ValueError, "could not parse XML"
        fh.close()
        root = dom.documentElement()
        node = root.firstChild().toElement()

        while node.toElement().tagName()!="ins":
            node=node.nextSibling()
        child=node.toElement().firstChild().toElement()
        while not child.isNull():
            x=child.attribute("x").toFloat()[0]
            y=child.attribute("y").toFloat()[0]
            d=child.attribute("d")
            input.append([x, y, d])
            child = child.nextSibling().toElement()

        while node.toElement().tagName()!="outs":
            node=node.nextSibling()
        child=node.toElement().firstChild().toElement()
        while not child.isNull():
            x=child.attribute("x").toFloat()[0]
            y=child.attribute("y").toFloat()[0]
            d=child.attribute("d")
            output.append([x, y, d])
            child = child.nextSibling().toElement()

        self.input=[]
        if input:
            for entrada in input:
                obj=QtGui.QGraphicsEllipseItem(self)
                obj.setRect(entrada[0]*self.boundingRect().width()-5, entrada[1]*self.boundingRect().height()-5, 10, 10)
                obj.direction=int(entrada[2])
                obj.setPen(QtGui.QColor(255, 255, 255))
                obj.setBrush(QtGui.QColor(Preferences.get("PFD", 'Color_Entrada')))
                self.input.append(obj)
        self.output=[]
        if output:
            for salida in output:
                obj=QtGui.QGraphicsEllipseItem(self)
                obj.setRect(salida[0]*self.boundingRect().width()-5, salida[1]*self.boundingRect().height()-5, 10, 10)
                obj.direction=int(salida[2])
                obj.setPen(QtGui.QColor(255, 255, 255))
                obj.setBrush(QtGui.QColor(Preferences.get("PFD", 'Color_Salida')))
                self.output.append(obj)
        self.showInput(False)

    @property
    def equipment(self):
        return self.scene().project.getItem(self.id)
    @property
    def entity(self):
        return self.equipment

    def mouseDoubleClickEvent(self, event=None):
        if self.dialogoId!=None:
            args=[self.equipment]
            #Añadimos argumentos adicionales en algunos equipos
            if isinstance(self.equipment, flux.Divider):         #divisor
                if not len(self.down):
                    return
                args.append(len(self.down))
            elif isinstance(self.equipment, flux.Mixer):         #mezclador
                if not len(self.up):
                    return
                args.append(len(self.up))
            elif isinstance(self.equipment, spreadsheet.Spreadsheet):         #spreadsheet
                self.equipment(project=self.scene().project)
                args.append(self.scene().project)

            dialog = self.dialogo(*args)
            if dialog.exec_():
                self.scene().project.setItem(self.id, dialog.Equipment)
#                self.up[0].setCorriente(dialog.Equipment.entrada)
                for i, corriente in enumerate(dialog.Equipment.salida):
                    self.down[i].setCorriente(corriente)

                self.itemChange(QtGui.QGraphicsItem.ItemPositionChange, 0)
        else:
            if self.output:
                self.down[0].mouseDoubleClickEvent()
            else:
                self.up[0].mouseDoubleClickEvent()

    def mousePressEvent(self, event):
        QtSvg.QGraphicsSvgItem.mousePressEvent(self, event)
        if self.scene().addObj:
            if self.scene().addType=="stream":
                if len(self.scene().Pos)==0:
                    punto=self.output
                    self.scene().up=self
                    x=self.down_used
                else:
                    punto=self.input
                    self.scene().down=self
                    x=self.up_used
                self.scene().Pos.append(self.mapToScene(punto[x].rect().center()))
                self.scene().points.append(punto[x])
            else:
                self.scene().Pos.append(event.pos())

    def mouseMoveEvent(self, event=None):
        if event:
            QtGui.QGraphicsPixmapItem.mouseMoveEvent(self, event)
        for i, corriente in enumerate(self.up):
            corriente.redraw(salida=self.mapToScene(self.input[i].rect().center()))
        for i, corriente in enumerate(self.down):
            corriente.redraw(entrada=self.mapToScene(self.output[i].rect().center()))

    def showInput(self, bool):
        for entrada in self.input:
            entrada.setVisible(bool)
        for salida in self.output:
            salida.setVisible(bool)

    def hoverEnterEvent(self, event):
        if self.scene().addObj and self.scene().addType=="stream":
            self.showInput(True)
        else:
            if self.dialogoId!=None:
                self.tabla=Table_Graphics(self.equipment, self.id, self.scene().parent().Preferences)
            else:
                if self.output:
                    self.tabla=Table_Graphics(self.down[0].corriente, self.down[0].id, self.scene().parent().Preferences)
                else:
                    self.tabla=Table_Graphics(self.up[0].corriente, self.up[0].id, self.scene().parent().Preferences)
            self.tabla.move(event.screenPos())
            self.tabla.show()

    def hoverLeaveEvent(self, event):
        self.showInput(False)
        if not self.scene().addObj and self.scene().addType=="stream":
            self.tabla.hide()
            self.tabla.deleteLater()

#    def hoverMoveEvent(self, event):
#        self.tabla.move(event.screenPos())
#        self.tabla.show()
#

    def keyPressEvent(self, event):
        if event.key()==QtCore.Qt.Key_Delete or event.key()==QtCore.Qt.Key_Backspace:
            self.delete()
        elif event.key()==QtCore.Qt.Key_Escape:
            self.setSelected(False)
        elif event.key()==QtCore.Qt.Key_Return or event.key()==QtCore.Qt.Key_Enter:
            self.mouseDoubleClickEvent()
        elif event.key()==QtCore.Qt.Key_Up:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()-factor))
            self.mouseMoveEvent()
        elif event.key()==QtCore.Qt.Key_Down:
            self.setPos(QtCore.QPointF(self.pos().x(), self.pos().y()+factor))
            self.mouseMoveEvent()
        elif event.key()==QtCore.Qt.Key_Left:
            self.setPos(QtCore.QPointF(self.pos().x()-factor, self.pos().y()))
            self.mouseMoveEvent()
        elif event.key()==QtCore.Qt.Key_Right:
            self.setPos(QtCore.QPointF(self.pos().x()+factor, self.pos().y()))
            self.mouseMoveEvent()

    def itemChange(self, change, variant):
        if self.scene():
            if change == QtGui.QGraphicsItem.ItemPositionChange:
                if self.scene().parent().dirty[self.scene().parent().idTab]==False:
                    self.scene().parent().dirty[self.scene().parent().idTab]=True
                    self.scene().parent().activeControl(True)
                    self.scene().parent().tabModified(self.scene().parent().idTab)
        return QtGui.QGraphicsItem.itemChange(self, change, variant)


    def contextMenu(self):
        if self.dialogoId!=None:
            ViewAction=createAction(QtGui.QApplication.translate("pychemqt", "View Properties"), slot=self.view, parent=self.scene())
            ViewAction.setEnabled(self.equipment.status)

            contextMenu= QtGui.QMenu("Equipment %i" %self.id, self.scene().parent())
            contextMenu.setIcon(QtGui.QIcon(QtGui.QPixmap(os.environ["pychemqt"]+"/images/equipment/%s.svg" % self.name)))
            contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Show/Hide Id Label"), self.idLabelVisibility)
            contextMenu.addSeparator()
            contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Edit"), self.mouseDoubleClickEvent)
            contextMenu.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/editDelete.png"), QtGui.QApplication.translate("pychemqt", "Delete"), self.delete)
            contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Export to spreadsheet"), self.exportExcel)
            contextMenu.addAction(QtGui.QApplication.translate("pychemqt", "Appearance"), self.format)
            contextMenu.addSeparator()
            contextMenu.addAction(ViewAction)
            contextMenu.addSeparator()
#            contextMenu.addAction("Run", self.mouseDoubleClickEvent)
        else:
            if self.output:
                contextMenu=self.down[0].contextMenu()
            else:
                contextMenu=self.up[0].contextMenu()

        contextMenu.addSeparator()
        self.menuTransform=QtGui.QMenu(QtGui.QApplication.translate("pychemqt", "Transform"))
        self.menuTransform.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/transform_rotate_90.png"), QtGui.QApplication.translate("pychemqt", "Rotate by 90"), partial(self.rotate, 90))
        self.menuTransform.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/transform_rotate_180.png"), QtGui.QApplication.translate("pychemqt", "Rotate by 180"), partial(self.rotate, 180))
        self.menuTransform.addAction(QtGui.QIcon(os.environ["pychemqt"]+"/images/button/transform_rotate_270.png"), QtGui.QApplication.translate("pychemqt", "Rotate by 270"), partial(self.rotate, 270))
        self.menuTransform.addSeparator()
        self.menuTransform.addAction(QtGui.QApplication.translate("pychemqt", "Mirror about X"), partial(self.rotate, 270))
        self.menuTransform.addAction(QtGui.QApplication.translate("pychemqt", "Mirror about Y"), partial(self.rotate, 270))
        contextMenu.addAction(self.menuTransform.menuAction())

        return contextMenu


    def delete(self):
        self.scene().delete(self)

    def format(self):
        pass

    def postDelete(self):
        while self.down:
            stream=self.down.pop()
            self.scene().delete(stream)
        while self.up:
            stream=self.up.pop()
            self.scene().delete(stream)

    def idLabelVisibility(self):
        self.idLabel.setVisible(not self.idLabel.isVisible())

    def rotate(self, angle):
        self.angle=angle
        transform=self.transform()
        transform.rotate(angle)
        self.setTransform(transform)
        self.mouseMoveEvent()
        for i, entrada in enumerate(self.up):
            new_angle=(self.input[i].direction+angle) % 360
            self.input[i].direction=new_angle
            entrada.Ang_salida=new_angle
            entrada.redraw()
        for i, salida in enumerate(self.down):
            new_angle=(self.output[i].direction+angle) % 360
            self.output[i].direction=new_angle
            salida.Ang_entrada=new_angle
            salida.redraw()



class GraphicsView(QtGui.QGraphicsView):
    mouseMove = QtCore.pyqtSignal(QtCore.QPointF)
#    mouseClick = QtCore.pyqtSignal("QtCore.QPointF")
    def __init__(self, PFD=True, parent=None):
        super(GraphicsView, self).__init__(parent)
        self.setDragMode(QtGui.QGraphicsView.RubberBandDrag)
        self.setRenderHint(QtGui.QPainter.Antialiasing)
        self.setRenderHint(QtGui.QPainter.TextAntialiasing)
        self.setBackgroundBrush(QtGui.QBrush(QtGui.QColor("#aaaaaa"),QtCore.Qt.Dense7Pattern))
        self.setMouseTracking(True)
        self.setAcceptDrops(True)
        self.PFD=PFD

#    def wheelEvent(self, event):
#        print event.delta()
#
#        factor = 1.41 ** (-event.delta() / 240.0)
#        self.zoom(factor)

    def mouseMoveEvent(self, event):
        QtGui.QGraphicsView.mouseMoveEvent(self, event)
        self.mouseMove.emit(event.posF())

    def mousePressEvent(self, event):
        QtGui.QGraphicsView.mousePressEvent(self, event)
        if not self.PFD:
            self.scene().views()[0].centerOn(self.mapToScene(event.pos()))

    def closeEvent(self, event):
        if self.PFD:
            event.ignore()

#    def dragEnterEvent(self, event):
#        if event.mimeData().hasFormat("application/x-equipment"):
#            event.accept()
#        else:
#            event.ignore()
#
#    def dragMoveEvent(self, event):
#        if event.mimeData().hasFormat("application/x-equipment"):
#            event.setDropAction(QtCore.Qt.CopyAction)
#            event.accept()
#        else:
#            event.ignore()
#
#    def dropEvent(self, event):
#        if event.mimeData().hasFormat("application/x-equipment"):
#            data = event.mimeData().data("application/x-equipment")
#            stream = QtCore.QDataStream(data, QtCore.QIODevice.ReadOnly)
#            icon = QtGui.QIcon()
#            stream >> icon
#            event.setDropAction(QtCore.Qt.CopyAction)
#            print event.pos()
#            event.accept()
#            self.updateGeometry()
#            self.update()
#        else:
#            event.ignore()

    def zoom(self, value):
        factor = value / 100.0
        self.resetMatrix()
        self.scale(factor, factor)




class GraphicsScene(QtGui.QGraphicsScene):
    copiedItem =QtCore.QByteArray()
    pasteOffset = 5
    points=[]
    addObj=False
    addType=""
    project=Project()
    objects={"txt": [], "square": [], "ellipse": [], "stream": {}, "in": {}, "out": {}, "equip": {}}

    def __init__(self, parent=None):
        super(GraphicsScene, self).__init__(parent)

    def mousePressEvent(self, event):
        QtGui.QGraphicsScene.mousePressEvent(self, event)
        if self.addObj and self.addType!="stream":
            self.Pos.append(event.scenePos())

    def addActions(self, menu, pos=None):
        menu.addAction(QtGui.QApplication.translate("pychemqt", "Redraw"), self.update)
        menu.addSeparator()
        menu.addAction(createAction(QtGui.QApplication.translate("pychemqt", "Select All"), slot=self.selectAll, shortcut=QtGui.QKeySequence.SelectAll, icon=os.environ["pychemqt"]+"/images/button/selectAll", parent=self))
        menu.addSeparator()
        actionCut=createAction(QtGui.QApplication.translate("pychemqt", "Cut"), slot=self.cut, shortcut=QtGui.QKeySequence.Cut, icon=os.environ["pychemqt"]+"/images/button/editCut", parent=self)
        menu.addAction(actionCut)
        actionCopy=createAction(QtGui.QApplication.translate("pychemqt", "Copy"), slot=self.copy, shortcut=QtGui.QKeySequence.Copy, icon=os.environ["pychemqt"]+"/images/button/editCopy", parent=self)
        menu.addAction(actionCopy)
        actionPaste=createAction(QtGui.QApplication.translate("pychemqt", "Paste"), slot=partial(self.paste, pos), shortcut=QtGui.QKeySequence.Paste, icon=os.environ["pychemqt"]+"/images/button/editPaste", parent=self)
        menu.addAction(actionPaste)
        actionDelete=createAction(QtGui.QApplication.translate("pychemqt", "Delete All"), slot=self.delete, shortcut=QtGui.QKeySequence.Delete, icon=os.environ["pychemqt"]+"/images/button/editDelete", parent=self)
        menu.addAction(actionDelete)
        menu.addSeparator()

        if self.copiedItem.isEmpty():
            actionPaste.setEnabled(False)
        items=self.selectedItems()
        if not items:
            actionCut.setEnabled(False)
            actionCopy.setEnabled(False)
            actionDelete.setEnabled(False)

        for item in items:
            menuEl=item.contextMenu()
            menu.addAction(menuEl.menuAction())
        return menu

    def contextMenuEvent(self, event):
        item=self.itemAt(event.scenePos())
        if item:
            item.setSelected(True)
        contextMenu= QtGui.QMenu()
        self.addActions(contextMenu, event.scenePos())
        contextMenu.exec_(event.screenPos())

    def selectAll(self):
        for item in self.items():
            item.setSelected(True)

    def copy(self, item=None):
        if not item:
            item=self.selectedItems()[0]
        self.copiedItem.clear()
        self.pasteOffset = 5
        stream = QtCore.QDataStream(self.copiedItem, QtCore.QIODevice.WriteOnly)
        self.writeItemToStream(stream, item)

    def cut(self):
        item=self.selectedItems()[0]
        self.copy(item)
        self.removeItem(item)
        del item

    def paste(self, pos=None):
        stream = QtCore.QDataStream(self.copiedItem, QtCore.QIODevice.ReadOnly)
        item=self.readItemFromStream(stream)
        if pos:
            item.setPos(pos)
        else:
            item.setPos(item.pos()+QtCore.QPointF(self.pasteOffset, self.pasteOffset))
            self.pasteOffset += 5
        self.addItem(item)

    def delete(self, items=None):
        if items:
            items=[items]
        else:
            items=self.selectedItems()
        for item in items:
            tipo=item.type
            if tipo in ["stream", "equip"]:
                item.postDelete()
                del self.objects[tipo][item.id]
            else:
                self.objects[tipo].remove(item)
            self.removeItem(item)
        self.update()
        self.parent().list.updateList(self.objects)


    def waitClick(self, numClick, type, object):
        self.object=object
        self.addType=type
        self.addObj=True
        self.views()[0].viewport().setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
        self.parent().statusbar.showMessage(QtGui.QApplication.translate("pychemqt", "Click in desire text position in screen"))
        self.Pos=[]
        self.clickCollector = WaitforClick(numClick, self)
        self.clickCollector.finished.connect(self.click)
        self.clickCollector.start()


    def click(self):
        if self.addType in ["equip", "in",  "out", "txt"]:
            self.object.setPos(self.Pos[0])
        elif self.addType=="stream":
            self.object.up=self.up
            self.object.down=self.down
            self.up.down=self.up.down+[self.object]
            self.down.up=self.down.up+[self.object]
            self.up.down_used+=1
            self.down.up_used+=1
            self.object.Ang_entrada=self.points[0].direction
            self.object.Ang_salida=self.points[1].direction
            self.object.redraw(self.Pos[0], self.Pos[1])
        elif self.addType in ["square", "ellipse"]:
            rect=QtCore.QRectF(self.Pos[0], self.Pos[1])
            self.object.setRect(rect)

        self.addItem(self.object)
        if self.addType=="equip":
            self.project.addItem("e%i" %(self.object.id), self.object.dialogo.Equipment)
        elif self.addType=="in":
            self.project.addItem("i%i" %(self.object.id), self.object.dialogo.corriente)
        elif self.addType=="out":
            self.project.addItem("o%i" %(self.object.id), self.object.dialogo.corriente)
        elif self.addType=="stream":
            self.project.addStream(self.object.id, "%s%i" %(self.up.tipo, self.up.id), "%s%i" %(self.down.tipo, self.down.id), Corriente(), self.up.down_used-1, self.down.up_used-1)

        self.parent().dirty[self.parent().idTab]=True
        self.parent().saveControl()
        self.update()
        self.object.setSelected(True)

        self.parent().statusbar.clearMessage()
        self.addObj=False
        if self.addType in ("txt", "square", "ellipse"):
            self.objects[self.addType].append(self.object)
        else:
            id=self.object.id
            self.objects[self.addType][id]=self.object
        self.parent().list.updateList(self.objects)
#        for item in self.items():
#            if isinstance(item, QtSvg.QGraphicsSvgItem):
#                item.hoverLeaveEvent(None)
#            item.setAcceptHoverEvents(False)




    def readItemFromStream(self, stream):
        type = QtCore.QString()
        matrix = QtGui.QTransform()
        stream >> type >> matrix
        if type == "txt":
            text = QtCore.QString()
            stream >> text
            item=TextItem(text)
        elif type == "square":
            rect=QtCore.QRectF()
            pen=QtGui.QPen()
            stream >> rect >> pen
            item=RectItem()
            item.setRect(rect)
            item.setPen(pen)
        elif type == "ellipse":
            rect=QtCore.QRectF()
            pen=QtGui.QPen()
            stream >> rect >> pen
            item=EllipseItem()
            item.setRect(rect)
            item.setPen(pen)
        elif type == "equip":
            name = QtCore.QString()
            stream >> name
            dialogoid=stream.readInt32()
            item=EquipmentItem(name, dialogoid)

        item.setTransform(matrix)
        return item


    def writeItemToStream(self, stream, item):
        stream << QtCore.QString(item.type) << item.transform()
        if isinstance(item, TextItem):
            stream << item.toHtml()
        elif isinstance(item, EllipseItem):
            stream  << item.rect()  << item.pen()
        elif isinstance(item, RectItem):
            stream << item.rect()  << item.pen()
        elif isinstance(item, EquipmentItem):
            stream << QtCore.QString(item.name)
            stream.writeInt32(item.dialogoId)


    def readFromFile(self, stream):
        n_txt = stream.readInt32()
        for obj in range(n_txt):
            txt=QtCore.QString()
            pos=QtCore.QPointF()
            stream >> txt
            stream >> pos
            s=TextItem(txt)
            s.setPos(pos)
            self.objects["txt"].append(s)
            self.addItem(s)

        n_square=stream.readInt32()
        for obj in range(n_square):
            rect=QtCore.QRectF()
            pen=QtGui.QPen()
            stream >> rect
            stream >> pen
            s=RectItem()
            s.setRect(rect)
            s.setPen(pen)
            self.objects["square"].append(s)
            self.addItem(s)

        n_circle=stream.readInt32()
        for obj in range(n_circle):
            rect=QtCore.QRectF()
            pen=QtGui.QPen()
            stream >> rect
            stream >> pen
            s=EllipseItem()
            s.setRect(rect)
            s.setPen(pen)
            self.objects["ellipse"].append(s)
            self.addItem(s)

        n_stream=stream.readInt32()
        id_stream=[]
        up_stream={}
        down_stream={}
        for obj in range(n_stream):
            entrada=QtCore.QPointF()
            salida=QtCore.QPointF()
            pen=QtGui.QPen()
            stream >> entrada
            stream >> salida
            stream >> pen
            id_stream.append(stream.readInt32())
            up_type=QtCore.QString()
            down_type=QtCore.QString()
            stream >> up_type
            up_id=stream.readInt32()
            stream >> down_type
            down_id=stream.readInt32()
            up_stream[id_stream[-1]]=up_type, up_id
            down_stream[id_stream[-1]]=down_type, down_id
            s=StreamItem()
            s.id=id_stream[-1]
            s.entrada=entrada
            s.salida=salida
            s.setPen(pen)
            s.Ang_entrada=stream.readInt32()
            s.Ang_salida=stream.readInt32()
            self.objects["stream"][id_stream[-1]]=s
            self.addItem(s)
            txt=QtCore.QString()
            pos=QtCore.QPointF()
            stream >> txt
            stream >> pos
            s.idLabel.setPos(pos)
            s.idLabel.setHtml(txt)

        n_in=stream.readInt32()
        angle_in={}
        for obj in range(n_in):
            pos=QtCore.QPointF()
            stream >> pos
            id=stream.readInt32()
            angle_in[id]=stream.readInt32()
            n_up=stream.readInt32()
            up=[]
            for i in range(n_up):
                up.append(self.objects["stream"][stream.readInt32()])
            n_down=stream.readInt32()
            down=[]
            for i in range(n_down):
                down.append(self.objects["stream"][stream.readInt32()])
            s=EquipmentItem("in", None)
            s.setPos(pos)
            s.down=down
            s.up=up
            self.objects["in"][id]=s
            self.addItem(s)

        n_out=stream.readInt32()
        angle_out={}
        for obj in range(n_out):
            pos=QtCore.QPointF()
            stream >> pos
            id=stream.readInt32()
            angle_out[id]=stream.readInt32()
            n_up=stream.readInt32()
            up=[]
            for i in range(n_up):
                up.append(self.objects["stream"][stream.readInt32()])
            n_down=stream.readInt32()
            down=[]
            for i in range(n_down):
                down.append(self.objects["stream"][stream.readInt32()])
            s=EquipmentItem("out", None)
            s.setPos(pos)
            s.down=down
            s.up=up
            self.objects["out"][id]=s
            self.addItem(s)

        n_equip=stream.readInt32()
        angle_equip={}
        for obj in range(n_equip):
            name=QtCore.QString()
            pos=QtCore.QPointF()
            stream >> name
            dialogoId=stream.readInt32()
            stream >> pos
            id=stream.readInt32()
            angle_equip[id]=stream.readInt32()
            n_up=stream.readInt32()
            up=[]
            for i in range(n_up):
                up.append(self.objects["stream"][stream.readInt32()])
            n_down=stream.readInt32()
            down=[]
            for i in range(n_down):
                down.append(self.objects["stream"][stream.readInt32()])
            s=EquipmentItem(name, dialogoId)
            s.setPos(pos)
            s.down=down
            s.up=up
            self.objects["equip"][id]=s
            self.addItem(s)
            txt=QtCore.QString()
            pos=QtCore.QPointF()
            stream >> txt
            stream >> pos
            s.idLabel.setPos(pos)
            s.idLabel.setHtml(txt)

        for id in id_stream:
            tipo, i=up_stream[id]
            self.objects["stream"][id].up=self.getObject(tipo, i)
            tipo, i=down_stream[id]
            self.objects["stream"][id].down=self.getObject(tipo, i)
            self.objects["stream"][id].redraw()

        for id, angle in angle_in.iteritems():
            self.objects["in"][id].rotate(angle)
        for id, angle in angle_out.iteritems():
            self.objects["out"][id].rotate(angle)
        for id, angle in angle_equip.iteritems():
            self.objects["equip"][id].rotate(angle)


    def saveToFile(self, stream):
        stream.writeInt32(len(self.objects["txt"]))
        for obj in self.objects["txt"]:
            stream << obj.toHtml()
            stream << obj.pos()

        stream.writeInt32(len(self.objects["square"]))
        for obj in self.objects["square"]:
            stream << obj.rect()
            stream << obj.pen()

        stream.writeInt32(len(self.objects["ellipse"]))
        for obj in self.objects["ellipse"]:
            stream << obj.rect()
            stream << obj.pen()

        stream.writeInt32(len(self.objects["stream"]))
        for id in self.objects["stream"]:
            stream << self.objects["stream"][id].entrada
            stream << self.objects["stream"][id].salida
            stream << self.objects["stream"][id].pen()
            stream.writeInt32(self.objects["stream"][id].id)
            stream << QtCore.QString(self.objects["stream"][id].up.tipo)
            stream.writeInt32(self.objects["stream"][id].up.id)
            stream << QtCore.QString(self.objects["stream"][id].down.tipo)
            stream.writeInt32(self.objects["stream"][id].down.id)
            stream.writeInt32(self.objects["stream"][id].Ang_entrada)
            stream.writeInt32(self.objects["stream"][id].Ang_salida)
            stream << self.objects["stream"][id].idLabel.toHtml()
            stream << self.objects["stream"][id].idLabel.pos()

        stream.writeInt32(len(self.objects["in"]))
        for id in self.objects["in"]:
            stream << self.objects["in"][id].pos()
            stream.writeInt32(self.objects["in"][id].id)
            stream.writeInt32(self.objects["in"][id].angle)
            stream.writeInt32(len(self.objects["in"][id].up))
            for up in self.objects["in"][id].up:
                stream.writeInt32(up.id)
            stream.writeInt32(len(self.objects["in"][id].down))
            for down in self.objects["in"][id].down:
                stream.writeInt32(down.id)

        stream.writeInt32(len(self.objects["out"]))
        for id in self.objects["out"]:
            stream << self.objects["out"][id].pos()
            stream.writeInt32(self.objects["out"][id].id)
            stream.writeInt32(self.objects["out"][id].angle)
            stream.writeInt32(len(self.objects["out"][id].up))
            for up in self.objects["out"][id].up:
                stream.writeInt32(up.id)
            stream.writeInt32(len(self.objects["out"][id].down))
            for down in self.objects["out"][id].down:
                stream.writeInt32(down.id)

        stream.writeInt32(len(self.objects["equip"]))
        for id in self.objects["equip"]:
            stream << QtCore.QString(self.objects["equip"][id].name)
            stream.writeInt32(self.objects["equip"][id].dialogoId)
            stream << self.objects["equip"][id].pos()
            stream.writeInt32(self.objects["equip"][id].id)
            stream.writeInt32(self.objects["equip"][id].angle)
            stream.writeInt32(len(self.objects["equip"][id].up))
            for up in self.objects["equip"][id].up:
                stream.writeInt32(up.id)
            stream.writeInt32(len(self.objects["equip"][id].down))
            for down in self.objects["equip"][id].down:
                stream.writeInt32(down.id)
            stream << self.objects["equip"][id].idLabel.toHtml()
            stream << self.objects["equip"][id].idLabel.pos()


    def getObject(self, tipo, id):
        if tipo=="e":
            lista=self.objects["equip"]
        elif tipo=="i":
            lista=self.objects["in"]
        elif tipo=="o":
            lista=self.objects["out"]
        elif tipo=="s":
            lista=self.objects["stream"]
        else:
            raise Exception
        return lista[id]

if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
#    dialogo = ConfLineDialog()
    dialogo = SelectStreamProject()
    dialogo.show()
    sys.exit(app.exec_())

