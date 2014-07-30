#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Módulo que define thread secundarios usados en el proyecto"""

from time import sleep

from PyQt4.QtCore import QThread, QMutex


class WaitforClick(QThread):
    """Thread usado en la scena para esperar el click que indica donde añadir
    un nuevo elemento gráfico"""
    def __init__(self, num, parent=None):
        super(WaitforClick, self).__init__(parent)
        self.setTerminationEnabled(True)
        self.num = num

    def run(self):
        while True:
            sleep(0.1)
            if len(self.parent().Pos) >= self.num:
                break


class Evaluate(QThread):
    """Thread usado para aislar el proceso de calculo de entity (corrientes,
    equipos, proyectos) de la aplicación principal"""
    def __init__(self, parent=None):
        super(Evaluate, self).__init__(parent)
        self.mutex = QMutex()

    def start(self, entity, kwargs):
        self.entity = entity
        self.kwargs = kwargs
        QThread.start(self)

    def run(self):
        self.mutex.lock()
        self.entity(**self.kwargs)
        self.mutex.unlock()
