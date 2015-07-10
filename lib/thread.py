#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Library for work with thread in pychemqt for improve UI response
#   - WaitforClick: Thread for draw stream in PFD
#   - Evaluate: Thread to insolate entity calculation from gui, used in streams,
#       equipment, and project
###############################################################################

from time import sleep

from PyQt5.QtCore import QThread, QMutex


class WaitforClick(QThread):
    """Thread used in PFD drawing to specified stream input and output or add
    equipment.
    TODO: Use to customize stream drawing"""
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
    """Thread used to insolate calculation process in entities (stream, project
    and equipment, so gui can response while calculation is in process"""
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
