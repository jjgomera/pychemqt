#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Tools to create a python shell with pychemqt libraries imported
# For now only work in linux with xterm as terminal
###############################################################################

import atexit

from PyQt4 import QtCore, QtGui

from lib.firstrun import which


class XTerm(QtGui.QX11EmbedContainer):
    """Gui container for terminal widget, QX11... name show only work in
    gnu/linux, it would great find any app equivalent for windows or mac"""
    def __init__(self, config, parent=None):
        super(XTerm, self).__init__(parent)
        self.config = config
        self.process = QtCore.QProcess(self)
        atexit.register(self.kill)
        self.show_term()
        if config.getboolean("Applications", "maximized"):
            self.showMaximized()

    def kill(self):
        self.process.kill()
        self.process.waitForFinished()

    def sizeHint(self):
        size = QtCore.QSize(400, 300)
        return size.expandedTo(QtGui.QApplication.globalStrut())

    def show_term(self):
        term = self.config.get("Applications", 'Shell')

        args = [
            "-into", str(self.winId()),
            "-bg", self.config.get("Applications", "backgroundColor"),
            "-fg", self.config.get("Applications", "foregroundColor"),
            # blink cursor
            "-bc"]

        if self.config.getboolean("Applications", 'ipython') and \
                which("ipython"):
            args.append("ipython")
        else:
            args.append("python")

        self.process.start(term, args)
        if self.process.error() == QtCore.QProcess.FailedToStart:
            print("xterm not installed")


if __name__ == "__main__":
    import sys
    from configparser import ConfigParser
    import os
    app = QtGui.QApplication(sys.argv)
    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    pychemqt_dir = os.environ["PWD"] + "/"
    preferences = ConfigParser()
    preferences.read(conf_dir+"pychemqtrc")

    terminal = XTerm(preferences)
    terminal.show()
    app.exec_()
