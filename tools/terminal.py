#!/usr/bin/python
# -*- coding: utf-8 -*-

###############################################################################
# Tools to create a python shell with pychemqt libraries imported
# For now only work in linux with xterm as terminal
###############################################################################

import atexit

from PyQt5 import QtCore, QtWidgets

from lib.firstrun import which


class XTerm(QtCore.QProcess):
    """Gui container for terminal widget"""
    def __init__(self, config, parent=None):
        super(XTerm, self).__init__(parent)
        self.config = config
        atexit.register(self.kill)
        self.show_term()

    def sizeHint(self):
        size = QtCore.QSize(400, 300)
        return size.expandedTo(QtWidgets.QApplication.globalStrut())

    def show_term(self):
        term = self.config.get("Applications", 'Shell')

        args = [
            "-bg", self.config.get("Applications", "backgroundColor"),
            "-fg", self.config.get("Applications", "foregroundColor"),
            # blink cursor
            "-bc", 
            # title
            "-T", QtWidgets.QApplication.translate(
                "pychemqt", "pychemqt python console")]
            
        if self.config.getboolean("Applications", "maximized"):
            args.append("-maximized")

        if self.config.getboolean("Applications", 'ipython') and \
                which("ipython"):
            args.append("ipython3")
        else:
            args.append("python3")

        self.start(term, args)
        if self.error() == QtCore.QProcess.FailedToStart:
            print("xterm not installed")


if __name__ == "__main__":
    import sys
    from configparser import ConfigParser
    import os
    app = QtWidgets.QApplication(sys.argv)
    conf_dir = os.path.expanduser('~') + "/.pychemqt/"
    pychemqt_dir = os.environ["PWD"] + "/"
    preferences = ConfigParser()
    preferences.read(conf_dir+"pychemqtrc")

    terminal = XTerm(preferences)
    app.exec_()
