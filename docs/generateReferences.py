#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Define the pychemqt environmentl
import os

os.environ["pychemqt"] = "/home/jjgomera/Programacion/pychemqt/"
os.environ["freesteam"] = "False"
os.environ["pybel"] = "False"
os.environ["CoolProp"] = "False"
os.environ["refprop"] = "False"
os.environ["ezodf"] = "False"
os.environ["openpyxl"] = "False"
os.environ["xlwt"] = "False"
os.environ["icu"] = "False"
os.environ["reportlab"] = "False"
os.environ["PyQt5.Qsci"] = "False"

# Generate the *-ref.rst files with list of references
import lib

all = lib.__all__ 
for mod in lib.EoS.__all__:
    all.append(".".join(mod.__name__.split(".")[1:]))

total = []
for library in all:
    if library in ["EoS", "mEoS"]:
        continue

    __import__("lib.%s" % library)
    if "EoS" in library:
        submodule, name = library.split(".")
        module = lib.__getattribute__(submodule).__getattribute__(name)
    else:
        module = lib.__getattribute__(library)
    if hasattr(module, "__doi__") and module.__doi__:
        with open("docs/lib.%s_ref.rst" % library, "w") as file:
            print("References", file=file)
            print("----------", file=file)
            for id, rf in module.__doi__.items():
                id = str(id)
                print(".. [%s] %s; %s. %s" % (
                    id, rf["autor"], rf["title"], rf["ref"]), file=file)
