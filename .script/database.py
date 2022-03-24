#!/usr/bin/python3
# -*- coding: utf-8 -*-


from math import log10, log
import csv
import sqlite3

connection = sqlite3.connect('/home/jjgomera/pychemqt/dat/databank.db')
# Add table for wagner vapor pressure equation
cursor = connection.cursor()

# Add table from caleb thermo for Antoine extended equation
with open('/home/jjgomera/Programacion/chemeng/thermo/thermo/Vapor Pressure/Wagner Original McGarry.tsv') as csvfile:
    reader = csv.reader(csvfile, delimiter="\t")
    header = True
    for CAS, name, A, B, C, D, Pc, Tc, Tmin in reader:
        if header:
            header = False
            continue

        A = float(A)
        B = float(B)
        C = float(C)
        D = float(D)

        query = "UPDATE compuestos SET wagner_A==%f, wagner_B==%f, " % (A, B)
        query += "wagner_C==%f, wagner_D==%f WHERE CAS=='%s'" % (C, D, CAS)
        print(query)
        cursor.execute(query)
        connection.commit()

