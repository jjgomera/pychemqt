#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Console tool for easy chemical elements translation configuration

import sqlite3

# Connection to database with element data
connection = sqlite3.connect('dat/elemental.db')
databank = connection.cursor()
connection2 = sqlite3.connect('dat/elemental.db')
data2 = connection2.cursor()


def create_table():
    # Create table if not exist
    databank.execute("""CREATE TABLE IF NOT EXISTS TRANSLATION \
                            (id INTEGER PRIMARY KEY, name TEXT)""")

    # Add name elements from ELEMENTS table
    databank.execute("SELECT id, name FROM ELEMENTS")
    for id, name in databank:
        data2.execute("insert into TRANSLATION (name) values (?)", (name, ))
    connection2.commit()


def addTranslation(lng):
    """Add schema for new language"""
    databank.execute("ALTER TABLE TRANSLATION ADD COLUMN name_%s TEXT" % lng)


def translate(lng):
    """Console procedure for edit translation"""
    column = "name_%s" % lng
    databank.execute("SELECT id, name, %s FROM TRANSLATION" % column)
    data = []
    for element in databank:
        data.append(element)

    for id, name, translation in data:
        if translation:
            print("%3i %-20s%s" % (id, name, translation))
        else:
            tr = input("%i: %s: " % (id, name))
            databank.execute("UPDATE TRANSLATION SET %s==? WHERE id==?" % column, (tr, id))
            connection.commit()

# addTranslation("ZH")
translate("ZH")


