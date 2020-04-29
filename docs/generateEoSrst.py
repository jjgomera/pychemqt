#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Generate the EoS-*.rst files with list of references

import os

import lib


# Generate index file
txt = "lib.EoS module" + os.linesep
txt += "==============" + os.linesep + os.linesep
txt += "Library to implement the equation of state supported in the main "
txt += "program, the available equation include:" + os.linesep + os.linesep
txt += ".. toctree::" + os.linesep
txt += "    :maxdepth: 3" + os.linesep + os.linesep

for mod in lib.EoS.__all__:
    txt += "    %s" % mod.__name__ + os.linesep

with open("docs/lib.EoS.rst", "w") as file:
    file.write(txt)


# Generate index file for Cubic submodule
txt = "lib.EoS.Cubic module" + os.linesep
txt += "====================" + os.linesep + os.linesep
txt += "Library with the implemented cubic equation of state" + os.linesep
txt += os.linesep + ".. toctree::" + os.linesep
txt += "    :maxdepth: 1" + os.linesep + os.linesep

for mod in lib.EoS.Cubic._all:
    txt += "    lib.EoS.Cubic.%s" % mod.__name__ + os.linesep

txt += os.linesep
txt += "Furthemore in a separated file it's the common functionality"
txt += os.linesep + os.linesep

txt += ".. automodule:: lib.EoS.cubic" + os.linesep
txt += "    :members:" + os.linesep
txt += "    :undoc-members:" + os.linesep
txt += "    :private-members:" + os.linesep
txt += "    :show-inheritance:" + os.linesep
txt += "    :member-order: bysource" + os.linesep + os.linesep

txt += os.linesep + ".. include:: EoSnotimplement.rst" + os.linesep

txt += "References" + os.linesep
txt += "----------" + os.linesep
for id, rf in lib.EoS.cubic.__doi__.items():
    id = str(id)
    txt += ".. [%s] %s; %s. %s" % (id, rf["autor"], rf["title"], rf["ref"])
    txt += os.linesep

with open("docs/lib.EoS.Cubic.rst", "w") as file:
    file.write(txt)

# Generate the cubic eos files
for mod in lib.EoS.Cubic._all:
    library = mod.__name__

    with open("docs/lib.EoS.Cubic.%s.rst" % library, "w") as file:
        print("lib.%s module" % library, file=file)
        print("="*(len(library)+4+7), file=file)
        print("", file=file)
        print(".. automodule:: lib.EoS.Cubic.%s" % library, file=file)
        print("    :members:", file=file)
        print("    :undoc-members:", file=file)
        print("    :private-members:", file=file)
        print("    :show-inheritance:", file=file)
        print("    :member-order: bysource", file=file)

        print("", file=file)
        print("References", file=file)
        print("----------", file=file)

        count = 1
        for lnk in mod.__doi__:
            ref = ".. [%i] %s; %s, %s" % (
                count, lnk["autor"], lnk["title"], lnk["ref"])
            if lnk["doi"]:
                ref += ", http://dx.doi.org/%s" % lnk["doi"]
            count += 1
            print(ref, file=file)

# Generate files for simple equation methods
for mod in lib.EoS.__all__:
    if len(mod._all) > 1:
        continue

    library = mod.__name__
    if library == "lib.EoS.virial":
        # Create only the references file
        with open("docs/%s_ref.rst" % library, "w") as file:
            print("", file=file)
            print("References", file=file)
            print("----------", file=file)

            for id, rf in mod.__doi__.items():
                id = str(id)
                print(".. [%s] %s; %s. %s" % (
                    id, rf["autor"], rf["title"], rf["ref"]), file=file)

    else:
        with open("docs/%s.rst" % library, "w") as file:
            print("%s module" % library, file=file)
            print("="*(len(library)+7), file=file)
            print("", file=file)
            print(".. automodule:: %s" % library, file=file)
            print("    :members:", file=file)
            print("    :undoc-members:", file=file)
            print("    :private-members:", file=file)
            print("    :show-inheritance:", file=file)
            print("    :member-order: bysource", file=file)

            print("", file=file)
            print("References", file=file)
            print("----------", file=file)

            count = 1
            for lnk in mod._all[0].__doi__:
                ref = ".. [%i] %s; %s, %s" % (
                    count, lnk["autor"], lnk["title"], lnk["ref"])
                if lnk["doi"]:
                    ref += ", http://dx.doi.org/%s" % lnk["doi"]
                count += 1
                print(ref, file=file)

