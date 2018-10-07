#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Generate the *-ref.rst files with list of references

import lib

all = lib.__all__ 
# for mod in lib.EoS.__all__:
    # all.append(".".join(mod.__name__.split(".")[1:]))

total = []
for library in all:
    if library in ["EoS", "mEoS"]:
        continue

    # Import module to intronspection
    __import__("lib.%s" % library)
    if "EoS" in library:
        submodule, name = library.split(".")
        module = lib.__getattribute__(submodule).__getattribute__(name)
    else:
        module = lib.__getattribute__(library)

    # Make lib.rst schemas
    with open("docs/lib.%s.rst" % library, "w") as file:
        print("lib.%s module" % library, file=file)
        print("="*(len(library)+4+7), file=file)
        print("", file=file)
        print(".. automodule:: lib.%s" % library, file=file)
        print("    :members:", file=file)
        print("    :undoc-members:", file=file)
        print("    :show-inheritance:", file=file)
        print("    :member-order: bysource", file=file)

        if hasattr(module, "__doi__") and module.__doi__:
            print("", file=file)
            print(".. include:: lib.%s_ref.rst" % library, file=file)

    if hasattr(module, "__doi__") and module.__doi__:
        with open("docs/lib.%s_ref.rst" % library, "w") as file:
            print("References", file=file)
            print("----------", file=file)
            for id, rf in module.__doi__.items():
                id = str(id)
                print(".. [%s] %s; %s. %s" % (
                    id, rf["autor"], rf["title"], rf["ref"]), file=file)
