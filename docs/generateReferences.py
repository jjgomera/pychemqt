#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Module autogeneration of documentation files. This script autogenerate the
lib.*.rst and equipment.*.rst files with the automodule schema

lib.? module
====================

.. automodule:: lib.?
    :members:
    :undoc-members:
    :show-inheritance:
    :member-order: bysource

.. include:: lib.*_ref.rst

Furthermore autogenerate the lib.*_ref.rst files with a list of bibliographic
references used in the module
"""

import glob
from os.path import basename, splitext

import lib
import equipment

total = []

all = lib.__all__
# for mod in lib.EoS.__all__:
    # all.append(".".join(mod.__name__.split(".")[1:]))

# library module autogeneration
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

        # Add the *-ref.rst file to lib.*.rst file if it's available
        if hasattr(module, "__doi__") and module.__doi__:
            print("", file=file)
            print(".. include:: lib.%s_ref.rst" % library, file=file)

    # Generate the *-ref.rst files with list of references
    if hasattr(module, "__doi__") and module.__doi__:
        with open("docs/lib.%s_ref.rst" % library, "w") as file:
            print("References", file=file)
            print("----------", file=file)
            for id, rf in module.__doi__.items():
                id = str(id)
                print(".. [%s] %s; %s. %s" % (
                    id, rf["autor"], rf["title"], rf["ref"]), file=file)

# equipment module autogeneration
lst = glob.glob("./equipment/[!^UI|__]*.py")
for mod in lst:
    library = splitext(basename(mod))[0]

    __import__("equipment.%s" % library)
    module = equipment.__getattribute__(library)

    # Make equipment.rst schemas
    with open("docs/equipment.%s.rst" % library, "w") as file:
        print("equipment.%s module" % library, file=file)
        print("="*(len(library)+10+7), file=file)
        print("", file=file)
        print(".. automodule:: equipment.%s" % library, file=file)
        print("    :members:", file=file)
        print("    :undoc-members:", file=file)
        print("    :show-inheritance:", file=file)
        print("    :member-order: bysource", file=file)
