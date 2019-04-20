#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Generate the *-ref.rst files with list of references

import lib

all = lib.__all__ 

total = []
for library in all:
    # Import module to intronspection
    __import__("lib.%s" % library)
    module = lib.__getattribute__(library)

    # Special case for metalibreries
    if library in ["mEoS", "EoS"]:
        for cmp in module.__doi__:
            for eq in module.__doi__[cmp]:
                rf = module.__doi__[cmp][eq]
                if rf not in total:
                    total.append(rf)
        continue

    # General case for simple libraries
    # Make lib.rst schemas
    with open("docs/lib.%s.rst" % library, "w") as file:
        print("lib.%s module" % library, file=file)
        print("="*(len(library)+4+7), file=file)
        print("", file=file)
        print(".. automodule:: lib.%s" % library, file=file)
        print("    :members:", file=file)
        print("    :undoc-members:", file=file)
        print("    :private-members:", file=file)
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
                if rf not in total:
                    total.append(rf)

# Global references file
with open("docs/references.rst", "w") as file:
    print("References", file=file)
    print("----------", file=file)

    id = 0
    for lnk in sorted(total, key=lambda lnk: str.upper(lnk["autor"])):
        if lnk["autor"] or lnk["title"] or lnk["ref"]:
            id += 1
            ref = ".. [%i] %s; %s, %s" % (
                id, lnk["autor"], lnk["title"], lnk["ref"])
            if lnk["doi"]:
                ref += ", http://dx.doi.org/%s" % lnk["doi"]
            print(ref, file=file)
