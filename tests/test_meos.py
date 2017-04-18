#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
from unittest import TestLoader, TestSuite

from lib.mEoS import __all__


loader = TestLoader()
path = os.path.dirname(os.path.realpath(sys.argv[0]))
tests = loader.discover(os.path.join(path, "lib", "mEoS"), pattern="*.py")
TestMEOS = TestSuite(tests)


# test_cmp = Test()
# TestMeos.addTest(test_cmp)

# runner = TextTestRunner(verbosity=2)
# runner.run(suite)



# from lib.mEoS import __all__



# class Test(unittest.TestCase):
    # def test_meos(self):
        # for module in __all__:
            # print(module.__module__)
            # if "Air" not in module.__module__:
                # continue




# for module in __all__:
    # if "D2O" not in module.__module__:
        # continue
    # print(module.__module__)
    # inst = module()
    # for eq in inst.eq[0:1]:
        # if "__test__" in eq:
            # inst.__doc__ += eq["__test__"]
    # if inst._viscosity is not None:
        # for eq in inst._viscosity:
            # if "__test__" in eq:
                # inst.__doc__ += eq["__test__"]
    # if inst._thermal is not None:
        # for eq in inst._thermal:
            # if "__test__" in eq:
                # inst.__doc__ += eq["__test__"]
    # doctest.run_docstring_examples(inst, globs={module.__name__: module})
# #    timeit.timeit("test()", setup="from __main__ import test", number=3)

