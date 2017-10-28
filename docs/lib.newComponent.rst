lib.newComponent module
=======================

.. automodule:: lib.newComponent
    :members:

.. exec::
    from lib.newComponent import __doi__
    print("References")
    print("----------")
    for id, rf in __doi__.items():
        print(".. [%i] %s; %s. %s" % (id, rf["autor"], rf["title"], rf["ref"]))
