lib.compuestos module
=====================

.. automodule:: lib.compuestos
    :members:
    :undoc-members:
    :show-inheritance:

.. exec::
    from lib import compuestos
    print("References")
    print("----------")
    for id, rf in compuestos.__doi__.items():
        print(".. [%i] %s; %s. %s" % (id, rf["autor"], rf["title"], rf["ref"]))
