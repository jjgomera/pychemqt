tools.UI_Tables package
=========================
Tool for use the multiparamter equation of state of several fluids to high accuracy property calculations.

The module can use several libraries as calculation procedures:

    * The internal mEoS library
    * CoolProp
    * RefProp

The fluids availables are:

.. include:: lib.mEoSlst.rst


All functionality are integrate in pychemqt main program, accesible in main menu Tools/MEoS. Furthermore these equations of state are available to use in process stream properties calculation.

.. image:: images/UITables.png 
    :alt: UITables

Configuration
-------------

**Fluid**

In a open project we must first define the fluid to calculate properties, Tools/MEoS/Fluid. When we have defined a fluid the name appear in this menu option.

.. image:: images/mEoSFluid.png 
    :alt: mEos Fluid

In this dialog we can choose the EoS to use, and if there several options, the method to calculate the viscosity or thermal conductivity.

**Reference**

We can define the reference state for entropy and enthalpy derived properties, Tools/MEoS/Reference. We can choose between any of standard reference states or define a custom. Now unimplemented

.. image:: images/mEoSReference.png 
    :alt: mEos Reference

**Properties**

It can be defined the properties to show, and the order too, Tools/MEoS/Properties.

.. image:: images/mEoSProperties.png 
    :alt: mEos Reference

The properties availables depend of backend used in library.

For configure the library to use and several parameters of generated plot, use the options Tools/MEoS/Configure, or use the tab in main program Preferences dialog.


