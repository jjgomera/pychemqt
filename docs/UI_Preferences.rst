Pychemqt has many option to customize for the user, in this dialog user can change default options as desired. The dialog has several tab to organized configuration

General
-------

General main program options

  * Hightlight color in input
  * Read only color in input
  * Number of recent files shows in menu
  * Load last sessions project at start
  * Show tray icon

.. image:: images/UI_Preferences_1_General.png 
    :alt: UI_Preferences_1_General

PFD
---

Proces Flow Diagram (PFD) options

  * Background brush color
  * Background brush style
  * Input stream color
  * Output stream color
  * Line format options (color, width, join, cap, dash)
  * PFD resolution

.. image:: images/UI_Preferences_2_PFD.png 
    :alt: UI_Preferences_2_PFD


Tooltips in PFD
---------------

When the mouse is over a entity in PFD, a popup dialog with calculated properties can be shown with the properties defined in this tab

.. image:: images/UI_Preferences_3_TooltipPFD.png
    :alt: UI_Preferences_3_TooltipPFD.png


Tooltips in units
-----------------

When the mouse is over a entry with units support, a popup dialog with other units values can be shown as defined in this tab

.. image:: images/UI_Preferences_4_TooltipUnits.png
    :alt: UI_Preferences_4_TooltipUnits.png


Numeric Format
--------------

Prefered numeric format when show a number in any text 

.. image:: images/UI_Preferences_5_NumericFormat.png
    :alt: UI_Preferences_5_NumericFormat.png

The options include the normal format options

  * Fixed decimal point defining total and decimal digits
  * Significant figures
  * Exponential format by default
  * Exponential format for big/small values
  * Show sign in positives values
  * Show thousand separator (for currency for example)

.. image:: images/UI_Preferences_5_NumericFormatConfig.png
    :alt: UI_Preferences_5_NumericFormatConfig.png


Pseudocomponents
----------------

Define preferred methods to calculate properties when define oil pseudocomponents

  * Molecular weight
  * Critic properties
  * Critic volume
  * Acentric factor
  * Zc
  * Boiling temperature
  * Specific gravity
  * Refractive index
  * PNA composition
  * Destilate curve conversion
  * Hydrogen composition %

.. image:: images/UI_Preferences_6_Pseudocomponent.png
    :alt: UI_Preferences_6_Pseudocomponent.png


Applications
------------

External applications call from pychemqt

  * Calculator
  * Text file viewer
  * Use external pdf viewer, pychemqt has a simple internal pdf viewer
  * Shell used along options (background, foreground color, use ipython)

.. image:: images/UI_Preferences_7_Applications.png
    :alt: UI_Preferences_7_Applications.png


Plotting
--------

Plot can be configurated using the complete set of `matplotlib config options <https://matplotlib.org/stable/users/explain/customizing.html>`_

  * Can be use the matlotlib predefined styled
  * Customizing any of option available in matplotlibrc

.. image:: images/UI_Preferences_8_Matplotlib.png
    :alt: UI_Preferences_8_Matplotlib.png


qtelemental
-----------

Configuration avilable for the periodic table tools

  * Properties to define color squeme
  * Color scale

.. image:: images/UI_Preferences_9_qtelemental.png
    :alt: UI_Preferences_9_qtelemental.png


MEoS
----

Advanced properties library calculation options

  * Use external libraries (Coolprop, RefProp) if it's availables
  * Saturation line style (width, dash, color, marker)
  * Isolines shows in plots, as range or a list of lines values including line style config
  * Plot definition using more point (more slow)
  * Draw grid
  * 3D plot options (draw mesh, mesh type, colormap, wireframe line config)

.. image:: images/UI_Preferences_10_MEoS.png
    :alt: UI_Preferences_10_MEoS.png


Psychrometric tools
-------------------

Configuration for the psychrometric tool

  * Chart type (ASHRAE, Mollier)
  * Use virial or external application if availables
  * Saturation line style (width, dash, color, marker)
  * Crux line style (width, dash, color, marker)
  * Isolines shows in plots, as range or a list of lines values including line style config, label and its config

.. image:: images/UI_Preferences_11_Psychrometric.png
    :alt: UI_Preferences_11_Psychrometric.png


Moody chart
-----------

Moody chart configuration 

  * Method
  * Calculate fanning friction factor (by default it calculates the darcy friction factor)
  * ε/D isolines
  * Reliative roughtness line style (width, dash, color, marker)
  * Crux line style (width, dash, color, marker)
  * Grid line style (width, dash, color, marker)

.. image:: images/UI_Preferences_12_Moody.png
    :alt: UI_Preferences_12_Moody.png


Drag sphere chart
-----------------

Drag sphere chart configuration

  * Method
  * Drag coefficient line style (width, dash, color, marker)
  * Crux line style (width, dash, color, marker)
  * Grid line style (width, dash, color, marker)

.. image:: images/UI_Preferences_13_DragSphere.png
    :alt: UI_Preferences_13_DragSphere.png


Standing-Katz chart
-------------------

Standing-Katz chart configuration

  * Method
  * Tr lines
  * Reduced temperature line style (width, dash, color, marker)
  * Crux line style (width, dash, color, marker)
  * Grid line style (width, dash, color, marker)

.. image:: images/UI_Preferences_14_StandingKatz.png
    :alt: UI_Preferences_14_StandingKatz.png


Openbabel
---------

Configuration of libray openbabel to adjust the generated image

  * Bond color
  * Background color
  * Heteroatoms in color
  * Atoms details (show all atoms, only the terminal atoms)
  * Thicker bond lines
  * Asymetric double bond
  * Show atom index

.. image:: images/UI_Preferences_15_Openbabel.png
    :alt: UI_Preferences_15_Openbabel.png


API reference
-------------
