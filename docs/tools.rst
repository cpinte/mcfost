Analysis tools
==============

For visualization and fitting tools, several packages are available so far for
MCFOST:

-  the python package `pymcfost <https://github.com/cpinte/pymcfost>`__ (lead author: Pinte).

   .. important:: pymcfost is the recommended and most up-to-date analysis tool for mcfost.

-  an alternative python package `mcfost-python <https://github.com/swolff9/mcfost-python>`__ lead author: S. Wolff).
-  an IDL GUI interface, MCRE, mostly designed to explore a large number of models (lead author: M. Perrin)
-  a more simple, command-line based IDL interface (lead author: G. Duchene)
-  a command-line based, `Yorick interface <https://github.com/cpinte/yomcfost>`__ (lead author: Pinte, deprecated and unsupported). Most functions (but not all) have now been ported to pymcfost.

Only the analysis of the results are done with IDL or Yorick or Python.
They are not required for the use of MCFOST, but make its use easier,
especially when dealing with large numbers of models.
