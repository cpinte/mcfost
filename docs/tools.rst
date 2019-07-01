Analysis tools
==============

.. important:: The recommended analysis tool is `pymcfost <https://github.com/cpinte/pymcfost>`__. Below are briefs descriptions of the alternative options, but support for them is now very limited.

For visualization and fitting tools, several packages are available so far for
MCFOST:

-  an IDL GUI interface (lead author: Perrin)
-  a more simple, command-line based IDL interface (lead author: Duchene)
-  a command-line based, Yorick interface (lead author: Pinte, deprecated).
-  the python package `pymcfost <https://github.com/cpinte/pymcfost>`__ (lead author: Pinte).
-  an alternative, unmaintained python package `mcfost-python <https://github.com/cpinte/pymcfost>`__ lead author: Perrin).

Only the analysis of the results are done with IDL or Yorick or Python.
They are not required for the use of MCFOST, but make its use easier,
especially when dealing with large numbers of models.


Tools for Creating Parameter Files
----------------------------------

**From IDL:** There is a script `crea_grid_mcfost` originally by M.
Perrin / G. Duchene and recently updated by Elodie Choquet.

1. Edit a parameter file by hand, to set all the parameters you want
   constant.

2. Run the IDL command ``crea_grid_mcfost``, with keywords setting the
   parameters you want to vary. e.g.::

     $ crea_grid_mcfost, mdust=[1e-9, 1e-8, 1e-7], amax=[1,10,100]

   would create a set of 9 models with varying dust mass and max
   grain size.

3. If you want to add additional elements onto an existing grid, call
   the same function again with the new desired parameters, and add
   the init_counter= parameter set to the starting number to use
   (one greater than the number of models output in the prior run)

**From Yorick:** This tool is maintained by C. Pinte. Usage is similar
to the IDL interface. It also includes grid fitting and Bayesian
inference tools, a genetic algorithm, an ant colony algorithm and
parallel MCMC sampler.

The IDL and yorick interfaces will be progressively deprecated and be
replaced by the python interface.

**From Python:** There are Python tools under development by Marshall
Perrin, Christophe Pinte, Schuyler Wolff et al., including both grid
creation, model display, and model fitting.

The code is available from `github <https://github.com/cpinte/mcfost-python>`__

Note that every time the parameter file format changes in a new version
of MCFOST, all of the above scripts must be updated to be compatible
with the changes..

Tools for Computing a Model Grid
--------------------------------

The right way to do this depends on what machine you're running on.

-  To run on a single computer, you can just use a for loop in a shell
   script, or the equivalent

-  On the FOSTINO cluster, use distrib_mcfost.py or
   distrib_addimages.py (ask Marshall) or the new yorick interface
   (ask Christophe)

-  On a Mac with xgrid set up, use mcfost_runxgrid.py (Marshall)

IDL Tools for Visualizing MCFOST Output
---------------------------------------

MCRE (MCFOST Results Explorer, author : Marshall) is a graphical
application for browsing and analyzing large (10k-100k+) model grids
computed using MCFOST. It has its own manual; please see that document.

There are a number of stand-alone simple tools that are bundled with
MCRE, and provide a simple IDL command line interface to plotting
various outputs.

-  ``mcseds`` Plot SED for all inclinations provided by a model

-  ``mcsed`` Plot SED for online one inclination

-  ``mcimgs`` Display images at all inclinations (for one wavelength)

-  ``mcimg`` Display image at one inclination and one wavelength

-  ``mc_rt`` Replaced mcimg for raytraced image - should merge code.

-  ``mcwaveimgs`` Display images at all wavelengths (broken?)

-  ``mcdiff`` Overplot 2 SEDs for comparison

-  ``mcrim`` Plot puffed-up inner rim graphically

-  ``mcgrid`` Compute coordinates for each grid cell in a computation
   (broken?)

-  ``mctemp`` Display temperature plot

Python Tools for Visualizing MCFOST Output
------------------------------------------

Marshall and Christophe have started developing a python module for
manipulating MCFOST output. Thus far it is very incomplete compared to
the IDL functionality (but does have nicer looking plots in general,
since it's Python). This is available from github, including its own
documentation. Very much a work in progress, with contributions
welcomed.

The python-code is available via `github <https://github.com/cpinte/mcfost-python>`__.

Yorick tools for Vizualizing MCFOST Output
------------------------------------------

The Yorick interface (author : Christophe) includes a number of command
line functions to display and manipulate various outputs. Detailed
documentation of the function is included in the yorick source code.
Here is a list of the most useful commands :

-  ``open_mcfost`` read a mcfost model and put the parameter file and
   resulating calculations in a yorick structure

-  ``read_params`` and ``write_params`` read a parameter file into a
   yorick structure (and write a yorick structure to parameter file)

-  ``plot_sed`` plot SED with extinction, can separate the various flux
   contribution

-  ``interfero (and calc_vis)`` plot (or output) visibility curves

-  ``plot_temp`` plot temperature structure

-  ``plot_dens`` plot density structure

-  ``surface_density`` compute and plot the surface density

-  ``plot_contrib`` plot the various flux contribution in an image

-  ``plot_pola`` compute and plot polarization maps from the Stokes
   maps

-  ``sum_mcfost`` merge several Monte-Carlo models

-  ``plot_line`` plot line data (integrated map and intergrated line
   profile)

-  ``spectral_index`` compute the spectral index of a model between 2
   wavelengths

-  ``partial_mass`` compute the model mass between 2 radii

-  ``ProDiMo_star_2_mcfost_star`` convert a ProDiMo ASCII input
   stellar spectrum to a MCFOST fits input stellar spectrum

-  ``create_grid`` create a uniform grid of models

-  ``sedfit`` fit and vis a SED or visibility data set with a grid of
   models and compute the Bayesian inference

-  ``start_grid`` set up and launch a grid on FOSTINO

-  ``stat_grid`` check the status of a running grid on FOSTINO

-  ``mcfost_grid2ascii`` output MCFOST grid fitting results to an ASCII
   file

-  ``mcfost_genetic`` run a genetic model fitting (to be executed on
   FOSTINO)

-  ``mcfost_eMCMC`` run a parallel MCMC (to be executed on FOSTINO)
