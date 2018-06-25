Input Data Files
================

All the data used by mcfost is generally read from the directory
``$MCFOST_UTILS``. An additional environment variable ``MY_MCFOST_UTILS`` can
be defined (giving the path to another directory) by the user to add
her/his own data file. This has to advantage to ensure that no personal
data file will be overwritten during an update.

When looking for data files, mcfost will first look in the local working
directory, then in ``MY_MCFOST_UTILS`` if it is defined and finally in
``MCFOST_UTILS``. If the data file used is in the local directory, it
will be copied into the output directory as well.

The data directory also has a version number which correponds to the
minimum MCFOST version which can use this data set (ie MCFOST version >
MCFOST_UTILS version). MCFOST will ask the user to update the data
directory if needed (note this also means that you won't be able to use
an old MCFOST version if the data directory is updated).

Dust Optical Constants
----------------------

Filenames that are underlined in the following list indicate
well-identified dust species (i.e., the origin of the file can be
tracked in the literature with reasonable certainty). Comments at the
top the files must start with a #.

The first line contains the material
density in g/cm\ :sup:`3` and the sublimation temperature. It is
followed by an empty line and then 3 columns containing the wavelengths
(in microns) and the real and imaginary parts of the optical indices.
Wavelengths can be in ascending or descending order but must be ordered.

-  ``CO2_50K-reduit.dat:`` CO2 ice at 50K (origin?)

-  ``Draine_Si.dat "``\ Astronomical silicates" from Bruce Draine
   *(very similar though not exactly the same as those shown in Fig
   3 of Draine & Lee 1984, ApJ, 285, 89 - use dlsi_opct.dat
   instead);* file ``Draine_Si_sUV.dat`` is the version which has
   been smoothed in the UV.

-  ``Graphite_para.dat`` graphite with a "parallel" orientation (from
   Draine & Lee 1984, ApJ, 285, 89?)

-  ``Graphite_perp.dat`` graphite with a "perpendicular" orientation
   (from Draine & Lee 1984, ApJ, 285, 89?)

-  ``MW89.dat`` Highly porous composite interstellar grains. From
   `Mathis and Whiffen
   1989 <http://adsabs.harvard.edu/cgi-bin/bib_query?1989ApJ...341..808M>`__.

-  ``Olivine_modif.dat`` Amorphous MgFeSiO4 olivine. From Dorschner,
   Begemann, et al. 1995. A&A 300 503, see their Fig 5). File
   ``Olivine.dat`` contains the same optical indices but has no
   information for lambda > 500 micron, forcing the code to
   extrapolate. File ``Olivine_no_SI.dat`` contains the same
   optical indices as ``Olivine_modif.dat``, except for an
   artificial linear interpolation across the 10 and 20 micron
   features to suppress them (no obvious reason why that would be
   useful...).

-  ``ac_opct.dat`` (amorphous) carbonaceous dust, Rouleau & Martin
   1991, ApJ, 377, 526

-  ``crstsi_opct.dat`` crystalline silicate dust (olivine), optical
   constants derived by Li, A., & Draine, B. T. 2001, ApJ, 550,
   L213, no corresponding optical constants figure in the literature

-  ``dlsi_opct.dat`` amorphous silicate dust ("astronomical
   silicates") from Draine, B. T., & Lee, H. M. 1984, ApJ, 285, 89

-  ``ice_opct.dat`` H2O-dominated "dirty" ice (mostly crystalline)
   compiled from several sources described in Li, A., & Greenberg,
   J. M. 1998, A&A, 331, 291. no corresponding optical constants
   figure in the literature

-  ``indices_acar.dat`` "amorphous carbon" optical indices file used
   by Sebastian Wolf with MC3D (origin?). *Judging from the shape of
   the plots, this material also contain some silicates, though.*

-  ``indices_acar2.dat`` "amorphous carbon" with high porosity (~85%)
   optical indices file used by Sebastian Wolf with MC3D (origin?).
   *Judging from the shape of the plots, this material also contain
   some silicates, though.*

-  ``indices_mc3d.dat`` "silicates" optical indices file used by
   Sebastian Wolf with MC3D; identical to the indices in file
   ``dlsi_opct.dat``, except that this file has no coverage for
   lambda < 0.2 micron. File ``indices_mc3d_sil.dat`` is
   identical.

-  ``lgsi_opct.dat`` amorphous olivine (MgFeSiO_4), based on lab
   data and synthesized by Li, A., & Greenberg, J. M. 1997, A&A,
   323, 566, optical constants displayed as solid line in Figure 1
   of Li, A., & Greenberg, J. M. 1997, A&A, 323, 566

-  ``porousice_james.dat`` Water ice with optical indices only
   provided at two wavelengths (very limited usefulness)

-  ``waterice210K.dat`` Water ice (at 210K?), as compiled in Pollack
   et al. 1994, ApJ, 421, 615

-  ``PAHion.dat`` ionized PAHs computed by B. T. Draine (compiled from
   Laor, A., & Draine, B.T. 1993, ApJ 402:441, Draine, B.T., & Lee,
   H.M. 1984, ApJ 285:89, and Li, A., & Draine, B.T. 2001, ApJ
   554:778)

-  ``PAHneu.dat`` neutral PAHs computed by B. T. Draine (compiled from
   Laor, A., & Draine, B.T. 1993, ApJ 402:441, Draine, B.T., & Lee,
   H.M. 1984, ApJ 285:89, and Li, A., & Draine, B.T. 2001, ApJ
   554:778)

Stellar Models
--------------

We should also add references for these.
