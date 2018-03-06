:orphan:


stsdas.analysis.fitting
=======================

Curve fitting tools.

Notes
-----

**For questions or comments please see** `our github
page <https://github.com/spacetelescope/stak>`__. **We encourage and
appreciate user feedback.**

**Most of these notebooks rely on basic knowledge of the Astropy FITS
I/O module. If you are unfamiliar with this module please see the**
`Astropy FITS I/O user
documentation <http://docs.astropy.org/en/stable/io/fits/>`__ **before
using this documentation**.

The stsdas.analysis.fitting package is used to do various types of 1D
and 2D fitting. Fitting has been well developed in Python with the
`Scipy <https://docs.scipy.org/doc/scipy/reference/>`__ and
`Astropy <http://docs.astropy.org/en/stable/>`__ libraries. We have
covered these tasks in other IRAF notebooks, and refer to those entries
below.

Contents:

-  `gfit1d-nfit1d-ngaussfit <#gfit1d-nfit1d-ngaussfit>`__
-  `function <#function>`__





gfit1d-nfit1d-ngaussfit
-----------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Gfit1d and nfit1d are used to interactively achieve a 1-d fit to an
image, tables or lists. For various interactive curve fitting on images
please see the `Python
imexam <http://imexam.readthedocs.io/en/v0.7.1/imexam/imexam_command.html#>`__.
For non-interactive fitting of a table or list, please see the tasks in
the **tables.ttools** notebook, and the **images.imfit.fit1d-lineclean**
tasks.



function
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The function task is used to apply functions to images, tables or lists.
For image function application see **images.imutil.imfunction-imexpr**.
Once the table or list is coverted to a numpy array, the method for
applying a fucntion is the same as **images.imutil.imfunction-imexpr**.
See **tables.ttools.taextract-tainsert** for table conversions.





Not Replacing
-------------

-  controlpars - Pset with algorithm control parameters. Deprecated.
-  errorpars - Pset with error-related parameters. Deprecated.
-  bbodypars - Pset with parameters for black-body function. Deprecated.
-  cgausspars - Pset with parameters for constrained Gaussians function.
   Deprecated.
-  comppars - Pset with parameters for composite (bb + powerlaw)
   function. Deprecated.
-  galprofpars - Pset with parameters for galaxy profile function.
   Deprecated.
-  gausspars - Pset with parameters for Gaussians function. Deprecated.
-  i2gaussfit - Iterative 2-d Gaussian fit to noisy images. Deprecated.
-  n2gaussfit - 2-d Gaussian fit to images. See
   **images.imfit.imsurfit**
-  powerpars - Pset with parameters for powerlaw function. Deprecated.
-  prfit - Print contents of fit tables created by fitting task.
   Deprecated.
-  samplepars - Pset with data sampling parameters. Deprecated.
-  tgausspars - Pset with parameters for two-dim Gaussian function.
   Deprecated.
-  twobbpars - Pset with parameters for two-black-body function.
   Deprecated.
-  userpars - Pset with parameters for user-defined function.
   Deprecated.
