:orphan:


utilities
=========

A miscellaneous utilities package.

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

The utilities package can be easily replicated with Python
functionality.

Contents:

-  `detab-entab <#detab-entab>`__
-  `lcase-ucase <#lcase-ucase>`__
-  `urand <#urand>`__





detab-entab-translit
--------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Detab and entab are used to replace tabs with blanks or vice versa.
Python contains a built-in `replace
method <https://docs.python.org/3.6/library/stdtypes.html#string-methods>`__
which can be used for this purpose.



lcase-ucase
-----------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

lcase and ucase are used to convert a file to lower case or uppercase.
Python contains `built in
functionality <https://docs.python.org/3.6/library/stdtypes.html#string-methods>`__
for string replacement, including ``lower`` and ``upper`` methods.



urand
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Urand provides a uniform random number generator. This utility is
contained in ``Numpy`` with the
`numpy.randome.uniform <https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.uniform.html>`__
task. ``Numpy`` also contains `other types of random
generation <https://docs.scipy.org/doc/numpy/reference/routines.random.html>`__.





Not Replacing
-------------

-  curfit - Fit data with Chebyshev, Legendre or spline curve. See
   **images.imfit.it1d-lineclean**
-  polyfit - Fit polynomial to list of X,Y data. See
   **images.imfit.fit1d-lineclean**
-  surfit - Fit a surface, z=f(x,y), to a set of x, y, z points. See
   **images.imfit.imsurfit**
-  split - Split a large file into smaller segments. Deprecated.
