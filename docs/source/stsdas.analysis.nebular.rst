:orphan:


stsdas.analysis.nebular
=======================

Tasks for analyzing nebular emission lines.

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

The nebular package has been replaced in Python by the PyNeb package,
documentation can be found `here <http://pythonhosted.org/PyNeb/>`__,
with a `homepage here <http://www.iac.es/proyecto/PyNeb/>`__. In this
notebook we will simply list the equivalent Python task for the original
IRAF task where possible.

Contents:

-  `abund <#abund>`__
-  `ionic <#ionic>`__
-  `ntcontour-ntplot <#ntcontour-ntplot>`__
-  `redcorr <#redcorr>`__
-  `temden <#temden>`__
-  `zones <#zones>`__





abund
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

abund - Derive ionic abundances relative to H+ in 3-zone nebula.

-  See the getIonAbundance function and section 1.9 of the PyNeb
   handbook.



ionic
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

ionic - Determine ionic abundance relative to H+.

-  see the Atom class in PyNeb



ntcontour-ntplot
----------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

ntcontour - Plot contours of N\_e- or T\_e-sensitive line ratios.

ntplot - Construct N\_e vs. T\_e plot for observed diagnostic ratios.

-  See section 2.4 of the PyNeb handbook.



redcorr
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

redcorr - Correct line flux for interstellar reddening.

-  See the RedCorr class in PyNeb.



temden
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

temden - Determine electron temperature or density from diagnostic
ratio.

-  See section 1.11 of the PyNeb handbook.



zones
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

zones - Determine electron temps & densities in 3-zone nebula.

-  See the Atom object in PyNeb.





Not Replacing
-------------

-  at\_data - Documentation on the atomic reference data. Derecated.
-  diagcols - Pset for zone temperature & density column names.
   Deprecated.
-  fluxcols - Pset for line flux column names. Deprecated.
-  nlevel - Documentation on the N-level atom approximation. Deprecated.
