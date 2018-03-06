:orphan:


noao.digiphot.daophot
=====================

The daophot package is used to preform Dao crowded-field photometry.

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

Large parts of the daophot functionality has been replaced in the
astropy.photutils package.

We recommend looking at the photutils documentation to explore the new
workflow. But if you are trying to achieve a specific subfunction of
daophot, this notebook maybe be useful. The tasks are grouped together
when they reference the same doumentation.

Contents:

-  `allstar-peak-phot-psf <#allstar-peak-phot-psf>`__
-  `daofind-nstar-substar <#daofind-nstar-substar>`__
-  `group-pstselect <#group-pstselect>`__
-  `grpselect-pcalc-pconcat-pdump-pselect-psort <#grpselect-pcalc-pconcat-pdump-pselect-psort>`__
-  `addstar <#addstar>`__
-  `pstselect <#pstselect>`__





allstar-peak-phot-psf
---------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

allstar - Group and fit psf to multiple stars simultaneously, with the
initial esitmates of positions already provided.

peak - Fit the psf to single stars.

phot - Compute sky values and initial magnitudes for a list of stars use
photometry call.

psf - Compute the point spread function.

-  see photutils `performing PSF photometry with fixed
   centriods <https://photutils.readthedocs.io/en/stable/photutils/psf.html#performing-psf-photometry-with-fixed-centroids>`__
   and
   `photutils.psf.BasicPSFPhotometry <https://photutils.readthedocs.io/en/stable/api/photutils.psf.BasicPSFPhotometry.html#photutils.psf.BasicPSFPhotometry>`__.

-  photutils `PSF
   photometry <http://photutils.readthedocs.io/en/stable/photutils/psf.html?highlight=psf%20fitting#psf-photometry>`__
   and
   `photutils.psf.DAOPhotPSFPhotometry <http://photutils.readthedocs.io/en/stable/api/photutils.psf.DAOPhotPSFPhotometry.html#photutils.psf.DAOPhotPSFPhotometry>`__.

-  `photutils.psf.IterativelySubtractedPSFPhotometry <http://photutils.readthedocs.io/en/stable/api/photutils.psf.IterativelySubtractedPSFPhotometry.html#photutils.psf.IterativelySubtractedPSFPhotometry>`__



daofind-nstar-substar
---------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

daofind - Find stars in an image using the DAO algorithm.

nstar - Fit the psf to predefined groups of stars.

substar - Subtract the fitted stars from the original image.

-  See photutils `source
   detection <https://photutils.readthedocs.io/en/stable/photutils/detection.html>`__
   and
   `photutils.DAOStarFinder <https://photutils.readthedocs.io/en/stable/api/photutils.DAOStarFinder.html#photutils.DAOStarFinder>`__.



group-pstselect
---------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

group - Group stars based on positional overlap and signal/noise.

-  See photutils
   `groups <https://photutils.readthedocs.io/en/stable/photutils/grouping.html>`__
   and
   `photutils.psf.DAOGroup <https://photutils.readthedocs.io/en/stable/api/photutils.DAOGroup.html#photutils.DAOGroup>`__.



grpselect-pcalc-pconcat-pdump-pselect-psort
-------------------------------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

grpselect - Select groups of a specified size from a daophot database.

pcalc - Do arithmetic operations on a list of daophot databases.

pconcat - Concatenate a list of daophot databases, astropy table.

pdump - Print selected fields from a list of daophot databases.

pselect - Select records from a daophot database.

psort - Sort a daophot database.

pfmerge - Merge a list of photometry databases.

-  The output of ``photutils.psf.DAOGroup``, linked in the
   `group <#notes>`__ entry, as well as most photutils outputs are an
   `Astropy
   Table <http://docs.astropy.org/en/stable/table/index.html>`__. This
   class can be `sorted, split, joined, and have a function
   applied <http://docs.astropy.org/en/stable/table/operations.html>`__.
   Basic arithmetic can be done automatically by simply calling the
   desired columns with the operator of your choice. See the
   **tables.ttools** notebook in STAK for the ttools mapping.



addstar
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

addstar - Add artificial stars to an image using the computed psf.

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress



pstselect
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

pstselect - Select candidate psf stars based on proximity.

-  see ``photutils.psf.DAOGroup`` ``find_group``
   `method <http://photutils.readthedocs.io/en/stable/api/photutils.psf.DAOGroup.html#photutils.psf.DAOGroup>`__.





Not Replacing
-------------

Build photutils doc pages, and put in link

-  centerpars - Edit the centering algorithm parameters. Deprecated.
-  daoedit - Review/edit algorithm parameters interactively. Deprecated.
-  daopars - Edit the daophot algorithms parameter set. See `photutils
   documentation <https://photutils.readthedocs.io/en/stable/>`__.
-  daotest - Run basic tests on the daophot package tasks. Deprecated.
-  datapars - Edit the image data dependent parameters. Deprecated
-  findpars - Edit the star detection parameters. Deprecated.
-  fitskypars - Edit the sky fitting algorithm parameters. Deprecated.
-  photpars - Edit the aperture photometry parameters. Deprecated.
-  setimpars - Save/restore parameter sets for a particular image.
   Depreciated.
-  pconvert - Convert a text database to a tables database, see `Astropy
   unified I/O <http://docs.astropy.org/en/stable/io/unified.html>`__
-  pexamine - Interactively examine and edit a daophot database.
   Deprecated, use Astropy Table tools
-  prenumber - Renumber stars in a daophot database, see
   `grpselect-pcalc-pconcat-pdump-pselect-psort <#grpselect-pcalc-pconcat-pdump-pselect-psort>`__
-  seepsf - convert a sampled PSF lookup table to a PSF image.
   Deprecated.
