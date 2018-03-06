:orphan:


images.tv
=========

The tv package contains interactive image display utilities.

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

Many of the images.tv functionality has been replaced by the Python
imexam package, which works in conjunction with both the
`ginga <https://ginga.readthedocs.io/en/latest/>`__ and
`ds9 <http://ds9.si.edu/site/Home.html>`__ viewers. In this notebook we
will simply provide the task name for imexam that covers the IRAF task.
We leave it to the user to reference the `imexam
documentation <http://imexam.readthedocs.io/en/latest/>`__ for further
information.

Contents:

-  `display <#display>`__
-  `imexamine <#imexamine>`__
-  `bpmedit-imedit-badpixcorr <#bpmedit-imedit-badpixcorr>`__
-  `tvmark <#tvmark>`__
-  `wcslab <#wcslab>`__



display
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

display - Load an image or image section into the display

-  see
   http://imexam.readthedocs.io/en/latest/imexam/examples.html#basic-usage
   and
   http://imexam.readthedocs.io/en/latest/imexam/imexam\_command.html#cutout-a-simple-fits-image



imexamine
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

imexamine - Examine images using image display, graphics, and text

-  see http://imexam.readthedocs.io/en/latest/imexam/iraf\_imexam.html



bpmedit-imedit-badpixcorr
-------------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

bpmedit-badpixcorr - examine and edit bad pixel masks associated with
images (badpixcorr in ginga)

-  see
   http://stginga.readthedocs.io/en/latest/stginga/plugins\_manual/badpixcorr.html

imedit - Examine and edit pixels in images

-  see
   http://imexam.readthedocs.io/en/v0.7.1/imexam/imexam\_command.html#pixel-coordinates-and-value,
   combine with
   `astropy.io.fits <http://docs.astropy.org/en/stable/io/fits/>`__



tvmark
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

tvmark - Mark objects on the image display

-  see
   http://ginga.readthedocs.io/en/latest/manual/plugins\_local/tvmark.html



wcslab
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

wcslab - Overlay a displayed image with a world coordinate grid.

-  see
   http://ginga.readthedocs.io/en/latest/manual/plugins\_local/wcsaxes.html?highlight=wcs



Not Replacing
-------------

-  iis subpackage - these utilities are replaced by
   `ginga <https://ginga.readthedocs.io/en/latest/index.html>`__ and
   `ds9 <http://ds9.si.edu/site/Home.html>`__.
