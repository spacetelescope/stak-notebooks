:orphan:


tables.fitsio
=============

The tables.fitsio package contains IO utilities for FITS and GEIS
images.

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

Many of the tasks in the this package are no longer in common usage and
are not covered here. If there is a task you would like to request
please contact the `STAK
team <http://stak.readthedocs.io/en/latest/>`__.

Contents:

-  `catfits <#catfits>`__
-  `stwfits <#stwfits>`__



catfits
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The catfits task was used to quickly produce a catalog of fits headers
from a file list. In the below example we provide the summary catalog
provided by ``astropy.io.fits``.

.. code:: ipython3

    # Standard Imports
    import glob
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change these values to your desired data files, glob will capture all wildcard matches
    test_data = glob.glob('/eng/ssb/iraf_transition/test_data/*.fits')
    
    for filename in test_data:
        fits.info(filename)


.. parsed-literal::

    Filename: /eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     266   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   float32   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/iczgs3ygq_newdtype_flt.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     266   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   int32   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     265   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   float32   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/imarith_out.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     266   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   float32   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/x31g0108t.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     112   ()      
      1  SCI         ImageHDU        28   (256, 256)   float32   
    Filename: /eng/ssb/iraf_transition/test_data/imfunction2_out.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     266   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   float64   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     270   ()      
      1  SCI         ImageHDU       200   (4096, 2048)   float32   
      2  ERR         ImageHDU        56   (4096, 2048)   float32   
      3  DQ          ImageHDU        48   (4096, 2048)   int16   
      4  SCI         ImageHDU       198   (4096, 2048)   float32   
      5  ERR         ImageHDU        56   (4096, 2048)   float32   
      6  DQ          ImageHDU        48   (4096, 2048)   int16   
      7  D2IMARR     ImageHDU        15   (64, 32)   float32   
      8  D2IMARR     ImageHDU        15   (64, 32)   float32   
      9  D2IMARR     ImageHDU        15   (64, 32)   float32   
     10  D2IMARR     ImageHDU        15   (64, 32)   float32   
     11  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     12  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     13  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     14  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     15  WCSCORR     BinTableHDU     59   14R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     270   ()      
      1  SCI         ImageHDU       200   (4096, 2048)   float32   
      2  ERR         ImageHDU        56   (4096, 2048)   float32   
      3  DQ          ImageHDU        48   (4096, 2048)   int16   
      4  SCI         ImageHDU       198   (4096, 2048)   float32   
      5  ERR         ImageHDU        56   (4096, 2048)   float32   
      6  DQ          ImageHDU        48   (4096, 2048)   int16   
      7  D2IMARR     ImageHDU        15   (64, 32)   float32   
      8  D2IMARR     ImageHDU        15   (64, 32)   float32   
      9  D2IMARR     ImageHDU        15   (64, 32)   float32   
     10  D2IMARR     ImageHDU        15   (64, 32)   float32   
     11  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     12  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     13  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     14  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     15  WCSCORR     BinTableHDU     59   14R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/imcopy_out.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     266   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   float32   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/imfunction_out.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     266   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   float32   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/imslice_out2.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  SCI         PrimaryHDU     199   (4096, 2048, 1)   float32   
    Filename: /eng/ssb/iraf_transition/test_data/hselect_test.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU       8   (6, 4)   float64   
      1              ImageHDU         9   (6, 4)   float64   
      2              ImageHDU         8   (6, 4)   float64   
    Filename: /eng/ssb/iraf_transition/test_data/imstack_out.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  SCI         PrimaryHDU     199   (4096, 2048, 2)   float32   
    Filename: /eng/ssb/iraf_transition/test_data/imslice_out1.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  SCI         PrimaryHDU     199   (4096, 2048, 1)   float32   
    Filename: /eng/ssb/iraf_transition/test_data/tester.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     270   ()      
      1  SCI         ImageHDU       200   (4096, 2048)   float32   
      2  ERR         ImageHDU        56   (4096, 2048)   float32   
      3  DQ          ImageHDU        48   (4096, 2048)   int16   
      4  SCI         ImageHDU       198   (4096, 2048)   float32   
      5  ERR         ImageHDU        56   (4096, 2048)   float32   
      6  DQ          ImageHDU        48   (4096, 2048)   int16   
      7  D2IMARR     ImageHDU        15   (64, 32)   float32   
      8  D2IMARR     ImageHDU        15   (64, 32)   float32   
      9  D2IMARR     ImageHDU        15   (64, 32)   float32   
     10  D2IMARR     ImageHDU        15   (64, 32)   float32   
     11  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     12  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     13  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     14  WCSDVARR    ImageHDU        15   (64, 32)   float32   
     15  WCSCORR     BinTableHDU     59   14R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/jczgx1ppq_rice.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     270   ()      
      1  SCI         CompImageHDU    200   (4096, 2048)   float32   
    Filename: /eng/ssb/iraf_transition/test_data/jczgx1ppq_gzip.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     270   ()      
      1  SCI         CompImageHDU    200   (4096, 2048)   float32   
    Filename: /eng/ssb/iraf_transition/test_data/fxsplit.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU       4   ()      
      1  DQ          ImageHDU        43   (1014, 1014)   int16   
    Filename: /eng/ssb/iraf_transition/test_data/fxinsert.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU       5   (100,)   float64   
      1  SCI         ImageHDU       200   (4096, 2048)   float32   
    Filename: /eng/ssb/iraf_transition/test_data/fxdelete.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU     265   ()      
      1  SCI         ImageHDU       140   (1014, 1014)   float32   
      2  ERR         ImageHDU        51   (1014, 1014)   float32   
      3  DQ          ImageHDU        43   (1014, 1014)   int16   
      4  SAMP        ImageHDU        37   (1014, 1014)   int16   
      5  TIME        ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR     BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   
    Filename: /eng/ssb/iraf_transition/test_data/empty.fits
    No.    Name         Type      Cards   Dimensions   Format
      0  PRIMARY     PrimaryHDU       4   ()      




stwfits
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

stwfits is used to translate a GEIS (Generic Edited Information Set),
STSDAS tables, or ascii file to an standard FITS(Flexible Image
Transport System) format. Here we will cover how to convert a GEIS file
to a FITS files using the ``stsci.tools.readgeis`` function. There are
two ways to use this function, through the command line, or through a
Python session or script. For instructions on running this task on the
command line see the ``stsci.tools`` `Conversion Utilities
documentation <http://ssb.stsci.edu/doc/stsci_python_dev/stsci.tools.doc/html/convert.html>`__.
Below we show an example of running this task in a python session. You
may or may not need to byteswap your image data depending on which
system it was originally written on.

.. code:: ipython3

    from stsci.tools import readgeis

.. code:: ipython3

    filename = '/eng/ssb/iraf_transition/test_data/x31g0108t.c0h'
    hdulist = readgeis.readgeis(filename)
    hdulist[1].data = hdulist[1].data.byteswap()
    del hdulist[1].header['CD1_1']
    del hdulist[1].header['CD2_2']
    hdulist.writeto('stwfits_out.fits', overwrite = True)


.. parsed-literal::

    ===================================
    = WARNING:                        =
    =  Input image:                   =
    /eng/ssb/iraf_transition/test_data/x31g0108t.c0h[1]
    =  had floating point data values =
    =  of NaN and/or Inf.             =
    ===================================
    ===================================
    =  This file may have been        =
    =  written out on a platform      =
    =  with a different byte-order.   =
    =                                 =
    =  Please verify that the values  =
    =  are correct or apply the       =
    =  '.byteswap()' method.          =
    ===================================
    




Not Replacing
-------------

-  fits\_example - used to provide more documentation for stwfits and
   strfits
-  fitscopy - used to produce a copy of a fits file, producing a copy of
   a fits file is straightforward in Python and the command line using
   exsisting libraries
-  geis - used to provide a description of GEIS file format
-  gftoxdim - GEIS conversion, no longer in common usage
-  strfits - converts FITS files to GEIS or STSDAS tables, no longer in
   common usage
-  xdimtogf - convert single group GEIS to multigroup GEIS, no longer in
   common usage
