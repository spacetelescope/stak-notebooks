:orphan:


noao.imred.crutil
=================

The noao.imred.crutil package contains various algorithms for finding
and replacing cosmic rays in single images or image sets.

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

Contents: \* `crgrow <#crgrow>`__ \* `crmedian <#crmedian>`__



crgrow
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The crgrow replacement uses the ``skimage.morphology`` package to grow
the values in any numpy array. The dilation task is a wrapper around
``scipy.ndimage.grey_dilation``. You can insert any kernal type where
``disk`` is called in this example. See the
`skimage.morphology.dilation <http://scikit-image.org/docs/dev/api/skimage.morphology.html#skimage.morphology.dilation>`__
for more information. More kernel shapes are also listed on this page.

.. code:: ipython3

    # Standard Imports
    from skimage.morphology import disk,dilation
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change this value to your desired data file
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    out_file = 'crgrow.fits'
    
    # Read in your fits file, when using some fits file, the byteswap call is required to
    # make sure your array data type is correct when the dilation function is used. This
    # may be due to a bug in the dilation funciton.
    # For this example we will work with the 3rd extensions, the DQ array
    hdu = fits.open(test_data,mode='update')
    hdu.info()
    dq1 = hdu[3].data.byteswap().newbyteorder('=')
    
    # Dilation used to grow the CR flags, here we use the disk radius 2 shape kernel
    grownDQ = dilation(dq1, disk(2))
    
    # Re-assign the changed array to our original fits file and close save the updated FITS to a new file.
    hdu[3].data = grownDQ
    hdu.writeto(out_file, overwrite=True)


.. parsed-literal::

    Filename: /eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     266   ()      
      1  SCI           1 ImageHDU       140   (1014, 1014)   float32   
      2  ERR           1 ImageHDU        51   (1014, 1014)   float32   
      3  DQ            1 ImageHDU        43   (1014, 1014)   int16   
      4  SAMP          1 ImageHDU        37   (1014, 1014)   int16   
      5  TIME          1 ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR       1 BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   




crmedian
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The crmedian task is a way to identify and replace cosmic rays in a
single image by detecting pixels that deviate a statistically
significant amount from the median by comparing to a median filtered
version of the image. The identified cosmic rays can then be replaced by
the median filtered value. A similar algorithm has been used in
`ccdproc.cosmicray\_median <http://ccdproc.readthedocs.io/en/latest/api/ccdproc.cosmicray_median.html#ccdproc.cosmicray_median>`__.
In ``ccdproc.cosmicray_median`` you also have the option of using an
error array. If none is provided the standard deviation of the data is
used. Ccdproc is an evolving package, please see `their
documentation <https://ccdproc.readthedocs.io/en/latest/>`__ for more
information on usage.

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy import units
    from ccdproc import cosmicray_median, fits_ccddata_reader

.. code:: ipython3

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    
    # First we need to pull out the science and error(uncertainty) array to 
    # create CCDData objects. Our acutal unit is electrons/sec, this is not
    # accepted by the current set of units
    image_data = fits_ccddata_reader(test_data, hdu=1, unit=units.electron/units.s, hdu_uncertainty=2)
    error_data = image_data.uncertainty.array
    
    # Now we run cosmicray_median, since we input a CCDData type, a CCDData type is returned
    # If a numpy.ndarray if the input data type, it will return a numpy.ndarray
    newdata = cosmicray_median(image_data, error_image=error_data, thresh=5, mbox=11, rbox=11, gbox=3)


.. parsed-literal::

    INFO: using the unit electron / s passed to the FITS reader instead of the unit ELECTRONS/S in the FITS file. [ccdproc.ccddata]




Not Replacing
-------------

-  cosmicrays - Remove cosmic rays using flux ratio algorithm.
-  craverage - Detect CRs against average and avoid objects.
-  crcombine - Combine multiple exposures to eliminate cosmic rays.
-  credit - Interactively edit cosmic rays using an image display.
-  crfix - Fix cosmic rays in images using cosmic ray masks.
-  crnebula - Detect and replace cosmic rays in nebular data.
