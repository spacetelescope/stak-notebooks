:orphan:


noao.imred.crutil
=================

The noao.imred.crutil package contains various algorithims for finding
and replacing comsic rays in single images or image sets.

Notes
-----

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress

Contents:

-  `cosmicrays <#cosmicrays>`__
-  `craverage <#craverage>`__
-  `crfix <#crfix>`__
-  `crgrow <#crgrow>`__
-  `crmedian <#crmedian>`__
-  `crnebula <#crnebula>`__



cosmicrays
----------

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress



craverage
---------

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress



crfix
-----

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress



crgrow
------

**Please review the `Notes <#notes>`__ section above before running any
examples in this notebook**

The crgrow replacement uses the ``skimage.morphology`` package to grow
the values in any numpy array. The dilation task is a wrapper around
``scipy.ndimage.grey_dilation``. You can insert any kernal type where
``disk`` is called in this example.

.. code:: ipython2

    # Standard Imports
    from skimage.morphology import disk,dilation
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change this value to your desired data file
    test_data = '/eng/ssb/iraf_transition/test_data/id0k16pdq_blv_tmp.fits'
    
    # Read in your fits file, when using a fits file, the bytesway call is required to
    # make sure your arry data type is correct.
    hdu = fits.open(test_data,mode='update')
    dq1 = hdu[3].data.byteswap().newbyteorder('=')
    
    # Dilation used to grow the CR flags
    grownDQ = dilation(dq1, disk(2))
    
    # Re-assign the changed array to our original fits file and close the file to save.
    hdu[3].data = grownDQ
    hdu.close()



crmedian
--------

**Please review the `Notes <#notes>`__ section above before running any
examples in this notebook**

The crmedian task is a way to indentify and replace cosmic rays in a
single image by detecting pixels that deviate a statistically
significant amount from the median by comparing to a median filtered
version of the image. The indentfied cosmic rays can then be replaced by
the median filtered value. A similar algorithim has been used in
`ccdproc.cosmicray\_median <http://ccdproc.readthedocs.io/en/latest/api/ccdproc.cosmicray_median.html#ccdproc.cosmicray_median>`__.
In ``ccdproc.cosmicray_median`` you also have the option of using an
error array. If none is provided the standard deviation of the data is
used.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy import units
    from ccdproc import cosmicray_median, fits_ccddata_reader

.. code:: ipython2

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    
    # First we need to pull out the science arrays to create CCDData objects
    # Our acutal unit is electrons/sec, this is not accepted by the current
    # set of units
    image_data = fits_ccddata_reader(test_data, hdu=1, unit=units.electron/units.s, hdu_uncertainty=2)
    error_data = image_data.uncertainty.array
    
    # Now we run cosmicray_median, since we input a CCDData type, a CCDData type is returned
    # If a numpy.ndarray if the input data type, it will return a numpy.ndarray
    newdata = cosmicray_median(image_data, error_image=error_data, thresh=5, mbox=11, rbox=11, gbox=3)


.. parsed-literal::

    INFO: using the unit electron / s passed to the FITS reader instead of the unit ELECTRONS/S in the FITS file. [ccdproc.ccddata]




crnebula
--------

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress



Not Replacing
-------------

-  crcombine - see **ctio.immatch.imcombine, work in progress**
-  credit - see **images.tv.imedit, work in progress**

For questions or comments please see `our github
page <https://github.com/spacetelescope/stak>`__. We encourage and
appreciate user feedback.
