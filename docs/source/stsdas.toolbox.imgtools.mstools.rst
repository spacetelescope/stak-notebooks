:orphan:


stsdas.toolbox.imgtools.mstools
===============================

Tasks to handle HST imsets.

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

The imgtools.mstools package contains tasks for working with STIS,
NICMOS, ACS, and WFC3 data. Some tasks are 'extensions" of existing
tasks in the STSDAS system, and support other instruments/file formats
as well.

Contents:

-  `ecdel-ecextract-extdel-msdel-msjoin-mssplit <#ecdel-ecextract-extdel-msdel-msjoin-mssplit>`__
-  `mscombine <#mscombine>`__
-  `msstatistics <#msstatistics>`__





ecdel-ecextract-extdel-msdel-msjoin-mssplit
-------------------------------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

These tasks contain various methods for deleting or moving extensions
around in FITS files. This can be easily done using the
``astropy.io.fits`` module in Astropy. Here is a `good
page <http://docs.astropy.org/en/stable/io/fits/>`__ to familarize
yourself with this package.



mscombine
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The original ``mscombine`` IRAF task performed image combination of
several ``SCI`` exetensions of HST data while allowing the user to
reject specified ``DQ`` bits. Additionally, the user could choose to
combine the stack using the average or the median.

This ``mscombine`` alternative uses the ``bitfield_to_boolean_mask``
task and ``numpy`` masked arrays to avoid using flagged pixels in the
``DQ`` array. In this simple example, we average-combine several
full-frame WFC3/UVIS images.

Tasks for image combination are currently being developed in the
``CCDPROC`` package, see the `CCDPROC doc
page <https://ccdproc.readthedocs.io/en/latest/#>`__ for more details or
the **images.imutil.imsum** task for a short usage example.

.. code:: ipython3

    # Standard Imports
    import glob
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits
    from stsci.tools.bitmask import bitfield_to_boolean_mask

.. code:: ipython3

    # Get the data
    test_data = glob.glob('/eng/ssb/iraf_transition/test_data/mscombine/*_blv_tmp.fits')

.. code:: ipython3

    # Create masked arrays
    masked_arrays_ext1, masked_arrays_ext2, masked_arrays_ext4, masked_arrays_ext5 = [], [], [], []
    for filename in test_data:
        with fits.open(filename) as hdulist:
            
                # For UVIS chip 2
                mask_ext3 = np.invert(bitfield_to_boolean_mask(hdulist[3].data, ignore_flags=0))
                masked_arrays_ext1.append(np.ma.masked_array(hdulist[1].data, mask=mask_ext3))
                masked_arrays_ext2.append(np.ma.masked_array(hdulist[2].data, mask=mask_ext3))
    
                # For UVIS chip 1            
                mask_ext6 = np.invert(bitfield_to_boolean_mask(hdulist[6].data, ignore_flags=0))
                masked_arrays_ext4.append(np.ma.masked_array(hdulist[4].data, mask=mask_ext6))
                masked_arrays_ext5.append(np.ma.masked_array(hdulist[5].data, mask=mask_ext6))

.. code:: ipython3

    # Average-combine SCI arrays
    comb_ext1 = np.ma.mean(masked_arrays_ext1, axis=0).data
    comb_ext4 = np.ma.mean(masked_arrays_ext4, axis=0).data

.. code:: ipython3

    # Propoagate uncertainties for ERR arrays, divide by zero expected
    weight_image_ext1 = np.zeros((2051, 4096))
    weight_image_ext4 = np.zeros((2051, 4096))
    for array in masked_arrays_ext1:
        mask = array.mask
        weight_image_ext1[np.where(mask == False)] += 1.0
    for array in masked_arrays_ext4:
        mask = array.mask
        weight_image_ext4[np.where(mask == False)] += 1.0
    masked_arrays_ext2_squared = [(item * (1/weight_image_ext1))**2 for item in masked_arrays_ext2]
    masked_arrays_ext5_squared = [(item * (1/weight_image_ext4))**2 for item in masked_arrays_ext5]
    comb_ext2 = np.sqrt(np.ma.sum(masked_arrays_ext2_squared, axis=0)).data
    comb_ext5 = np.sqrt(np.ma.sum(masked_arrays_ext5_squared, axis=0)).data


.. parsed-literal::

    /Users/ogaz/miniconda2/envs/irafdev3/lib/python3.6/site-packages/ipykernel_launcher.py:10: RuntimeWarning: divide by zero encountered in true_divide
      # Remove the CWD from sys.path while we load stuff.
    /Users/ogaz/miniconda2/envs/irafdev3/lib/python3.6/site-packages/ipykernel_launcher.py:11: RuntimeWarning: divide by zero encountered in true_divide
      # This is added back by InteractiveShellApp.init_path()


.. code:: ipython3

    # Create empty DQ arrays
    comb_ext3 = np.zeros((2051, 4096))
    comb_ext6 = np.zeros((2051, 4096))

.. code:: ipython3

    # Build and save the combined file, using the first final for the header
    hdu0 = fits.PrimaryHDU(header=fits.getheader(test_data[0], 0))
    hdu1 = fits.ImageHDU(comb_ext1, header=fits.getheader(test_data[0], 0))
    hdu2 = fits.ImageHDU(comb_ext2, header=fits.getheader(test_data[0], 1))
    hdu3 = fits.ImageHDU(comb_ext3, header=fits.getheader(test_data[0], 2))
    hdu4 = fits.ImageHDU(comb_ext4, header=fits.getheader(test_data[0], 3))
    hdu5 = fits.ImageHDU(comb_ext5, header=fits.getheader(test_data[0], 4))
    hdu6 = fits.ImageHDU(comb_ext6, header=fits.getheader(test_data[0], 5))
    hdulist = fits.HDUList([hdu0, hdu1, hdu2, hdu3, hdu4, hdu5, hdu6])
    hdulist.writeto('mscombine_test.fits', overwrite=True)



msstatistics
------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The msstatictics task is similiar to images.imutil.imstatistics, but
with the added capability to mask using HST bit masking in the Data
Quality extensions. Here we show an example of the ``stsci.tools``
`bitfield\_to\_boolean\_mask <https://github.com/spacetelescope/stsci.tools/blob/master/lib/stsci/tools/bitmask.py>`__
function.

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy import stats
    from stsci.tools.bitmask import bitfield_to_boolean_mask

.. code:: ipython3

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    # multiple reads file
    #test_data = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits
    hdulist = fits.open(test_data)
    
    # Make mask, using bit flags 32 and 4
    boolean_mask = bitfield_to_boolean_mask(hdulist[3].data,"~4,128")
    
    # The sigma_clipped_stats function returns the mean, median, and stddev respectively
    mean, median, std = stats.sigma_clipped_stats(hdulist[1].data, mask=boolean_mask, sigma=2.0, iters=3)
    print("mean: {}".format(mean))
    print("median: {}".format(median))
    print("standard deviation: {}".format(std))
    
    # Close fits file
    hdulist.close()


.. parsed-literal::

    mean: 2.05237880697
    median: 0.847149193287
    standard deviation: 2.83607972172






Not Replacing
-------------

-  msarith - Image arithmetic with NICMOS and STIS files. See
   **images.imutil.imarith**.
-  mscopy - Copy image sets of a multi-extension FITS file. See
   **images.imutil.imcopy**
-  mssort - Sort a FITS file to get all extensions of like version
   number. Deprecated.
