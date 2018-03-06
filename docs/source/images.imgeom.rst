:orphan:


images.imgeom
=============

The images.imgeom package contains various image spatial manipulation
and interpolation tasks.

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

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress

**General Note about this package:**

The tasks in this IRAF package include support for WCS updates after the
image manipulations are preformed. We are currently working on the
replacements for these WCS capabilities. For the moment we have included
the array manipulation part of these tasks in this notebook.

**Boundary Condition and Interpolation Options:**

In the IRAF package the interpolation options are as follows: nearest,
linear, poly3, poly5, spline3, sinc/lsinc. In the ``scipy.ndimage``
functions, the interpolation option are spline degrees 0-5, where
spline-0 is nearest and spline-1 is linear.

The boundary condition options for IRAF and ``scipy`` are the same:
nearest, wrap, reflect, and constant.

Contents:

-  `blkavg <#blkavg>`__
-  `blkrep <#blkrep>`__
-  `im3dtran-imtranspose <#im3dtran-imtranspose>`__
-  `imshift-shiftlines <#imshift-shiftlines>`__
-  `magnify <#magnify>`__
-  `rotate <#rotate>`__



blkavg
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The blkavg task takes in an arbitrary dimensioned image and performs a
block average across the requested box size. We can preform the same
task with the ``astropy.nddata.utils.block_reduce`` function. In fact,
this function is more generalized as you can apply any function (not
just averaging) that accepts an ``ndarray`` and axis keyword to the
block filter.

Below we show an example that reads in a FITS file, runs
``block_reduce`` on the first image extension, and saves out the updated
file.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy.nddata.utils import block_reduce

.. code:: ipython3

    # Change this value to your desired data file, here were creating a filename
    # for our new changed data
    orig_data = '/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits'
    new_fits = 'iczgs3ygq_newdtype_flt.fits'
    
    # Read in your fits file
    hdu = fits.open(orig_data)
    
    # Print FITS summary
    hdu.info()
    
    # Run block reduce on first sci extension
    red_data = block_reduce(hdu[1].data, [2,2], func=np.mean)
    
    # Re-insert the data into the HDUList
    hdu[1].data = red_data
    
    # Save out fits to a new file
    hdu.writeto(new_fits, overwrite=True)
    hdu.close()


.. parsed-literal::

    Filename: /eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     270   ()      
      1  SCI           1 ImageHDU       200   (4096, 2048)   float32   
      2  ERR           1 ImageHDU        56   (4096, 2048)   float32   
      3  DQ            1 ImageHDU        48   (4096, 2048)   int16   
      4  SCI           2 ImageHDU       198   (4096, 2048)   float32   
      5  ERR           2 ImageHDU        56   (4096, 2048)   float32   
      6  DQ            2 ImageHDU        48   (4096, 2048)   int16   
      7  D2IMARR       1 ImageHDU        15   (64, 32)   float32   
      8  D2IMARR       2 ImageHDU        15   (64, 32)   float32   
      9  D2IMARR       3 ImageHDU        15   (64, 32)   float32   
     10  D2IMARR       4 ImageHDU        15   (64, 32)   float32   
     11  WCSDVARR      1 ImageHDU        15   (64, 32)   float32   
     12  WCSDVARR      2 ImageHDU        15   (64, 32)   float32   
     13  WCSDVARR      3 ImageHDU        15   (64, 32)   float32   
     14  WCSDVARR      4 ImageHDU        15   (64, 32)   float32   
     15  WCSCORR       1 BinTableHDU     59   14R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   




blkrep
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The task ``blkrep`` is used to block replicate an n-dimensional image.
Astropy has the equivalent function ``block_replicate``.

For an example of how to read in a FITS extension, edit the image array,
and save out the updated file, see `blkavg <#blkavg>`__ above.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.nddata.utils import block_replicate

.. code:: ipython3

    # test data
    my_arr = np.array(([[0., 1.], [2., 3.]]))
    print("input array")
    print(my_arr)
    
    # conservation of the array sum is the default
    out = block_replicate(my_arr, 3)
    print("sum conservation")
    print(out)
    
    # you can changes this using conserve_sum=False
    out = block_replicate(my_arr, 3, conserve_sum=False)
    print("no sum conservation")
    print(out)


.. parsed-literal::

    input array
    [[ 0.  1.]
     [ 2.  3.]]
    sum conservation
    [[ 0.          0.          0.          0.11111111  0.11111111  0.11111111]
     [ 0.          0.          0.          0.11111111  0.11111111  0.11111111]
     [ 0.          0.          0.          0.11111111  0.11111111  0.11111111]
     [ 0.22222222  0.22222222  0.22222222  0.33333333  0.33333333  0.33333333]
     [ 0.22222222  0.22222222  0.22222222  0.33333333  0.33333333  0.33333333]
     [ 0.22222222  0.22222222  0.22222222  0.33333333  0.33333333  0.33333333]]
    no sum conservation
    [[ 0.  0.  0.  1.  1.  1.]
     [ 0.  0.  0.  1.  1.  1.]
     [ 0.  0.  0.  1.  1.  1.]
     [ 2.  2.  2.  3.  3.  3.]
     [ 2.  2.  2.  3.  3.  3.]
     [ 2.  2.  2.  3.  3.  3.]]




im3dtran-imtranspose
--------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Tasks used to transpose images.
`numpy.transpose <https://docs.scipy.org/doc/numpy/reference/generated/numpy.transpose.html>`__
can handle any number of dimensions.

For an example of how to read in a FITS extension, edit the image array,
and save out the updated file, see `blkavg <#blkavg>`__ above.

.. code:: ipython3

    # Standard Imports
    import numpy as np

.. code:: ipython3

    # in_array constructs a 3 x 5 array of integer values from 0 to 14
    in_array = np.arange(15).reshape(5,3)
    # we then transpose it it to a 5 x 3 array
    out_array = np.transpose(in_array)
    
    print('Original array:')
    print(in_array)
    print('Transpose of original array')
    print(out_array)


.. parsed-literal::

    Original array:
    [[ 0  1  2]
     [ 3  4  5]
     [ 6  7  8]
     [ 9 10 11]
     [12 13 14]]
    Transpose of original array
    [[ 0  3  6  9 12]
     [ 1  4  7 10 13]
     [ 2  5  8 11 14]]




imshift-shiftlines
------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The task imshift can shift an image in x and y by float values and will
use interpolation to create the output image. Shiftlines preformed
similar functionality but We will be using
`scipy.ndimage.shift <https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.ndimage.shift.html#scipy.ndimage.shift>`__,
where you can shift in any axis of your image. See the
`Notes <#notes>`__ at the top of the notebook for fitting and boundary
options.

For an example of how to read in a FITS extension, edit the image array,
and save out the updated file, see `blkavg <#blkavg>`__ above.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    from scipy.ndimage import shift

.. code:: ipython3

    # Don't forget that Python uses (y,x) format when specifiying shifts
    in_array = np.arange(25).reshape(5,5)
    out_array = shift(x, (0.8,0.8), order=3, mode='constant', cval=2)
    
    print('Original array:')
    print(in_array)
    print('A zoom of 0.5 in y and 2 in x with nearest')
    print(out_array)


.. parsed-literal::

    Original array:
    [[ 0  1  2  3  4]
     [ 5  6  7  8  9]
     [10 11 12 13 14]
     [15 16 17 18 19]
     [20 21 22 23 24]]
    A zoom of 0.5 in y and 2 in x with nearest
    [[ 2  2  2  2  2]
     [ 2  0  2  2  4]
     [ 2  6  7  8  9]
     [ 2 11 12 13 14]
     [ 2 16 18 19 20]]




magnify
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The task magnify takes an image and magnifies the image by the desired
amount, using a chosen iterpolation. The interpolation options
avaialable for the magnify task are nearest, linear, poly3, poly5,
spine3, sinc, lsinc, and drizzle. We will be using
`scipy.ndimage.zoom <https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.ndimage.zoom.html#scipy.ndimage.zoom>`__
as a python equivalent. For this task, the available interpolation
options are nearest, and spline0-5 fits.

For an example of how to read in a FITS extension, edit the image array,
and save out the updated file, see `blkavg <#blkavg>`__ above.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    from scipy.ndimage import zoom

.. code:: ipython3

    # Don't forget that Python uses (y,x) format when specifiying magnification
    in_array = np.arange(25).reshape(5,5)
    out_array = zoom(in_array, (0.5,2.5), order=0)
    
    print('Original array:')
    print(in_array)
    print('A zoom of 0.5 in y and 2.5 in x with nearest')
    print(out_array)


.. parsed-literal::

    Original array:
    [[ 0  1  2  3  4]
     [ 5  6  7  8  9]
     [10 11 12 13 14]
     [15 16 17 18 19]
     [20 21 22 23 24]]
    A zoom of 0.5 in y and 2.5 in x with nearest
    [[ 0  0  1  1  1  2  2  2  3  3  3  4  4]
     [10 10 11 11 11 12 12 12 13 13 13 14 14]
     [20 20 21 21 21 22 22 22 23 23 23 24 24]]




rotate
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The task rotate is used to rotate and shift images. We will only cover
rotation here, for shifting please see `shiftlines <#shiftlines>`__. We
will be using
`scipy.ndimage.rotate <https://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.ndimage.interpolation.rotate.html>`__
for rotation using interpolation. For a simple 90 degree unit rotation
we will use
`numpy.rot90 <https://docs.scipy.org/doc/numpy/reference/generated/numpy.rot90.html#numpy.rot90>`__
(you can do a 90 degree rotation using ``scipy.ndimage.rotate``).

For an example of how to read in a FITS extension, edit the image array,
and save out the updated file, see `blkavg <#blkavg>`__ above.

Rotation using interpolation:

.. code:: ipython3

    # Standard Imports
    import numpy as np
    from scipy.ndimage import rotate

.. code:: ipython3

    in_array = np.arange(25).reshape(5,5)
    # Rotate by 60 degrees
    out_array = rotate(in_array, 60, axes=(1,0))
    
    print('Original array:')
    print(in_array)
    print('A rotation of 60 degrees')
    print(out_array)


.. parsed-literal::

    Original array:
    [[ 0  1  2  3  4]
     [ 5  6  7  8  9]
     [10 11 12 13 14]
     [15 16 17 18 19]
     [20 21 22 23 24]]
    A rotation of 60 degrees
    [[ 0  0  0  0  0  0  0]
     [ 0  0  3  9  0  0  0]
     [ 0  0  5 11 15 21  0]
     [ 0  2  7 12 17 22  0]
     [ 0  3  9 13 19  0  0]
     [ 0  0  0 15 21  0  0]
     [ 0  0  0  0  0  0  0]]


Rotation in increments of 90 degrees:

.. code:: ipython3

    # Standard Imports
    import numpy as np

.. code:: ipython3

    in_array = np.arange(25).reshape(5,5)
    # Rotate by 270 degrees
    out_array = np.rot90(in_array, 3)
    
    print('Original array:')
    print(in_array)
    print('A rotation of 270 degrees')
    print(out_array)


.. parsed-literal::

    Original array:
    [[ 0  1  2  3  4]
     [ 5  6  7  8  9]
     [10 11 12 13 14]
     [15 16 17 18 19]
     [20 21 22 23 24]]
    A rotation of 270 degrees
    [[20 15 10  5  0]
     [21 16 11  6  1]
     [22 17 12  7  2]
     [23 18 13  8  3]
     [24 19 14  9  4]]






Not Replacing
-------------

-  imlintran - see `images.imgeom.magnify <#magnify>`__,
   `images.imgeom.rotate <#rotate>`__, and
   `images.imgeom.imshift <#imshift>`__
