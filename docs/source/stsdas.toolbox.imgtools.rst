:orphan:


stsdas.toolbox.imgtools
=======================

Tasks for performing operations on images and masks.

Notes
-----

**For questions or comments please see** `our github
page <https://github.com/spacetelescope/stak>`__. **We encourage and
appreciate user feedback.**

The various image tasks found in the stsdas.toolbox.imgtools package
have been replaced in the `Numpy <https://docs.scipy.org/doc/numpy/>`__
and `Astropy <http://docs.astropy.org/en/stable/>`__ libraries.

Contents:

-  `addmasks <#addmasks>`__
-  `iminsert <#iminsert>`__
-  `improject <#improject>`__
-  `mkgauss <#mkgauss>`__
-  `pixlocate <#pixlocate>`__
-  `rd2xy-xy2rd <#rd2xy-xy2rd>`__





addmasks
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Addmasks is used to combine several masks or bad pixel lists. We can do
this using the ``numpy`` bitwise tasks:
`bitwise\_or <https://docs.scipy.org/doc/numpy/reference/generated/numpy.bitwise_or.html>`__,
`bitwise\_and <https://docs.scipy.org/doc/numpy/reference/generated/numpy.bitwise_and.html>`__,
and
`invert <https://docs.scipy.org/doc/numpy/reference/generated/numpy.invert.html>`__,
along with a slew of `other numpy bit
functions <https://docs.scipy.org/doc/numpy/reference/routines.bitwise.html>`__.
Below we show examples of ``bitwise_and`` and ``bitwise_or``.

.. code:: ipython2

    # Standard Imports
    import numpy as np

.. code:: ipython2

    a = np.array([1,4,10])
    b = np.array([1,0,8])
    
    # OR
    print(np.bitwise_or(a,b))
    
    # AND
    print(np.bitwise_and(a,b))


.. parsed-literal::

    [ 1  4 10]
    [1 0 8]




iminsert
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Iminsert is used to insert a small image into a larger image. This is
easy to do with the Numpy array indexing after you've read in your
images with ``Astropy.io.fits``. Below we'll show a quick array example.

.. code:: ipython2

    # Standard Imports
    import numpy as np

.. code:: ipython2

    # generate test arrays
    my_array = np.random.rand(7,7)
    ones = np.array(([1,1,1],[1,1,1],[1,1,1]))
    
    # replace middle 3x3 square with ones
    my_array[2:5,2:5] = ones
    
    # Print result
    print(my_array)


.. parsed-literal::

    [[ 0.23716783  0.49467955  0.17823965  0.17247774  0.39261896  0.63899393
       0.1590178 ]
     [ 0.29311263  0.89693848  0.85860996  0.78426103  0.79732931  0.42044396
       0.18684822]
     [ 0.48983825  0.34300706  1.          1.          1.          0.31296056
       0.2697361 ]
     [ 0.30018196  0.00991694  1.          1.          1.          0.17958182
       0.90966563]
     [ 0.89446453  0.74249067  1.          1.          1.          0.23503997
       0.73016154]
     [ 0.2939247   0.71674139  0.99701222  0.24283608  0.14886197  0.61169257
       0.25279554]
     [ 0.65994208  0.7224029   0.76381896  0.68119012  0.77388175  0.22266843
       0.77680789]]




improject
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Improject is used to sum or average an image along one axis. This can be
accomplised using the
`numpy.average <https://docs.scipy.org/doc/numpy/reference/generated/numpy.average.html>`__
or the
`numpy.sum <https://docs.scipy.org/doc/numpy/reference/generated/numpy.sum.html>`__
functions and choosing which dimensions you wish to collapse. Below we
show an example using ``numpy.average``.

.. code:: ipython2

    # Standard Imports
    import numpy as np

.. code:: ipython2

    # build random test array
    my_array = np.random.rand(10,10,3)
    
    # reduce third dimension down
    new_array = np.average(my_array, axis=2)
    print(new_array.shape)


.. parsed-literal::

    (10, 10)




mkgauss
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The mkgauss funtionality has been replicated in the Photutils package
with
`photutils.datasets.make\_gaussian\_sources <http://photutils.readthedocs.io/en/stable/api/photutils.datasets.make_gaussian_sources.html>`__.



pixlocate
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Pixlocate is used to print positions matching a certain value condition.
This is replicated with the
`numpy.where <https://docs.scipy.org/doc/numpy/reference/generated/numpy.where.html>`__
function.



rd2xy-xy2rd
-----------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Rd2xy and xy2rd are used to translate RA/Dec to the pixel coordinate and
vice-versa. This capability is well covered in the ``Astropy.wcs``
package. Please see the
`documentation <http://docs.astropy.org/en/stable/wcs/>`__ for more
details on usage.





Not Replacing
-------------

-  boxinterp - Fill areas with smoothed values from surrounding area.
   See **images.imfit** notebook.
-  countfiles - Count how many files are in the input file template.
   Deprecated.
-  gcombine - Combine a set of GEIS images into one image. Deprecated,
   for FITS see **stsdas.toolbox.imgtools.mstools.mscombine**
-  gcopy - Generic multi-group copy utility. GEIS, deprecated.
-  gstatistics - Compute and print image pixel statistics for all
   groups. GEIS, deprecated. For FITS see **images.imutil.imstatistics**
-  imcalc - Perform general arithmetic operations on images. See
   **images.imtuil.imarith**.
-  imfill - Set fill value in image according to a mask. See
   **images.imutil.imreplace**.
-  listarea - Print an area of an image. See `numpy basics
   documentation <https://docs.scipy.org/doc/numpy-dev/user/quickstart.html>`__.
-  moveheader - Combine the header and pixels from two images. GEIS,
   deprecated.
-  pickfile - Get the file name picked from the input file template.
   Deprecated.
-  pixedit - Screen editor for image pixels. See **images.tv.imedit**
-  rbinary - Create an image from a binary file. Deprecated.
-  stack - Stack images to form a new image with one more dimension. See
   **images.imutil.imstack**
-  xyztable - Interpolate table values, writing results to a table. See
   **images.imfit.imsurfit** and **tables.ttools.tcopy-tdump**
-  xyztoim - Interpolate table values, writing results to an image. See
   **images.imfit.imsurfit**, `Astropy Tables
   documentation <http://docs.astropy.org/en/stable/table/>`__, and
   **tables.ttools.tcopy-tdump**.
