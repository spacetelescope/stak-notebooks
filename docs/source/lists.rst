:orphan:


lists
=====

List processing package.

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

The simple tasks in this package are already covered by ``Astropy`` or
other Python built-ins

Contents:

-  `lintran <#lintran>`__
-  `raverage <#raverage>`__
-  `table-unique <#table-unique>`__



lintran
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The lintran task will linear transform a list of coordinates. There are
various tasks in the `numpy linear algebra
package <https://docs.scipy.org/doc/numpy/reference/routines.linalg.html>`__
that can be used to achieve this in Python. Below is a simple example of
a 90 degree rotation.

.. code:: ipython3

    # Standard Imports
    import numpy as np

.. code:: ipython3

    # 90 degree rotation counterclockwise
    x_points = np.array([1,4,5,7,3])
    y_points = np.array([4,5,2,2,4])
    points = np.vstack([x_points,y_points])
    transform_matrix = a = np.array([[0, -1],
                                     [1, 0]])
    
    new_points = np.linalg.inv(a).dot(points)
    
    print("new x values: {}".format(new_points[0]))
    print("new y values: {}".format(new_points[1]))


.. parsed-literal::

    new x values: [ 4.  5.  2.  2.  4.]
    new y values: [-1. -4. -5. -7. -3.]




raverage
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

raverage computes the running average and standard deviation for a
1-dimensional array. We can efficiently do this in Python using numpy's
`convolve
function <https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.convolve.html>`__.
For a more complex but more efficient implementatin of a rolling average
that will also give you an output standard deviation you can use strides
in ``numpy`` as seen in `this example
here <http://www.rigtorp.se/2011/01/01/rolling-statistics-numpy.html>`__.

.. code:: ipython3

    # Standard Imports
    import numpy as np

.. code:: ipython3

    # setup test data
    data = np.arange(20)
    # kernel size
    N = 4
    
    out_data = np.convolve(data, np.ones((N,))/N, mode='valid')
    
    print("original array: {}".format(data))
    print("out_data: {}".format(out_data))


.. parsed-literal::

    original array: [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19]
    out_data: [  1.5   2.5   3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5  11.5  12.5
      13.5  14.5  15.5  16.5  17.5]




table-unique
------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The table task a list of text and transfers it to a table. We will show
an example of this using `Astropy
Tables <http://docs.astropy.org/en/stable/table/>`__ and a text file
with return seperated values. There are various parameters for reading
in ascii files, documentation `found
here <http://docs.astropy.org/en/stable/io/ascii/read.html#io-ascii-read-parameters>`__.
We will continue this example to show how to extract that list and then
apply a unique requirement, and save the list back out to a file, as in
the IRAF task unique. More information on going to and from ``Astropy``
``Tables`` see:

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy.table import Table, unique

.. code:: ipython3

    # Read text file into table
    text_file = '/eng/ssb/iraf_transition/test_data/table.txt'
    tab = Table.read(text_file, format='ascii.no_header')
    tab.pprint()


.. parsed-literal::

     col1
    -----
    star1
    star2
    star3
    star3
    star4
    star5
    star5
    star6


.. code:: ipython3

    # Run unique and print table
    out_table = unique(tab)
    out_table.pprint()
    
    # Save results to new text file
    out_table.write("out_table.txt",format='ascii')


.. parsed-literal::

     col1
    -----
    star1
    star2
    star3
    star4
    star5
    star6




Not Replacing
-------------

-  average - see **images.imutil.imstatistics** or `numpy
   tools <https://docs.scipy.org/doc/numpy/reference/routines.statistics.html>`__
-  columns - used to convert multicolumn data into CL lists, deprecated
-  rgcursor - Read the graphics cursor, deprecated
-  rimcursor - Read the image display cursor, see **images.tv** or
   `Python imexam <http://imexam.readthedocs.io/en/v0.7.1/>`__
-  tokens - deprecated
-  words - deprecated, see Python's built in `file
   reader <https://docs.python.org/3.6/tutorial/inputoutput.html>`__
