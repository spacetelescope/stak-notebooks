:orphan:


Python Introduction
===================

This notebook will introduce new Python and Astropy users to various
aspects of the Python astronomy workflow that may help them to take
better advantage of the rest of the STAK notebooks. If you are a new
Python or Astropy user, we highly recommend reading through all sections
of this notebook. We will discuss:

-  `Workflow Philosophy, the difference in workflow between IRAF and
   Python <#workflow>`__.
-  `FITS File I/O, opening and updating FITS files in Python <#fits>`__.
-  `IRAF/IDL to Python gotchas <#gotchas>`__.
-  `Helpful links and additional resources <#links>`__.

Workflow Philosophy
-------------------

The traditional IRAF workflow consists of individual tasks that take a
FITS file as input, run a data processing algorithm on the input file,
and return an updated output file for the user. The workflow in Python
is very different. Like most interpreted languages, data manipulation is
more a flow of operations expressed in regular code than individual
tasks that operate on files.

For example, instead of providing an input FITS file to a task and
getting an updated output file, for many astronomy Python functions the
user will open the FITS file themselves, access the file object, and
continue with any analysis or data processing from there. Once the data
manipulation is done (which could be a long sequence of manipulations),
if the user would like to update their FITS files they actively save the
changes to the file. This could include saving out intermediate steps to
new files, if the user so desires. While this may take a few more lines
of code to open and close files, you gain the major advantage of
flexibility, customization and freedom, as well as potentially more
efficient I/O. Oftentimes the increase in the number of lines of code
allows the user to be more explicit about what the code is doing,
instead of relying on complex parameter files and black box code. If a
user finds themselves running the same code numerous time and would find
a one line call useful, they can always put their workflow into a
function in a Python file and it will be ready to import and use in any
local Python session.

Particularly useful for this kind of scientific workflow is the `Jupyter
Notebook <https://jupyter-notebook.readthedocs.io/en/stable/>`__. This
tool allows users to run fully interactive Python sessions with the
ability to save the code that was run, its screen output, the order it
was run in, and other helpful tools. If you are new to Python, we highly
recommend taking advantage of this tool. All STAK notebooks were written
in Jupyter Notebooks. This means the user can downloads these notebooks
and run all code examples locally, adding additional code to use their
data in the provided examples. For a more in-depth example of how to
effectively use Python and the Jupyter Notebooks for data analysis,
there are some good video tutorials by `Jake VanderPlas
here <https://jakevdp.github.io/blog/2017/03/03/reproducible-data-analysis-in-jupyter/>`__.

FITS File I/O
-------------

Utilities to open, edit, and save FITS files are contained in the
`Astropy package <http://docs.astropy.org/en/stable/>`__. Although there
are many documentation resources for this package (which we will
reference below in the `links and resources section <#links>`__), it can
be a bit overwhelming for new users. To help new users adjust to a more
Pythonic workflow, some of the examples shown in the STAK notebooks do
not show the opening and closing of the FITS file. Instead, we show the
basic calls to do this below.

Before we show any code, we will go over a basic explanation of the
Astropy classes and object types for FITS files. The highest level class
for a FITS file is the
`HDUList <http://docs.astropy.org/en/stable/io/fits/api/hdulists.html>`__.
This is the object type you get when you first open a FITS file in
Astropy. The HDUList contains a sequence of HDU (Header Data Unit)
objects. There are three common flavors of HDU object, a
`PrimaryHDU <http://docs.astropy.org/en/stable/io/fits/api/hdus.html#>`__,
an
`ImageHDU <http://docs.astropy.org/en/stable/io/fits/api/images.html#astropy.io.fits.ImageHDU>`__
and a
`BinTableHDU <http://docs.astropy.org/en/stable/io/fits/api/tables.html#astropy.io.fits.BinTableHDU>`__.
As you would expect, the PrimaryHDU object contains the primary header,
the ImageHDU object contains an image array and its associated header,
and a BinTableHDU object contains a binary table and its associated
header.

The next level of objects are the individual header and data objects.
Headers are stored in a
`Header <http://docs.astropy.org/en/stable/io/fits/api/headers.html>`__
object. Images are stored in a `Numpy
array <https://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html>`__
object. Binary tables are stored in a `record
array <http://docs.astropy.org/en/stable/io/fits/usage/table.html>`__.

Now that we've covered the object types, we will show with code snippets
how to extract and work with these objects using a FITS file containing
image data. For more information on table data see the `Astropy FITS
Table
tutorial <http://www.astropy.org/astropy-tutorials/FITS-tables.html>`__.
For more detailed coverage of the available tools in the Astropy FITS
module see the `narrative documentation
page <http://docs.astropy.org/en/stable/io/fits/>`__.

.. code:: ipython3

    from astropy.io import fits

.. code:: ipython3

    # To pull our sample data file we will use an astroquery call
    from astroquery.mast import Observations
    obsid = '2004615006'
    Observations.download_products(obsid,productFilename="iczgs3ygq_flt.fits")


.. parsed-literal::

    INFO: Found cached file ./mastDownload/HST/iczgs3ygq/iczgs3ygq_flt.fits with expected size 16534080. [astroquery.query]




.. raw:: html

    <i>Table length=1</i>
    <table id="table90371172504" class="table-striped table-bordered table-condensed">
    <thead><tr><th>Local Path</th><th>Status</th><th>Message</th><th>URL</th></tr></thead>
    <thead><tr><th>str47</th><th>str8</th><th>object</th><th>object</th></tr></thead>
    <tr><td>./mastDownload/HST/iczgs3ygq/iczgs3ygq_flt.fits</td><td>COMPLETE</td><td>None</td><td>None</td></tr>
    </table>



.. code:: ipython3

    # Assign filename to the variable test_file
    test_file = './mastDownload/HST/ICZGS3YGQ/iczgs3ygq_flt.fits'

Below we will open the FITS file. You can open the file in various
modes, for this example we will open in update mode. The default mode is
read only.

.. code:: ipython3

    # Open the FITS file with Astropy
    HDUList_object = fits.open(test_file, mode='update')

Next we will show the info print out for this HDUList object using the
info() method. Notice the No. and Type columns. These will be useful for
indexing the HDUList.

.. code:: ipython3

    # HDUList info call
    HDUList_object.info()


.. parsed-literal::

    Filename: ./mastDownload/HST/ICZGS3YGQ/iczgs3ygq_flt.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     265   ()      
      1  SCI           1 ImageHDU       140   (1014, 1014)   float32   
      2  ERR           1 ImageHDU        51   (1014, 1014)   float32   
      3  DQ            1 ImageHDU        43   (1014, 1014)   int16   
      4  SAMP          1 ImageHDU        37   (1014, 1014)   int16   
      5  TIME          1 ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR       1 BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   


Now we will extract the primary header into the variable
``primary_header``

.. code:: ipython3

    # Extract primary header
    primary_header = HDUList_object[0].header
    
    # Index header object with keyword name and print value
    print(primary_header['FILENAME'])


.. parsed-literal::

    iczgs3ygq_flt.fits


Next we extract the image data into a variable called ``image_data``
from the first image extension. We will index this using the index
number from the No. column returned by ``info()``. This variable is a
numpy array object and the object that allows you to directly interact
with the image data. For more information on indexing here is a useful
`Numpy documentation
page <https://docs.scipy.org/doc/numpy/user/basics.indexing.html>`__.

.. code:: ipython3

    # Extract image data from the first extension
    image_data = HDUList_object[1].data
    print(image_data)


.. parsed-literal::

    [[  0.88747919   0.83535039   0.80967814 ...,   3.19881892   4.66315889
       12.94333744]
     [  0.94745165   0.80834782   0.76161045 ...,   0.91167408   3.91721344
        2.38371158]
     [  0.86024958   0.86270761   0.85969168 ...,   2.71301699  -4.11855459
        2.52296972]
     ..., 
     [ 33.32045746  23.79335022   4.87152386 ...,  22.54588509  21.88571739
       23.2428627 ]
     [ 47.97618103   1.16626728  13.08955574 ...,  12.46915627  21.59257698
       16.61116219]
     [ 30.99951744  29.15618515  46.40042877 ...,   0.           9.47169876
       20.67056084]]


We now have two options for saving out the FITS information.

We can save it out to the original file by using our ``HDUList`` file
object and the ``close`` argument. If the file was opened using the
update mode, this will flush (write) the file changes. If the file was
opened in the default readonly mode, it will **not** be updated when
closed.

We can also use the ``writeto`` method to save the ``HDUObject`` to a
new file. ``writeto`` will close the new file for you.

The ``flush`` method will also save to the original file. In this case,
the original file handling object will still need to be closed at some
point in the session.

**No matter which mode you used to open a FITS file, you should still
call the close method to close any open FITS file.**

.. code:: ipython3

    # Save using the writeto method to a new file, writeto will close the new file for you
    HDUList_object.writeto("wfc3data_new.fits")
    
    # Save using the writeto method, overwriting the original file
    HDUList_object.flush()

.. code:: ipython3

    # Save to same file using close
    # We show this last because we need to close the original copy of the file we opened, even after using a writeto
    HDUList_object.close()

IRAF/IDL to Python Gotchas
--------------------------

There are some important differences in syntax between IRAF, IDL, and
Python that new users should keep in mind. For more in depth information
about indexing and slicing ``Numpy`` arrays see `their indexing
documentation
here <https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.indexing.html>`__.

x versus y
~~~~~~~~~~

When working with images (2-dimensional) arrays, IRAF and IDL both have
the index order ``[x, y]``. In Python's ``Numpy`` package, the order is
reversed, ``[y, x]``.

index 0
~~~~~~~

IRAF indexes begin at 1 whereas Python and IDL both index arrays
starting at zero. So to pull out the first element of a 1-dimensional
array you would use ``array[0]``. To pull out the lower left corner of a
2-dimensional array you would use ``array[0,0]``.

slicing
~~~~~~~

Slicing in IRAF and IDL is inclusive for the right side of the slice. In
Python the right side of the slice is exclusive. For example, if you end
a slice with the 4th index, ``array[0:4]``, the fourth index element
(actually the 5th element in the array since index begins at 0) will
**not** be included in the slice.

matplotlib origin
~~~~~~~~~~~~~~~~~

The default origin location for ``matplotlib`` plots (a common Python
plotting library) will be in the upper-left. To change this to the lower
left (common for images) you can use the ``origin=lower`` parameter in
the ``imshow`` call as follows: ``plt.imshow(..., origin='lower')``.

close your files!
~~~~~~~~~~~~~~~~~

Almost any file handling object you open in Python (and this included a
FITS file opened with the Astropy open function!) will need to be closed
in your Python session with the appropriate close command. See above
section for examples.

Links and Resources
-------------------

Astropy
~~~~~~~

**Main user documentation page:** http://docs.astropy.org/en/stable/

**Main FITS page:** http://docs.astropy.org/en/stable/io/fits/index.html

**Tutorials:** http://www.astropy.org/astropy-tutorials/

Scipy and Numpy
~~~~~~~~~~~~~~~

**Main documentation pages:** https://docs.scipy.org/doc/

**Numpy indexing guide:**
https://docs.scipy.org/doc/numpy-1.13.0/reference/arrays.indexing.html

CCDProc
~~~~~~~

**Main documentation page:** http://ccdproc.readthedocs.io/en/latest/

Matplotlib
~~~~~~~~~~

**Main documentation page:** https://matplotlib.org/

Ginga
~~~~~

**Main documentation page:** http://ginga.readthedocs.io/en/latest/

Astroconda
~~~~~~~~~~

**Main documentation page:** http://astroconda.readthedocs.io/en/latest/
