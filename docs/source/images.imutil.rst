:orphan:


images.imutil
=============

The images.imutil package provides general FITS image tools such as
header editing and image arithmetic.

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

Contents:

-  `chpixtype <#chpixtype>`__
-  `hedit <#hedit>`__
-  `hselect <#hselect>`__
-  `imarith-imdivide <#imarith-imdivide>`__
-  `imcopy <#imcopy>`__
-  `imfunction-imexpr <#imfunction-imexpr>`__
-  `imheader <#imheader>`__
-  `imhistogram <#imhistogram>`__
-  `imreplace <#imreplace>`__
-  `imslice <#imslice>`__
-  `imstack <#imstack>`__
-  `imstatistics <#imstatistics>`__
-  `imsum <#imsum>`__
-  `listpixels <#listpixels>`__



chpixtype
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Chpixtype is a task that allows you to change the pixel type of a FITS
image. There is built in functionality in ``astropy.io.fits`` to preform
this task with the ``scale`` method. Below you will find a table that
translates the chpixtype newpixtype options into their equivalent
`numpy/astropy
type <http://docs.scipy.org/doc/numpy/user/basics.types.html>`__.

**Type Conversions**

+--------------+----------------------+
| Chpixtype    | Numpy/Astropy Type   |
+==============+======================+
| ``ushort``   | ``uint16``           |
+--------------+----------------------+
| ``short``    | ``int16``            |
+--------------+----------------------+
| ``int``      | ``int32``            |
+--------------+----------------------+
| ``long``     | ``int64``            |
+--------------+----------------------+
| ``real``     | ``float32``          |
+--------------+----------------------+
| ``double``   | ``float64``          |
+--------------+----------------------+

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change this value to your desired data file, here were creating a filename
    # for our new changed data
    orig_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    new_data = 'iczgs3ygq_newdtype_flt.fits'
    
    # Read in your FITS file
    hdu = fits.open(orig_data)
    
    # Print info about FITS file
    hdu.info()
    
    # Edit the datatype for the first sci extension
    hdu[1].scale(type='int32')
    
    # Save changed hdu object to new file
    # The overwrite argument tells the writeto method to overwrite if file already exists
    hdu.writeto(new_data, overwrite=True)
    hdu.close()


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
    <class 'astropy.io.fits.hdu.hdulist.HDUList'>




hedit
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The hedit task allows users to edit an image header. This functioanlity
is covered in ``astropy.io.fits``. Take note that to make changes to a
FITS file, you must use the ``mode='update'`` keyword in the
``fits.open`` call. The default mode for ``fits.open`` is ``readonly``.
Below you'll find examples of editing a keyword if it does/doesn't
exist, and how to delete keywords from the header. Also provided is an
example of updating multiple files at once using the `convience function
setval <http://docs.astropy.org/en/stable/io/fits/api/files.html#setval>`__.

For examples on printing/viewing header keywords please see
`hselect <#hselect>`__

.. code:: ipython3

    # Standard Imports
    from glob import glob
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change this value to your desired data file
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    
    # Open FITS file, include the mode='update' keyword
    hdu = fits.open(test_data, mode='update')
    
    # Simple header change, will add keyword if it doesn't exist
    hdu[0].header['MYKEY1'] = 'Editing this keyword'
    
    # Only add keyword if it does not already exist:
    if 'MYKEY2' not in hdu[0].header:
        hdu[0].header['MYKEY2'] = 'Also editing this'
    
    # To delete keywords, first check if they exist:
    if 'MYKEY2' in hdu[0].header:
        del hdu[0].header['MYKEY2']
        
    # Close FITS file, this will save your changes
    hdu.close()

Below we will show an example of how to update a keyword in multiple
FITS files using the Astropy convenience function
`astropy.io.fits.setval <http://docs.astropy.org/en/stable/io/fits/api/files.html#setval>`__
and the `glob <https://docs.python.org/3/library/glob.html>`__ function.
``Astropy.io.fits.setval`` will add the keyword if it does not already
exist.

.. code:: ipython3

    # Change this value to your desired search
    data_list = glob('/eng/ssb/iraf_transition/test_data/hedit/*.fits')
    
    # Now we loop over the list of file and use the setval function to update keywords
    # Here we update the keyword MYKEY1 value to the integer 5.
    for filename in data_list:
        fits.setval(filename, 'MYKEY1', value=5)



hselect
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The hselect task allows users to search for keyword values in the FITS
headers. This functionality has been replaced by the `CCDProc
ImageFileCollection
class <http://ccdproc.readthedocs.io/en/stable/api/ccdproc.ImageFileCollection.html>`__.
This class stores the header keyword values in an `Astropy Table
object <http://docs.astropy.org/en/stable/table/index.html#module-astropy.table>`__.
There is also an executable script provided by Astropy called
`fitsheader <http://docs.astropy.org/en/stable/io/fits/usage/scripts.html#module-astropy.io.fits.scripts.fitsheader>`__.
You'll find examples of both below.

If you wish to save your output to a text file, please see the `Astropy
Table Documentation <http://docs.astropy.org/en/stable/table/io.html>`__
and the `Astropy Unified I/O
page <http://docs.astropy.org/en/stable/io/unified.html>`__.

.. code:: ipython3

    # Astronomy Specific Imports
    from ccdproc import ImageFileCollection

.. code:: ipython3

    # first we make the ImageFileCollection object
    collec = ImageFileCollection('/eng/ssb/iraf_transition/test_data', 
                                 keywords=["filetype","date","exptime","filter"],
                                 glob_include="icz*.fits", ext=0)
    
    # header keywords values are stored in an Astropy Table in the summary attribute 
    out_table = collec.summary
    out_table




.. raw:: html

    &lt;Table masked=True length=3&gt;
    <table id="table4542238616" class="table-striped table-bordered table-condensed">
    <thead><tr><th>file</th><th>filetype</th><th>date</th><th>exptime</th><th>filter</th></tr></thead>
    <thead><tr><th>str27</th><th>str3</th><th>str10</th><th>float64</th><th>str5</th></tr></thead>
    <tr><td>iczgs3y5q_flt.fits</td><td>SCI</td><td>2016-06-02</td><td>652.937744</td><td>F125W</td></tr>
    <tr><td>iczgs3ygq_flt.fits</td><td>SCI</td><td>2016-06-02</td><td>602.937317</td><td>F140W</td></tr>
    <tr><td>iczgs3ygq_newdtype_flt.fits</td><td>SCI</td><td>2016-06-02</td><td>602.937317</td><td>F140W</td></tr>
    </table>



.. code:: ipython3

    # Now we can filter our table based on keyword values using Python bitwise operators
    filtered_table = out_table[(out_table['exptime'] > 602) & (out_table['filter'] == 'F140W')]
    filtered_table




.. raw:: html

    &lt;Table masked=True length=2&gt;
    <table id="table4542100368" class="table-striped table-bordered table-condensed">
    <thead><tr><th>file</th><th>filetype</th><th>date</th><th>exptime</th><th>filter</th></tr></thead>
    <thead><tr><th>str27</th><th>str3</th><th>str10</th><th>float64</th><th>str5</th></tr></thead>
    <tr><td>iczgs3ygq_flt.fits</td><td>SCI</td><td>2016-06-02</td><td>602.937317</td><td>F140W</td></tr>
    <tr><td>iczgs3ygq_newdtype_flt.fits</td><td>SCI</td><td>2016-06-02</td><td>602.937317</td><td>F140W</td></tr>
    </table>



.. code:: ipython3

    # Now let's extract the filename list from our filtered table into a python List object
    filelist = filtered_table['file'].data
    print(filelist)
    
    for filename in filelist:
        print(filename)
        # Do your analysis here


.. parsed-literal::

    ['iczgs3ygq_flt.fits' 'iczgs3ygq_newdtype_flt.fits']
    iczgs3ygq_flt.fits
    iczgs3ygq_newdtype_flt.fits




Also available is the Astropy executable script fitsheader. Fitsheader
can be run from the command line.

.. code:: ipython3

    # the "!" character tells the notebook to run this command as if it were in a terminal window
    !fitsheader --help


.. parsed-literal::

    usage: fitsheader [-h] [-e HDU] [-k KEYWORD] [-t [FORMAT]] [-c]
                      filename [filename ...]
    
    Print the header(s) of a FITS file. Optional arguments allow the desired
    extension(s), keyword(s), and output format to be specified. Note that in the
    case of a compressed image, the decompressed header is shown by default.
    
    positional arguments:
      filename              path to one or more files; wildcards are supported
    
    optional arguments:
      -h, --help            show this help message and exit
      -e HDU, --extension HDU
                            specify the extension by name or number; this argument
                            can be repeated to select multiple extensions
      -k KEYWORD, --keyword KEYWORD
                            specify a keyword; this argument can be repeated to
                            select multiple keywords; also supports wildcards
      -t [FORMAT], --table [FORMAT]
                            print the header(s) in machine-readable table format;
                            the default format is "ascii.fixed_width" (can be
                            "ascii.csv", "ascii.html", "ascii.latex", "fits", etc)
      -c, --compressed      for compressed image data, show the true header which
                            describes the compression rather than the data


.. code:: ipython3

    # print out only the keyword names that match FILE* or NAXIS*
    !fitsheader --keyword FILE* --keyword NAXIS* /eng/ssb/iraf_transition/test_data/hedit/*.fits


.. parsed-literal::

    # HDU 0 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    FILENAME= 'jczgx1ppq_flc.fits' / name of file                                   
    FILETYPE= 'SCI      '          / type of data found in data file                
    NAXIS   =                    0                                                  
    
    # HDU 1 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 2 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 3 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 4 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 5 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 6 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 7 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 8 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 9 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 10 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 11 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 12 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 13 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 14 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 15 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                  455 / length of dimension 1                          
    NAXIS2  =                   14 / length of dimension 2                          
    # HDU 0 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    FILENAME= 'jczgx1q1q_flc.fits' / name of file                                   
    FILETYPE= 'SCI      '          / type of data found in data file                
    NAXIS   =                    0                                                  
    
    # HDU 1 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 2 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 3 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 4 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 5 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 6 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2                                                  
    NAXIS1  =                 4096                                                  
    NAXIS2  =                 2048                                                  
    
    # HDU 7 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 8 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 9 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 10 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 11 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 12 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 13 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 14 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                   64                                                  
    NAXIS2  =                   32                                                  
    
    # HDU 15 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    NAXIS   =                    2 / number of array dimensions                     
    NAXIS1  =                  455 / length of dimension 1                          
    NAXIS2  =                   14 / length of dimension 2                          


.. code:: ipython3

    # print out only the first extension and keyword names that match FILE* or NAXIS*
    !fitsheader --extension 0 --keyword FILE* --keyword NAXIS* /eng/ssb/iraf_transition/test_data/hedit/*.fits


.. parsed-literal::

    # HDU 0 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1ppq_flc.fits:
    FILENAME= 'jczgx1ppq_flc.fits' / name of file                                   
    FILETYPE= 'SCI      '          / type of data found in data file                
    NAXIS   =                    0                                                  
    # HDU 0 in /eng/ssb/iraf_transition/test_data/hedit/jczgx1q1q_flc.fits:
    FILENAME= 'jczgx1q1q_flc.fits' / name of file                                   
    FILETYPE= 'SCI      '          / type of data found in data file                
    NAXIS   =                    0                                                  




imarith-imdivide
----------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Imarith and imdivide both provide functionality to apply basic operators
to whole image arrays. This task can be achieved with basic
``astropy.io.fits`` functionality along with ``numpy`` array
functionality. We show a few examples below. In the first code cell we
adding and dividing two image arrays together. In the second code cell
we show how to use a data quality array to decide which image array
values to replace with zero.

The basic operands (``+``,\ ``-``,\ ``/``,\ ``*``) can all be used with
an assignment operator in python (``+=``,\ ``-=``,\ ``/=``,\ ``*=``).
See http://www.tutorialspoint.com/python/python\_basic\_operators.htm
for more details

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Basic operands (+,-,/,*)
    # Change these values to your desired data files
    test_data1 = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    test_data2 = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    output_data = 'imarith_out.fits'
    output_data2 = 'imarith_new.fits'
    
    
    # Open FITS file
    hdu1 = fits.open(test_data1)
    hdu2 = fits.open(test_data2)
    
    # Print information about the FITS file we opened
    hdu1.info()
    hdu2.info()
    
    # Here we add hdu2-ext1 to hdu1-ext1 by using the shortcut += operator
    hdu1[1].data += hdu2[1].data
    
    # If you are dividing and need to avoid zeros in the image use indexing
    indx_zeros = hdu2[1].data == 0
    indx_nonzeros = hdu2[1].data != 0
    
    # Set this value as you would the divzero parameter in imarith
    # Here we're working with the error arrays of the image
    set_zeros = 999.9
    hdu1[2].data[indx_nonzeros] /= hdu2[2].data[indx_nonzeros]
    hdu1[2].data[indx_zeros] = 999.9
    
    # Save your new file
    # The overwrite argument tells the writeto method to overwrite if file already exists
    hdu1.writeto(output_data, overwrite=True)
    
    # If you want to save you updated array to a new file with just the updated image array 
    # we can repackage the extension into a new HDUList
    image_array = hdu1[1].data
    new_hdu = fits.PrimaryHDU(image_array)
    new_hdu.writeto(output_data2, overwrite=True)
    
    # Close hdu files
    hdu1.close()
    hdu2.close()


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
    Filename: /eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     265   ()      
      1  SCI           1 ImageHDU       140   (1014, 1014)   float32   
      2  ERR           1 ImageHDU        51   (1014, 1014)   float32   
      3  DQ            1 ImageHDU        43   (1014, 1014)   int16   
      4  SAMP          1 ImageHDU        37   (1014, 1014)   int16   
      5  TIME          1 ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR       1 BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   


.. code:: ipython3

    # Here we show an example of using an HST DQ array to
    # replace only certain values with zero in an image array
    
    # Change these values to your desired data files
    test_data1 = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    output_file = 'iczgs3ygq_updated.fits'
    
    # Open FITS file
    hdulist = fits.open(test_data1)
    
    # First we should use the DQ array to make a boolean mask
    DQ_mask = hdulist[3].data > 16384
    
    # Now we can use the mask to replace values in the image array
    # with 0.
    hdulist[1].data[DQ_mask] = 0
    
    # Now we can save out the edited FITS to a new file
    hdulist.writeto(output_file)
    
    # And finally, close the original FITS file
    # The orignially file will not be updated since we did not
    # open the file in 'update' mode
    hdulist.close()



imcopy
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Imcopy allows users to copy a FITS image to a new file. We can
accomplish this using ``astropy.io.fits`` by saving our FITS file to a
new filename.

Imcopy will also make a cutout of an image and save the cutout to a new
file with an updated WCS. We show an exampe of this in Python using the
`Cutout2D <http://docs.astropy.org/en/stable/api/astropy.nddata.Cutout2D.html>`__
tool in ``Astropy``. For more information on how to use ``Cutout2D``
please see `this tutorial
page <http://docs.astropy.org/en/stable/nddata/utils.html#cutout-images>`__.

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy import wcs
    from astropy.io import fits
    from astropy.nddata import Cutout2D

Simple example of a file copy

.. code:: ipython3

    # Change these values to your desired filenames
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    output_data = 'imcopy_out.fits'
    
    hdulist = fits.open(test_data)
    # The overwrite argument tells the writeto method to overwrite if file already exists
    hdulist.writeto(output_data, overwrite=True)
    hdulist.close()

Example using a new cutout, here we will take a 50x50 pixel cutout from
all image extensions centered at x:200, y:300

.. code:: ipython3

    # Change these values to your desired filenames
    test_data = '/eng/ssb/iraf_transition/test_data/jcw505010_drz.fits'
    output_data = 'imcopy_cutout_out.fits'
    
    hdulist = fits.open(test_data)
    
    # Create iterable list of tuples to feed into Cutout2D, 
    # seperate list for extensions with wcs, as feeding the wcs 
    # back into the FITS file takes more work.
    ext_list = [1,2]
    for ext in ext_list:
        orig_wcs = wcs.WCS(hdulist[ext].header)
        cutout = Cutout2D(hdulist[ext].data, (200,300), (50,50), wcs=orig_wcs)
        hdulist[ext].data = cutout.data
        hdulist[ext].header.update(cutout.wcs.to_header())
        
    hdulist.writeto(output_data, overwrite=True)
    
    hdulist.close()



imfunction-imexpr
-----------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Imfunction will apply a function to the image pixel values in an image
array. Imexpr gives you similiar functionality with the added capability
to combine different images using a user created expression. We can
accomplish this using the built in funcitonality of the `numpy
library <http://docs.scipy.org/doc/numpy/reference/routines.math.html>`__.

If there is a particular function you would like to apply to your image
array that you cannot find in the ``numpy`` library you can use the
``np.vectorize`` function, which can make any python function apply to
each element of your array. But keep in mind that
`np.vectorize <http://docs.scipy.org/doc/numpy/reference/generated/numpy.vectorize.html>`__
is esentially looping over the array, and may not be the most efficient
method.

Example using exsisting numpy function:

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    output_data = 'imfunction_out.fits'
    
    # Here we use the cosine function as an example
    hdu = fits.open(test_data)
    sci = hdu[1].data
    
    # When you call your new function, make sure to reassign the array to
    # the new values if the original function is not changing values in place
    hdu[1].data = np.cos(hdu[1].data)
    
    # Now save out to a new file, and close the original file, changes will
    # not be applied to the oiginal FITS file.
    hdu.writeto(output_data, overwrite=True)
    hdu.close()

Example using user defined function and ``np.vectorize``:

.. code:: ipython3

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    output_data = 'imfunction2_out.fits'
    
    # Here we use the following custom function as an example
    def my_func(x):
        return (x**2)+(x**3)
    
    # Now we open our file, and vectorize our function
    hdu = fits.open(test_data)
    sci = hdu[1].data
    vector_func = np.vectorize(my_func)
    
    # When you call your new function, make sure to reassign the array to
    # the new values if the original function is not changing values in place
    hdu[1].data = vector_func(hdu[1].data)
    
    # Now save out to a new file, and close the original file, changes will
    # not be applied to the oiginal FITS file.
    hdu.writeto(output_data, overwrite=True)
    hdu.close()



imheader
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The imheader task allows the user to list header parameters for a list
of images. Here we can use the ``astropy`` convenience function,
``fits.getheader()``. We also show in this example how to save a header
to a text file, see the `Python file I/O
documentation <https://docs.python.org/3/tutorial/inputoutput.html>`__
for more details.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    import glob
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change these values to your desired data files, glob will capture all wildcard matches
    test_data = glob.glob('/eng/ssb/iraf_transition/test_data/iczgs3y*')
    out_text = 'imheader_out.txt'
    
    for filename in test_data:
        # Pull the header from extension 1 using FITS convenience function.
        # To access multiple header it's better to use the fits.open() function.
        head = fits.getheader(filename, ext=1)
        
        # Using repr function to format output
        print(repr(head))
        
        # Save header to text file
        with open(out_text, mode='a') as out_file:
            out_file.write(repr(head))
            out_file.write('\n\n')


.. parsed-literal::

    XTENSION= 'IMAGE   '           / IMAGE extension                                
    BITPIX  =                  -32                                                  
    NAXIS   =                    2                                                  
    NAXIS1  =                 1014                                                  
    NAXIS2  =                 1014                                                  
    PCOUNT  =                    0 / required keyword; must = 0                     
    GCOUNT  =                    1 / required keyword; must = 1                     
    ORIGIN  = 'HSTIO/CFITSIO March 2010'                                            
    DATE    = '2016-06-02' / date this file was written (yyyy-mm-dd)                
    INHERIT =                    T / inherit the primary header                     
    EXTNAME = 'SCI     '           / extension name                                 
    EXTVER  =                    1 / extension version number                       
    ROOTNAME= 'iczgs3ygq                         ' / rootname of the observation set
    EXPNAME = 'iczgs3ygq                ' / exposure identifier                     
    BUNIT   = 'ELECTRONS/S'        / brightness units                               
                                                                                    
                  / World Coordinate System and Related Parameters                  
                                                                                    
    WCSAXES =                    2 / number of World Coordinate System axes         
    CRPIX1  =                507.0 / x-coordinate of reference pixel                
    CRPIX2  =                507.0 / y-coordinate of reference pixel                
    CRVAL1  =       36.85374208875 / first axis value at reference pixel            
    CRVAL2  =       48.92264646942 / second axis value at reference pixel           
    CTYPE1  = 'RA---TAN-SIP'       / the coordinate type for the first axis         
    CTYPE2  = 'DEC--TAN-SIP'       / the coordinate type for the second axis        
    CD1_1   = -3.1758778512629E-05 / partial of first axis coordinate w.r.t. x      
    CD1_2   = -1.8099259044494E-05 / partial of first axis coordinate w.r.t. y      
    CD2_1   = -2.0157648752092E-05 / partial of second axis coordinate w.r.t. x     
    CD2_2   = 2.83052387051731E-05 / partial of second axis coordinate w.r.t. y     
    LTV1    =        0.0000000E+00 / offset in X to subsection start                
    LTV2    =        0.0000000E+00 / offset in Y to subsection start                
    LTM1_1  =                  1.0 / reciprocal of sampling rate in X               
    LTM2_2  =                  1.0 / reciprocal of sampling rate in Y               
    PA_APER =              -32.556 / Position Angle of reference aperture center (de
    VAFACTOR=   9.999085821139E-01 / velocity aberration plate scale factor         
    ORIENTAT=              -32.556 / position angle of image y axis (deg. e of n)   
    RA_APER =   3.685374208875E+01 / RA of aperture reference position              
    DEC_APER=   4.892264646942E+01 / Declination of aperture reference position     
                                                                                    
                  / REPEATED EXPOSURES INFORMATION                                  
                                                                                    
    NCOMBINE=                    1 / number of image sets combined during CR rejecti
                                                                                    
                  / READOUT DEFINITION PARAMETERS                                   
                                                                                    
    CENTERA1=                  513 / subarray axis1 center pt in unbinned dect. pix 
    CENTERA2=                  513 / subarray axis2 center pt in unbinned dect. pix 
    SIZAXIS1=                 1024 / subarray axis1 size in unbinned detector pixels
    SIZAXIS2=                 1024 / subarray axis2 size in unbinned detector pixels
    BINAXIS1=                    1 / axis1 data bin size in unbinned detector pixels
    BINAXIS2=                    1 / axis2 data bin size in unbinned detector pixels
                                                                                    
                  / READOUT PARAMETERS                                              
                                                                                    
    SAMPNUM =                   13 / MULTIACCUM sample number                       
    SAMPTIME=           602.937317 / total integration time (sec)                   
    DELTATIM=            50.000412 / integration time of this sample (sec)          
    ROUTTIME=   5.740229030181E+04 / UT time of array readout (MJD)                 
    TDFTRANS=                    0 / number of TDF transitions during current sample
                                                                                    
                  / DATA PACKET INFORMATION                                         
                                                                                    
    FILLCNT =                    0 / number of segments containing fill             
    ERRCNT  =                    0 / number of segments containing errors           
    PODPSFF =                    F / podps fill present (T/F)                       
    STDCFFF =                    F / science telemetry fill data present (T=1/F=0)  
    STDCFFP = '0x5569'             / science telemetry fill pattern (hex)           
                                                                                    
                  / IMAGE STATISTICS AND DATA QUALITY FLAGS                         
                                                                                    
    NGOODPIX=               990475 / number of good pixels                          
    SDQFLAGS=                31743 / serious data quality flags                     
    GOODMIN =       -2.8782272E+00 / minimum value of good pixels                   
    GOODMAX =        1.1788658E+04 / maximum value of good pixels                   
    GOODMEAN=        9.9831134E-01 / mean value of good pixels                      
    SNRMIN  =        1.8871337E-02 / minimum signal to noise of good pixels         
    SNRMAX  =        6.3982178E+01 / maximum signal to noise of good pixels         
    SNRMEAN =        5.3425826E-02 / mean value of signal to noise of good pixels   
    SOFTERRS=                    0 / number of soft error pixels (DQF=1)            
    MEANDARK=        1.2191877E+01 / average of the dark values subtracted          
    MEANBLEV=        1.4332316E+04 / average of all bias levels subtracted          
    RADESYS = 'ICRS    '                                                            
    OCX10   = 0.000786257500294596                                                  
    OCX11   =   0.1354287266731262                                                  
    OCY10   =   0.1209582984447479                                                  
    OCY11   = -0.00042557646520435                                                  
    IDCSCALE=   0.1282500028610229                                                  
    IDCTHETA=                 45.0                                                  
    IDCXREF =                507.0                                                  
    IDCYREF =                507.0                                                  
    IDCV2REF=    1.019000053405762                                                  
    IDCV3REF=  -0.5070000290870667                                                  
    WCSNAMEO= 'OPUS    '                                                            
    WCSAXESO=                    2                                                  
    CRPIX1O =                507.0                                                  
    CRPIX2O =                507.0                                                  
    CDELT1O =                  1.0                                                  
    CDELT2O =                  1.0                                                  
    CUNIT1O = 'deg     '                                                            
    CUNIT2O = 'deg     '                                                            
    CTYPE1O = 'RA---TAN'                                                            
    CTYPE2O = 'DEC--TAN'                                                            
    CRVAL1O =       36.85374208875                                                  
    CRVAL2O =       48.92264646942                                                  
    LONPOLEO=                180.0                                                  
    LATPOLEO=       48.92264646942                                                  
    RADESYSO= 'ICRS    '                                                            
    CD1_1O  =         -3.17711E-05                                                  
    CD1_2O  =         -1.80786E-05                                                  
    CD2_1O  =         -2.01487E-05                                                  
    CD2_2O  =          2.83166E-05                                                  
    IDCTAB  = 'iref$w3m18525i_idc.fits'                                             
    B_1_3   = 1.69983940010457E-13                                                  
    B_0_3   = -2.2777970488111E-10                                                  
    A_2_2   = 1.11275247848408E-13                                                  
    B_0_4   = 1.03978470894974E-12                                                  
    A_0_4   = -2.0083179974495E-13                                                  
    B_3_1   = 3.81044199963010E-13                                                  
    A_3_0   = -1.9851733613323E-10                                                  
    B_4_0   = -5.7352409055905E-13                                                  
    B_0_2   = 2.98815054868485E-05                                                  
    A_1_3   = 6.08832045645843E-13                                                  
    A_4_0   = -3.2156784473326E-13                                                  
    B_ORDER =                    4                                                  
    A_0_2   = 2.77482030873749E-08                                                  
    A_2_1   = 1.22255499299390E-10                                                  
    B_2_0   = 6.92276069494587E-06                                                  
    A_2_0   = -2.0701735553551E-07                                                  
    A_3_1   = 4.13947711822547E-13                                                  
    A_1_2   = 3.11477338242516E-11                                                  
    A_ORDER =                    4                                                  
    B_1_2   = 7.47270961118588E-11                                                  
    B_2_2   = 1.38557115814168E-13                                                  
    A_0_3   = 4.55691839657869E-11                                                  
    B_2_1   = -2.3836656728517E-10                                                  
    B_3_0   = 5.14014553890418E-11                                                  
    B_1_1   = -2.8538202053351E-07                                                  
    A_1_1   = 2.44176437155426E-05                                                  
    WCSNAME = 'IDC_w3m18525i'                                                       
    MDRIZSKY=   0.8125642368041847 / Sky value computed by AstroDrizzle             
    XTENSION= 'IMAGE   '           / IMAGE extension                                
    BITPIX  =                   32                                                  
    NAXIS   =                    2                                                  
    NAXIS1  =                 1014                                                  
    NAXIS2  =                 1014                                                  
    PCOUNT  =                    0 / required keyword; must = 0                     
    GCOUNT  =                    1 / required keyword; must = 1                     
    ORIGIN  = 'HSTIO/CFITSIO March 2010'                                            
    DATE    = '2016-06-02' / date this file was written (yyyy-mm-dd)                
    INHERIT =                    T / inherit the primary header                     
    EXTNAME = 'SCI     '           / extension name                                 
    EXTVER  =                    1 / extension version number                       
    ROOTNAME= 'iczgs3ygq                         ' / rootname of the observation set
    EXPNAME = 'iczgs3ygq                ' / exposure identifier                     
    BUNIT   = 'ELECTRONS/S'        / brightness units                               
                                                                                    
                  / World Coordinate System and Related Parameters                  
                                                                                    
    WCSAXES =                    2 / number of World Coordinate System axes         
    CRPIX1  =                507.0 / x-coordinate of reference pixel                
    CRPIX2  =                507.0 / y-coordinate of reference pixel                
    CRVAL1  =       36.85374208875 / first axis value at reference pixel            
    CRVAL2  =       48.92264646942 / second axis value at reference pixel           
    CTYPE1  = 'RA---TAN-SIP'       / the coordinate type for the first axis         
    CTYPE2  = 'DEC--TAN-SIP'       / the coordinate type for the second axis        
    CD1_1   = -3.1758778512629E-05 / partial of first axis coordinate w.r.t. x      
    CD1_2   = -1.8099259044494E-05 / partial of first axis coordinate w.r.t. y      
    CD2_1   = -2.0157648752092E-05 / partial of second axis coordinate w.r.t. x     
    CD2_2   = 2.83052387051731E-05 / partial of second axis coordinate w.r.t. y     
    LTV1    =        0.0000000E+00 / offset in X to subsection start                
    LTV2    =        0.0000000E+00 / offset in Y to subsection start                
    LTM1_1  =                  1.0 / reciprocal of sampling rate in X               
    LTM2_2  =                  1.0 / reciprocal of sampling rate in Y               
    PA_APER =              -32.556 / Position Angle of reference aperture center (de
    VAFACTOR=   9.999085821139E-01 / velocity aberration plate scale factor         
    ORIENTAT=              -32.556 / position angle of image y axis (deg. e of n)   
    RA_APER =   3.685374208875E+01 / RA of aperture reference position              
    DEC_APER=   4.892264646942E+01 / Declination of aperture reference position     
                                                                                    
                  / REPEATED EXPOSURES INFORMATION                                  
                                                                                    
    NCOMBINE=                    1 / number of image sets combined during CR rejecti
                                                                                    
                  / READOUT DEFINITION PARAMETERS                                   
                                                                                    
    CENTERA1=                  513 / subarray axis1 center pt in unbinned dect. pix 
    CENTERA2=                  513 / subarray axis2 center pt in unbinned dect. pix 
    SIZAXIS1=                 1024 / subarray axis1 size in unbinned detector pixels
    SIZAXIS2=                 1024 / subarray axis2 size in unbinned detector pixels
    BINAXIS1=                    1 / axis1 data bin size in unbinned detector pixels
    BINAXIS2=                    1 / axis2 data bin size in unbinned detector pixels
                                                                                    
                  / READOUT PARAMETERS                                              
                                                                                    
    SAMPNUM =                   13 / MULTIACCUM sample number                       
    SAMPTIME=           602.937317 / total integration time (sec)                   
    DELTATIM=            50.000412 / integration time of this sample (sec)          
    ROUTTIME=   5.740229030181E+04 / UT time of array readout (MJD)                 
    TDFTRANS=                    0 / number of TDF transitions during current sample
                                                                                    
                  / DATA PACKET INFORMATION                                         
                                                                                    
    FILLCNT =                    0 / number of segments containing fill             
    ERRCNT  =                    0 / number of segments containing errors           
    PODPSFF =                    F / podps fill present (T/F)                       
    STDCFFF =                    F / science telemetry fill data present (T=1/F=0)  
    STDCFFP = '0x5569'             / science telemetry fill pattern (hex)           
                                                                                    
                  / IMAGE STATISTICS AND DATA QUALITY FLAGS                         
                                                                                    
    NGOODPIX=               990475 / number of good pixels                          
    SDQFLAGS=                31743 / serious data quality flags                     
    GOODMIN =       -2.8782272E+00 / minimum value of good pixels                   
    GOODMAX =        1.1788658E+04 / maximum value of good pixels                   
    GOODMEAN=        9.9831134E-01 / mean value of good pixels                      
    SNRMIN  =        1.8871337E-02 / minimum signal to noise of good pixels         
    SNRMAX  =        6.3982178E+01 / maximum signal to noise of good pixels         
    SNRMEAN =        5.3425826E-02 / mean value of signal to noise of good pixels   
    SOFTERRS=                    0 / number of soft error pixels (DQF=1)            
    MEANDARK=        1.2191877E+01 / average of the dark values subtracted          
    MEANBLEV=        1.4332316E+04 / average of all bias levels subtracted          
    RADESYS = 'ICRS    '                                                            
    OCX10   = 0.000786257500294596                                                  
    OCX11   =   0.1354287266731262                                                  
    OCY10   =   0.1209582984447479                                                  
    OCY11   = -0.00042557646520435                                                  
    IDCSCALE=   0.1282500028610229                                                  
    IDCTHETA=                 45.0                                                  
    IDCXREF =                507.0                                                  
    IDCYREF =                507.0                                                  
    IDCV2REF=    1.019000053405762                                                  
    IDCV3REF=  -0.5070000290870667                                                  
    WCSNAMEO= 'OPUS    '                                                            
    WCSAXESO=                    2                                                  
    CRPIX1O =                507.0                                                  
    CRPIX2O =                507.0                                                  
    CDELT1O =                  1.0                                                  
    CDELT2O =                  1.0                                                  
    CUNIT1O = 'deg     '                                                            
    CUNIT2O = 'deg     '                                                            
    CTYPE1O = 'RA---TAN'                                                            
    CTYPE2O = 'DEC--TAN'                                                            
    CRVAL1O =       36.85374208875                                                  
    CRVAL2O =       48.92264646942                                                  
    LONPOLEO=                180.0                                                  
    LATPOLEO=       48.92264646942                                                  
    RADESYSO= 'ICRS    '                                                            
    CD1_1O  =         -3.17711E-05                                                  
    CD1_2O  =         -1.80786E-05                                                  
    CD2_1O  =         -2.01487E-05                                                  
    CD2_2O  =          2.83166E-05                                                  
    IDCTAB  = 'iref$w3m18525i_idc.fits'                                             
    B_1_3   = 1.69983940010457E-13                                                  
    B_0_3   = -2.2777970488111E-10                                                  
    A_2_2   = 1.11275247848408E-13                                                  
    B_0_4   = 1.03978470894974E-12                                                  
    A_0_4   = -2.0083179974495E-13                                                  
    B_3_1   = 3.81044199963010E-13                                                  
    A_3_0   = -1.9851733613323E-10                                                  
    B_4_0   = -5.7352409055905E-13                                                  
    B_0_2   = 2.98815054868485E-05                                                  
    A_1_3   = 6.08832045645843E-13                                                  
    A_4_0   = -3.2156784473326E-13                                                  
    B_ORDER =                    4                                                  
    A_0_2   = 2.77482030873749E-08                                                  
    A_2_1   = 1.22255499299390E-10                                                  
    B_2_0   = 6.92276069494587E-06                                                  
    A_2_0   = -2.0701735553551E-07                                                  
    A_3_1   = 4.13947711822547E-13                                                  
    A_1_2   = 3.11477338242516E-11                                                  
    A_ORDER =                    4                                                  
    B_1_2   = 7.47270961118588E-11                                                  
    B_2_2   = 1.38557115814168E-13                                                  
    A_0_3   = 4.55691839657869E-11                                                  
    B_2_1   = -2.3836656728517E-10                                                  
    B_3_0   = 5.14014553890418E-11                                                  
    B_1_1   = -2.8538202053351E-07                                                  
    A_1_1   = 2.44176437155426E-05                                                  
    WCSNAME = 'IDC_w3m18525i'                                                       
    MDRIZSKY=   0.8125642368041847 / Sky value computed by AstroDrizzle             
    XTENSION= 'IMAGE   '           / IMAGE extension                                
    BITPIX  =                  -32                                                  
    NAXIS   =                    2                                                  
    NAXIS1  =                 1014                                                  
    NAXIS2  =                 1014                                                  
    PCOUNT  =                    0 / required keyword; must = 0                     
    GCOUNT  =                    1 / required keyword; must = 1                     
    ORIGIN  = 'HSTIO/CFITSIO March 2010'                                            
    DATE    = '2016-06-02' / date this file was written (yyyy-mm-dd)                
    INHERIT =                    T / inherit the primary header                     
    EXTNAME = 'SCI     '           / extension name                                 
    EXTVER  =                    1 / extension version number                       
    ROOTNAME= 'iczgs3y5q                         ' / rootname of the observation set
    EXPNAME = 'iczgs3y5q                ' / exposure identifier                     
    BUNIT   = 'ELECTRONS/S'        / brightness units                               
                                                                                    
                  / World Coordinate System and Related Parameters                  
                                                                                    
    WCSAXES =                    2 / number of World Coordinate System axes         
    CRPIX1  =                507.0 / x-coordinate of reference pixel                
    CRPIX2  =                507.0 / y-coordinate of reference pixel                
    CRVAL1  =       36.85747964213 / first axis value at reference pixel            
    CRVAL2  =       48.92227663477 / second axis value at reference pixel           
    CTYPE1  = 'RA---TAN-SIP'       / the coordinate type for the first axis         
    CTYPE2  = 'DEC--TAN-SIP'       / the coordinate type for the second axis        
    CD1_1   = -3.1760811272930E-05 / partial of first axis coordinate w.r.t. x      
    CD1_2   = -1.8097365221752E-05 / partial of first axis coordinate w.r.t. y      
    CD2_1   = -2.0155198493371E-05 / partial of second axis coordinate w.r.t. x     
    CD2_2   = 2.83091348126201E-05 / partial of second axis coordinate w.r.t. y     
    LTV1    =        0.0000000E+00 / offset in X to subsection start                
    LTV2    =        0.0000000E+00 / offset in Y to subsection start                
    LTM1_1  =                  1.0 / reciprocal of sampling rate in X               
    LTM2_2  =                  1.0 / reciprocal of sampling rate in Y               
    PA_APER =             -32.5531 / Position Angle of reference aperture center (de
    VAFACTOR=   9.999381116940E-01 / velocity aberration plate scale factor         
    ORIENTAT=             -32.5531 / position angle of image y axis (deg. e of n)   
    RA_APER =   3.685747964213E+01 / RA of aperture reference position              
    DEC_APER=   4.892227663477E+01 / Declination of aperture reference position     
                                                                                    
                  / REPEATED EXPOSURES INFORMATION                                  
                                                                                    
    NCOMBINE=                    1 / number of image sets combined during CR rejecti
                                                                                    
                  / READOUT DEFINITION PARAMETERS                                   
                                                                                    
    CENTERA1=                  513 / subarray axis1 center pt in unbinned dect. pix 
    CENTERA2=                  513 / subarray axis2 center pt in unbinned dect. pix 
    SIZAXIS1=                 1024 / subarray axis1 size in unbinned detector pixels
    SIZAXIS2=                 1024 / subarray axis2 size in unbinned detector pixels
    BINAXIS1=                    1 / axis1 data bin size in unbinned detector pixels
    BINAXIS2=                    1 / axis2 data bin size in unbinned detector pixels
                                                                                    
                  / READOUT PARAMETERS                                              
                                                                                    
    SAMPNUM =                   14 / MULTIACCUM sample number                       
    SAMPTIME=           652.937744 / total integration time (sec)                   
    DELTATIM=            50.000412 / integration time of this sample (sec)          
    ROUTTIME=   5.740226431774E+04 / UT time of array readout (MJD)                 
    TDFTRANS=                    0 / number of TDF transitions during current sample
                                                                                    
                  / DATA PACKET INFORMATION                                         
                                                                                    
    FILLCNT =                    0 / number of segments containing fill             
    ERRCNT  =                    0 / number of segments containing errors           
    PODPSFF =                    F / podps fill present (T/F)                       
    STDCFFF =                    F / science telemetry fill data present (T=1/F=0)  
    STDCFFP = '0x5569'             / science telemetry fill pattern (hex)           
                                                                                    
                  / IMAGE STATISTICS AND DATA QUALITY FLAGS                         
                                                                                    
    NGOODPIX=               990476 / number of good pixels                          
    SDQFLAGS=                31743 / serious data quality flags                     
    GOODMIN =       -2.9155195E+00 / minimum value of good pixels                   
    GOODMAX =        2.6231844E+04 / maximum value of good pixels                   
    GOODMEAN=        9.3451303E-01 / mean value of good pixels                      
    SNRMIN  =        1.1295157E-02 / minimum signal to noise of good pixels         
    SNRMAX  =        9.8745354E+01 / maximum signal to noise of good pixels         
    SNRMEAN =        4.9034115E-02 / mean value of signal to noise of good pixels   
    SOFTERRS=                    0 / number of soft error pixels (DQF=1)            
    MEANDARK=        1.3298962E+01 / average of the dark values subtracted          
    MEANBLEV=        1.4334856E+04 / average of all bias levels subtracted          
    RADESYS = 'ICRS    '                                                            
    OCX10   = 0.000779107213020324                                                  
    OCX11   =   0.1354261934757233                                                  
    OCY10   =    0.120962917804718                                                  
    OCY11   = -0.00042105099419131                                                  
    IDCSCALE=   0.1282500028610229                                                  
    IDCTHETA=                 45.0                                                  
    IDCXREF =                507.0                                                  
    IDCYREF =                507.0                                                  
    IDCV2REF=    1.019000053405762                                                  
    IDCV3REF=  -0.5070000290870667                                                  
    WCSNAMEO= 'OPUS    '                                                            
    WCSAXESO=                    2                                                  
    CRPIX1O =                507.0                                                  
    CRPIX2O =                507.0                                                  
    CDELT1O =                  1.0                                                  
    CDELT2O =                  1.0                                                  
    CUNIT1O = 'deg     '                                                            
    CUNIT2O = 'deg     '                                                            
    CTYPE1O = 'RA---TAN'                                                            
    CTYPE2O = 'DEC--TAN'                                                            
    CRVAL1O =       36.85747964213                                                  
    CRVAL2O =       48.92227663477                                                  
    LONPOLEO=                180.0                                                  
    LATPOLEO=       48.92227663477                                                  
    RADESYSO= 'ICRS    '                                                            
    CD1_1O  =         -3.17721E-05                                                  
    CD1_2O  =         -1.80771E-05                                                  
    CD2_1O  =         -2.01471E-05                                                  
    CD2_2O  =          2.83175E-05                                                  
    IDCTAB  = 'iref$w3m18525i_idc.fits'                                             
    B_1_2   = 2.35150691092754E-11                                                  
    A_3_0   = -1.8769691205859E-10                                                  
    B_ORDER =                    4                                                  
    A_2_1   = 9.33802326056672E-11                                                  
    A_1_1   = 2.44489619913889E-05                                                  
    A_2_2   = 5.99856272799014E-15                                                  
    B_0_3   = -2.0092851573342E-10                                                  
    B_3_1   = 1.00607112230593E-13                                                  
    B_3_0   = 3.66824943640799E-11                                                  
    A_2_0   = -1.8678411786277E-07                                                  
    B_1_3   = -6.9677270201133E-15                                                  
    A_0_2   = 4.73630640333079E-08                                                  
    A_1_3   = 5.55221560333543E-13                                                  
    B_0_4   = 7.52827599670567E-13                                                  
    B_2_2   = -1.1683621160870E-13                                                  
    A_0_4   = -2.0852050771470E-13                                                  
    B_0_2   = 2.99875048026693E-05                                                  
    A_4_0   = -3.1314754837293E-13                                                  
    B_4_0   = -6.4384058620497E-13                                                  
    A_ORDER =                    4                                                  
    A_0_3   = 2.65011000430244E-11                                                  
    B_2_1   = -2.8558390691514E-10                                                  
    A_1_2   = 5.07616164062598E-11                                                  
    B_1_1   = -2.0379403931148E-07                                                  
    A_3_1   = 5.25748787891111E-13                                                  
    B_2_0   = 6.97816138011029E-06                                                  
    WCSNAME = 'IDC_w3m18525i'                                                       
    MDRIZSKY=   0.7757664823972165 / Sky value computed by AstroDrizzle             




imhistogram
-----------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Imhistogram will plot a customized histogram of the provided image data.
To make a histogram in Python we are going to use Matplotlib's ``hist``
function. See the ``hist``
`documentation <http://matplotlib.org/api/pyplot_api.html>`__ for
options to change the histogram type, scaling, bin sizes, and more.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits
    
    # Plotting Imports/Setup
    import matplotlib.pyplot as plt
    %matplotlib inline

.. code:: ipython3

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    
    # Pull out the first science array, we also need to flatten the data to a 
    # 1D array before sending it to hist
    sci1 = fits.getdata(test_data,ext=1)
    sci1f = sci1.flatten()
    
    # Now we can plot our histogram, using some of the optional keywords in hist
    # The hist function returns the values of the histogram bins (n), the edges
    # of the bins (obins), and the patches used to create the histogram
    fig = plt.figure()
    n, obins, patches = plt.hist(sci1f,bins=100,range=(0,2))
    
    # Save resulting figure to png file
    fig.savefig('hist.png')



.. image:: images.imutil_files/images.imutil_59_0.png




imreplace
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Imreplace is used to replace array sections with a constant. We can use
simple ``numpy`` array manipulation to replicate imreplace. For details
on how to grow the boolean array for replacement see crgrow, or the
`skimage.dilation
documentation <http://scikit-image.org/docs/0.12.x/api/skimage.morphology.html?highlight=dilation#skimage.morphology.dilation>`__.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    out_file = 'imreplace_out.fits'
    
    # Pull out the first science array
    hdu = fits.open(test_data)
    sci1 = hdu[1].data
    
    print("cutout of array before replacements:")
    print(sci1[50:55, 50:55])
    
    # Make boolean mask with your requirements, here we produce a boolean mask 
    # where all array elements with values >0.5 and <0.6 are set to True.
    mask1 = np.logical_and(sci1>0.8, sci1<0.82)
    
    # Use mask to replace values
    sci1[mask1] = 99
    
    print("\ncoutout of array after replacements:")
    print(sci1[50:55, 50:55])
    
    # Take updated array and write out new FITS file
    hdu[1].data = sci1
    hdu.writeto(out_file, overwrite=True)
    
    # Close FITS file
    hdu.close()


.. parsed-literal::

    cutout of array before replacements:
    [[ 0.89118606  0.87640154  0.81239933  0.77495182  0.80048275]
     [ 0.83939391  0.79715788  0.71130604  0.83452195  0.74553812]
     [ 0.82984501  0.82536161  0.82937354  0.82661521  0.80760878]
     [ 0.88277584  0.78050691  0.85906219  0.80846858  0.8092978 ]
     [ 0.85532236  0.73028219  0.81455106  0.76300722  0.85437953]]
    
    coutout of array after replacements:
    [[  0.89118606   0.87640154  99.           0.77495182  99.        ]
     [  0.83939391   0.79715788   0.71130604   0.83452195   0.74553812]
     [  0.82984501   0.82536161   0.82937354   0.82661521  99.        ]
     [  0.88277584   0.78050691   0.85906219  99.          99.        ]
     [  0.85532236   0.73028219  99.           0.76300722   0.85437953]]


.. code:: ipython3

    # We can also use numpy where to pull out index numbers
    mask2 = np.where(sci1 > 1000)
    print("Index values where sci1 is > 1,000")
    print(mask2)


.. parsed-literal::

    Index values where sci1 is > 1,000
    (array([ 474,  474,  606,  607,  607,  607,  608,  608,  608,  608,  609,
            609,  609,  609,  610,  610,  610,  804,  804,  809,  809,  810,
            883,  883, 1002, 1013]), array([455, 456, 285, 284, 285, 286, 284, 285, 286, 287, 284, 285, 286,
           287, 284, 285, 286, 349, 350,  53, 575,  53, 161, 162, 104, 460]))




imslice
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Imslice can take a 3-D datacube FITS image and return multiple 2D images
sliced through the chosen dimension. Keep in mind for the python
equivalent workflow that the header file from the original input image
will be used for all output images, including WCS information. We will
be using
`numpy.split <https://docs.scipy.org/doc/numpy/reference/generated/numpy.split.html#numpy.split>`__.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Pull image data array and image header
    orig_hdu = fits.open('/eng/ssb/iraf_transition/test_data/imstack_out.fits')
    
    print("Here's the extensions in our input file:")
    orig_hdu.info()
    
    header1 = orig_hdu[0].header
    image1 = orig_hdu[0].data
    orig_hdu.close()
    
    print("\noriginal array - the dimension order is listed " +
          "in reverse order \nnow that we have read the array into a numpy array:")
    print(image1.shape)
    
    # Slice images easily by using numpy.split, which returns a list of the output arrays
    # THen numpy.squeeze is used to remove the extra length one dimensions left over from
    # numpy.split.
    arr_list = np.split(image1, 2)
    arr_list = np.squeeze(arr_list)
    print("\nfinal shape of a slice is:")
    print(arr_list[0].shape)
    
    # Now we can write this new array into a new FITS files by packing it back into an HDU object
    hdu1 = fits.PrimaryHDU(arr_list[0],header1)
    hdu1.writeto('imslice_out1.fits', overwrite=True)
    hdu2 = fits.PrimaryHDU(arr_list[1],header1)
    hdu2.writeto('imslice_out2.fits', overwrite=True)


.. parsed-literal::

    Here's the extensions in our input file:
    Filename: /eng/ssb/iraf_transition/test_data/imstack_out.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  SCI           1 PrimaryHDU     199   (4096, 2048, 2)   float32   
    
    original array - the dimension order is listed in reverse order 
    now that we have read the array into a numpy array:
    (2, 2048, 4096)
    
    final shape of a slice is:
    (2048, 4096)




imstack
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

imstack can take multiple FITS images and stack the data, writing out a
new file where the FITS data is 1-dimension higher then the input
images. Here we show that manipulation using the ``astropy`` library and
`numpy.stack <https://docs.scipy.org/doc/numpy/reference/generated/numpy.stack.html#numpy.stack>`__.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Pull two image data arrays and an image header
    header1 = fits.getheader('/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits',ext=1)
    image1 = fits.getdata('/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits')
    image2 = fits.getdata('/eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits')
    
    # Stack arrays, the new dimension will be put first, unless otherwise specified with the axis keyword
    outstack = np.stack((image1,image2))
    print("final shape is:")
    print(outstack.shape)
    
    # Now we can write this new array into a new FITS file by packing it back into an HDU object
    hdu = fits.PrimaryHDU(outstack,header1)
    hdu.writeto('imstack_out.fits', overwrite=True)


.. parsed-literal::

    final shape is:
    (2, 2048, 4096)




imstatistics
------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

We will use the ``astropy.stats.sigma_clipped_stats`` function here,
which has some wider capabilites then the imstatistics function. Please
see the ``stats`` `package
documentation <http://docs.astropy.org/en/stable/api/astropy.stats.sigma_clipped_stats.html>`__
for details on the advanced usage. We also use some Numpy functions for
additional statistics.

.. code:: ipython3

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy import stats

.. code:: ipython3

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    sci1 = fits.getdata(test_data, ext=1)
    
    # The sigma_clipped_stats function returns the mean, median, and stddev respectively
    # To more closely replicate the IRAF version that is using n-1 in it's calculations
    # we use the std_ddof parameter
    output = stats.sigma_clipped_stats(sci1, sigma=3.0, iters=3, std_ddof=1)
    print("mean, median, standard deviation:")
    print(output)
    
    # To see the min and max of an array we can use numpy.min and numpy.max
    array_min = np.min(sci1)
    array_max = np.max(sci1)
    print("\nmin, max")
    print("{}, {}".format(array_min, array_max))
    
    # To find out how many pixels are greater then a particular value we can use numpy.where
    where_result = np.where(sci1 > 1000)
    count = len(where_result[0])
    print("\nNumber of pixels above 1,000:")
    print(count)


.. parsed-literal::

    mean, median, standard deviation:
    (0.82595410841884809, 0.81768394, 0.074634554991261454)
    
    min, max
    -4007.712890625, 27569.6015625
    
    Number of pixels above 1,000:
    26




imsum
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Imsum is used to compute the sum, average, or mean of a set of images.
We will be using the ``ccdproc`` ``Combiner`` class here. Keep in mind
that the original FITS header is not retained in the ``CCDData`` object.
Please see the `ccdproc
documentation <http://ccdproc.readthedocs.io/en/latest/ccdproc/image_combination.html>`__
for more details.

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy import units
    from ccdproc import CCDData, Combiner

.. code:: ipython3

    # Change these values to your desired data files
    test_data1 = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    test_data2 = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    
    # First we need to pull out the science arrays to create CCDData objects
    # Our actual unit is electrons/sec, this is not accepted by the current
    # set of units
    cdata1 = CCDData.read(test_data1, hdu=1, unit=units.electron/units.s)
    cdata2 = cdata1.copy()
    cdata3 = CCDData.read(test_data2, hdu=1, unit=units.electron/units.s)
    cdata4 = cdata3.copy()
    combiner = Combiner([cdata1, cdata2, cdata3, cdata4])
    
    # Now we can make our mask for extrema clipping
    # The equivalent of low_reject, high_reject parameter
    combiner.clip_extrema(nlow=1, nhigh=1)
    
    # And finally to combine...
    final_combine = combiner.average_combine()
    print(final_combine.data)


.. parsed-literal::

    INFO: using the unit electron / s passed to the FITS reader instead of the unit ELECTRONS/S in the FITS file. [astropy.nddata.ccddata]
    INFO: using the unit electron / s passed to the FITS reader instead of the unit ELECTRONS/S in the FITS file. [astropy.nddata.ccddata]
    [[  0.87720111   0.82106587   0.79521415 ...,   3.87308204   7.41545987
        9.01969481]
     [  0.89028609   0.7884455    0.8240625  ...,   0.86163342   4.53510189
        0.99109203]
     [  0.81683022   0.83273572   0.82175627 ...,   3.60699821  -7.82266164
        2.95994186]
     ..., 
     [ 40.72796059  15.36561799  -8.79329443 ...,  22.68277168  25.31048012
       28.829813  ]
     [ 46.28870392  -4.50218874   1.74757147 ...,  13.24364138  25.70440292
       11.0971849 ]
     [ 42.8106432   29.66250706  63.18441772 ...,   0.           9.80057049
       22.66858006]]




listpixels
----------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Listpixels was used to list an indexed section of a FITS data array.
This is easy to do using ``astropy``, but **keep in mind that Python
indexes from zero, and with the y-axis leading, i.e. [y,x]**. You also
want to end the cut with the pixel *after* the end pixel. So to get 1-10
in x and 5-15 in y, you will index like so: array[4:15,0:10]. To see
listpixels results for more then one file, you will need to loop over a
list of files, see information about Python loops
`here <http://www.pythonforbeginners.com/loops/for-while-and-nested-loops-in-python>`__.

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # Change this value to your desired data files
    test_data1 = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    
    # To quickly pull out the data array you can use the astropy convenience function
    data_arr = fits.getdata(test_data1,ext=1)
    
    # Now we can index the array as desired
    # We're cutting out 5 in y, and 2 in x
    print(data_arr[0:5,0:2])


.. parsed-literal::

    [[ 0.86692303  0.80678135]
     [ 0.83312052  0.76854318]
     [ 0.77341086  0.80276382]
     [ 0.80539584  0.78261763]
     [ 0.78274417  0.82206035]]




Not Replacing
-------------

-  imrename - can use command line utilities or the Python ``os``
   package for this functionality.
-  imdelete - can use command line utilities or the Python ``os``
   package for this functionality.
-  imtile - **may** replace infuture
-  sections - IRAF utility function
-  imgets - see `images.imutil.hselect <#hselect>`__
-  minmax - see `images.imutil.imstatistics <#imstatistics>`__
