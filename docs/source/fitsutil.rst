:orphan:


fitsutil
========

General fits file utilities.

Notes
-----

**For questions or comments please see** `our github
page <https://github.com/spacetelescope/stak>`__. **We encourage and
appreciate user feedback.**

For the compression tasks included in fitsutil, astropy has replaced
this functionality with the
`CompImageHDU <http://docs.astropy.org/en/stable/io/fits/api/images.html#astropy.io.fits.CompImageHDU>`__
class. We list both compression tasks together in this notebook with a
few examples to show the usage of ``CompImageHDU``. the
`astropy.io.fits <http://docs.astropy.org/en/stable/io/fits/index.html>`__
can natively open compressed file with a standard ``fits.open`` command.
To uncompress the file, you can then save the fits file object out to a
new file.

``astropy.io.fits`` is the library responsible for opening and closing
fits files. It opens the file into a ``HDUList`` object, which contains
multiple ``HDU`` objects that can be indexed with integers. Each ``HDU``
object contains array data object, and a header object.

Contents:

-  `fpack-ricepack <#fpack-ricepack>`__
-  `fxcopy-fxinsert <#fxcopy-fxinsert>`__
-  `fxdelete-fxsplit-fxextract <#fxdelete-fxsplit-fxextract>`__
-  `fxdummyh <#fxdummyh>`__
-  `fxheader <#fxheader>`__
-  `fxplf <#fxplf>`__





fpack-ricepack
--------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

We can compress an image HDU using the
`fits.CompImageHDU <http://docs.astropy.org/en/stable/io/fits/api/images.html#astropy.io.fits.CompImageHDU>`__
class in ``astropy``. This class has several compression options (RICE,
PLIO, GZIP, and HCOMPRESS). Here we show one example using RICE
compression and another using GZIP compression

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # RICE example
    
    # test files
    test_data = '/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits'
    outfile_rice = '/eng/ssb/iraf_transition/test_data/jczgx1ppq_rice.fits'
    
    hdulist = fits.open(test_data)
    hdu_rice = fits.CompImageHDU(data=hdulist[1].data, header=hdulist[1].header, compression_type='RICE_1')
    hdulist_rice = fits.HDUList([hdulist[0],hdu_rice])
    hdulist_rice.writeto(outfile_rice, overwrite=True)
    hdulist.close()

.. code:: ipython2

    # Gzip example
    
    # test files
    test_data = '/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits'
    outfile_gzip = '/eng/ssb/iraf_transition/test_data/jczgx1ppq_gzip.fits'
    
    hdulist = fits.open(test_data)
    hdu_gzip = fits.CompImageHDU(data=hdulist[1].data, header=hdulist[1].header, compression_type='GZIP_1')
    hdulist_gzip = fits.HDUList([hdulist[0],hdu_rice])
    hdulist_gzip.writeto(outfile_gzip, overwrite=True)
    hdulist.close()



fxcopy-fxinsert
---------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Here we show how to copy out and add new HDU objects, the ``astropy``
equivalent of fxcopy and fxinsert.

.. code:: ipython2

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # test files
    test_data = '/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits'
    outfile = '/eng/ssb/iraf_transition/test_data/fxinsert.fits'
    
    # open fits file, this outputs an hdulist object
    hdulist = fits.open(test_data)
    
    # now let's pull out a reference (copy) of an HDU object from this HDUList
    my_hdu = hdulist[1]
    
    # Now let's create a new array to make a new HDU object, this will be the primary HDU
    new = np.arange(100.0)
    new_hdu = fits.PrimaryHDU(new)
    
    # Now we can create a new HDUList object to put our HDU objects into
    my_hdulist = fits.HDUList([new_hdu,my_hdu])
    
    # Now we close write our new HDUList to a file, and close our test_data file
    my_hdulist.writeto(outfile, overwrite=True)
    hdulist.close()



fxdelete-fxsplit-fxextract
--------------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

fxdelete will delete a FITS extension in place, and fxsplit and
fxextract will take a multiple extension FITS file and break them out
into single FITS files. Both these tasks can be done using
`astropy.io.fits <http://docs.astropy.org/en/stable/io/fits/index.html>`__.
Below we show some a short example. We will pull out the 3rd extension
from the test file, save it to a new fits file, and delete that
extension from the original ``HDUList``

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # FITS filenames
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    outfile_1 = '/eng/ssb/iraf_transition/test_data/fxsplit.fits'
    outfile_2 = '/eng/ssb/iraf_transition/test_data/fxdelete.fits'
    
    # Print out some stats for this file
    fits.info(test_data)
    
    # Open FITS file
    hdulist = fits.open(test_data)
    
    # Pull out single HDU extension and put into new FITS file
    single_HDU = hdulist[3]
    primary_HDU = fits.PrimaryHDU()
    new_hdulist = fits.HDUList([primary_HDU,single_HDU])
    new_hdulist.writeto(outfile_1, overwrite=True)
    
    
    # Now save a new copy of the original file without that thrid extension
    #cut_hdulist = hdulist[0:2] + hdulist[]
    hdulist.pop(3)
    hdulist.writeto(outfile_2, overwrite=True)
    
    # Close original file
    hdulist.close()


.. parsed-literal::

    Filename: /eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU     265   ()      
      1  SCI           1 ImageHDU       140   (1014, 1014)   float32   
      2  ERR           1 ImageHDU        51   (1014, 1014)   float32   
      3  DQ            1 ImageHDU        43   (1014, 1014)   int16   
      4  SAMP          1 ImageHDU        37   (1014, 1014)   int16   
      5  TIME          1 ImageHDU        37   (1014, 1014)   float32   
      6  WCSCORR       1 BinTableHDU     59   7R x 24C   [40A, I, A, 24A, 24A, 24A, 24A, D, D, D, D, D, D, D, D, 24A, 24A, D, D, D, D, J, 40A, 128A]   




fxdummyh
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Fxdummyh will create an empty fits file.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Write empty file
    hdu = fits.PrimaryHDU()
    hdu.writeto('/eng/ssb/iraf_transition/test_data/empty.fits', overwrite=True)
    
    # Let's look at the file we made
    fits.info('/eng/ssb/iraf_transition/test_data/empty.fits')


.. parsed-literal::

    Filename: /eng/ssb/iraf_transition/test_data/empty.fits
    No.    Name      Ver    Type      Cards   Dimensions   Format
      0  PRIMARY       1 PrimaryHDU       4   ()      




fxheader
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

Fxheader lists one line of header description per FITS unit. This
funtionality has been replaced in a convience function in ``astropy``,
`astropy.io.fits.info <http://docs.astropy.org/en/stable/io/fits/#convenience-functions>`__.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # run fits.info
    fits.info('/eng/ssb/iraf_transition/test_data/tester.fits')


.. parsed-literal::

    Filename: /eng/ssb/iraf_transition/test_data/tester.fits
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




fxplf
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

fxplf is used to convert a pixel list file into a BINTABLE extension. We
show a simple example below, see the `Astropy unified read/write
documentation <http://docs.astropy.org/en/stable/io/unified.html#fits>`__
for more details.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy.table import Table

.. code:: ipython2

    # Define input and output files
    infile = '/eng/ssb/iraf_transition/test_data/table3.txt'
    outfile = '/eng/ssb/iraf_transition/test_data/table3.fits'
    
    # read txt, write to fits
    t = Table.read(infile, format='ascii')
    print(t)
    t.write(outfile, overwrite=True)


.. parsed-literal::

    col1 col2
    ---- ----
     200   45
      34  222
       3    4
     100  200
       8   88
      23  123






Not Replacing
-------------

-  funpack - Uncompress FITS file, can be done by opening and resaving
   file with
   `astropy.io.fits <http://docs.astropy.org/en/stable/io/fits/index.html>`__
-  fxconvert - Convert between IRAF image types. See
   **images.imutil.imcopy**
-  fgread - Read a MEF file with FOREIGN extensions. Deprecated.
-  fgwrite - Create a MEF file with FOREIGN extensions. Deprecated.
