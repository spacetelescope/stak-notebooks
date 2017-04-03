:orphan:


images.imutil
=============

The images.imutil package provides general fits image tools such as
header editing and image arithimetic.

Notes
-----

Contents:

-  `chpixtype <#chpixtype>`__
-  `hedit <#hedit>`__
-  `hselect <#hselect>`__
-  `imarith/imdivide <#imarith>`__
-  `imcopy <#imcopy>`__
-  `imfunction/imexpr <#imfunction>`__
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

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

Chpixtype is a task that allows you to change the pixel type of a fits
image. There is built in functionality in ``astropy.io.fits`` to preform
this task with the ``scale`` method. Below you will find a table that
translates the chpixtype newpixtype options into their equivalent
``numpy``/``astropy`` type
(http://docs.scipy.org/doc/numpy/user/basics.types.html).

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

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change this value to your desired data file, here were creating a filename
    # for our new changed data
    orig_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    new_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_newdtype_flt.fits'
    
    # Read in your fits file
    hdu = fits.open(orig_data)
    
    # Edit the datatype
    hdu[1].scale(type='int32')
    
    # Save changed hdu object to new file
    # The clobber argument tells the writeto method to overwrite if file already exists
    hdu.writeto(new_data, clobber=True)
    hdu.close()



hedit
-----

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

The hedit task allows users to edit an image header. This functioanlity
is covered in ``astropy.io.fits``. Take note that to make changes to a
fits file, you must use the ``mode='update'`` keyword in the
``fits.open`` call. Below you'll find examples of editing a keyword if
it does/doesn't exist, and how to delete keywords from the header.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change this value to your desired data file
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    
    # Open fits file, include the mode='update' keyword
    hdu = fits.open(test_data, mode='update')
    
    # Simple header change, will add keyword if it doesn't not exist
    hdu[0].header['MYKEY1'] = 'Editing this keyword'
    
    # Only add keyword if it does not already exist:
    if 'MYKEY2' not in hdu[0].header:
        hdu[0].header['MYKEY2'] = 'Also editing this'
    
    # To delete keywords, first check if they exist:
    if 'MYKEY2' in hdu[0].header:
        del hdu[0].header['MYKEY2']
        
    # Close fits file, this will save your changes
    hdu.close()



hselect
-------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

hselect is used to pull out specific header keywords. You can provide
any filename string as you would in IRAF and it will be exapanded
(wildcards are accepted). You can also use specific keyword values to
filter files. We will be using the ``stak`` package ``Hselect`` class.
The output table is an ``astropy.table`` object and stored in the
``table`` attribute.

.. code:: ipython2

    # Astronomy Specific Imports
    from stak import Hselect

.. code:: ipython2

    # Create Hselect object
    myList = Hselect("/eng/ssb/iraf_transition/test_data/jcz*", "BUNIT,TIME-OBS", extension="0,1,2,3")
    # Display output astropy table object in nice notebook formatting
    myList.table.show_in_notebook()




.. raw:: html

    &lt;Table masked=True length=8&gt;
    <table id="table4497701264-588354" class="table-striped table-bordered table-condensed">
    <thead><tr><th>idx</th><th>Filename</th><th>ExtNumber</th><th>BUNIT</th><th>TIME-OBS</th></tr></thead>
    <tr><td>0</td><td>/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits</td><td>0</td><td>--</td><td>01:04:51</td></tr>
    <tr><td>1</td><td>/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits</td><td>1</td><td>ELECTRONS</td><td>--</td></tr>
    <tr><td>2</td><td>/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits</td><td>2</td><td>ELECTRONS</td><td>--</td></tr>
    <tr><td>3</td><td>/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits</td><td>3</td><td>UNITLESS</td><td>--</td></tr>
    <tr><td>4</td><td>/eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits</td><td>2</td><td>ELECTRONS</td><td>--</td></tr>
    <tr><td>5</td><td>/eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits</td><td>3</td><td>UNITLESS</td><td>--</td></tr>
    <tr><td>6</td><td>/eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits</td><td>0</td><td>--</td><td>02:16:10</td></tr>
    <tr><td>7</td><td>/eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits</td><td>1</td><td>ELECTRONS</td><td>--</td></tr>
    </table><style>table.dataTable {clear: both; width: auto !important; margin: 0 !important;}
    .dataTables_info, .dataTables_length, .dataTables_filter, .dataTables_paginate{
    display: inline-block; margin-right: 1em; }
    .paginate_button { margin-right: 5px; }
    </style>
    <script>
    require.config({paths: {
        datatables: 'https://cdn.datatables.net/1.10.9/js/jquery.dataTables.min'
    }});
    require(["datatables"], function(){
        console.log("$('#table4497701264-588354').dataTable()");
        $('#table4497701264-588354').dataTable({
            "order": [],
            "iDisplayLength": 50,
            "aLengthMenu": [[10, 25, 50, 100, 500, 1000, -1], [10, 25, 50, 100, 500, 1000, 'All']],
            "pagingType": "full_numbers"
        });
    });
    </script>




.. code:: ipython2

    # Create Hselect object using expression parsing
    myList2 = Hselect("/eng/ssb/iraf_transition/test_data/jcz*", "BUNIT", extension="0,1,2,3",
                     expr="BUNIT='ELECTRONS'")
    # Display output astropy table object with a standard print
    print(myList2.table)


.. parsed-literal::

                           Filename                       ExtNumber   BUNIT  
    ----------------------------------------------------- --------- ---------
    /eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits         2 ELECTRONS
    /eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits         1 ELECTRONS
    /eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits         2 ELECTRONS
    /eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits         1 ELECTRONS




imarith - imdivide
------------------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

Imarith and imdivide both provide functionality to apply basic operators
to whole image arrays. This task can be achieved with basic
``astropy.io.fits`` functionality along with ``numpy`` array
functionality.

The basic operands (``+``,\ ``-``,\ ``/``,\ ``*``) can all be used with
an assignment operator in python (``+=``,\ ``-=``,\ ``/=``,\ ``*=``).
See http://www.tutorialspoint.com/python/python\_basic\_operators.htm
for more details

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Basic operands (+,-,/,*)
    # Change these values to your desired data files
    test_data1 = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    test_data2 = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    output_data = '/eng/ssb/iraf_transition/test_data/imarith_out.fits'
    
    # Open fits file
    hdu1 = fits.open(test_data1)
    hdu2 = fits.open(test_data2)
    
    # Here we add hdu2-ext1 to hdu1-ext1 by using the shortcute += operator
    hdu1[1].data += hdu2[1].data
    
    # If you are dividing and need to avoid zeros in the image use indexing
    indx_zeros = [hdu2[4].data == 0]
    indx_nonzeros = [hdu2[4].data != 0]
    # Set this value as you would the divzero parameter in imarith
    set_zeros = 999.9
    hdu1[4].data[indx_nonzeros] /= hdu2[4].data[indx_nonzeros]
    hdu1[4].data[indx_zeros] = 999.9
    
    # Save your new file
    # The clobber argument tells the writeto method to overwrite if file already exists
    hdu1.writeto(output_data, clobber=True)
    
    # Close hdu files
    hdu1.close()
    hdu2.close()



imcopy
------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

Imcopy allows users to copy a fits image to a new file. We can
accomplish this using ``astropy.io.fits`` by saving our fits file to a
new filename.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change these values to your desired filenames
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    output_data = '/eng/ssb/iraf_transition/test_data/imcopy_out.fits'
    
    hdu = fits.open(test_data)
    # The clobber argument tells the writeto method to overwrite if file already exists
    hdu.writeto(output_data, clobber=True)
    hdu.close()



imfunction - imexpr
-------------------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

Imfunction will apply a function to the image pixel values in an image
array. Imexpr gives you similiar functionality with the added capability
to combine different images using a user created expression. We can
accomplish this using the built in funcitonality of the ``numpy``
library (http://docs.scipy.org/doc/numpy/reference/routines.math.html)

If there is a particular function you would like to apply to your image
array that you cannot find in the ``numpy`` library you can use the
``np.vectorize`` function, which can make any python function apply to
each element of your array. But keep in mind that ``np.vectorize`` is
esentially looping over the array, and may not be the most efficient
method
(http://docs.scipy.org/doc/numpy/reference/generated/numpy.vectorize.html).

Example using exsisting numpy function:

.. code:: ipython2

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    output_data = '/eng/ssb/iraf_transition/test_data/imfunction_out.fits'
    
    # Here we use the cosine function as an example
    hdu = fits.open(test_data)
    sci = hdu[1].data
    
    # When you call your new function, make sure to reassign the array to
    # the new values if the original function is not changing values in place
    hdu[1].data = np.cos(hdu[1].data)
    
    # Now save out to a new file, and close the original file, changes will
    # not be applied to the oiginal fits file.
    hdu.writeto(output_data, clobber=True)
    hdu.close()

Example using user defined function and ``np.vectorize``:

.. code:: ipython2

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    output_data = '/eng/ssb/iraf_transition/test_data/imfunction2_out.fits'
    
    # Here we use the following custom function as an example
    def my_func(x):
        return (x**2)+(x**3)
    
    # Now we open our file, and vectorize our function
    hdu = fits.open(test_data)
    sci = hdu[1].data
    vcos = np.vectorize(my_func)
    
    # When you call your new function, make sure to reassign the array to
    # the new values if the original function is not changing values in place
    hdu[1].data = vcos(hdu[1].data)
    
    # Now save out to a new file, and close the original file, changes will
    # not be applied to the oiginal fits file.
    hdu.writeto(output_data)
    hdu.close()



imheader
--------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

The imheader task allows the user to list header parameters for a list
of images. Here we can use the ``astropy`` convenience function,
``fits.getheader()``

.. code:: ipython2

    # Standard Imports
    import numpy as np
    import glob
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change these values to your desired data files, glob will capture all wildcard matches
    test_data = glob.glob('/eng/ssb/iraf_transition/test_data/iczgs3y*')
    
    for filename in test_data:
        # Pull the header from extension 1
        head = fits.getheader(filename, ext=1)
        print repr(head)


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

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

Imhistogram will plot a customized histogram of the provided image data.
To make a histogram in Python we are going to use matplotlibs ``hist``
function. See the ``hist`` documentation for options to change the
histogram type, scaling, bin sizes, and more
(http://matplotlib.org/api/pyplot\_api.html)

.. code:: ipython2

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits
    
    # Plotting Imports/Setup
    import matplotlib.pyplot as plt
    %matplotlib inline

.. code:: ipython2

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    
    # Pull out the first science array, we also need to flatten the data before sending it to hist
    sci1 = fits.getdata(test_data,ext=1)
    sci1f = sci1.flatten()
    
    # Now we can plot our histogram, using some of the optional keywords in hist
    # The hist function returns the values of the histogram bins (n), the edges
    # of the bins (obins), and the patches used to create the histogram
    n, obins, patches = plt.hist(sci1f,bins=100,range=(0,2))



.. image:: images.imutil_files/images.imutil_48_0.png




imreplace
---------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

We can use simple ``numpy`` array manipulation to replicate imreplace.
For details on how to grow the boolean array for replacement see crgrow,
or the ```skimage.dilation``
documentation <http://scikit-image.org/docs/0.12.x/api/skimage.morphology.html?highlight=dilation#skimage.morphology.dilation>`__.

.. code:: ipython2

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    
    # Pull out the first science array, make boolean mask with your requirements
    hdu = fits.open(test_data)
    sci1 = hdu[1].data
    hdu.close()
    mask1 = np.logical_and(sci1>0.5, sci1<0.6)
    
    # Use mask to replace values
    sci1[mask1] = 999
    
    # We can also use numpy where to pull out index numbers
    mask2 = np.where(sci1 > 1000)
    print mask2


.. parsed-literal::

    (array([ 474,  474,  606,  607,  607,  607,  608,  608,  608,  608,  609,
            609,  609,  609,  610,  610,  610,  804,  804,  809,  809,  810,
            883,  883, 1002, 1013]), array([455, 456, 285, 284, 285, 286, 284, 285, 286, 287, 284, 285, 286,
           287, 284, 285, 286, 349, 350,  53, 575,  53, 161, 162, 104, 460]))




imslice
-------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

Imslice can take a 3-D datacube fits image and return multiple 2D images
sliced through the chosen dimension. Keep in mind for the python
equivalent workflow that the header file from the original input image
will be used for all output images, including WCS information. We will
be using
```numpy.split`` <https://docs.scipy.org/doc/numpy/reference/generated/numpy.split.html#numpy.split>`__.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Pull image data array and image header
    orig_hdu = fits.open('/eng/ssb/iraf_transition/test_data/imstack_out.fits')
    header1 = orig_hdu[0].header
    image1 = orig_hdu[0].data
    orig_hdu.close()
    
    # Slice images easily by using numpy.split, which returns a list of the output arrays
    arr_list = np.split(image1, 2)
    print("final shape of a slice is:")
    print(arr_list[0].shape)
    
    # Now we can write this new array into a new fits files by packing it back into an HDU object
    hdu1 = fits.PrimaryHDU(arr_list[0],header1)
    hdu1.writeto('/eng/ssb/iraf_transition/test_data/imslice_out1.fits', clobber=True)
    hdu2 = fits.PrimaryHDU(arr_list[1],header1)
    hdu2.writeto('/eng/ssb/iraf_transition/test_data/imslice_out2.fits', clobber=True)


.. parsed-literal::

    final shape of a slice is:
    (1, 2048, 4096)




imstack
-------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

imstack can take multiple fits images and stack the data, writing out a
new file where the fits data is 1-dimension higher then the input
images. Here we show that manipulation using the ``astropy`` library and
```numpy.stack`` <https://docs.scipy.org/doc/numpy/reference/generated/numpy.stack.html#numpy.stack>`__.

.. code:: ipython2

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Pull two image data arrays and image header
    header1 = fits.getheader('/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits',ext=1)
    image1 = fits.getdata('/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits')
    image2 = fits.getdata('/eng/ssb/iraf_transition/test_data/jczgx1q1q_flc.fits')
    
    # Stack arrays, the new dimension will be put first, unless otherwise specified with the axis keyword
    outstack = np.stack((image1,image2))
    print("final shape is:")
    print(outstack.shape)
    
    # Now we can write this new array into a new fits file by packing it back into an HDU object
    hdu = fits.PrimaryHDU(outstack,header1)
    hdu.writeto('/eng/ssb/iraf_transition/test_data/imstack_out.fits', clobber=True)


.. parsed-literal::

    final shape is:
    (2, 2048, 4096)




imstatistics
------------

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

We will use the ``astropy.stats.sigma_clipped_stats`` function here,
which has some wider capabilites then the imstatistics function. Please
see the ``stats`` `package
documentation <http://docs.astropy.org/en/stable/api/astropy.stats.sigma_clipped_stats.html>`__
for details on the advanced usage .

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy import stats

.. code:: ipython2

    # Change these values to your desired data files
    test_data = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    sci1 = fits.getdata(test_data,ext=1)
    
    # The sigma_clipped_stats function returns the mean, median, and stddev respectively
    output = stats.sigma_clipped_stats(sci1,sigma=2.0,iters=3)
    print output


.. parsed-literal::

    (0.82121155347072006, 0.81694626808166504, 0.058198063937460652)




imsum
-----

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

We will be using the ``ccdproc`` ``Combiner`` class here. Keep in mind
that the original fits header is not retained in the ``CCDData`` object.
Please see the documentation for more details
(http://ccdproc.readthedocs.io/en/latest/ccdproc/image\_combination.html).

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy import units
    from ccdproc import CCDData, Combiner

.. code:: ipython2

    # Change these values to your desired data files
    test_data1 = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    test_data2 = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    
    # First we need to pull out the science arrays to create CCDData objects
    # Our acutal unit is electrons/sec, this is not accepted by the current
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
    print final_combine.data


.. parsed-literal::

    INFO: using the unit electron / s passed to the FITS reader instead of the unit ELECTRONS/S in the FITS file. [ccdproc.ccddata]
    INFO: using the unit electron / s passed to the FITS reader instead of the unit ELECTRONS/S in the FITS file. [ccdproc.ccddata]
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

\*\* Please review the `Notes <#notes>`__ section above before running
any examples in this notebook \*\*

listpixels was used to list an indexed section of a FITs data array.
This is easy to do using ``astropy``, but keep in mind that Python
indexs from zero, and with the y-axis leading, i.e. [y,x]. You also want
to end the cut with the pixel *after* the end pixel. So to get 1-10 in x
and 5-15 in y, you will index like so: array[4:15,0:10]

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # Change these values to your desired data files
    test_data1 = '/eng/ssb/iraf_transition/test_data/iczgs3y5q_flt.fits'
    
    # To quickly pull out the data array you can use the astropy convience fucntion
    data_arr = fits.getdata(test_data1,ext=1)
    
    # Now we can index the array as desired, we're cutting out 5 in y, and 2 in x
    print data_arr[0:5,0:2]


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
-  imgets - see `**images.imutil.hselect** <#hselect>`__
-  minmax - see `**images.imutil.imstat** <#imstat>`__

For questions or comments please see `our github
page <https://github.com/spacetelescope/stak>`__. We encourage and
appreciate user feedback.
