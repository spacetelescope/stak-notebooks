:orphan:


stsdas.toolbox.headers
======================

This package provides utilities for comparing and editing image headers.

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

Many of the headers tasks can be replaced with utilities in ``astropy``,
as seen below.

Contents:

-  `eheader <#eheader>`__
-  `hdiff <#hdiff>`__
-  `stfhistory <#stfhistory>`__
-  `upreffile <#upreffile>`__





eheader
-------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress



hdiff
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The hdiff task will take two FITS headers and report the differences
between them. This functionality has been replaced and improved upon in
``astropy`` with the ``astropy.io.fits.Differs`` class, which can be
easily called with the `printdiff convenience
function <http://docs.astropy.org/en/stable/io/fits/api/files.html#printdiff>`__.
For more details on a more advanced differ result using
``astropy.io.fits.Differs`` directly, see the `API
doc <http://docs.astropy.org/en/stable/io/fits/api/diff.html>`__.

.. code:: ipython3

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy.io.fits import printdiff

.. code:: ipython3

    file1 = '/eng/ssb/iraf_transition/test_data/iczgs3ygq_flt.fits'
    file2 = '/eng/ssb/iraf_transition/test_data/jczgx1ppq_flc.fits'
    
    # printdiff example ignoring HISTORY and COMMENT cards, and only extension 0
    printdiff(file1, file2, ext=0, ignore_keywords=('HISTORY', 'COMMENT'))


.. parsed-literal::

    
     Headers contain differences:
       Headers have different number of cards:
        a: 225
        b: 233
       Extra keyword 'ANG_SIDE' in a: 0.0
       Extra keyword 'BIACFILE' in a: 'N/A'
       Extra keyword 'CCDOFSAB' in a: 190
       Extra keyword 'CCDOFSCD' in a: 190
       Extra keyword 'CSMID'  in a: 'IR'
       Extra keyword 'DWELL_LN' in a: 0
       Extra keyword 'DWELL_TM' in a: 0.0
       Extra keyword 'FILTER' in a: 'F140W'
       Extra keyword 'MYKEY1' in a: 'Editing this keyword'
       Extra keyword 'NLINCORR' in a: 'COMPLETE'
       Extra keyword 'NLINFILE' in a: 'iref$u1k1727mi_lin.fits'
       Extra keyword 'NO_LINES' in a: 0
       Extra keyword 'NSAMP'  in a: 14
       Extra keyword 'PHOTBW' in a: 1132.39
       Extra keyword 'PHOTFLAM' in a: 1.4737148e-20
       Extra keyword 'PHOTFNU' in a: 9.5291135e-08
       Extra keyword 'PHOTMODE' in a: 'WFC3 IR F140W'
       Extra keyword 'PHOTPLAM' in a: 13922.907
       Extra keyword 'PHOTZPT' in a: -21.1
       Extra keyword 'SAACRMAP' in a: 'N/A'
       Extra keyword 'SAA_DARK' in a: 'N/A'
       Extra keyword 'SAA_EXIT' in a: '2016.015:05:44:00'
       Extra keyword 'SAA_TIME' in a: 3808
       Extra keyword 'SAMPZERO' in a: 2.911756
       Extra keyword 'SAMP_SEQ' in a: 'SPARS50'
       Extra keyword 'SCAN_ANG' in a: 0.0
       Extra keyword 'SCAN_COR' in a: 'C'
       Extra keyword 'SCAN_LEN' in a: 0.0
       Extra keyword 'SCAN_RAT' in a: 0.0
       Extra keyword 'SCAN_TYP' in a: 'N'
       Extra keyword 'SCAN_WID' in a: 0.0
       Extra keyword 'SNKCFILE' in a: 'N/A'
       Extra keyword 'SUBTYPE' in a: 'FULLIMAG'
       Extra keyword 'UNITCORR' in a: 'COMPLETE'
       Extra keyword 'ZOFFCORR' in a: 'COMPLETE'
       Extra keyword 'ZSIGCORR' in a: 'COMPLETE'
       Extra keyword 'ATODCORR' in b: 'OMIT'
       Extra keyword 'BIASCORR' in b: 'COMPLETE'
       Extra keyword 'CCDOFSTA' in b: 1
       Extra keyword 'CCDOFSTB' in b: 1
       Extra keyword 'CCDOFSTC' in b: 1
       Extra keyword 'CCDOFSTD' in b: 1
       Extra keyword 'CFLTFILE' in b: 'N/A'
       Extra keyword 'CRSPLIT' in b: 1
       Extra keyword 'CTEDIR' in b: 'NONE'
       Extra keyword 'CTEIMAGE' in b: 'NONE'
       Extra keyword 'CTE_NAME' in b: 'PixelCTE 2012'
       Extra keyword 'CTE_VER' in b: '3.3'
       Extra keyword 'DARKTIME' in b: 581.247202
       Extra keyword 'EXPSCORR' in b: 'COMPLETE'
       Extra keyword 'FILTER1' in b: 'CLEAR1L'
       Extra keyword 'FILTER2' in b: 'F814W'
       Extra keyword 'FLASHCUR' in b: 'OFF'
       Extra keyword 'FLASHDUR' in b: 0.0
       Extra keyword 'FLASHSTA' in b: 'NOT PERFORMED'
       Extra keyword 'FLSHCORR' in b: 'OMIT'
       Extra keyword 'FW1ERROR' in b: False
       Extra keyword 'FW1OFFST' in b: 0
       Extra keyword 'FW2ERROR' in b: False
       Extra keyword 'FW2OFFST' in b: 0
       Extra keyword 'FWSERROR' in b: False
       Extra keyword 'FWSOFFST' in b: 0
       Extra keyword 'JWROTYPE' in b: 'DS_int'
       Extra keyword 'LRFWAVE' in b: 0.0
       Extra keyword 'MLINTAB' in b: 'N/A'
       Extra keyword 'PCTECORR' in b: 'COMPLETE'
       Extra keyword 'PCTEFRAC' in b: 2.011591423149
       Extra keyword 'PCTENSMD' in b: 1
       Extra keyword 'PCTERNCL' in b: 5.0
       Extra keyword 'PCTESHFT' in b: 7
       Extra keyword 'PCTESMIT' in b: 7
       Extra keyword 'PCTETRSH' in b: -10.0
       Extra keyword 'PHOTTAB' in b: 'N/A'
       Extra keyword 'SHADCORR' in b: 'OMIT'
       Extra keyword 'SHADFILE' in b: 'N/A'
       Extra keyword 'SHUTRPOS' in b: 'A'
       Extra keyword 'SPOTTAB' in b: 'N/A'
       Extra keyword 'STATFLAG' in b: False
       Extra keyword 'WRTERR' in b: True
       Inconsistent duplicates of keyword ''      :
        Occurs 19 time(s) in a, 17 times in (b)
       Keyword         [8] has different values:
          a>       / INSTRUMENT CONFIGURATION INFORMATION
           ?                                 ------------
          b>       / SCIENCE INSTRUMENT CONFIGURATION
           ?        ++++++++
       Keyword         [9] has different values:
          a>       / POST-SAA DARK KEYWORDS
          b>       / CALIBRATION SWITCHES: PERFORM, OMIT, COMPLETE
       Keyword         [10] has different values:
          a>       / SCAN KEYWORDS
          b>       / CALIBRATION REFERENCE FILES
       Keyword         [11] has different values:
          a>       / CALIBRATION SWITCHES: PERFORM, OMIT, COMPLETE, SKIPPED
          b>       / COSMIC RAY REJECTION ALGORITHM PARAMETERS
       Keyword         [12] has different values:
          a>       / CALIBRATION REFERENCE FILES
          b>       / OTFR KEYWORDS
       Keyword         [13] has different values:
          a>       / COSMIC RAY REJECTION ALGORITHM PARAMETERS
          b>       / PATTERN KEYWORDS
       Keyword         [14] has different values:
          a>       / PHOTOMETRY KEYWORDS
          b>       / POST FLASH  PARAMETERS
       Keyword         [15] has different values:
          a>       / OTFR KEYWORDS
          b>       / ENGINEERING PARAMETERS
       Keyword         [16] has different values:
          a>       / PATTERN KEYWORDS
          b>       / CALIBRATED ENGINEERING PARAMETERS
       Keyword         [17] has different values:
          a>       / ENGINEERING PARAMETERS
          b>       / ASSOCIATION KEYWORDS
       Keyword APERTURE has different values:
          a> IR-FIX
          b> WFCENTER
       Keyword ASN_ID   has different values:
          a> NONE
          b> JCZGX1020
       Keyword ASN_MTYP has different values:
          b> EXP-DTH
       Keyword ASN_TAB  has different values:
          a> NONE
          b> jczgx1020_asn.fits
       Keyword ATODGNA  has different values:
          a> 2.3399999
          b> 2.02
       Keyword ATODGNB  has different values:
          a> 2.3699999
          b> 1.886
       Keyword ATODGNC  has different values:
          a> 2.3099999
          b> 2.017
       Keyword ATODGND  has different values:
          a> 2.3800001
          b> 2.0109999
       Keyword ATODTAB  has different comments:
          b> analog to digital correction file
       Keyword BIASFILE has different values:
          a> N/A
          b> jref$0a415460j_bia.fits
       Keyword BIASFILE has different comments:
          b> bias image file name
       Keyword BIASLEVA has different values:
          a> 0.0
          b> 4221.167
       Keyword BIASLEVB has different values:
          a> 0.0
          b> 4029.7476
       Keyword BIASLEVC has different values:
          a> 0.0
          b> 4441.6982
       Keyword BIASLEVD has different values:
          a> 0.0
          b> 4631.4839
       Keyword BLEVCORR has different comments:
          a> subtract bias level computed from ref pixels
           ?                                    ^^^^ ^^^^
          b> subtract bias level computed from overscan img
           ?                                   +++ ^^^^^ ^^
       Keyword BPIXTAB  has different values:
          a> iref$y711520di_bpx.fits
          b> jref$t3n1116nj_bpx.fits
       Keyword CAL_VER  has different values:
          a> 3.3(28-Jan-2016)
          b> 8.3.4 (07-Jul-2016)
       Keyword CAL_VER  has different comments:
          a> CALWF3 code version
           ?    ^^^
          b> CALACS code version
           ?    ^^^
       Keyword CCDGAIN  has different values:
          a> 2.5
          b> 2.0
       Keyword CCDTAB   has different values:
          a> iref$t2c16200i_ccd.fits
          b> jref$xa81715gj_ccd.fits
       Keyword CCDTAB   has different comments:
          a> detector calibration parameters
           ? ^^^^^^^^
          b> CCD calibration parameters
           ? ^^^
       Keyword CRCORR   has different values:
          a> COMPLETE
          b> OMIT
       Keyword CRCORR   has different comments:
          a> identify cosmic ray hits
          b> combine observations to reject cosmic rays
       Keyword CRDS_CTX has different values:
          a> hst_0453.pmap
           ?       ^^
          b> hst_0484.pmap
           ?       ^^
       Keyword CRDS_VER has different values:
          a> 1.9.0, opus_2015_3dms, 2767
          b> 7.0.1, opus_2016.1-universal, af27872
       Keyword CRREJTAB has different values:
          a> iref$u6a1748ri_crr.fits
          b> N/A
       Keyword CSYS_VER has different values:
          a> opus_2015_3a_dms
          b> hstdp-2016.1
       Keyword D2IMFILE has different values:
          a> N/A
          b> jref$02c1450oj_d2i.fits
       Keyword D2IMFILE has different comments:
          b> Column Correction Reference File
       Keyword DARKFILE has different values:
          a> iref$xag19296i_drk.fits
          b> jref$0a41546sj_drk.fits
       Keyword DATE     has different values:
          a> 2016-06-02
          b> 2016-10-16
       Keyword DATE-OBS has different values:
          a> 2016-01-15
           ?       -  ^
          b> 2016-10-16
           ?      +   ^
       Keyword DEC_TARG has different values:
          a> 48.92264646942
          b> 65.84194444444
       Keyword DETECTOR has different values:
          a> IR
          b> WFC
       Keyword DETECTOR has different comments:
          a> detector in use: UVIS or IR
          b> detector in use: WFC, HRC, or SBC
       Keyword DGEOFILE has different values:
          a> N/A
          b> jref$qbu16429j_dxy.fits
       Keyword DISTNAME has different values:
          a> iczgs3ygq_w3m18525i-NOMODEL-NOMODEL
          b> jczgx1ppq_0461802dj-02c1450rj-02c1450oj
       Keyword DRIZCORR has different values:
          a> COMPLETE
          b> PERFORM
       Keyword DRKCFILE has different values:
          a> N/A
          b> jref$0a41546rj_dkc.fits
       Keyword DRKCFILE has different comments:
          b> De-trailed Dark Reference File
       Keyword EXPEND   has different values:
          a> 57402.29030181
          b> 57677.05173856
       Keyword EXPSTART has different values:
          a> 57402.28332292
          b> 57677.04503644
       Keyword EXPTIME  has different values:
          a> 602.937317
          b> 578.0
       Keyword FILENAME has different values:
          a> iczgs3ygq_flt.fits
          b> jczgx1ppq_flc.fits
       Keyword FLSHFILE has different comments:
          b> post flash correction file name
       Keyword IDCTAB   has different values:
          a> iref$w3m18525i_idc.fits
          b> jref$0461802dj_idc.fits
       Keyword IMPHTTAB has different values:
          a> iref$wbj1825ri_imp.fits
          b> jref$08b18470j_imp.fits
       Keyword INSTRUME has different values:
          a> WFC3
          b> ACS
       Keyword LINENUM  has different values:
          a> S3.008
          b> X1.009
       Keyword MDRIZTAB has different values:
          a> iref$ubi1853pi_mdz.fits
           ? ^      ^ -  ^^
          b> jref$ub21537aj_mdz.fits
           ? ^      ^   ^^^
       Keyword MOONANGL has different values:
          a> 57.153374
          b> 92.141869
       Keyword NEXTEND  has different values:
          a> 6
          b> 15
       Keyword NPOLFILE has different values:
          a> N/A
          b> jref$02c1450rj_npl.fits
       Keyword NPOLFILE has different comments:
          b> Non-polynomial Offsets Reference File
       Keyword OBSMODE  has different values:
          a> MULTIACCUM
          b> ACCUM
       Keyword OPUS_VER has different values:
          a> HSTDP 2015_3a
           ?          ^ ^
          b> HSTDP 2016_1a
           ?          ^ ^
       Keyword OSCNTAB  has different values:
          a> iref$q911321mi_osc.fits
          b> jref$lch1459bj_osc.fits
       Keyword OSCNTAB  has different comments:
          a> detector overscan table
          b> CCD overscan table
       Keyword PA_V3    has different values:
          a> 282.776093
          b> 88.003448
       Keyword PCTETAB  has different values:
          a> N/A
          b> jref$xa81724cj_cte.fits
       Keyword PCTETAB  has different comments:
          b> CTE Correction Table
       Keyword PFLTFILE has different values:
          a> iref$uc721143i_pfl.fits
          b> jref$qb12257pj_pfl.fits
       Keyword PROCTIME has different values:
          a> 57541.85384259
          b> 57677.3616088
       Keyword PROPAPER has different values:
          b> WFCENTER
       Keyword PYWCSVER has different values:
          a> 1.1.2
          b> 1.2.1
       Keyword RA_TARG  has different values:
          a> 36.85374208875
          b> 127.7389583333
       Keyword READNSEA has different values:
          a> 20.200001
          b> 4.3499999
       Keyword READNSEB has different values:
          a> 19.799999
          b> 3.75
       Keyword READNSEC has different values:
          a> 19.9
          b> 4.0500002
       Keyword READNSED has different values:
          a> 20.1
          b> 5.0500002
       Keyword ROOTNAME has different values:
          a> iczgs3ygq
          b> jczgx1ppq
       Keyword RPTCORR  has different comments:
          a> combine individual repeat observations
           ? ^^^^^^^
          b> add individual repeat observations
           ? ^^^
       Keyword SIPNAME  has different values:
          a> iczgs3ygq_w3m18525i
          b> jczgx1ppq_0461802dj
       Keyword SUNANGLE has different values:
          a> 112.720184
          b> 91.557938
       Keyword SUN_ALT  has different values:
          a> 3.227515
          b> 54.863163
       Keyword TARGNAME has different values:
          a> ANY
          b> ACO-665
       Keyword TIME-OBS has different values:
          a> 06:47:59
          b> 01:04:51
       Keyword T_SGSTAR has different values:
          b> N/A
    




stfhistory
----------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The stfhistory task will read history information from a text file and
add it to an image header. Here we will show how to do this with a FITS
file using Python's built in i/o functionality and the
``astropy.io.fits`` package.

.. code:: ipython3

    # Standard Imports
    import shutil
    
    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython3

    # open our text file and fits file objects, we're going to make a copy of a fits file, and edit the copy
    my_file = open('/eng/ssb/iraf_transition/test_data/history_info.txt', 'r')
    shutil.copyfile('/eng/ssb/iraf_transition/test_data/stfhist.fits','stfhist_copy.fits')
    test_data = fits.open('stfhist_copy.fits', mode='update')
    
    # loop through lines in text file and write to fits file
    # here we add the HISTORY lines to the zeroth header
    for line in my_file:
        test_data[0].header.add_history(line.strip('\n'))
        
    # make sure to close your files after the edits are done
    test_data.close()
    my_file.close()



upreffile
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress





Not Replacing
-------------

-  groupmod - GEIS header editing. Deprecated, for FITS header editing
   see **images.imutil.hedit**
-  hcheck - see **images.imutil.hselect**
-  iminfo - see **images.imutil.imheader**
