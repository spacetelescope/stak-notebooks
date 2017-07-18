:orphan:


stsdas.toolbox.headers
======================

This package provides utilities for comparing and editing image headers.

Notes
-----

**For questions or comments please see** `our github
page <https://github.com/spacetelescope/stak>`__. **We encourage and
appreciate user feedback.**

Many of the headers tasks can be replaced with utilities in ``astropy``.

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

.. code:: ipython2

    # Standard Imports
    import numpy as np
    
    # Astronomy Specific Imports
    from astropy.io import fits
    
    # Plotting Imports/Setup
    import matplotlib.pyplot as plt
    %matplotlib inline

.. code:: ipython2

    # code goes here



hdiff
-----

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The hdiff task will take two FITS headers. This functionality has been
replaced and improved upon in ``astropy`` with the
``astropy.io.fits.Differs`` class, which can be easily called with the
``printdiff`` convenience function. **when is this going in?**

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits
    from astropy.io.fits import printdiff

.. code:: ipython2

    # plain header example
    # example ignoring HISTORY and COMMENT cards
    # example with image data!



stfhistory
----------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The stfhistory task will read history information from a text file and
adds it to an image header. Here we will show how to do this with a FITS
file using Python's build in i/o functionality and the
``astropy.io.fits`` package.

**check, do we need to trim newlines?**

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.io import fits

.. code:: ipython2

    # open our text file and fits file objects
    my_file = open('history_info.txt', 'r')
    my_fits = fits.open('my_image.fits', mode='append')
    
    # loop through lines in text file and write to fits file
    # here we add the HISTORY lines to the zeroth header
    for line in my_file:
        my_fits[0].header.add_history(line)
        
    # make sure to close your fits file after the edits are done
    my_fits.close()



upreffile
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

this is probably in crds

.. figure:: static/150pxblueconstuc.png
   :alt: Work in progress




Not Replacing
-------------

-  groupmod - GEIS header editing. Deprecated, for FITS header editing
   see **images.imutil.hedit**
-  hcheck - see **images.imutil.hselect**
-  iminfo - see **images.imutil.imheader**

