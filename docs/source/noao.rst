:orphan:


noao
====

Top level of the noao package, only contains the task observatory.

Notes
-----

**For questions or comments please see** `our github
page <https://github.com/spacetelescope/stak>`__. **We encourage and
appreciate user feedback.**

Observatory is the only task on the top noao module level. See the
Python example below.



observatory
-----------

The observatory task will return various information about an
observatory location. This task is included in ``Astropy`` as the
``astropy.coordinates.EarthLocation`` class. See the doc pages
`here <http://docs.astropy.org/en/stable/coordinates/index.html#convenience-methods>`__
and
`here <http://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html#astropy.coordinates.EarthLocation>`__
for more information.

.. code:: ipython2

    # Astronomy Specific Imports
    from astropy.coordinates import EarthLocation

.. code:: ipython2

    EarthLocation.of_site('Apache Point Observatory')




.. math::

    (-1463969.3, -5166673.3, 3434985.7) \; \mathrm{m}




