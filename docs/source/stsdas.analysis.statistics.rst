:orphan:


stsdas.analysis.statistics
==========================

The statictics package contains statistical analysis tasks.

Notes
-----

**For questions or comments please see** `our github
page <https://github.com/spacetelescope/stak>`__. **We encourage and
appreciate user feedback.**

Contents:

-  `bhkmethod <#bhkmethod>`__
-  `buckleyjames-kmestimate <#buckleyjames-kmestimate>`__
-  `coxharzard <#coxharzard>`__
-  `kolmov <#kolmov>`__
-  `spearman <#spearman>`__
-  `twosampt <#twosampt>`__





bhkmethod
---------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The bhkmethod task is used to compute the generalized Kendall's tau
correlation coefficient. We show a short example here taken from the
`scipy.stats.kendalltau <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kendalltau.html>`__
documentation.

.. code:: ipython2

    # Standard Imports
    from scipy import stats

.. code:: ipython2

    x1 = [12, 2, 1, 12, 2]
    x2 = [1, 4, 7, 1, 0]
    tau, p_value = stats.kendalltau(x1, x2)
    print("tau: {}".format(tau))
    print("p_value: {}".format(tau))


.. parsed-literal::

    tau: -0.471404520791
    p_value: -0.471404520791




buckleyjames-kmestimate
-----------------------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The buckleyjames and kestimate tasks compute linear regression
coefficients and esitmators with the Kaplan-Meier estimator. There is
currently a Python package called ``lifelines`` that `have this
fitter <http://lifelines.readthedocs.io/en/latest/Quickstart.html#kaplan-meier-and-nelson-aalen>`__.



coxharzard
----------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The coxhazard task is used to compute the correlation probability by
Cox's proportional hazard model. The ``lifelines`` package contains
`this fitter
here <https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html#cox-s-proportional-hazard-model>`__.



kolmov
------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The kolmov task uses the Kolmogorov-Smirnov test for goodness of fit.
You can find both the one-sided and two-sided test in ``scipy``:

-  `one-sided <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ksone.html#scipy.stats.ksone>`__
-  `two-sided <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.kstwobign.html#scipy.stats.kstwobign>`__



spearman
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The spearman task is used to compute regression coefficients by Scmitt's
method. ``Scipy`` contains a version of `this
task <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.spearmanr.html#scipy.stats.spearmanr>`__.

.. code:: ipython2

    # Standard Imports
    from scipy import stats

.. code:: ipython2

    rho, pvalue = stats.spearmanr([1,2,3,4,5],[5,6,7,8,7])
    print("rho: {}".format(rho))
    print("p-value: {}".format(pvalue))


.. parsed-literal::

    rho: 0.820782681668
    p-value: 0.0885870053135




twosampt
--------

**Please review the** `Notes <#notes>`__ **section above before running
any examples in this notebook**

The twosampt task is used to determine if two sets of data are from the
same population. It provided the following types of two sample test:
geham-permute, gehan-hyper, logrank, peto-peto, and peto-prentice. These
tests do not currently have an equivalent in Scipy, but the following
two sample tests are availalbe:

-  `Ranksums <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ranksums.html>`__
-  `Wilcoxon <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.wilcoxon.html>`__
-  `Man-Whitney <https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.mannwhitneyu.html#scipy.stats.mannwhitneyu>`__





Not Replacing
-------------

-  censor - Information about the censoring indicator in survival
   analysis. Deprecated.
-  emmethod - Compute linear regression for censored data by EM method.
   Deprecated.
-  schmittbin - Compute regression coefficients by Schmitt's method.
   Deprecated.
-  survival - Provide background & overview of survival analysis.
   Deprecated.
