=======================
Old Fisher matrix codes
=======================

This is an old (mostly historical) set of tools for Fisher matrices used in the FoMSWG (2009), WFIRST concept studies, and up through the Science Requirements Document (2018). It is no longer maintained, but is being provided here, lightly edited, for comparison with those older studies.

Compilation
===========

You can compile the BAO tools here by running the ``make`` command, producing the file ``test-b.x``.


Usage
=====

I recommend you make a ``data/`` subdirectory to store your datafiles. A BAO/RSD datafile is of the form:

.. code-block:: rst

    46 0.065449785 0.001
    0.550	0.050	1.385	3.568E-04
    0.600	0.050	1.420	4.720E-04
    0.650	0.050	1.455	5.673E-04
    ...

Here the first line contains the number of lines of data (number of redshift bins); the fractional sky coverage (e.g., ``0.5`` for half of the sky or 20626 deg^2); and the redshift uncertainty sigma(z)/(1+z).

Then each subsequent line contains the redshift z; the redshift width delta z of the slice; the galaxy bias b; and the number density n in units of Mpc^-3 (**no** *h*!!).

Then you can call the BAO code (based on the Seo & Eisenstein 2007 forecast algorithm):

.. code-block:: rst

    ./test-b.x data/my_infile.txt /dev/null

and get output of the form:

.. code-block:: rst

    0.550  5.97879e+08  2.58092e+00 0.041320 0.080340  0.40918 0.029041
    0.600  6.71677e+08  3.41558e+00 0.036002 0.070802  0.40880 0.025392
    0.650  7.44072e+08  4.10446e+00 0.032551 0.064354  0.40859 0.022996
    ...

In this output, the columns are redshift; comoving volume (Mpc^3, again no *h*); nP at 0.2 *h*/Mpc; fractional distance uncertainty; fractional Hubble uncertainty; correlation coefficient; and fractional total scale uncertainty.

The RSD code may be called via

.. code-block:: rst

    perl red_driver.pl < data/my_infile.txt

and you will get output of the form:

.. code-block:: rst

    0.550	0.050	1.385	3.568E-04	0.04291
    0.600	0.050	1.420	4.720E-04	0.03859
    0.650	0.050	1.455	5.673E-04	0.03575
    ...

The first 4 columns are the inputs: z, delta z, bias, and number density n in units of Mpc^-3 (**no** *h*!!). The last column is the fractional uncertainty in f sigma_8 with all the cosmological parameters fixed (but galaxy bias floating) and scale cut kmax = 0.2 *h*/Mpc (the latter is currently hard-coded in ``rsdmat.py``). This was designed for consistency with the DESI forecasts (DESI Collaboration, https://arxiv.org/abs/1611.00036, Table 2.3, last column, although note that the z binning this code outputs is whatever you give it).

Notes
=====

The ``nr_utils.c`` and ``nr_utils.h`` are slight modifications of public domain Numerical Recipes code. See: https://numerical.recipes
