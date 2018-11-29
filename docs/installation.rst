
.. _installation:

Installation
============

`Katuali` is a `Snakemake <https://snakemake.readthedocs.io>`_ pipeline comprising a `Snakefile <https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#>`_ and `config <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`_.

As such, all that is required to run the pipeline is `Snakemake`. 

A Makefile is provided to create a `virtual environment
<https://docs.python.org/3/tutorial/venv.html>`_ into which to install snakemake as well the katuali convenience wrapper. 

Katuali has been tested on Linux, specifically Ubuntu 16.

To setup the environment run:

.. code-block:: bash

    git clone https://github.com/nanoporetech/katuali.git
    cd katuali
    make install
    source ./venv/bin/activate


.. _dependencies:

Dependencies
------------

Katuali makes use of a number of tools to do basecalling, assembly and
polishing that you will need to install.  

You only need to install what you will use.

Tools are inputs to the targets; if something can't be found, Snakemake will tell you it is missing.

The following config default parameters can be changed in the config or on the command line to point to your installations of these tools: 

    * `SCRAPPIE <https://github.com/nanoporetech/scrappie>`_: "~/git/scrappie"
    * `FLAPPIE <https://github.com/nanoporetech/flappie>`_: "~/git/github/flappie"
    * GUPPY: "/usr/bin/guppy_basecaller"
    * `IN_POMOXIS <https://github.com/nanoporetech/pomoxis>`_: "~/git/pomoxis/venv/bin/activate"
    * `CANU_EXEC <https://github.com/marbl/canu>`_: "~/git/canu-1.7.1/Linux-amd64/bin/canu"
    * `NANOPOLISH <https://github.com/jts/nanopolish>`_: "~/git/nanopolish"
    * `IN_MEDAKA <https://github.com/nanoporetech/medaka>`_: "~/git/medaka/venv/bin/activate"
    * IN_MIYAGI: "~/git/miyagi/venv/bin/activate"

Please refer to the documentation of each of these tools for installation instructions.


.. _tests:

Tests
-----

The easiest way to test the pipeline is to run the tests, which will basecall,
assemble and polish a small dataset that comes bundled with `Katuali`. 
The tests require scrappie, pomoxis, canu and nanopolish to be installed, and can be run with:

.. code-block:: bash

    make test

