
.. _installation:

Installation
============

`Katuali` is a `Snakemake <https://snakemake.readthedocs.io>`_ pipeline comprising a `Snakefile <https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#>`_ and `config <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`_.

As such, all that is required to run the pipeline is `Snakemake`. 

A Makefile is provided to create a `virtual environment
<https://docs.python.org/3/tutorial/venv.html>`_ into which all necessary dependencies will be installed. 

Katuali has been tested on Linux (specifically Ubuntu 16) and MacOS.

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
polishing that need to be installed to perform such steps of a pipeline. 
Katuali does not enforce that all possible tools are present; only those
required to run a user's pipeline need be installed. If Katuali cannot
find a tool, it will alert the user. 

The list below indicates possible options for pipeline steps and the default
filesystem location where Katuali expects to find them. These locations can be
changed by the user, see the
'documentation' for more information.
Please refer to the documentation of each of these tools for installation
instructions.
The default config parameters can be changed in the config or on the command 
line to point to your installations of these tools. 

Basecalling Tools
^^^^^^^^^^^^^^^^^

    * `GUPPY <https://community.nanoporetech.com/downloads>`_: "/usr/bin/guppy_basecaller"
    * `SCRAPPIE <https://github.com/nanoporetech/scrappie>`_: "~/git/scrappie"
    * `FLAPPIE <https://github.com/nanoporetech/flappie>`_: "~/git/github/flappie"

    Guppy is recommended, both scrappie and flappie are research prototypes

Assembly Tools
^^^^^^^^^^^^^^

    * `IN_POMOXIS <https://github.com/nanoporetech/pomoxis>`_: "~/git/pomoxis/venv/bin/activate"
    * `CANU_EXEC <https://github.com/marbl/canu>`_: "~/git/canu-1.7.1/Linux-amd64/bin/canu"

Consensus Tools
^^^^^^^^^^^^^^^

    * `NANOPOLISH <https://github.com/jts/nanopolish>`_: "~/git/nanopolish"
    * `IN_MEDAKA <https://github.com/nanoporetech/medaka>`_: "~/git/medaka/venv/bin/activate"

Medeka is recommended in combination witht he lastest "flip-flop" algorithm in the guppy basecaller, 
and for rapid consensus. Nanopolish may still by preferred if using older basecallers.  


