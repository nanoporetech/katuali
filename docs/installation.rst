
.. _installation:

Installation
============

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

Katuali makes use of a number of tools to perform basecalling, assembly, and
polishing that need to be installed in order to perform such steps of a pipeline.
Katuali does not enforce that all possible tools are present; only those
required to run a user's pipeline need be installed. If Katuali cannot find a
tool, it will alert the user. 

The list below indicates possible options for pipeline steps and the default
filesystem location where Katuali expects to find them. Each item is shown as
must be specified in the pipeline configuration, the locations can be
changed by the user as appropriate.

Refer to the documentation of each of these tools for installation instructions.

Basecalling Tools
^^^^^^^^^^^^^^^^^

.. parsed-literal::

    `GUPPY <https://community.nanoporetech.com/downloads>`_: "/usr/bin/guppy_basecaller"
    `SCRAPPIE <https://github.com/nanoporetech/scrappie>`_: "~/git/scrappie"
    `FLAPPIE <https://github.com/nanoporetech/flappie>`_: "~/git/github/flappie"

Guppy is recommended, both scrappie and flappie are research prototypes.

Assembly Tools
^^^^^^^^^^^^^^

.. parsed-literal::

    `IN_POMOXIS <https://github.com/nanoporetech/pomoxis>`_: "~/git/pomoxis/venv/bin/activate"
    `CANU_EXEC <https://github.com/marbl/canu>`_: "~/git/canu-1.8/Linux-amd64/bin/canu"

In the context of Katuali, pomoxis functions as a wrapper to
`miniasm <https://github.com/lh3/miniasm>`_. Canu is utilised in the standard
workflow, while pomoxis is used in the fast workflow.

Consensus Tools
^^^^^^^^^^^^^^^

.. parsed-literal::

    `NANOPOLISH <https://github.com/jts/nanopolish>`_: "~/git/nanopolish"
    `IN_MEDAKA <https://github.com/nanoporetech/medaka>`_: "~/git/medaka/venv/bin/activate"

Medeka is recommended in combination with he lastest "flip-flop" algorithm in
the guppy basecaller, and for rapid consensus. Nanopolish may still by
preferred if using older basecallers.  


