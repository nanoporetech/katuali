
.. _installation:

Installation
============

A Makefile is provided to create a `virtual environment
<https://docs.python.org/3/tutorial/venv.html>`_ into which ``katuali`` will be installed. 

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

General Tools
^^^^^^^^^^^^^

.. parsed-literal::

    `pomoxis <https://github.com/nanoporetech/pomoxis>`_: "~/git/pomoxis/venv/bin/activate"

Pomoxis is used as a general utility throughout much of katuali, its
installation is required for most functionality. 

Pomoxis >= 0.3.2 is required for the medaka variant-calling training pipelines. 


Basecalling Tools
^^^^^^^^^^^^^^^^^

.. parsed-literal::

    `guppy <https://community.nanoporetech.com/downloads>`_: "/usr/bin/guppy_basecaller"


Assembly Tools
^^^^^^^^^^^^^^

.. parsed-literal::

    `canu <https://github.com/marbl/canu>`_: "~/git/canu-1.8/Linux-amd64/bin/canu"
    `flye <https://github.com/fenderglass/Flye>`_: "~/git/Flye/bin/flye"
    `shasta <https://github.com/chanzuckerberg/shasta>`_: "~/git/shasta-Linux-0.3.0"
    `pomoxis <https://github.com/nanoporetech/pomoxis>`_: "~/git/pomoxis/venv/bin/activate"

In the context of Katuali, pomoxis functions as a wrapper to
`miniasm <https://github.com/lh3/miniasm>`_. Canu is utilised in the standard
workflow, while pomoxis is used in the fast workflow.


Consensus Tools
^^^^^^^^^^^^^^^

.. parsed-literal::

    `medaka <https://github.com/nanoporetech/medaka>`_: "~/git/medaka/venv/bin/activate"

.. note:: Note that ``tensorflow-gpu`` must be installed in your medaka
    environment if you wish to run medaka using a GPU (required for medaka
    training, optional for running medaka consensus)

Katuali can also use `racon <https://github.com/lbcb-sci/racon>`_ which is included with
pomoxis.


Other tools
^^^^^^^^^^^

.. parsed-literal::

    `liftover <https://anaconda.org/bioconda/ucsc-liftover>`_: "~/miniconda3/envs/liftover/bin/liftOver"
    `vcfcreatemulti <https://github.com/vcflib/vcflib>`_: "~/git/vcflib/bin/vcfcreatemulti"

LiftOver and vcfcreatemulti are used in the medaka variant-calling training pipelines and only
needs to be installed if you intend to use these pipelines.
