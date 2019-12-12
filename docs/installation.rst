
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

General Tools
^^^^^^^^^^^^^

.. parsed-literal::

    `pomoxis <https://github.com/nanoporetech/pomoxis>`_: "~/git/pomoxis/venv/bin/activate"

Pomoxis is used as a general utility throughout much of katuali, its
installation is required for most functionality.

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

Medaka is recommended in combination with the latest "flip-flop" algorithm in
the guppy basecaller, and for rapid consensus. Nanopolish may still by
preferred if using older basecallers. 


