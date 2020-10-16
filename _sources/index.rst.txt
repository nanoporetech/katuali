Katuali
=======

`Katuali` is a flexible consensus pipeline implemented in 
`Snakemake <https://snakemake.readthedocs.io>`_ to basecall, assemble, and
polish Oxford Nanopore Technologies' sequencing data.

Features
--------

  * fast5 to high quality consensus in a single command.
  * Recommended fixed `standard` and `fast` pipelines.
  * Interchange basecaller, assembler, and consensus components of the
    pipelines simply by changing the target filepath. 
  * Medaka training pipeline including generation of training data, model
    training and model evaluation. 
  * Seemless distribution of tasks over local or distributed compute.
  * Open source (Mozilla Public License 2.0).

.. admonition:: Research Release

    Research releases are provided as technology demonstrators to provide early
    access to features or stimulate Community development of tools. Support for
    this software will be minimal and is only provided directly by the developers.
    Feature requests, improvements, and discussions are welcome and can be
    implemented by forking and pull requests. However much as we would
    like to rectify every issue and piece of feedback users may have, the 
    developers may have limited resource for support of this software. Research
    releases may be unstable and subject to rapid iteration by Oxford Nanopore
    Technologies.


.. _quickstart:

Quickstart
----------

The `Katuali` :ref:`tests` contain examples of how to basecall,
assemble, and polish a small dataset that comes bundled with `Katuali`. 

To run with other data, start by creating a directory of reads (which could
contain subdirectories of reads) within a run directory (`run1` in this example):

.. code-block:: bash

    mkdir -p run1
    cd run1
    ln -s /path/to/fast5 reads  # create a softlink to the fast5 data
    cd ..
    
Then make a copy of the katuali config into your working directory;

.. code-block:: bash

    katuali_config my_config.yaml

and update the katuali config to reflect your data:

.. code-block:: yaml

    DATA:
        'run1':
            'GENOME_SIZE': '4.0M'  # for canu/flye we need to specify genome size

There are two predefined pipelines that can be used starting from fast5 input:

1. To basecall the reads, assemble them with miniasm, and polish the assembly with
   racon and medaka simply run: 

   .. code-block:: bash
  
       katuali fast_assm_polish


2. To basecall, assemble with canu then polish with racon and medaka run: 

   .. code-block:: bash
  
       katuali standard_assm_polish


See :ref:`introduction` for details on creating flexible multistep pipelines.


Table of contents
-----------------

.. toctree::
   :maxdepth: 1

   installation
   tests
   examples
   medaka_train
   medaka_train_variant
   configuration
   faq


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
