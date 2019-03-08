Katuali
=======

`Katuali` is a flexible consensus pipeline implemented in 
`Snakemake <https://snakemake.readthedocs.io>`_ to basecall, assemble, and
polish Oxford Nanopore Technologies' sequencing data.

Features
--------

  * Run a pipeline processing fast5s from multiple runs into multiple consensuses in a single command.
  * Recommended fixed `standard` and `fast` pipelines.
  * Interchange basecaller, assembler, and consensus components of the
    pipelines simply by changing the target filepath. 
  * Medaka training pipeline including generation of training data, model training and model evaluation. 
  * Seemless distribution of tasks over local or distributed compute.
  * Highly configurable.  
  * Open source (Mozilla Public License 2.0).


.. _quickstart:

Quickstart
----------

The `Katuali` :ref:`tests` contain examples of how to basecall,
assemble, and polish a small dataset that comes bundled with `Katuali`. 

To run with other data, start by creating a directory of reads (which could
contain subdirectories of reads) within a run directory (`run1` in this example):

.. code-block:: bash

    mkdir -p run1 && cd run1 && ln -s /path/to/fast5 reads && cd ../
    
Then make a copy of the katuali config into your working directory;

.. code-block:: bash

    cp ~/git/katuali/config.yaml .

and update the katuali config to reflect your data:

.. code-block:: yaml

    DATA:
        'run1':
            'GENOME_SIZE': '4.0M'  # for canu we need to specify genome size

Then calculate any of the outputs the pipeline knows how to make by running e.g.:

.. code-block:: bash

    katuali fast_assm_polish

This will basecall the reads, assemble them with miniasm, and polish the
assembly with racon and medaka. 

.. code-block:: bash

    katuali standard_assm_polish

will instead basecall (flipflop), assemble with canu, then polish with racon and medaka (flipflop). 


Table of contents
-----------------

.. toctree::
   :maxdepth: 2

   installation
   examples
   medaka_train
   configuration
   faq


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
