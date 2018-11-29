Katuali
=======

`Katuali` is a flexible consensus pipeline implemented in `Snakemake <https://snakemake.readthedocs.io>`_ to basecall, assemble and polish 
Oxford Nanopore Technologies' sequencing data.

Features
--------

  * Run a pipeline processing fast5s to a consensus in a single command.
  * Recommended fixed `standard` and `fast` pipelines.
  * Interchange basecaller, assembler and consensus components of the pipelines simply by changing the target filepath. 
  * Seemless distribution of tasks over local or distributed compute.
  * Highly configurable.  
  * Open source (Mozilla Public License 2.0).


.. _quickstart:

Quickstart
----------

The `Katuali` :ref:`tests` contain examples of how to basecall,
assemble and polish a small dataset that comes bundled with `Katuali`. 

To run with your own data, start by creating a dictory of reads (which could contain subdirectories of reads):

.. code-block:: bash

    ln -s /path/to/fast5 reads
    
Then you can then make any number of targets the pipline knows how to make by running e.g.:

.. code-block:: bash

    katuali fast_assm_polish

This will basecall the reads, then assemble them with miniasm, and polish with racon and medaka. 

Running

.. code-block:: bash

    katuali standard_assm_polish

will instead basecall, assemble with canu and the polish with nanopolish. 


Table of contents
-----------------

.. toctree::
   :maxdepth: 2

   installation
   examples


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
