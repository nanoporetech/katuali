
.. _tests:

Basic Usage and Tests
=====================

The easiest way to test the pipeline is to run the tests, which will basecall,
assemble and polish a small dataset that comes bundled with `Katuali`.  The
tests require suppy, pomoxis, canu, flye and medaka to be installed, and can
be run with:

.. code-block:: bash

    make test

The output is placed in nested directories under the folder `test/`.  For
example, if the test has run correctly a `.bam` alignment file will have been
produced containing 25X coverage of reads basecalled with guppy aligning to and
E.coli plasmid.

.. code-block:: bash

    test/                                  
      MinIonRun1/              # dataset name as defined in config. 
        guppy/            # guppy basecaller with default options
          align/             # alignment of bases
            all_contigs/     # extraction of all_contigs in alignment  
              25X/           # subsampling alignments
                sub_sample_25X_ecoli_SCS110_plasmid2.calls2ref.bam


Predefined Workflows
--------------------

`Katuali` comes with a number of predefined workflows. To use these with your
own data, start by creating a directory of reads (which could contain
subdirectories of reads) within a run directory:

.. code-block:: bash

    mkdir -p run1
    cd run1
    ln -s /path/to/fast5 reads  # create a softlink to the fast5 data
    cd ..
    
Then make a copy of the `katuali` config into your working directory,

.. code-block:: bash

    katuali_config my_config.yaml

and update this file to reflect your data:
    
.. code-block:: bash

    DATA:
        'run1':
            'GENOME_SIZE': '4.0M'  # for canu/flye we need to specify genome size

There are two standard workflows available:

1. To perform basecalling, a quick assembly with miniasm, and consensus with racon
   and medaka it is sufficient to run:
   
   .. code-block:: bash
  
       katuali --configfile my_config.yaml all_fast_assm_polish

2. Alternatively to assemble with canu/flye (depending on ASSEMBLER option in my_config.yaml) run:

   .. code-block:: bash
  
       katuali --configfile my_config.yaml all_standard_assm_polish


The :ref:`introduction` section describes how to create a pipeline with any
combination of basecallers, assemblers and polishers. 
