
.. _introduction:

Custom pipelines
================

One of the main aims of `Katuali` is to enable the bolting-together of any
number of basecallers, assemblers, and polishers into flexible (i.e. not
predefined) pipelines. Basecallers, assemblers, and polishers can be
interchanged because their Snakefile rules have standardised inputs and
outputs.

`Katuali` makes strong use of `Snakemake`'s powerful `wildcards
<https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards>`_
to extract pipeline parameters from the target file path. For example
which assembler to use as well as parameters for these tools.

For example, while

.. code-block:: bash

    katuali run1/guppy/canu/consensus.fasta

will basecall with guppy then assemble with canu

.. code-block:: bash

    katuali run1/guppy/shasta/consensus.fasta

will basecall with guppy then perform assembly with shasta (starting from `./run1/reads`), 

Each step of a multi-step pipeline stores its outputs in a subdirectory of the
previous stage. In this way, we can define a multistep-pipeline to basecall
with guppy, assemble with canu, perform consensus with racon, and finally polish
with medaka:

.. code-block:: bash

    katuali run1/guppy/canu/racon/medaka/consensus.fasta

This nested working directory stores the data in such a way that is it obvious
what went into and out of each stage of the pipeline.

It also enables forks in the pipelines which might share basecalls (or other
intermediate results), but differ in assembly or consensus methods.

For example, if we wish to basecall with guppy, assemble with canu, run
racon followed by medaka on the canu assembly, and seperately run medaka directly on
the canu assembly, we could use the targets: 

.. code-block:: bash

    katuali run1/guppy/canu/racon/medaka/consensus.fasta \
            run1/guppy/canu/medaka/consensus.fasta

`Snakemake` will create a graph of tasks to perform the common basecall
and assembly tasks, then run separately two indepdendent medaka tasks from the same
inputs (and in parallel, given enough resource).


Basecalling
-----------

Supported basecallers include guppy.

.. code-block:: bash

    katuali run1/guppy/basecalls.fasta


Assembly
--------

Reads can be assembled in four ways at present:

.. code-block:: bash

    # assemble with canu
    katuali run1/guppy/canu/consensus.fasta  

    # use pomoxis mini_assemble to assemble with miniasm, then form consensus
    # with racon
    katuali run1/guppy/miniasm_racon/consensus.fasta  

    # assemble with flye 
    katuali run1/guppy/flye/consensus.fasta

    # assemble with shasta
    katuali run1/guppy/shasta/consensus.fasta


Polishing
---------

Sequences can be polished with Racon or medaka to create higher
accuracy consensus sequences. Consensus methods can also be combined (e.g.
racon/medaka) meaning that the input to medaka will be the racon consensus. 
The last example requests two rounds of medaka (something not generally
required or encouraged).

.. code-block:: bash

    katuali run1/guppy/canu/racon/consensus.fasta
    katuali run1/guppy/canu/racon/medaka/consensus.fasta
    katuali run1/guppy/canu/racon/medaka/medaka/consensus.fasta


Pipeline restrictions
---------------------

`Katuali` aims to be as flexible as possible, but there are some obvious
restrictions:

    * basecalling must be performed before assembly.
    * assembly must come before polishing (use of polishing targets to
      error correct reads is not supported).


Automatic generation of custom pipeline targets
-----------------------------------------------

If your pipeline involves the creation of many targets by looping over some
variable(s), for example datasets, regions, assemblers, you can get
katuali to automatically generate all the targets for you by creating a
template of the target containing named
placeholders of the config variable(s) that will be looped over. 

The fast_assm_polish workflow is implemented with the following target template:

.. code-block:: yaml

    PIPELINES:
        all_fast_assm_polish: [
            "{DATA}/guppy/miniasm/racon/medaka/consensus.fasta"
    ]


Running

.. code-block:: bash

    katuali all_fast_assm_polish

will expand all the variables in the target template. In this example ``{DATA}`` will be expanded
to all the datsets defined in ``config[DATA]``. 

You can use any config parameter as a placeholder, however there are some rules
concerning variables which are dataset-specific:


1. Dataset-specific variables are defined within the config section for that
   dataset (e.g. the ``MEDAKA_EVAL_REGIONS`` for dataset ``MinIonRun1`` are
   defined in ``config[DATA][MinIonRun1][MEDAKA_EVAL_REGIONS]``), so that pipelines can be 
   customised in a data-set specific way. 

2. The ``{GENOME_SIZE}`` placeholder, used to provide some assemblers an
   estimate of genome size,  is treated differently from other placholders. If
   the ``GENOME_SIZE`` variable is present in the config section of a dataset,
   this value will be used. However,
   if you have a reference and wish to assemble contigs independently (as is
   done in e.g. the medaka training pipeline), if ``config[DATA][MinIonRun1][GENOME_SIZE]``
   is not present, but ``config[DATA][MinIonRun1][REFERENCE]`` is present,
   ``{GENOME_SIZE}`` will be automatically calculated from the reference
   sequence for each of the contigs/regions defined for that dataset. Any
   placeholder containing the string ``REGION`` will be used in this way to
   calculate genome/region sizes.  The region definitions can be contig names
   or full samtools region strings with start and end. 

Config pipeline entries are lists so that multiple target templates can be used in a single pipeline. 

As an example, the ``all_consensus_eval`` pipeline contains two target templates, to
evaluate both the pre- and post-medaka consensus accuracy, in this case over a range of
datasets, regions, depths, and medaka models, generating hundreds of targets in the process. 

.. code-block:: yaml

    PIPELINES:
        all_consensus_eval: [
            "{DATA}/guppy/align/{MEDAKA_EVAL_REGIONS}/{DEPTHS}X/canu/racon/medaka{MEDAKA_EVAL_SUFFIXES}/assess{ASSESS_ASSM_SUFFIXES}/consensus_summ.txt",
            "{DATA}/guppy/align/{MEDAKA_EVAL_REGIONS}/{DEPTHS}X/canu/racon/assess{ASSESS_ASSM_SUFFIXES}/consensus_summ.txt"
        ]

The final step of each pipeline is to create an empty file with the name of the
pipeline (e.g. ``all_standard_assm_polish``) which indicates the pipeline has
finished.  If you rerun the pipeline this file will be automatically deleted
and recreated upon pipeline completion. 


.. _starting_from_basecalls:

Starting from existing basecalls
--------------------------------

If you have already basecalled your data, mocking out the working space as if
katuali had basecalled allows any derived targets to be created.

.. code-block:: bash
   
    # Input files
    BASECALLS=/path/to/basecalls.fastq
    SUMMARY=/path/to/sequencing_summary.txt

    # These should be set as in the config.yaml file used for running the
    workflow. RUN is # the top level key of the DATA section
    RUN=run1
    BASECALLER=guppy_flipflop
    IN_POMOXIS=~/git/pomoxis/venv/bin/activate

    # ...no need to edit below here
    BCDIR=${RUN}/${BASECALLER}/
    mkdir -p ${BCDIR}
    mkdir ${RUN}/reads
    ln -s ${SUMMARY} ${BCDIR}/sequencing_summary.txt

    source ${IN_POMOXIS}
    seqkit fq2fa ${BASECALLS} > ${BCDIR}/basecalls.fasta

Now katuali can be run as normal, for example:

.. code-block:: bash

    katuali --configfile my_config.yaml all_standard_assm_polish


Calculating read coverage depth
-------------------------------

It is often useful to know the read coverage depth of a dataset. 
This requires a reference.fasta to be specified in the config to which the reads will be aligned. 

.. code-block:: yaml

    DATA:
        'run1':
            'REFERENCE':/path/to/ref.fasta

The read coverage depth can then be calculated as follows: 

.. code-block:: bash

    katuali run1/guppy/align/depth

The depth directory will contain a text file per reference contig with coverage
vs genomic coordinate, as well as a file containing summary statistics for all
contigs.


Creating subsampled datasets
----------------------------

Katuali also supports the generation of datasets with even coverage over a
reference at a given depth. 
This requires a reference.fasta to be specified in the config to which the reads will be aligned. 

.. code-block:: yaml

    DATA:
        'run1':
            'REFERENCE':/path/to/ref.fasta

Once the reference is the config, running:

.. code-block:: bash

    katuali run1/guppy/align/all_contigs/25X/miniasm_racon/consensus.fasta

will perform the following steps:

    * basecall the reads to create:
      `run1/guppy/basecalls.fasta`
    * align the basecalls to the reference to create:
      `run1/guppy/align/calls2ref.bam`
    * subsample all contigs in the .bam file to 25X to create (in one step):
      `run1/guppy/align/all_contigs/25X/basecalls.fasta`
    * perform a ref-guided assembly and racon consensus to create:
      `run1/guppy/align/all_contigs/25X/miniasm_racon/consensus.fasta`


.. note:: The rule to create subsampled datasets differs from other rules in
    that it creates two levels of nested directories in a single step (in this case
    `all_contigs/25X`). The extraction of specific regions/contigs without
    subsampling to a specific depth is not currently supported.


Subsampling a single reference contig
-------------------------------------

It is also possible to subsample just one of the contigs in your reference by
specifying targets such as:

.. code-block:: bash

    katuali run1/guppy/align/ecoli_SCS110_plasmid2/25X/miniasm_racon/consensus.fasta 

which will just process the reference sequence `ecoli_SCS110_plasmid2`.


Subsampling a specified region
------------------------------

It is also possible to subsample a region specifed as a samtools string:

.. code-block:: bash

    katuali run1/guppy/align/ecoli_SCS110_chromosome:50000-150000/25X/miniasm_racon/consensus.fasta


.. _train_medaka:

Medaka training pipelines
-------------------------

It is possible to train both medaka consensus models and medaka variant-calling models starting from folders of fast5s. 
See :ref:`medaka_train` and :ref:`medaka_train_variant`. 
