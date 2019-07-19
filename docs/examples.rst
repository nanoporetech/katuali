
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
which basecaller or assembler to use as well as parameters for these tools.

For example, while

.. code-block:: bash

    katuali run1/basecall/guppy/basecalls.fasta

will basecall with guppy. 

.. code-block:: bash

    katuali run1/basecall/scrappie/basecalls.fasta

will perform basecalling with scrappie (starting from `./run1/reads`), 

Each step of a multi-step pipeline stores its outputs in a subdirectory of the
previous stage. In this way, we can define a multistep-pipeline to basecall
with guppy, assemble with canu, then polish with medaka, and finally run
nanopolish by contructing the target: 

.. code-block:: bash

    katuali run1/basecall/guppy/canu_gsz_4.0M/medaka/nanopolish/consensus.fasta

This nested working directory stores the data in such a way that is it obvious
what went into and out of each stage of the pipeline.

It also enables forks in the pipelines which might share basecalls (or other
intermediate results), but differ in assembly or consensus methods.

For example, if we wish to basecall with guppy, assemble with canu, run
medaka on the canu assembly, and seperately run nanopolish on the canu assembly,
we could use the targets: 

.. code-block:: bash

    katuali run1/basecall/guppy/canu_gsz_4.0M/nanopolish/consensus.fasta run1/basecall/guppy/canu_gsz_4.0M/medaka/consensus.fasta

`Snakemake` will create a graph of tasks to perform the common basecall
and assembly tasks, then run separately nanopolish and medaka from the same
inputs (and in parallel, given enough resource).


Basecalling
-----------

Supported basecallers include guppy, scrappie, flappie:

.. code-block:: bash

    katuali run1/basecall/guppy/basecalls.fasta
    katuali run1/basecall/scrappie/basecalls.fasta
    katuali run1/basecall/flappie/basecalls.fasta


Assembly
--------

Reads can be assembled in three ways at present:

.. code-block:: bash

    # assemble with canu, specifying the genome size in the target name. 
    katuali run1/basecall/scrappie/canu_gsz_4.0M/consensus.fasta  

    # use pomoxis mini_assemble to assemble with miniasm, then form consensus
    # with racon
    katuali run1/basecall/scrappie/miniasm_racon/consensus.fasta  

    # assemble with flye, specifying the genome size in the target name. 
    katuali run1/basecall/scrappie/flye_gsz_4.0M/consensus.fasta


Polishing
---------

Sequences can be polished with Racon, Nanopolish and medaka to create higher
accuracy consensus sequences. Consensus methods can also be combined (e.g.
racon/medaka) meaning that the input to medaka will be the racon consensus. 
The last example requests two rounds of medaka. 

.. code-block:: bash

    katuali run1/basecall/guppy_flipflop/canu_gsz_4.0M/racon/nanopolish/consensus.fasta
    katuali run1/basecall/guppy_flipflop/canu_gsz_4.0M/racon/medaka_flipflop/consensus.fasta
    katuali run1/basecall/guppy_flipflop/canu_gsz_4.0M/racon/medaka_flipflop/medaka_flipflop/consensus.fasta


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
variable(s), for example datasets, regions, basecallers, assemblers, you can get
katuali to automatically generate all the targets for you by creating a
template of the target containing named
placeholders of the config variable(s) that will be looped over. 

The fast_assm_polish workflow is implemented with the following target template:

.. code-block:: yaml

    PIPELINES:
        fast_assm_polish: [
            "{DATA}/basecall/{BASECALLER}{BASECALLER_SUFFIX}/miniasm_racon/medaka{BASECALLER_SUFFIX}/consensus.fasta"
        ]

Running

.. code-block:: bash

    katuali fast_assm_polish

will expand all the variables in the target template. ``{DATA}`` will be expanded
to all the datsets defined in ``config[DATA]``. As ``{BASECALLER}`` and
``{BASECALLER_SUFFIX}`` are single-entries in the config (rather than being a
list of strings), their placeholders are simply replaced with their values.

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

As an example, the ``medaka_eval`` pipeline contains two target templates, to
evaluate both the pre- and post-medaka consensus accuracy, in this case over a range of
datasets, regions, depths, and medaka models, generating hundreds of targets in the process. 

.. code-block:: yaml

    PIPELINES:
        medaka_eval: [
            "{DATA}/basecall/{BASECALLER}{BASECALLER_SUFFIX}/align/{MEDAKA_EVAL_REGIONS}/{DEPTHS}X/canu_gsz_{GENOME_SIZE}/racon/medaka{MEDAKA_EVAL_SUFFIXES}/consensus_to_truth_summ.txt",
            "{DATA}/basecall/{BASECALLER}{BASECALLER_SUFFIX}/align/{MEDAKA_EVAL_REGIONS}/{DEPTHS}X/canu_gsz_{GENOME_SIZE}/racon/consensus_to_truth_summ.txt"
        ]

The final step of each pipeline is to create an empty file with the name of the
pipeline (e.g. ``standard_assm_polish``) which indicates the pipeline has
finished.  If you wish to rerun the pipeline after changing config variables
which affect the pipeline targets, the empty file needs to
be deleted before rerunning the pipeline; without deleting it, `katuali` will
not create the new targets.


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
    BCDIR=${RUN}/basecall/${BASECALLER}/
    mkdir -p ${BCDIR}
    mkdir ${RUN}/reads
    ln -s ${SUMMARY} ${BCDIR}/sequencing_summary.txt

    source ${IN_POMOXIS}
    seqkit fq2fa ${BASECALLS} > ${BCDIR}/basecalls.fasta

Now katuali can be run as normal, for example:

.. code-block:: bash

    katuali --configfile my_config.yaml standard_assm_polish


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

    katuali run1/basecall/scrappie/align/depth

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

    katuali run1/basecall/scrappie/align/all_contigs/25X/miniasm_racon/consensus.fasta

will perform the following steps:

    * basecall the reads to create:
      `run1/basecall/scrappie/basecalls.fasta`
    * align the basecalls to the reference to create:
      `run1/basecall/scrappie/align/calls2ref.bam`
    * subsample all contigs in the .bam file to 25X to create (in one step):
      `run1/basecall/scrappie/align/all_contigs/25X/basecalls.fasta`
    * perform a ref-guided assembly and racon consensus to create:
      `run1/basecall/scrappie/align/all_contigs/25X/miniasm_racon/consensus.fasta`


.. note:: The rule to create subsampled datasets differs from other rules in
    that it creates two levels of nested directories in a single step (in this case
    `all_contigs/25X`). The extraction of specific regions/contigs without
    subsampling to a specific depth is not currently supported.


Subsampling a single reference contig
-------------------------------------

It is also possible to subsample just one of the contigs in your reference by
specifying targets such as:

.. code-block:: bash

    katuali run1/basecall/scrappie/align/ecoli_SCS110_plasmid2/25X/miniasm_racon/consensus.fasta 

which will just process the reference sequence `ecoli_SCS110_plasmid2`.


Subsampling specified regions
-----------------------------

It is also possible to subsample only specified regions specifed as samtools
strings:

.. code-block:: bash

    REGIONS="ecoli_SCS110_chromosome:50000-150000 ecoli_SCS110_chromosome:200000-250000"
    katuali run1/basecall/scrappie/align/my_regions/25X/miniasm_racon/consensus.fasta --config REGIONS="$REGIONS"


.. _train_medaka:

Medaka training pipeline
------------------------

It is possible to train medaka models starting from folders of fast5s in a
single command once the config has been modified to reflect your input data
(fast5s and genomes for each run as well as training and evaluation region
definitions).

`MEDAKA_TRAIN_REGIONS` and `MEDAKA_EVAL_REGIONS` define regions for training
and evaluation.  In the example below we train from the `minion` run using
`ecoli` and `yeast` contigs in the reference and evaluate on the `gridion` run
using the contigs `ecoli`, `yeast` and `na12878_chr21` in the reference.

.. code-block:: yaml

    DATA:
        'MinIonRun1': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': ['ecoli', 'yeast']
            'MEDAKA_EVAL_REGIONS': []
        'MinIonRun2': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': ['ecoli', 'yeast']
            'MEDAKA_EVAL_REGIONS': []
        'GridIonRun1': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': []
            'MEDAKA_EVAL_REGIONS': ['ecoli', 'yeast', 'na12878_chr21']
        'GridIonRun2': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': []
            'MEDAKA_EVAL_REGIONS': ['ecoli', 'yeast', 'na12878_chr21']

Running:

.. code-block:: bash

    katuali medaka_train_features --keep-going

will:

* basecall all the runs
* align each run to its reference
* create subsampled sets of basecalls over the desired regions and depths
* assemble those sets of basecalls
* create medaka training features for all those sets

Running:

.. code-block:: bash

    katuali medaka_train_replicates --keep-going

will do all the tasks of ``medaka_train_features`` and additionally launch
multiple medaka model-training replicates.

If some of your input runs have insufficient coverage-depth for some of the
training regions, some of the training feature files will not be made. In this
case the config flag ``USE_ONLY_EXISTING_MEDAKA_FEAT`` can be set to true to allow
katuali to train using only those features which exist already:

.. code-block:: yaml

    USE_ONLY_EXISTING_MEDAKA_FEAT: true 

Refer to comments in the config (katuali/data/config.yaml) to see how this process
can be controlled. 
