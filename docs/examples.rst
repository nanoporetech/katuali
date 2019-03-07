
.. _introduction:

Examples
========

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

    katuali run1/basecall/scrappie/basecalls.fasta

will perform basecalling with scrappie (starting from `./run1/reads`), 

.. code-block:: bash

    katuali run1/basecall/guppy/basecalls.fasta

will basecall with guppy. 

Each step of a multi-step pipeline stores its outputs in a subdirectory of the
previous stage. In this way, we can define a multistep-pipeline to basecall
with scrappie, assemble with canu, then polish with medaka, and finally run
nanopolish by contructing the target: 

.. code-block:: bash

    katuali run1/basecall/scrappie/canu_gsz_4.0M/medaka/nanopolish/consensus.fasta

This nested working directory stores the data in such a way that is it obvious
what went into and out of each stage of the pipeline.

It also enables forks in the pipelines which might share basecalls (or other
intermediate results), but differ in assembly or consensus methods.

For example, if we wish to basecall with scappie, assemble with canu, run
medaka on the canu assembly, and seperately run nanopolish on the canu assembly,
we could use the targets: 

.. code-block:: bash

    katuali run1/basecall/scrappie/canu_gsz_4.0M/nanopolish/consensus.fasta run1/basecall/scrappie/canu_gsz_4.0M/medaka/consensus.fasta

`Snakemake` will create a graph of tasks to perform the common basecall
and assembly tasks, then run separately nanopolish and medaka from the same
inputs (and in parallel, given enough resource).


Basecalling
-----------

Supported basecallers include scrappie, flappie and guppy:

.. code-block:: bash

    katuali run1/basecall/scrappie/basecalls.fasta
    katuali run1/basecall/flappie/basecalls.fasta
    katuali run1/basecall/guppy/basecalls.fasta


Assembly
--------

Reads can be assembled in two ways at present:

.. code-block:: bash

    # assemble with canu, specifying the genome size in the target name. 
    katuali run1/basecall/scrappie/canu_gsz_4.0M/consensus.fasta  

    # use pomoxis mini_assemble to assemble with miniasm, then form consensus
    # with racon
    katuali run1/basecall/scrappie/miniasm_racon/consensus.fasta  



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

.. _starting_from_basecalls:

Starting from existing basecalls
--------------------------------

If you have already basecalled your data, mocking out the working space as if
katuali had basecalled allows any derived targets to be created. 

.. code-block:: bash
    
    START=${PWD}
    BCDIR=run1/basecall/mybasecalls/
    mkdir -p ${BCDIR}
    cd ${BCDIR}
    source ~/git/pomoxis/venv/bin/activate
    seqkit fq2fa /path/to/basecalls.fastq > basecalls.fasta
    # these next two steps are only required if you wish to use nanopolish.  
    ln -s /path/to/sequencing_summary.txt sequencing_summary.txt
    ln -s /path/to/fast5 reads
    cd ${START}
    # now we can run katuali to assemble and polish
    katuali ${BCDIR}/miniasm_racon/consensus.fasta

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

It is possible to train medaka models starting from
folders of fast5s in a single command once the config has been modified to
reflect your input data (fast5s and genomes for each run as well as training
and evaluation region definitions).

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

    katuali all_medaka_train_features --keep-going

will:

* basecall all the runs
* align each run to its reference
* create subsampled sets of basecalls over the desired regions and depths
* assemble those sets of basecalls
* create medaka training features for all those sets


Running:

    katuali medaka_train_replicates --keep-going

will do all the tasks of `all_medaka_train_features` and additionally launch multiple medaka model-training replicates.

If some of your input runs have insufficient coverage-depth for some of the
training regions, some of the training feature files will not be made. In this
case the config flag USE_ONLY_EXISTING_MEDAKA_FEAT can be set to true to allow katuali to train using only those features which exist already:

    USE_ONLY_EXISTING_MEDAKA_FEAT: true 

Refer to comments in the config (katuali/config.yaml) to see how this process can be controlled. 
