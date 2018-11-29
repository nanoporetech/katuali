
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

    katuali basecall/scrappie/basecalls.fasta

will perform basecalling with scrappie (starting from `./reads`), 

.. code-block:: bash

    katuali basecall/guppy/basecalls.fasta

will basecall with guppy. 

Each step of a multi-step pipeline stores its outputs in a subdirectory of the
previous stage. In this way, we can define a multistep-pipeline to basecall
with scrappie, assemble with canu, then polish with medaka, and finally run
nanopolish by contructing the target: 

.. code-block:: bash

    katuali basecall/scrappie/canu/medaka/nanopolish/consensus.fasta

This nested working directory stores the data in such a way that is it obvious
what went into and out of each stage of the pipeline.

It also enables forks in the pipelines which might share basecalls (or other
intermediate results), but differ in assembly or consensus methods.

For example, if we wish to basecall with scappie, assemble with canu, run
medaka on the canu assembly, and seperately run nanopolish on the canu assembly,
we could use the targets: 

.. code-block:: bash

    katuali basecall/scrappie/canu/nanopolish/consensus.fasta basecall/scrappie/canu/medaka/consensus.fasta

`Snakemake` will create a graph of tasks to perform the common basecall
and assembly tasks, then run separately nanopolish and medaka from the same
inputs (and in parallel, given enough resource).


Basecalling
-----------

Supported basecallers include scrappie, flappie and guppy:

.. code-block:: bash

    katuali basecall/scrappie/basecalls.fasta
    katuali basecall/flappie/basecalls.fasta
    katuali basecall/guppy/basecalls.fasta


Assembly
--------

Reads can be assembled in three ways at present:

.. code-block:: bash

    # assemble with canu katuali basecall/scrappie/canu/consensus.fasta  

    # use pomoxis mini_assemble to assemble with miniasm, then form consensus
    # with racon
    katuali basecall/scrappie/miniasm_racon/consensus.fasta  

    # use pomoxis mini_assemble to align reads to a reference then form
    # consensus with racon 
    ln -s /path/to/ref.fasta ref.fasta
    katuali basecall/scrappie/ref_guided_racon/consensus.fasta  


Polishing
---------

Sequences can be polished with either Nanopolish or medaka to create higher
accuracy consensus sequences:

.. code-block:: bash

    katuali basecall/scrappie/canu/nanopolish/consensus.fasta
    katuali basecall/scrappie/canu/medaka/consensus.fasta


Pipeline restrictions
---------------------

`Katuali` aims to be as flexible as possible, but there are some obvious
restrictions:

    * basecalling must be performed before assembly.
    * assembly must come before polishing (use of polishing targets to
      error correct reads is not supports).


Starting from existing basecalls
================================

If you have already basecalled your data, mocking out the working space as if katuali had basecalled allows any derived targets to be created. 

.. code-block:: bash
    
    START=${PWD}
    BCDIR=basecall/mybasecalls/
    mkdir -p ${BCDIR}
    cd ${BCDIR}
    source ~/git/pomoxis/venv/bin/activate
    seqkit fq2fa /path/to/basecalls.fastq > basecalls.fasta
    # these next two steps are only required if you wish to use signal-level polishers such as nanopolish. 
    ln -s /path/to/sequencing_summary.txt sequencing_summary.txt
    ln -s /path/to/fast5 reads
    cd ${START}
    # now we can run katuali to assemble and polish
    katuali ${BCDIR}/miniasm_racon/consensus.fasta

Read coverage depth and subsampling
===================================

Calculating read coverage depth
-------------------------------

It is often useful to know the read coverage depth of a dataset. 
This can be calculated as follows:

.. code-block:: bash

    ln -s /path/to/ref.fasta ref.fasta
    katuali basecall/scrappie/align/depth

The depth directory will contain a text file per reference contig with coverage
vs genomic coordinate, as well as a file containing summary statistics for all
contigs.


Creating subsampled datasets
----------------------------

Katuali also supports the generation of datasets with even coverage over a
reference at a given depth.

.. code-block:: bash

    ln -s /path/to/ref.fasta ref.fasta
    katuali basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta

will perform the following steps:

    * basecall the reads to create:
      `basecall/scrappie/basecalls.fasta`
    * align the basecalls to the reference to create:
      `basecall/scrappie/align/calls2ref.bam`
    * subsample all contigs in the .bam file to 25X to create (in one step):
      `basecall/scrappie/align/all_contigs/25X/basecalls.fasta`
    * perform a ref-guided assembly and racon consensus to create:
      `basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta`


.. note:: The rule to create subsampled datasets differs from other rules in
    that it creates two levels of nested directories in a single step (in this case
    `all_contigs/25X`). The extraction of specific regions/contigs without
    subsampling to a specific depth is not currently supported.


Subsampling a single reference contig
-------------------------------------

It is also possible to subsample just one of the contigs in your reference by
specifying targets such as:

.. code-block:: bash

    katuali basecall/scrappie/align/ecoli_SCS110_plasmid2/25X/ref_guided_racon/consensus.fasta 

which will just process the reference sequence `ecoli_SCS110_plasmid2`.


Subsampling specified regions
-----------------------------

It is also possible to subsample only specified regions specifed as samtools
strings:

.. code-block:: bash

    REGIONS="ecoli_SCS110_chromosome:50000-150000 ecoli_SCS110_chromosome:200000-250000"
    katuali basecall/scrappie/align/my_regions/25X/ref_guided_racon/consensus.fasta --config REGIONS="$REGIONS"

