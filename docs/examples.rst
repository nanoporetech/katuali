

.. _introduction:

Introduction
============

One of the main aims of `Katuali` is to enable the bolting-together of any number of basecallers, assemblers and polishers into flexible (i.e. not predefined) pipelines.

Basecallers, assemblers and polishers can be interchanged because their Snakefile rules have standardised inputs and outputs.

`Katuali` makes strong use of `Snakemake`'s powerful `wildcards <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#wildcards>`_
to extract pipeline parameters (e.g. which basecaller, assembler as well as their parameters) from the target file path. 

For example, while

.. code-block:: bash

    katuali basecall/scrappie/basecalls.fasta

will basecall with scrappie (starting from `./reads`), 

.. code-block:: bash

    katuali basecall/guppy/basecalls.fasta

will basecall with guppy. 

Each step of a multi-step pipelines stores its outputs in a subdirectory of the previous stage. 

In this way, we can define a multistep-pipeline to basecall with scrappie,
assemble with canu, then polish with medaka followed by nanopolish with the
following target: 

.. code-block:: bash

    katuali basecall/scrappie/canu/medaka/nanopolish/consensus.fasta

This nested working directory stores the data in such a way that is it obvious what went into and out of each stage of the pipeline.

It also enables forks in the pipelines which might share e.g. basecalls, but differ in assembly or consensus methods. 

For example, if we want to basecall with scappie, assemble with canu, then run medaka on the canu assembly and seperately run nanopolish on the canu assembly, we could use the targets: 

.. code-block:: bash

    katuali basecall/scrappie/canu/nanopolish/consensus.fasta  basecall/scrappie/canu/medaka/consensus.fasta

`Snakemake` will create DAG of tasks to basecall then assemble, then fork the pipeline to run nanopolish and medaka seperately (and in parallel, given enough resource). 


Basecalling
===========

Supported basecallers include scrappie, flappie and guppy:

.. code-block:: bash

    katuali basecall/scrappie/basecalls.fasta
    katuali basecall/flappie/basecalls.fasta
    katuali basecall/guppy/basecalls.fasta

Assembling
==========

Reads can be assembled in three ways at present:

.. code-block:: bash

    # assemble with canu
    katuali basecall/scrappie/canu/consensus.fasta  

    # use pomoxis mini_assemble to assemble with miniasm, then form consensus with racon
    katuali basecall/scrappie/miniasm_racon/consensus.fasta  

    # use pomoxis mini_assemble to align reads to a reference then form consensus with racon
    katuali basecall/scrappie/ref_guided_racon/consensus.fasta  


Polishing
=========

Nanopolish, medaka and miyagi can be used to polish consensuses:

.. code-block:: bash

    katuali basecall/scrappie/canu/nanopolish/consensus.fasta 
    katuali basecall/scrappie/canu/medaka/consensus.fasta 
    katuali basecall/scrappie/canu/miyagi/consensus.fasta 


Pipeline restrictions
=====================

`Katuali` aims to be as flexible as possible, but there are some obvious restrictions:
    * you must basecall before you assemble.
    * you must assemble before you polish.


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


Creating subsampled datasets
============================

Katuali also supports the generation of datasets with even coverage at a given depth.

.. code-block:: bash

    katuali basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta


will perform the following steps:

    * basecall the reads to create target:

        basecall/scrappie/basecalls.fasta

    * align the basecalls to the reference to create target:

        basecall/scrappie/align/calls2ref.bam

    * subsample all contigs in the bam to 25X to create target (in one step):

        basecall/scrappie/align/all_contigs/25X/basecalls.fasta

    * perform a ref-guided assembly and racon consensus to create target: 

        basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta

.. note:: The rule to create subsampled datasets differs from other rules in
    that it creates two levels of nested directories in a single step (in this case
    `all_contigs/25X`). 
    The extraction of specific regions/contigs without
    subsampling to a specific depth is not currently supported.  


Subsampling only one reference contig
=====================================

It is also possible to subsample just one of the contigs in your reference by specifying targets such as:

.. code-block:: bash

    katuali basecall/scrappie/align/ecoli_SCS110_plasmid2/25X/ref_guided_racon/consensus.fasta 

which will just process the ecoli_SCS110_plasmid2.


Subsampling only specified regions
==================================

It is also possible to subsample only specified regions specifed as samtools strings:

.. code-block:: bash

    REGIONS="ecoli_SCS110_chromosome:50000-150000 ecoli_SCS110_chromosome:200000-250000"
    katuali basecall/scrappie/align/my_regions/25X/ref_guided_racon/consensus.fasta --config REGIONS="$REGIONS"


Pipeline configuration
======================

`Snakemake` allows pipeline parameters to be provided in a `config file, or on the command line <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`_ .

If you use the `katuali` wrapper script (rather than running `Snakemake` directly), by default your pipeline will use the yaml config provided with `katuali`. 

The default config file can be overridden using the `--configfile` option, and individual config parameters can be overwriddden with the `--config` option:

.. code-block:: bash

    # use a custom config
    katuali basecall/scrappie/miniasm_racon/consensus.fasta --configfile myconfig.yaml

.. code-block:: bash

    # override MINI_ASSEMBLE_OPTS config on the command line
    katuali basecall/scrappie/miniasm_racon/consensus.fasta --config MINI_ASSEMBLE_OPTS="-c"


Nested configs
==============

Nested configs allow you to access specific settings using a target suffix.
The nested config entry below defines different mini_assemble options to be used with different suffixes:
    
.. code-block:: yaml

    MINI_ASSEMBLE_OPTS: 
        "": ""  # use the mini_assemble defaults 
        "_c": "-c"  # run basecalls through pore-chop before assembly
        "_ce": "-c -e 10"  # run basecalls through pore-chop and error correct longest 10% of reads prior to assembly

The following katuali targets with then run with either the defails, or `_ce` options:

.. code-block:: bash

    # use default MINI_ASSEMBLE_OPTS (suffix is empty string "")
    katuali basecall/scrappie/miniasm_racon/consensus.fasta
    # use MINI_ASSEMBLE_OPTS specified by suffix "_ce"
    katuali basecall/scrappie/miniasm_racon_ce/consensus.fasta

You can add a suffix to most targets to specify options. If the suffix does not exist in the nested config, the default parameters will be used as if the suffix were empty.
This can be useful if you want to run the same target twice (to e.g. sample any random error or e.g. different versions of a code) and want the output files to be distinct. 

Furthermore, settings in the config file can be overridden on the command line:

.. code-block:: bash

    katuali fast_assm_polish --config MINI_ASSEMBLE_OPTS="-c -e 5"

However, this only works if you use the katuali wrapper, not if you run Snakemake directly (you can't change a nested config entry on the command line using `Snakemake`).

The wrapper achieves this merging command line `--config` options with the input `--configfile` and saving the merged YAML config before running snakemake with the merged config. 


Automatic saving of logs and configs
====================================

If you use the convenience wrapper `katuali` rather than calling snakemake
directly, the `katuali` wrapper will write a copy of all logs to the directory
`./logs` and all katuali configs to `./configs`. 


Processing and resource
=======================

The pipeline can be used on the local machine, or submitted to a queuing system. 

There are two parameters which control cpu usage:

* the --jobs `N` (or -j for short) option to Snakemake control the total number of threads which can be used at a time, and hence the number of tasks which can be run simultaneously. 

* the --config THREADS_PER_JOB=`n` config parameter, which determines the number of threads that a single multi-threaded job can use.

Note that --jobs will control the total number of threads used; i.e. if 
THREADS_PER_JOB is set to 4 and --jobs is set to 8, up to two multi-threaded jobs can run at a time.

When submitting to a queuing system, the --jobs will limit the number of queue slots used simultaneously.

The `katuali` wrapper has an `--sge` option which can handle submission to a default sge queue:
    
.. code-block:: bash

    NSLOTS=100
    target=fast_assm_polish
    katuali -j ${NSLOTS} --sge ${target}

which is equivalent to running: 

.. code-block:: bash

    NSLOTS=100
    target=fast_assm_polish
    katuali -j ${NSLOTS} --cluster-sync "${qsub_path} -V -cwd -P research -l {params.sge} -sync yes" --latency-wait 300 ${target}

The local snakemake task will then submit all tasks to the queue for execution. The --latency-wait parameter is useful for ensuring that pipelines don't crash due to output files not appearing on the node where snakemake is run due to latencies on networked file systems. 

