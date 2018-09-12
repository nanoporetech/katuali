
Katuali
=======

`Katuali` is a consensus pipeline implemented in SnakeMake to combine consensus tools 
focussed on reducing the errors in assemblies using Oxford Nanopore Technologies' data. 

Installation
============

The pipeline assumes the existance of various other installations including:

    * scrappie executable (default ~/git/scrappie/build/scrappie)
    * nanopolish source directory (default ~/git/nanopolish)
    * pomoxis venv (default ~/git/pomoxis/venv/bin/activate) 
    * medaka venv (default ~/git/medaka/venv/bin/activate) 
    * miyagi venv (default ~/git/miyagi/venv/bin/activate)

the paths to which can be specified in the config or on the command line. 

To install create a venv for running Snakemake run:
    
    make install

Quickstart
==========

The easiest way to test the pipeline is to run the tests, which will basecall,
assemble and polish a small dataset that comes bundled with `Katuali`.

    make test

If you want to run with your own data, start by creating a dictory of reads and reference

    ln -s /path/to/fast5/directories(s) reads
    ln -s /path/to/reference/ ref.fasta 
    
Then you can then make any number of targets the pipline knows how to make by running e.g.:

    snakemake -s ~/git/katuali/Snakefile analysis/basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta --printshellcmds

which will perform the following steps:

    * basecall the reads to create target:

        analysis/basecall/scrappie/basecalls.fasta

    * align the basecalls to the reference to create target:

        analysis/basecall/scrappie/align/calls2ref.bam

    * subsample all contigs in the bam to 25X to create target:

        analysis/basecall/scrappie/align/all_contigs/25X/basecalls.fasta

    * perform a ref-guided assembly and racon consensus to create target: 

        analysis/basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta

Pipeline Flexibility
====================

The pipeline is setup to allow for multiple rounds of polishing with different consensus methods, so making eg.

    analysis/basecall/scrappie/align/all_contigs/25X/ref_guided_racon/medaka/nanopolish/miyagi/consensus.fasta

will:
    * form a ref-guided racon consensus
    * polish the racon consensus with medaka
    * polish the medaka consensus with nanopolish
    * polish the nanopolish consensus with miyagi

The order of polishing tools does not matter here, except that for the moment polishing must start from a ref-guided racon consensus. 

Also, if you wish for example to test the impact of various scrappie basecalling settings
or models on racon consensus accuracy, you could launch several snakemake jobs, specifying the scrappy
options on the command line, and adding a suffix to each scrappie directory:
    
    set1="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 0.65 --temperature2 1.7"
    set2="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 1.6 --temperature2 1.7"
    target1=analysis/basecall/scrappie_t1_065_t2_170/align/all_contigs/25X/ref_guided_racon/consensus.fasta
    target2=analysis/basecall/scrappie_t1_160_t2_170/align/all_contigs/25X/ref_guided_racon/consensus.fasta

    snakemake -s ~/git/katuali/Snakefile $target1 --config SCRAPPIE_OPTS="${set1}" --printshellcmds
    snakemake -s ~/git/katuali/Snakefile $target2 --config SCRAPPIE_OPTS="${set2}" --printshellcmds

You can also add a suffix to the ref_guided_racon targets to allow e.g. a grid scan over different mini_assemble settings.

Processing just one reference contig/chromosome
===============================================

It is also possible to process just one of the contigs in your reference by specifying targets such as:

   analysis/basecall/scrappie/align/ecoli_SCS110_chromosome/25X/ref_guided_racon/consensus.fasta 

which will just process the ecoli_SCS110_chromosome.

Processing and resource
=======================

For the moment the pipeline is designed to be used on the local machine.
 
Basecalling with scrappie, forming a racon consensus, and running medaka tasks use threads to speed things up,
while running nanopolish and miyagi creates a queue of jobs to run over regions/batches. 

Though this is not currently set up, SnakeMake has integrated SGE support, so that these jobs could be sent to a queue.

Pipeline config
===============

Things like scrappie or nanopolish options can be set in the YAML config, or over-run on the command line using e.g. --config NSLOTS=40
You can also use a different config to the default (~/git/katuali/config.yaml) by specifying a --configfile path/to/config.yml

TODO list
=========
* Add guppy basecaller support
* Add suffix support (already setup for scrappie and ref_guided_racon) to other tools
* Setup SGE config so pipelines can be submitted with each target being a SGE job
* Add post-processing tools such as assess_assembly
* Add convenience analysis of e.g. distribution of depth, speed
* Add other assembly options
