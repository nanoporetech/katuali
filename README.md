
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

   snakemake -s ~/git/katuali/Snakefile fast_assm_polish

This will basecall the reads, then assemble them with miniasm, and polish with racon and medaka. 

Running

   snakemake -s ~/git/katuali/Snakefile standard_assm_polish

will instead basecall, assemble with canu and the polish with nanopolish. 

Creating subsampled datasets
============================

Katuali also supports the generation of datasets with even coverage at a given depth.
The command 

    snakemake -s ~/git/katuali/Snakefile basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta

will perform the following steps:

    * basecall the reads to create target:

        basecall/scrappie/basecalls.fasta

    * align the basecalls to the reference to create target:

        basecall/scrappie/align/calls2ref.bam

    * subsample all contigs in the bam to 25X to create target:

        basecall/scrappie/align/all_contigs/25X/basecalls.fasta

    * perform a ref-guided assembly and racon consensus to create target: 

        basecall/scrappie/align/all_contigs/25X/ref_guided_racon/consensus.fasta


Processing just one reference contig/chromosome
===============================================

It is also possible to process just one of the contigs in your reference by specifying targets such as:

   basecall/scrappie/align/ecoli_SCS110_plasmid2/25X/ref_guided_racon/consensus.fasta 

which will just process the ecoli_SCS110_plasmid2.


Pipeline Flexibility
====================

The pipeline is setup to allow for multiple rounds of polishing with different consensus methods, so making eg.

    basecall/scrappie/align/all_contigs/25X/ref_guided_racon/medaka/nanopolish/consensus_to_truth_summ.txt

will:
    * form a ref-guided racon consensus
    * polish the racon consensus with medaka
    * polish the medaka consensus with nanopolish
    * assess consensus accuracy using `pomoxis assess_assembly`

The order of polishing tools does not matter here, except that for the moment
polishing (medaka/nanopolish) must start from an existing consensus
(ref_guided_racon/minias_racon/canu). 

Also, if you wish for example to test the impact of various scrappie basecalling settings
or models on racon consensus accuracy, you could launch several snakemake jobs, specifying the scrappy
options on the command line, and adding a suffix to each scrappie directory:
    
    set1="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 0.65 --temperature2 1.7"
    set2="raw -H mean --model rgrgr_r94 --local 10.0 --uuid --temperature1 1.6 --temperature2 1.7"
    target1=basecall/scrappie_t1_065_t2_170/align/all_contigs/25X/ref_guided_racon/consensus.fasta
    target2=basecall/scrappie_t1_160_t2_170/align/all_contigs/25X/ref_guided_racon/consensus.fasta

    snakemake -s ~/git/katuali/Snakefile $target1 --config SCRAPPIE_OPTS="${set1}" --printshellcmds
    snakemake -s ~/git/katuali/Snakefile $target2 --config SCRAPPIE_OPTS="${set2}" --printshellcmds

You can also add a suffix to the ref_guided_racon targets to allow e.g. a grid scan over different mini_assemble settings.


Providing mini_assemble options on the command line
===================================================

The following will chop off adapters, error-correct the longest 20% of reads, then assemble with miniasm, polish with racon and medaka. 

    snakemake -s ~/git/katuali/Snakefile fast_assm_polish --config MINI_ASSEMBLE_OPTS="-c -e 20"


Processing and resource
=======================

The pipeline can be used on the local machine, or submitted to a queuing system. 

There are two parameters which control cpu usage:

* the --jobs `N` (or -j for short) option to Snakemake control the total number of threads which can be used at a time, and hence the number of tasks which can be run simultaneously. 

* the --config THREADS_PER_JOB=`n` config parameter, which determines the number of threads that a single multi-threaded job can use.

Note that --jobs will control the total number of threads used; i.e. if 
THREADS_PER_JOB is set to 4 and --jobs is set to 8, up to two multi-threaded jobs can run at a time.

When submitting to a queuing system, the --jobs will limit the number of queue slots used simultaneously. 
To submit to an SGE/UGE grid engine using use: 

    source ~/git/katuali/venv/bin/activate
    qsub_path=$(which qsub)
    NSLOTS=100
    target=fast_assm_polish
    snakemake -s ~/git/katuali/Snakefile -j ${NSLOTS} --cluster-sync "${qsub_path} -V -cwd -P research -l {params.sge} -sync yes" --latency-wait 300 ${target}

The local snakemake task will then submit all tasks to the queue for execution. The --latency-wait parameter is useful for ensuring that pipelines don't crash due to output files not appearing on the node where snakemake is run due to latencies on networked file systems. 


Pipeline config
===============

Things like scrappie or nanopolish options can be set in the YAML config, or over-run on the command line using e.g. --config THREADS_PER_JOB=40
You can also use a different config to the default (~/git/katuali/config.yaml) by specifying a --configfile path/to/config.yml

