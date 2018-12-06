![Oxford Nanopore Technologies logo](images/ONT_logo_590x106.png)


Katuali
=======

`Katuali` is a flexible consensus pipeline implemented in
[Snakemake](https://snakemake.readthedocs.io) to basecall, assemble, and polish
Oxford Nanopore Technologies' sequencing data.

Features
--------

  * Run a pipeline processing fast5s to a consensus in a single command.
  * Recommended fixed "standard" and "fast" pipelines.
  * Interchange basecaller, assembler, and consensus components of the
    pipelines simply by changing the target filepath. 
  * Seemless distribution of tasks over local or distributed compute.
  * Highly configurable.  
  * Open source (Mozilla Public License 2.0).


Documentation can be found at https://nanoporetech.github.io/katuali/.

Installation
------------

`Katuali` is a [Snakemake](https://snakemake.readthedocs.io) pipeline comprising a
[Snakefile](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#)
and a [configuration](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).
As such, all that is required to run the pipeline is `Snakemake`. 

A Makefile is provided to create a
[virtual environment](https://docs.python.org/3/tutorial/venv.html) into which
to install snakemake as well the katuali convenience wrapper. Katuali has been
tested on Linux, specifically Ubuntu 16.

To setup the environment run:

    git clone https://github.com/nanoporetech/katuali.git
    cd katuali
    make install
    source ./venv/bin/activate


Dependencies
------------

Katuali makes use of a number of tools to do basecalling, assembly, and
polishing that you will need to install. You need only install those tools
required for the analyses you intend to run. These tools are specified as
requirements to analysis outputs (pipeline targets); if something cannot be
found, Snakemake will tell you that it is missing.

The following default parameters can be changed in the configuration file or
on the command line to point to your installations of these tools: 

* [SCRAPPIE](https://github.com/nanoporetech/scrappie): "~/git/scrappie"
* [FLAPPIE](https://github.com/nanoporetech/flappie): "~/git/github/flappie"
* [IN_POMOXIS](https://github.com/nanoporetech/pomoxis): "~/git/pomoxis/venv/bin/activate"
* [CANU_EXEC](https://github.com/marbl/canu): "~/git/canu-1.7.1/Linux-amd64/bin/canu"
* [NANOPOLISH](https://github.com/jts/nanopolish): "~/git/nanopolish"
* [IN_MEDAKA](https://github.com/nanoporetech/medaka): "~/git/medaka/venv/bin/activate"
* GUPPY: "/usr/bin/guppy_basecaller"
* IN_RAY: "~/git/ray/venv/bin/activate"

Please refer to the documentation of each of these tools for installation
instructions.


Usage
-----

The `Katuali` tests contain examples of how to basecall, assemble, and polish
a small dataset that comes bundled with `Katuali`.

To run with other data, start by creating a directory of reads (which could
contain subdirectories of reads):

    ln -s /path/to/fast5 reads
    
Then calculate any of the outputs the pipeline knows how to make by running e.g.:

    katuali fast_assm_polish

This will basecall the reads, assemble them with miniasm, and polish the
assembly with racon and medaka. Running

    katuali standard_assm_polish

will instead basecall, assemble with canu, and the polish with nanopolish.

Medaka training branch
----------------------

Using this branch, it is now possible to train medaka models starting from folders of fast5s in a single command. 

The config is currently configured for 3 runs of yeast and 3 runs of lbrevis each with median 250X coverage.

The fast5 files should be located in ./runid/reads/. Here `reads` can be symbolic links to fast5 directories. 

Running:

    katuali medaka_train_replicates

will:
    * basecall all the runs
    * align each run to its reference
    * create subsampled sets of basecalls over the desired regions and depths
    * assemble those sets of basecalls
    * create medaka training features for all those sets
    * launch multiple medaka training replicates using all those features

The following config items can be set to control this process:

    # reference to be used for each region set
    REFERENCES: 
        "": "ref.fasta"
        "ecoli": "/nfs/groups_ech/res_data/active/refs/git_references/ecoli/ecoli_SCS110.fasta"
        "yeast": "/nfs/groups_ech/res_data/active/refs/git_references/yeast/yeast_S288C.fasta"
        "lbrevis": "/nfs/groups_ech/res_data/stored/datasets/apps_atcc_genomes/analysis/demultiplexed/refs_for_assembly/NC_008497.1.fa"

    # regions used for subsampling
    REGION_DEFINITIONS:
        "ecoli": "ecoli_SCS110_chromosome"
        "yeast": "yeast_S288C_chromosomeI yeast_S288C_chromosomeII yeast_S288C_chromosomeIII yeast_S288C_chromosomeIV yeast_S288C_chromosomeV yeast_S288C_chromosomeVI yeast_S288C_chromosomeVII yeast_S288C_chromosomeVIII yeast_S288C_chromosomeIX yeast_S288C_chromosomeX yeast_S288C_chromosomeXI yeast_S288C_chromosomeXII yeast_S288C_chromosomeXIII yeast_S288C_chromosomeXIV yeast_S288C_chromosomeXV yeast_S288C_chromosomeXVI"
        "lbrevis": "NC_008497.1"

    # train on these region names (these must be keys for REGION_DEFINITIONS, REFERENCES and RUNIDS)
    MEDAKA_TRAIN_REGIONS: ["yeast", "lbrevis"]

    # Run a training replicate for each of these suffixes
    MEDAKA_TRAIN_REPLICATES:
        ["rep_1", "rep_2", "rep_3"]

    # Define runs - fast5s should be under ./runid/reads/ and can be in subdirectories.
    RUNIDS:
        "yeast": ["80a94504", "f81e9478", "fad5950a"]
        "lbrevis": ["95730659", "15e93eaf", "c50fb4ef"]

    # create sets of basecalls, assemblies and training features for these depths. 
    # medaka models do not perform well (and can make consensus accuracy worse than
    # the draft) for depths they were not exposed to during training. 
    DEPTHS:
        [25, 50, 70, 100, 125, 150, 175, 200, 225, 250]
