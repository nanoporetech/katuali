![Oxford Nanopore Technologies logo](images/ONT_logo_590x106.png)


Katuali
=======

`Katuali` is a flexible consensus pipeline implemented in
[Snakemake](https://snakemake.readthedocs.io) to basecall, assemble, and polish
Oxford Nanopore Technologies' sequencing data.

Features
--------

  * Run a pipeline processing fast5s from multiple runs into multiple consensuses in a single command.
  * Recommended fixed "standard" and "fast" pipelines.
  * Interchange basecaller, assembler, and consensus components of the
    pipelines simply by changing the target filepath.
  * Medaka training pipeline including generation of training data, model training and model evaluation. 
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
* [IN_MEDAKA](https://github.com/nanoporetech/medaka): "~/git/medaka/venv/bin/activate"
* [NANOPOLISH](https://github.com/jts/nanopolish): "~/git/nanopolish"
* GUPPY: "/usr/bin/guppy_basecaller"

Please refer to the documentation of each of these tools for installation
instructions.


Usage
-----

The `Katuali` tests contain examples of how to basecall, assemble, and polish
a small dataset that comes bundled with `Katuali`.

To run with other data, start by creating a directory of reads (which could
contain subdirectories of reads) within a run directory:

    mkdir -p run1 && cd run1 && ln -s /path/to/fast5 reads && cd ../
    
Then make a copy of the katuali config into your working directory;

    cp ~/git/katuali/config.yaml .

and update the katuali config to reflect your data:
    
    DATA:
        'run1':
            'GENOME_SIZE': '4.0M'  # for canu we need to specify genome size

Then calculate any of the outputs the pipeline knows how to make by running e.g.:

    katuali fast_assm_polish

This will basecall the reads, assemble them with miniasm, and polish the
assembly with racon and medaka. 

Running

    katuali standard_assm_polish

will instead basecall, assemble with canu, then polish with racon and medaka. 

Alternatively,  

    katuali standard_assm_nanopolish

will instead basecall, assemble with canu, then polish with racon and nanopolish. 


Medaka training pipeline
------------------------

It is possible to train medaka models starting from
folders of fast5s in a single command once the config has been modified to
reflect your input data (fast5s and genomes for each run as well as training
and evaluation region definitions).

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
