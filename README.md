![Oxford Nanopore Technologies logo](images/ONT_logo_590x106.png)


Katuali
=======

`Katuali` is a flexible consensus pipeline implemented in
[Snakemake](https://snakemake.readthedocs.io) to basecall, assemble, and polish
Oxford Nanopore Technologies' sequencing data.

Features
--------

  * fast5 to high quality assembly in a single command.
  * Recommended fixed "standard" and "fast" pipelines.
  * Interchange basecaller, assembler, and consensus components of the
    pipelines simply by changing the target filepath.
  * [Medaka](https://github.com/nanoporetech/medaka) training pipeline
    including generation of training data, model training, and model evaluation.
  * Seemless distribution of tasks over local or distributed compute.
  * Open source (Mozilla Public License 2.0).

Documentation can be found at https://nanoporetech.github.io/katuali/.

Installation
------------

A Makefile is provided to create a
[virtual environment](https://docs.python.org/3/tutorial/venv.html) into which
all necessary dependencies will be installed. Katuali has been tested on Linux
(specifically Ubuntu 16) and MacOS.

To setup the environment run:

    git clone https://github.com/nanoporetech/katuali.git
    cd katuali
    make install
    source ./venv/bin/activate

A conda package is forthcoming.

Dependencies
------------

Katuali makes use of a number of tools to do basecalling, assembly, and
polishing that need to be installed to perform such steps of a pipeline.
Katuali does not enforce that all possible tools are present; only those
required to run a user's pipeline need be installed. If Katuali cannot
find a tool, it will alert the user.

The list below indicates possible options for pipeline steps and the default
filesystem location where Katuali expects to find them. These locations can be
changed by the user, see the
[documentation](https://nanoporetech.github.io/katuali/) for more information.
Please refer to the documentation of each of these tools for installation
instructions.

**Basecalling Tools**

* [scrappie](https://github.com/nanoporetech/scrappie): "~/git/scrappie"
* [flappie](https://github.com/nanoporetech/flappie): "~/git/github/flappie"
* [guppy](https://community.nanoporetech.com/downloads): "/usr/bin/guppy_basecaller"

Guppy is recommended, both scrappie and flappie are research prototypes.

**Assembly Tools**

* [canu](https://github.com/marbl/canu): "~/git/canu-1.7.1/Linux-amd64/bin/canu"
* [pomoxis](https://github.com/nanoporetech/pomoxis): "~/git/pomoxis/venv/bin/activate"

In the context of Katuali, pomoxis functions as a wrapper to
[miniasm](https://github.com/lh3/miniasm). Canu is utilised in the standard
workflow, while pomoxis is used in the fast workflow.

**Consensus Tools**

* [medaka](https://github.com/nanoporetech/medaka): "~/git/medaka/venv/bin/activate"
* [nanopolish](https://github.com/jts/nanopolish): "~/git/nanopolish"

Medaka is recommended in combination with the latest "flip-flop" algorithm
in the guppy basecaller, and for rapid consensus. Nanopolish may still be
preferred if using older basecallers.

Usage
-----

The `Katuali` tests contain examples of how to basecall, assemble, and polish
a small dataset that comes bundled with `Katuali`. These tests can be run with:

    make test

from the `Katuali` source directory. Running the tests requires installation of
`scrappie`, `pomoxis`, `medaka`, and `nanopolish`. See the
[documentation](https://nanoporetech.github.io/katuali/installation.html#tests)
for more information on the outputs of running the above.

To run with other data, start by creating a directory of reads (which could
contain subdirectories of reads) within a run directory:

    mkdir -p run1
    cd run1
    ln -s /path/to/fast5 reads  # create a softlink to the fast5 data
    cd ../
    
Then make a copy of the katuali config into your working directory,

    katuali_config my_config.yaml

and update the katuali config to reflect your data:
    
    DATA:
        'run1':
            'GENOME_SIZE': '4.0M'  # for canu we need to specify genome size

To perform basecalling, a quick assembly with miniasm, and consensus with racon
and medaka it is sufficient to run:

    katuali --configfile my_config.yaml fast_assm_polish

Alternatively to assemble with canu run:

    katuali --configfile my_config.yaml standard_assm_polish

To polish the assembly with nanopolish (not recommended in tandem with
flip-flop basecalling, where medaka is preferred):

    katuali --configfile my_config.yaml standard_assm_nanopolish


Medaka training pipeline
------------------------

It is possible to train medaka models starting from folders of fast5s in a
single command once the config has been modified to reflect your input data
(fast5s and genomes for each run as well as training and evaluation region
definitions).

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

will do all the tasks of `all_medaka_train_features` and additionally launch
multiple medaka model-training replicates.

If some of your input runs have insufficient coverage-depth for some of the
training regions, some of the training feature files will not be made. In this
case the config flag `USE_ONLY_EXISTING_MEDAKA_FEAT` can be set to true to allow
katuali to train using only those features which exist already:

    USE_ONLY_EXISTING_MEDAKA_FEAT: true 

Refer to comments in the config (katuali/config.yaml) to see how this process
can be controlled. 
