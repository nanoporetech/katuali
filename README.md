
Katuali
=======

`Katuali` is a flexible consensus pipeline implemented in [Snakemake](https://snakemake.readthedocs.io) to basecall, assemble and polish 
Oxford Nanopore Technologies' sequencing data.

Features
--------

  * Run a pipeline processing fast5s to a consensus in a single command.
  * Recommended fixed "standard" and "fast" pipelines.
  * Interchange basecaller, assembler and consensus components of the pipelines simply by changing the target filepath. 
  * Seemless distribution of tasks over local or distributed compute.
  * Highly configurable.  
  * Open source (Mozilla Public License 2.0).


Documentation can be found at https://nanoporetech.github.io/katuali/.

Installation
------------

`Katuali` is a [Snakemake](https://snakemake.readthedocs.io) pipeline comprising a [Snakefile](https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html#) and [config](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html).

As such, all that is required to run the pipeline is `Snakemake`. 

A Makefile is provided to create a [virtual environment](https://docs.python.org/3/tutorial/venv.html) into which to install snakemake as well the katuali convenience wrapper. 

Katuali has been tested on Linux, specifically Ubuntu 16.

To setup the environment run:

    git clone https://github.com/nanoporetech/katuali.git
    cd katuali
    make install
    source ./venv/bin/activate


Dependencies
------------

Katuali makes use of a number of tools to do basecalling, assembly and
polishing that you will need to install.  

You only need to install what you will use.

Tools are inputs to the targets; if something can't be found, Snakemake will tell you it is missing.

The following config default parameters can be changed in the config or on the command line to point to your installations of these tools: 

* [SCRAPPIE](https://github.com/nanoporetech/scrappie): "~/git/scrappie"
* [FLAPPIE](https://github.com/nanoporetech/flappie): "~/git/github/flappie"
* [IN_POMOXIS](https://github.com/nanoporetech/pomoxis): "~/git/pomoxis/venv/bin/activate"
* [CANU_EXEC](https://github.com/marbl/canu): "~/git/canu-1.7.1/Linux-amd64/bin/canu"
* [NANOPOLISH](https://github.com/jts/nanopolish): "~/git/nanopolish"
* [IN_MEDAKA](https://github.com/nanoporetech/medaka): "~/git/medaka/venv/bin/activate"
* GUPPY: "/usr/bin/guppy_basecaller"
* IN_MIYAGI: "~/git/miyagi/venv/bin/activate"
* IN_RAY: "~/git/ray/venv/bin/activate"

Please refer to the documentation of each of these tools for installation instructions.


Usage
-----

The `Katuali` tests contain examples of how to basecall,
assemble and polish a small dataset that comes bundled with `Katuali`. 

To run with your own data, start by creating a dictory of reads (which could contain subdirectories of reads):

    ln -s /path/to/fast5 reads
    
Then you can then make any number of targets the pipline knows how to make by running e.g.:

    katuali fast_assm_polish

This will basecall the reads, then assemble them with miniasm, and polish with racon and medaka. 

Running

    katuali standard_assm_polish

will instead basecall, assemble with canu and the polish with nanopolish. 
