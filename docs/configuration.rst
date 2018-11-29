

.. _configuration:

Pipeline configuration
======================

`Snakemake` allows pipeline parameters to be provided in a `config file, or on
the command line
<https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`_ .

If you use the `katuali` wrapper script (rather than running `Snakemake`
directly), by default your pipeline will use the yaml config provided with
`katuali`.

The default config file can be overridden using the `--configfile` option, and
individual config parameters can be overwriddden with the `--config` option:

.. code-block:: bash

    # use a custom config katuali
    basecall/scrappie/miniasm_racon/consensus.fasta --configfile myconfig.yaml

    # override MINI_ASSEMBLE_OPTS config on the command line katuali
    basecall/scrappie/miniasm_racon/consensus.fasta --config
    MINI_ASSEMBLE_OPTS="-c"


Nested configuration
--------------------

Nested configs allow access of specific settings using a target suffix.
The nested config entry below defines different mini_assemble options to be
used with different suffixes:
    
.. code-block:: yaml

    MINI_ASSEMBLE_OPTS:
        "": ""             # use the mini_assemble defaults
        "_c": "-c"         # run basecalls through pore-chop before assembly
        "_ce": "-c -e 10"  # run basecalls through pore-chop and error correct longest 10% of reads prior to assembly

The following katuali targets will then run with either the defaults, or `_ce`
options:

.. code-block:: bash

    # use default MINI_ASSEMBLE_OPTS (suffix is empty string "")
    katuali basecall/scrappie/miniasm_racon/consensus.fasta
    # use MINI_ASSEMBLE_OPTS specified by suffix "_ce"
    katuali basecall/scrappie/miniasm_racon_ce/consensus.fasta

A suffix can be added to most targets to specify options. If the suffix does
not exist in the nested config, the default parameters will be used as if the
suffix were empty. This can be useful if you want to run the same target twice
(maybe to sample any random error or to test different versions of a code) and
want the output files to be distinct. 

Further, settings in the config file can be overridden on the command line:

.. code-block:: bash

    katuali fast_assm_polish --config MINI_ASSEMBLE_OPTS="-c -e 5"

However, this only works if you use the katuali wrapper, not if you run
Snakemake directly; a nested config entry cannot be changed on the command line
using `Snakemake`.

The wrapper achieves this merging command line `--config` options with the
input `--configfile` and saving the merged YAML config before running snakemake
with the merged config. 


Automatic saving of logs and configuration
------------------------------------------

If you use the convenience wrapper `katuali` rather than calling snakemake
directly, the `katuali` wrapper will write a copy of all logs to the directory
`./logs` and all katuali configs to `./configs`. 


Processing and resource
-----------------------

The pipeline can be used on the local machine, or submitted to a queuing
system. 

There are two parameters which control cpu usage:

    * the `--jobs N` (or `-j` for short) option to Snakemake controls the total CPU
      resource which can be used at a time, and hence the number of tasks which
      can be run simultaneously. 
    
    * the `--config THREADS_PER_JOB=n` config parameter, which determines the
      number of threads that a single multi-threaded task can use.

Note that `--jobs` will control the total number of threads used; i.e. if
`THREADS_PER_JOB` is set to 4 and `--jobs` is set to 8, up to two multi-threaded
jobs can run at a time.

When submitting to a queuing system, the `--jobs` option will limit the number
of queue slots used simultaneously.

The `katuali` wrapper has an `--sge` option which can handle submission to a
default sge queue:
    
.. code-block:: bash

    NSLOTS=100
    target=fast_assm_polish
    katuali -j ${NSLOTS} --sge ${target}

which is equivalent to running: 

.. code-block:: bash

    NSLOTS=100
    target=fast_assm_polish
    katuali -j ${NSLOTS} --cluster-sync "${qsub_path} -V -cwd -l {params.sge} -sync yes" --latency-wait 300 ${target}

The local snakemake task will then submit all tasks to the queue for execution.
The `--latency-wait` parameter is useful for ensuring that pipelines don't crash
due to output files not appearing on the node where snakemake is run due to
latencies on networked file systems. 

