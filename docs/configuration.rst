
.. _configuration:

Pipeline configuration
======================

`Katuali` uses `Snakemake` which allows pipeline 
`parameters <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`_
to be provided in a file, or on the command line.

If you use the `katuali` wrapper script (rather than running `Snakemake`
directly), by default your pipeline will use the yaml config provided with
`katuali`.

The default config file can be overridden using the ``--configfile`` option, and
individual config parameters can be overwriddden with the ``--config`` option:

.. code-block:: bash

    # use a custom config
    katuali basecall/scrappie/miniasm_racon/consensus.fasta --configfile myconfig.yaml

    # override MINI_ASSEMBLE_OPTS config on the command line
    katuali basecall/scrappie/miniasm_racon/consensus.fasta --config MINI_ASSEMBLE_OPTS="-c"


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

    katuali all_fast_assm_polish --config MINI_ASSEMBLE_OPTS="-c -e 5"

The ``katuali`` program achieves this merging command line ``--config`` options
with the input ``--configfile`` and saving the merged YAML config before running
snakemake. 


Processing and resource
-----------------------

The pipeline can be used on the local machine, or submitted to a cluster.

There are two parameters which control CPU usage:

    * the ``--cores N`` option, which limits the totol number of threads which can be simultaneously used by all Snakemake tasks.
    
    * the ``--config THREADS_PER_JOB=n`` config parameter, determines the maximum
      number of `threads
      <https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-1-specifying-the-number-of-used-threads>`_
      that a single multi-threaded rule will use. When fewer cores than threads
      are provided, the number of threads a task uses will be reduced to the
      number of given cores.

As an example, if ``THREADS_PER_JOB`` is set to 4 and ``--cores`` is set to 8, up to two multi-threaded
tasks can run at a time.


Running on the local machine
----------------------------

When running on a local machine using GPUs (e.g. while basecalling with guppy
or training medaka models), `katuali` can limit the number of concurrent GPU
tasks scheduled so as not to saturate GPU resource by informing katuali how
many GPUs are present on the machine:

.. code-block:: bash

    NCPUS=$(nproc)  # how many cores available on the machine
    NGPUS=$(nvidia-smi --list-gpus | wc -l)  # how many GPUs available on the machine
    katuali --cores ${NCPUS} --resources gpu=${NGPUS} ${targets}

here ``--resources gpu=${NGPUS}`` specifies the maximum number of GPUs which can be used
simultaneously by concurrent tasks.

.. note:: Note that if ``--cores`` is not specified, it defaults to 1, while if
    ``--resources`` it defaults to 0 (unlimited) and that Snakemake manages
    `threads/cores
    <https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-1-specifying-the-number-of-used-threads>`_
    separately from other `resources
    <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#snakefiles-resources>`_. 


.. _using_cluster:

Submitting tasks to a cluster
-----------------------------

When submitting to a queuing system, the ``--cores`` option will limit the number
of queue slots used simultaneously.

The `katuali` wrapper has an ``--autocluster`` option which can handle submission to a
default cluster using DRMAA:
    
.. code-block:: bash

    NSLOTS=100
    target=all_fast_assm_polish
    katuali --cores ${NSLOTS} --autocluster ${target}

The ``--autocluster`` option makes us of the default `katuali` `cluster config
<https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration>`_ 
to submit jobs to an SGE cluster. The use of cluster configs allows
us to abstract away details specific to a given cluster, and easily switch
between clusters simply by changing the cluster config. See the `Snakemake documentation
on cluster configs for futher details
<https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration>`_. 

Using the default `katuali` cluster config in conjuction with the ``--autocluster`` option is equivalent to running:

.. code-block:: bash

    NSLOTS=100
    target=all_fast_assm_polish
    katuali --cores ${NSLOTS} --latency-wait 300 --drmaa "-V -cwd -l gpu={resources.gpu} -pe mt {threads} -o logs -j y"

Here, ``"-V -cwd -l gpu={resources.gpu} -pe mt {threads} -o logs -j y"`` are the
options specific to the SGE scheduler informing it what resources a task
requires.  Note that the resource requirements are expressed in brackets
(``{resources.gpu}`` and ``{threads}``) and will be replaced with actual values
depending on the rule generating the task being submitted.

`katuali` abstracts away these SGE-specific details by using its default cluster config:

.. code-block:: yaml

    __default__:
        n_cpu: "-pe mt "
        n_gpu: "-l gpu="
        export_env: "-V"
        cwd: "-cwd"
        logdir: "-o "
        misc: "-j y"


Using this cluster config, the `katuali` ``--autocluster`` option can support
any DRMAA-enabled cluster using an appropriate cluster-config as the command
line call to Snakemake is expressed in terms of cluster config entries. 
The ``--autocluster`` option implements:

.. code-block:: bash

    NSLOTS=100
    target=all_fast_assm_polish
    cluster_config=$(katuali_datafile cluster_config.yaml)
    katuali --cores ${NSLOTS} --latency-wait 300 --drmaa " {cluster.export_env} {cluster.cwd} {cluster.n_gpu}{resources.gpu} {cluster.n_cpu}{threads} {cluster.logdir}logs {cluster.misc}" --cluster-config ${cluster_config} ${target}

Here all ``{cluster.<variable_name>}`` templates are replaced by values from the cluster config. 

Hence running on another DRMAA cluster should be as simple as creating a new
cluster config with terms equivalent to those in the default katuali
cluster-config, then running:

.. code-block:: bash

    NSLOTS=100
    target=all_fast_assm_polish
    katuali --cores ${NSLOTS} --latency-wait 300 --autocluster --cluster-config my_cluster_config.yaml ${target}


When running on a cluster, the local snakemake task will submit all tasks to
the queue for execution.  The ``--latency-wait`` parameter is useful for ensuring
that pipelines don't crash due to output files not appearing on the node where
snakemake is run due to latencies on networked file systems.
