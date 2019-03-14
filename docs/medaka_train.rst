
.. _medaka_train:

Medaka training pipeline
========================

It is possible to train and evaluate medaka models starting from
folders of fast5s in a single command once the config has been modified to
reflect your input data (fast5s and genomes for each run as well as training
and evaluation region definitions).


Input data specification
------------------------

Read .fast5 files can be under top-level folders named by a `RUNID` (these need
not be actual run UUIDs). Fast5s should be in `RUNID`/reads.  The keys within
this data dictionary are `RUNIDs`.

`MEDAKA_TRAIN_REGIONS` and `MEDAKA_EVAL_REGIONS` define regions for training
and evaluation.

In the example below we train from the `minion` run using
`ecoli` and `yeast` contigs in the reference and evaluate on the `gridion` run
using the contigs `ecoli`, `yeast` and `na12878_chr21` in the reference.

.. code-block:: yaml

    DATA:
        'MinIonRun1': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': ['ecoli', 'yeast']
            'MEDAKA_EVAL_REGIONS': []
        'MinIonRun2': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': ['ecoli', 'yeast']
            'MEDAKA_EVAL_REGIONS': []
        'GridIonRun1': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': []
            'MEDAKA_EVAL_REGIONS': ['ecoli', 'yeast', 'na12878_chr21']
        'GridIonRun2': 
            'REFERENCE': '/path/to/references.fasta'   
            'MEDAKA_TRAIN_REGIONS': []
            'MEDAKA_EVAL_REGIONS': ['ecoli', 'yeast', 'na12878_chr21']



Coverage depths specification
-----------------------------

Read depths at which to create assemblies for training, this should
span the range of depths at which the model is to be used
Use the --keep-going option of Snakemake if you are happy the relax
the constraint of requiring all depths for all reference sequences. 

.. code-block:: yaml

    DEPTHS:
        [25, 50, 75, 100, 125, 150, 175, 200]


Creating training features
--------------------------

Running:

.. code-block:: bash

    katuali all_medaka_train_features --keep-going  --restart-times 3

will:

* basecall all the runs
* align each run to its reference
* create subsampled sets of basecalls over the desired regions and depths
* assemble those sets of basecalls
* create medaka training features for all those sets

The `--keep-going` flag tells Snakemake to continue with remaining tasks even
if some tasks fail (e.g. to insufficient coverage depth) while the
`--restart-times` option tells Snakemake to attempt each target a maximum of
three times. When building hundreds of targets on a cluster, we have observed
that some targets error out on the first attempt, only to succeed on later
attempts. 

Training models
---------------

Running:

.. code-block:: bash

    katuali medaka_train_replicates --keep-going

will do all the tasks of `all_medaka_train_features` and additionally launch
multiple medaka model-training replicates.

Coping with missing feature files
---------------------------------

If some of your input runs have insufficient coverage-depth for some of the
training regions, some of the training feature files will not be made. In this
case the config flag `USE_ONLY_EXISTING_MEDAKA_FEAT` can be set to true to
allow katuali to train using only those features which exist already.

.. code-block:: yaml

    USE_ONLY_EXISTING_MEDAKA_FEAT: true 

    
.. note:: Note that you need to first attempt to create all features with the
    ``all_medaka_train_features`` rule with ``USE_ONLY_EXISTING_MEDAKA_FEAT`` set to
    false, and then run ``medaka_train_replicates`` with the flag set to true. 

Refer to comments in the config (katuali/config.yaml) for further details. 
