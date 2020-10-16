
.. _medaka_train:

Medaka consensus training pipeline
==================================

It is possible to train and evaluate medaka consensus models starting from folders of
``.fast5`` or ``.fasta/q`` files in a single command.


Input data specification
------------------------

Read ``.fast5`` files should be placed under top-level folders. Multiple top-level
folders can be used, perhaps corresponding to multiple runs.

Within the ``DATA`` section of a katuali configuration file, these top-level
folders should be listed along with details of reference sequence files and sequences.
``MEDAKA_TRAIN_REGIONS`` and ``MEDAKA_EVAL_REGIONS`` define genomic regions for training
and evaluation.

In the example below data from the first two top-level directories
(``MinIonRun1`` and ``MinIonRun2``) will be used for training using ``ecoli``
and ``yeast`` reference sequences. Evaluation of the trained models will be
performed using the third and fourth top-level directories using the ``ecoli``,
``yeast``, and ``na12878_chr21`` sequences.

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

Read depths at which to create assemblies for training are specified by the
``DEPTH`` key of the katuali configuration. This list should span the range of
depths at which the model is to be used.

.. code-block:: yaml

    DEPTHS:
        [25, 50, 75, 100, 125, 150, 175, 200]


For some datasets it my not be possible to create assemblies for all reference
sequences at all depths. To avoid katuali exiting early when such trivial failures
occur the ``--keep-going`` option can be used. This allows tasks to continue
unaffected by the failure of unrelated tasks.


Creating training features
--------------------------

To create training data ("features") for medaka, ``katuali`` must:

* basecall data from all top-level directories (if ``.fast5`` s are provided),
* align all basecalls to the specified reference sequences,
* create subsampled sets of basecalls over the desired regions and depths,
* form draft assemblies from these read sets, and finally
* create medaka training features data and labels.

There is a single medaka target to perform the above tasks:

.. code-block:: bash

    katuali all_medaka_feat

``Katuali`` uses the ``Snakemake`` ``--keep-going`` flag instructs to continue processing tasks when
unrelated tasks fail. 

Having run the ``all_medaka_feat`` target, two files will be produced for
every valid combination of dataset (top-level folder), coverage depth, and reference
sequence. For example the files:

.. code-block:: bash

    4bf50792/guppy/align/senterica1/25X_prop/canu_gsz_4.8m/racon/medaka_train/medaka_train.hdf
    4bf50792/guppy/align/senterica1/25X_prop/canu_gsz_4.8m/racon/medaka_train/medaka_train_rc.hdf

will be produced for a top-level folder named ``4bf5079``, a reference sequence
``senterica1`` at coverage of ``25``-fold.


.. _training_models:

Training models
---------------

When the production of all the training data is complete, training can be commenced
by running:

.. code-block:: bash

    katuali all_medaka_train --keep-going

This step requires the the use of GPUs.  

.. note:: Note that to ``tensorflow-gpu`` must be installed in your medaka environment for medaka training. 


.. _missing_feat:

Coping with missing feature files
---------------------------------

If input datasets have insufficient coverage-depth for some of the
training regions, some training feature files will not be produced. In this
case the config flag ``USE_ONLY_EXISTING_MEDAKA_FEAT`` can be set to ``true`` to
allow katuali to train using only those features which exist already.

.. code-block:: yaml

    USE_ONLY_EXISTING_MEDAKA_FEAT: true 

    
.. note:: Note that you need to first attempt to create all features with the
    ``medaka_train_feat`` rule with ``USE_ONLY_EXISTING_MEDAKA_FEAT`` set to
    false, and then run ``all_medaka_train`` with the flag set to true. 

Refer to comments in the katuali configuration file for further details. 
