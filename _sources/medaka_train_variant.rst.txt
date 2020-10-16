
.. _medaka_train_variant:

Medaka variant-calling training pipeline
========================================

Katuali supports the training of medaka variant-calling models from ``.fast5``
files, or preprocessed data such as basecalls or alignments if they are
available.

Models can be trained both for the initial SNP calling from mixed reads prior
to phasing, referred to below as SNP-calling models, and for final variant
(SNP and indel) calling from phased reads, referred to below as variant-calling
models. See the medaka `docs
<https://nanoporetech.github.io/medaka/snp.html#>`_ for more information on the
medaka variant calling pipeline.

A specialised config file (``config_variant.yaml``) for medaka variant calling
is provided by running 
 
 .. code-block:: bash

    katuali_config --variant my_config_variant.yaml. 


Input data specification
------------------------

The data required for variant calling training is separated into two 
directory structures:

* training data 
* variant calling resources: references and truth data.


Training data
-----------------------------

As for the :ref:`medaka_train` read ``.fast5`` files should be placed 
under top-level folders. Multiple top-levelfolders can be used.

Within the ``DATA`` section of a katuali configuration file, these top-level
folders should be listed along with details of the reference sequence.
``MEDAKA_TRAIN_REGIONS``  lists the chromosomes included in training

In the example below data from the run  (``dna_prom``) will be used for training for 
chromosomes 15 -19 of the GIAB sample NA24385.


.. code-block:: yaml

    DATA:
        'dna_prom': 
            'REFERENCE': '/medaka_variant_resources/grch38/grch38.fna.gz.fasta'   
            'MEDAKA_TRAIN_REGIONS': ['chr15', 'chr16', 'chr17', 'chr18', 'chr19']
            'MEDAKA_EVAL_REGIONS': []
            'SAMPLE: "NA24385"


Variant calling resources
-----------------------------
            
The protocol for training for variant calling require various reference and
truth data. This can be configured within the ``VARIANT_RESOURCES`` section of
a katuali config file.  The default config is set up for the GIAB sample
NA24385 (HG002), with configurable URLs which will be used to download the GIAB
high-confidence variants (vcf file) and regions (bed file).

Medaka training features for variant calling (but not SNP-calling) are formed
from pileups of reads aligned to the GRCh38 reference scuffed up with variants
from the 1000 Genomes Project.  The ``VARIANTS_TO_ADD`` config entry defines
the template filepath of the per-chromosome VCF file which will be used to
scuff up the reference. If the file is not present, it will automatically be
downloaded. It is possible to use alternative variants to scuff the reference
by changing ``VARIANTS_TO_ADD`` to point at your own VCF. Similarly, if the
GRCh38 reference is not not found at the expected filepath
(``medaka_variant_resources/grch38/grch38.fna.gz``), it will automatically be
downloaded. If you already have this file downloaded and wish to avoid
downloading it again, you can place a symlink at the expected filepath. 

It is possible to train from a different, or multiple samples by extending the
``SAMPLES`` config section:


.. code-block:: yaml

    VARIANT_RESOURCES:
        "SAMPLES":
            "NA24385":
                "TRUTH_HIGH_CONF_BED_URL": "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed" 
                "TRUTH_HIGH_CONF_VCF_URL": "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz" 

        "VARIANTS_TO_ADD": "medaka_variant_resources/grch38/1kgenomes/ALL.{chr}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"


Based on the information provided in the config, ``Katuali`` will download and
process reference and truth data into various intermediate files required for
creating training features. This data required is stored under a directory
structure ``medaka_training_resources``

There should be no need for the user to delve into the contents of the
directory, but its contents are nevertheless documented here. The following
example shows the folder structure for the grch38 reference and the NA24385 sample.


.. code-block:: bash

    medaka_variant_resources/         
        grch38/                         # grch38 reference directory (currently only the grch38 is supported).
            grch38.fna.gz               # grch38 reference file which will be automatically downloaded. 
            grch38.fna.gz_per_chr/      # per chromosome references stored as chrX.fasta, automatically created. 
            1kgenomes/                  # VARIANTS_TO_ADD: vcf file per chromosome, automatically downloaded. 
            NA24385/                    # SAMPLE Truth reference {SAMPLE.fasta}: NA24385.fasta, created automatically from the truth vcf and reference. 
                NA24385_per_chr/        # Phased reference per chromosome (e.g. chr18.fasta
                                        # chr18_NA24385_maternal.fasta and chr18_NA24385_paternal.fasta, created automatically).
        scuffed_ref/
            NA24385/ 
                chrN/                   # Scuffed up reference ref_edited.fasta, ref_edited_rc.fasta, automatically created. 
                                                

Creating training features
--------------------------

SNP and variant-calling models can be trained in a single command using the
``all_medaka_train`` pipeline, see :ref:`medaka_train_variant_models`.

As with the :ref:`medaka_train`, it is also possible to just create the
features for training outside of katuali. 

Features for training SNP-calling models can be created with the
``all_medaka_snp_feat`` pipeline: 

.. code-block:: bash

    katuali all_medaka_snp_feat 


This will create a SNP-calling feature file for every valid combination of
dataset (top-level folder), region and coverage depth . For example the file

.. code-block:: bash

    dna_prom/guppy/align/chr15/10X_prop/medaka_diploid_snp_features/dna_prom_chr15_10X_medaka_train.hdf

will be produced for a top-level folder named ``dna_prom``, chromosome 15, at
coverage of ``10``-fold.


Similarly, features for training variant-calling models can be created with the
``all_medaka_variant_feat`` pipeline:

.. code-block:: bash

    katuali all_medaka_variant 

This will create two variant-calling feature files for every valid combination
of dataset (top-level folder), region and coverage depth. For example the files

.. code-block:: bash

    dna_prom/guppy/align/chr16/calls2scuffed_ref/20X_prop/medaka_features/medaka_train_variant_dna_prom_chr16_20X.hdf
    dna_prom/guppy/align/chr16/calls2scuffed_ref/20X_prop/medaka_features/medaka_train_variant_dna_prom_chr16_20X_rc.hdf


will be produced for a top-level folder named ``dna_prom``, chromosome 16, at
coverage of ``20``-fold.


.. _medaka_train_variant_models:

Training models
---------------

SNP and variant-calling models can be trained in a single command using the
same ``all_medaka_train`` pipeline as used for training consensus models, see :ref:`training_models` 

The features used for model training are fully configurable in the config,
allowing a single ``all_medaka_train`` pipeline to train either a medaka
consensus model, or medaka SNP and/or variant-calling models. 

This can be controled via the ``MEDAKA_TRAIN_REPLICATES`` config entry.

This example from the default config trains three medaka-consensus model replicates using the
``all_medaka_feat`` pipeline to create medaka-consensus training features.  

.. code-block:: yaml

    # Run multiple training replicates - output will be in medaka_train_{key},
    # values should be a key of the PIPELINES dictionary in this file. In simple
    # cases this allows running technical replicates of the training, but also
    # allows the pipeline to be changed to for example create features in a
    # different manner. For the latter change the value component to an alternative
    # user defined pipeline.
    MEDAKA_TRAIN_REPLICATES:
        "_rep_1": "all_medaka_feat" 
        "_rep_2": "all_medaka_feat"
        "_rep_3": "all_medaka_feat"


In contrast, this example from the variant-calling config trains three medaka
SNP-calling replicates and three medaka variant-calling replicates, each using their
respective feature-generation pipelines. 

.. code-block:: yaml

    MEDAKA_TRAIN_REPLICATES:
        "_snp_rep_1": "all_medaka_snp_feat" 
        "_snp_rep_2": "all_medaka_snp_feat"
        "_snp_rep_3": "all_medaka_snp_feat"
        "_variant_rep_1": "all_medaka_variant_feat" 
        "_variant_rep_2": "all_medaka_variant_feat"
        "_variant_rep_3": "all_medaka_variant_feat"


Missing feature files
---------------------------------

As for the :ref:`medaka_train`, if input datasets have
insufficient coverage-depth for some of the training chromosomes, some training
feature files will not be produced. See :ref:`missing_feat` for how to cope
with this. 
