
.. _tests:

Tests
=====

The easiest way to test the pipeline is to run the tests, which will basecall,
assemble and polish a small dataset that comes bundled with `Katuali`.
The tests require scrappie, pomoxis, canu and nanopolish to be installed, and can be run with:

.. code-block:: bash

    make test

The output is placed in nested directories under the folder test/.
For example, if the test has run correctly the following is one of the filepaths in the 
tree that should be produced: 

.. code-block:: bash

    test/                                   
        MinIonRun1/                        # as defined in config. 
            basecall/                       
                scrappie/                  # scrappie basecaller with default options
                    align/                 # alignment of bases
                        all_contigs/       
                            25X/           # subsampling alignments
                                sub_sample_25X_ecoli_SCS110_plasmid2.calls2ref.bam
                                           # specific contig


