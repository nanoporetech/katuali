
.. _faq:

Frequently asked questions
==========================

.. existing:

How do I start a pipeline using existing basecalls?
  Information on how to do this can be found at
  :ref:`starting_from_basecalls`.

.. medakatrain:

How do I use Katuali to train a model for medaka?
  The :ref:`train_medaka` section has details on how to accomplish this. The
  standard config. file template should also act as a guide.


.. cluster:

How do I use Katuali with a cluster?
  The :ref:`using_cluster` section has details of how to use Snakemake's
  cluster capabilities.

.. libraries:

I get an error about not finding libraries. What do I do?
  This commonly occurs when using ``medaka`` with a GPU e.g.:

  .. code-block:: bash

      ImportError: libcublas.so.9.0: cannot open shared object file: No such file or directory
  
  Setting ``LD_LIBRARY_PATH`` in the
  `config <https://github.com/nanoporetech/katuali/blob/master/katuali/data/config.yaml#L78>`_
  should resolve the issue.

.. features:

Can you implement using X program for performing task Y?
  Please feel free to make a feature request on
  `github <https://github.com/nanoporetech/katuali/issues>`_.
