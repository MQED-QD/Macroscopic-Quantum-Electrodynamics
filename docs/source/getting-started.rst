Getting started
===============

Installation
------------

.. code-block:: bash

   conda env create -f environment.yaml
   conda activate mqed
   pip install -e .

Core workflow
-------------

1. Generate a dyadic Green's function file.
2. Run Lindblad or NHSE dynamics using that Green's function.
3. Plot MSD/RMSD/IPR/PR from the dynamics outputs.

Run from repository root:

.. code-block:: bash

   mqed_GF_Sommerfeld
   mqed_nhse
   mqed_plot_msd

Configuration notes
-------------------

All workflows are configured with Hydra YAML files under ``configs/``.
You can either modify YAML directly or override keys from the CLI.

.. code-block:: bash

   mqed_nhse simulation.Nmol=50 simulation.d_nm=4.0

For detailed examples and beta-test references, see ``README.md``.
