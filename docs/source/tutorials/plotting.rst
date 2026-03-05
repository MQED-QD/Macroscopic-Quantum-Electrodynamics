.. _tutorial-plotting:

==================
Plotting & Analysis
==================

Goal
----

In this tutorial, you will use different command lines to visualiza simulation results — field enhancement (same as previous tutorial)
displacement (MSD), inverse participation ratio (IPR), and participation ratio (PR).

Prerequisites
-------------

Make sure you have installed the package and activated the environment as
described in :doc:`/installation`.

Depending on which plots you want to generate, you may need cached output
from earlier pipeline stages:

- Field enhancement results — see :ref:`tutorial-field-enhancement`
- Quantum dynamics results — see :ref:`tutorial-quantum-dynamics`

Quick start
-----------
To plot the MSD, run the following command:
.. code-block:: bash

   mqed_plot_msd

The similar commands for square root MSD, IPR, PR and compare field enhancement from different inputs are as following:
.. code-block:: bash

   mqed_plot_sqrt_msd
   mqed_plot_IPR
   mqed_plot_PR
   mqed_compare_enhancement

Available plot commands
-----------------------

.. TODO: Fill in descriptions for each command:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Command
     - Description
   * - ``mqed_plot_msd``
     - Plot mean squared displacement over time, user can choose single curve or multiple curves in configs.
   * - ``mqed_plot_sqrt_msd``
     - Plot square-root MSD over time, similar with MSD.
   * - ``mqed_plot_IPR``
     - Plot inverse participation ratio, similar with MSD.
   * - ``mqed_plot_PR``
     - Plot participation ratio, similar with MSD.
   * - ``mqed_compare_enhancement``
     - Compare field enhancement across different data inputs, similar with MSD.

Configuration walkthrough
-------------------------

.. TODO: Paste and annotate the relevant YAML block from configs/, e.g.:
   .. code-block:: yaml

      # your annotated config here

Customising plots
-----------------

.. TODO: Describe how to change figure size, labels, colour maps, save
   formats, etc.  Reference matplotlib settings in the YAML if applicable.

Example figures
---------------

.. TODO: Add example output figures, e.g.:
   .. figure:: /_static/example_msd.png
      :width: 500
      :align: center

      MSD vs time for a 10-emitter chain above a silver surface.

What's next?
------------

- :doc:`/getting-started` — return to the workflow overview.
- :ref:`tutorial-gf-sommerfeld` — re-run with different spectral parameters
  to explore new regimes.
