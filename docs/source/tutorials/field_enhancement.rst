.. _tutorial-field-enhancement:

=================
Field Enhancement
=================

Goal
----

In this tutorial you will compute the **electric field enhancement** —
the dipole-dipole interaction (DDI) ratio
:math:`V_{\alpha\beta}/V_{0,\alpha\beta}` and the generalised dissipation
ratio :math:`\Gamma_{\alpha\beta}/\Gamma_{0,\alpha\beta}` — from a
previously cached dyadic Green's function.

By the end you will know how to:

- run the ``mqed_FE`` command with default and custom parameters,
- choose donor and acceptor orientations,
- select which components (real, imaginary, or both) to plot,
- locate and interpret the output figure.

.. seealso::

   :doc:`/theory/MacroscopicQED` for the relation between the electric field
   enhancement and the dyadic Green's function.


Prerequisites
-------------

Make sure you have installed the package and activated the environment as
described in :doc:`/installation`.

This tutorial uses example data bundled under ``data/example/GF_data/``.
To generate your own Green's function cache, see
:ref:`tutorial-gf-sommerfeld`.


Quick start
-----------

Run from the repository root with all defaults:

.. code-block:: bash

   mqed_FE

This uses the YAML configuration at ``configs/analysis/FE.yaml``.
The default settings compute the enhancement from the bundled example file
``data/example/GF_data/Fresnel_GF_planar_Ag_height_8nm_665nm.hdf5`` with
donor and acceptor oriented at the magic angle in the XY plane.
Both :math:`V_{\alpha\beta}/V_{0,\alpha\beta}` (``real``) and
:math:`\Gamma_{\alpha\beta}/\Gamma_{0,\alpha\beta}` (``imag``) are
plotted by default.

.. tip::

   The example data ships with the repository under ``data/example/``.
   No prior simulation run is required to follow this tutorial.

.. tip::

   Override any key on the command line using
   `Hydra <https://hydra.cc/>`_ syntax, for example:

   .. code-block:: bash

      mqed_FE orientations.donor.phi_deg=0.0 orientations.acceptor.phi_deg=0.0


Configuration reference
-----------------------

The full default configuration file
(``configs/analysis/FE.yaml``) is reproduced below with
annotations.

.. code-block:: yaml

   # ── Input ───────────────────────────────────────────────
   input_file: ${oc.env:MQED_ROOT,${hydra:runtime.cwd}}/data/example/GF_data/Fresnel_GF_planar_Ag_height_8nm_665nm.hdf5

   # ── Orientations ────────────────────────────────────────
   orientations:
     donor:    { theta_deg: 90.0, phi_deg: "magic" }
     acceptor: { theta_deg: 90.0, phi_deg: "magic" }

   # ── Plot settings ───────────────────────────────────────
   plot_settings:
     save_plot: true
     dpi: 400

     # choose x-range by value or by indices, not both
     x_range_nm: [0.0, 100.0]       # mask points where 0 ≤ Rx_nm ≤ 100
     # x_index_range: [0, 9]        # alternative: mask first 10 points

     components: ["real", "imag"]    # options: ["real"], ["imag"], ["real","imag"]

     # Text & style
     xlabel: "Donor–Acceptor Distance (nm)"
     ylabel: "Enhancement"
     title_template: "Emitter energy = {energy:.3f} eV"
     legend:
       real_label: "$V_{\\alpha\\beta}/ V_{0,\\alpha\\beta}$"
       imag_label: "$\\Gamma_{\\alpha\\beta}/ \\Gamma_{0,\\alpha\\beta}$"

     lw: 1
     real_style: "r--"
     imag_style: "b--"

     xscale: linear                  # or "log"
     yscale: linear                  # or "log"

     # Optional axis limits
     xlim: [0, 50]                   # e.g., [0, 500]
     ylim: null                      # e.g., [0, 10]
     filename_prefix: "enhancement_magic_angle"
     grid: true

.. list-table:: Key parameters at a glance
   :header-rows: 1
   :widths: 30 50 20

   * - Parameter
     - Description
     - Default
   * - ``input_file``
     - Path to the cached dyadic Green's function HDF5 file.
     - ``data/example/GF_data/...``
   * - ``orientations.donor``
     - Polar and azimuthal angles of the donor dipole.
     - ``{ theta_deg: 90.0, phi_deg: "magic" }``
   * - ``orientations.acceptor``
     - Polar and azimuthal angles of the acceptor dipole.
     - ``{ theta_deg: 90.0, phi_deg: "magic" }``
   * - ``plot_settings.components``
     - Which ratios to plot: ``"real"`` for :math:`V/V_0`, ``"imag"`` for :math:`\Gamma/\Gamma_0`.
     - ``["real", "imag"]``
   * - ``plot_settings.x_range_nm``
     - Horizontal distance window for the plot (nm).
     - ``[0.0, 100.0]``
   * - ``plot_settings.xlim``
     - Matplotlib axis limits for the x-axis.
     - ``[0, 50]``


Expected output
---------------

After the simulation finishes, a success message is printed to the terminal
(and to the Hydra log file
``outputs/FE/.../FE.log``):

.. code-block:: text

   2025-10-24 11:40:02.818 | SUCCESS | mqed.analysis.FE:plot_field_enhancement:149
   - Simulation complete. Output saved to:
     /.../MacroscopicQED/outputs/FE/Y-M-D/H-M-S/enhancement_magic_angle_1.864eV.png

The output is a **PNG** figure showing the selected enhancement components as
a function of donor-acceptor distance.

.. figure:: /_static/FE_result/example_FE.png
   :width: 500
   :align: center

   Field enhancement vs donor-acceptor distance.

What's next?
------------

- :ref:`tutorial-quantum-dynamics` — use the Green's function and field
  enhancement results to drive open quantum dynamics.
- :ref:`tutorial-plotting` — visualise the enhancement spectra and other
  transport observables.
