.. _tutorial-field-enhancement:

=================
Field Enhancement
=================

Goal
----
In this tutorial, you will compute the dipole-dipole interaction (DDI) raio and generalized dissipation ratio,
or **electric field enhancement**, :math:`V_{\alpha\beta}/V_{0,\alpha\beta}` and :math:`\Gamma_{\alpha\beta}/\Gamma_{0,\alpha\beta}`
from dyadic Green's function. The definition is given in our **theory/MacroscopicQED**

.. seealso::

   :ref:`theory/MacroscopicQED` for relation of electric field enhancement 
   with dyadic Green's function.

Prerequisites
-------------

Make sure you have installed the package and activated the environment as
described in :doc:`/installation`.

You also need a cached Green's function HDF5 file.
See :ref:`tutorial-gf-sommerfeld` for how to generate one.

Quick start
-----------

.. code-block:: bash

   mqed_FE

This uses the YAML configuration at ``configs/analysis/FE.yaml``.
The default settings computes from input file ``/data/GF_cache/Fresnel_GF_planar_Ag_height_8nm_665nm.hdf5`` 
and the orientation of donor-acceptor are at magic angle in XY plane. The plot settings can be modified by user.
The default settings plot :math:`V_{\alpha\beta}/V_{0,\alpha\beta}` and :math:`\Gamma_{\alpha\beta}/\Gamma_{0,\alpha\beta}` together 
and denoted as 'real' and 'imag' in **plot_settings.components** which can be modified by user. 

Configuration walkthrough
-------------------------
The full default configuration file
(``configs/analysis/FE.yaml``) is reproduced below with
annotations.

.. code-block:: yaml
   input_file: ${oc.env:MQED_ROOT,${hydra:runtime.cwd}}/data/GF_cache/Fresnel_GF_planar_Ag_height_8nm_665nm.hdf5 # 

   orientations:
   donor:    { theta_deg: 90.0, phi_deg: "magic" }
   acceptor: { theta_deg: 90.0, phi_deg: "magic" }

   plot_settings:
   save_plot: true
   dpi: 400

   # choose x-range by value or by indices, not both  
   x_range_nm: [0.0,100.0] #-> masks points where 1 ≤ Rx_nm ≤ 100
   # x_index_range: [0, 9] #-> masks first 10 points (index 0 to 9)

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

   xscale: linear      # or "log"
   yscale: linear      # or "log"

   # Optional axis limits (override selection if you want)
   xlim: [0, 50]  # e.g., [0, 500]
   ylim: null  # e.g., [0, 10]
   filename_prefix: "enhancement_magic_angle"
   grid: true

Key parameters
--------------

.. TODO: Fill in the parameter table:

.. list-table::
   :header-rows: 1
   :widths: 25 50 25

   * - Parameter
     - Description
     - Default
   * - ``input_file``
     - Source of the dyadic Green's function hdf5 file.
     - ``${oc.env:MQED_ROOT,${hydra:runtime.cwd}}/data/GF_cache/Fresnel_GF_planar_Ag_height_8nm_665nm.hdf5``
   * - ``donor``
     - Orientation of the donor.
     - ``{ theta_deg: 90.0, phi_deg: "magic" }``
   * - ``acceptor``
     - Orientation of the acceptor
     - ``{ theta_deg: 90.0, phi_deg: "magic" }``
   * - ``components``
     - The plot components of :math:`V_{\alpha\beta}/V_{0,\alpha\beta}` or :math:`\Gamma_{\alpha\beta}/\Gamma_{0,\alpha\beta}` chose by real and imag.
     - ``["real", "imag"]``

Expected output
---------------

After the simulation finishes, a success message is printed to the terminal
(and to the Hydra log file
``outputs/FE/.../FE.log``):

.. code-block:: text

   2025-10-24 11:40:02.818 | SUCCESS | mqed.analysis.FE:plot_field_enhancement:149
   - Simulation complete. Output saved to:
     /.../MacroscopicQED/outputs/FE/Y-M-D/H-M-S/enhancement_magic_angle_1.864eV.png

What's next?
------------

- :ref:`tutorial-quantum-dynamics` — use the Green's function and field
  enhancement results to drive open quantum dynamics.
- :ref:`tutorial-plotting` — visualise the enhancement spectra.
