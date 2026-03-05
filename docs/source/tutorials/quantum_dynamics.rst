.. _tutorial-quantum-dynamics:

================
Quantum Dynamics
================

Goal
----

In this tutorial you will compute the **quantum dynamics for open quantum system** 
using the Lindblad master equation and non-Hermitian Schrödinger equation solvers implemented in MQED-QD.
The orginal implementation of the Lindblad solver is based on the `QuTiP <https://qutip.org/>`_ package.
The input parameters are specified in a YAML configuration file, and the output is written to an HDF5 file.

By the end you will know how to:
- run the ``mqed_lindblad`` and ``mqed_nhse`` commands with default and custom parameters.
- locate and interpret the HDF5 output file.

Prerequisites
-------------

Make sure you have installed the package and activated the environment as
described in :doc:`/installation`.

You also need a cached Green's function HDF5 file.
See :ref:`tutorial-gf-sommerfeld` for how to generate one.

Quick start
-----------

Run from the repository root with all defaults:
.. code-block:: bash

   mqed_lindblad

Or for the non-Hermitian approach (recommended in our simulations):
.. code-block:: bash

   mqed_nhse

Lindblad vs NHSE
----------------


See (put citation later) for the theoretical background.

Configuration walkthrough
-------------------------
The non-Hermitian Schrödinger equation YAML file looks like:
.. code-block:: yaml
   #quantum_dynamics_nhse.yaml
   greens:
      h5_path: ${oc.env:MQED_ROOT,${hydra:runtime.cwd}}/data/BEM_cache/BEM_GF_planar_Ag_665nm_height_2nm.hdf5 # update as needed

   height: 2 # nm, for naming purposes only, the actual height is read from the GF file
   #Material information
   material:
      name: silver
      geometry: planar # or 'sphere', 'planar', match with your simulation to avoid confusion.

   observables:
   - name: root_MSD
      kind: derived
      enabled: true

   - name: X_shift
      kind: operator

   - name: X_shift2
      kind: operator
   
   # conditional versions (recommended for NHSE transport comparison)
   - name: X_shift_cond
      type: callable

   - name: X_shift2_cond
      type: callable

   - name: IPR_site
      kind: callable
      params:
         Nmol: ${simulation.Nmol}    # pass through from your sim config

   # # Optionally: one or many site populations
   # - name: pop_site
   #   kind: operator
   #   params:
   #     site: 1                     # 1-based site index (your code uses +1 internally)

   # Add more site populations by repeating the item:
   # - { name: pop_site, kind: operator, params: { site: 10 } }
   # Simulation controls
   simulation:
   # time grid in picoseconds
      t_ps:
         start: 0.0
         stop: 150.0
         output_step: 5e-3

      coupling_limit:
         enable: false
         V_hop_radius: 1
         keep_V_on_site: false
         Gamma_rule: leave # "leave" | "same_as_V" | "diagonal_only" | "limit_by_hops"
         Gamma_hop_radius: None
         keep_Gamma_on_site: true

      # system
      Nmol: 100 # number of sites
      d_nm: 3 # lattice spacing (nm)

      #emitter frequency
      lambda_nm: 665        # wavelength of emitter in nm

      # Green's function data simulation method
      gf_method: BEM  # 'BEM' or 'Fresnel'


      # dipoles / orientation
      mu_D_debye: 3.8
      mu_A_debye: 3.8
      theta_deg: 90.0
      phi_deg: 'magic' # or a number
      mode: stationary # or 'disorder'
      # disorder_sigma_phi_deg: 15.0 # uncomment for disordered orientations


   # Initial condition
   initial_state:
   site_index: 1 # start the exciton at site |1>


   # Solver selection
   solver:
   method: NonHermitian # 'Lindblad' or 'NonHermitian'


   # Output file (goes under Hydra's run dir)
   output:
   filename: ${simulation.gf_method}_${material.name}_${material.geometry}_${simulation.lambda_nm}nm_N${simulation.Nmol}_height_${height}nm_inter_${simulation.d_nm}nm.hdf5

The Lindblad solver configs are inside ``configs/Lindblad/quantum_dynamics.yaml``. We separate these two configs for convenience.
Key parameters
--------------


.. list-table::
   :header-rows: 1
   :widths: 25 50 25

   * - Parameter
     - Description
     - Default
   * - ``h5_path``
     - Input path of dyadic Green's function hdf5 file. 
     - ``${oc.env:MQED_ROOT,${hydra:runtime.cwd}}/data/BEM_cache/BEM_GF_planar_Ag_665nm_height_2nm.hdf5 # update as needed``
   * - ``Nmol``
     - Number of emitters in simulation.
     - ``100``
   * - ``mu_D_debye``
     - dipole moment intensity of donor in Debye unit, same apply to ``mu_A_debye`` for acceptor
     - ``3.8``
   * - ``d_nm``
     - intermolecular distance for the simulation
     - ``3``
   * - ``coupling_limit.enabled, coupling_limit.V_hop_radius``
     - This is a tool to help user truncate the dipole-dipole interaction (DDI). If enabled, user can choose ``V_hop_radius`` to **1 for nearest-neighbor truncation**
   * - ''false``
Expected output
---------------

.. TODO: Describe what the command prints and what files it creates, e.g.:
   The solver writes an HDF5 file to the Hydra output directory containing
   the time-evolved density matrix / state vector ...
After the simulation finishes, a success message is printed to the terminal
(and to the Hydra log file
``outputs/NonHermitian/.../NonHermitian.log``):

.. code-block:: text

   2025-10-24 11:40:02.818 | SUCCESS | mqed.Lindblad.run_quantum_dynamics:app_run:229 
   - Simulation complete. Output saved to:
     /.../MacroscopicQED/outputs/NonHermitian/.../NAME.hdf5
What's next?
------------

- :ref:`tutorial-plotting` — visualise MSD, participation ratio, and other
  transport observables from the dynamics output.
