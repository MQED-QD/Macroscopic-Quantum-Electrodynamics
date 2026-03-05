.. _tutorial-quantum-dynamics:

================
Quantum Dynamics
================

Goal
----

In this tutorial you will compute the **quantum dynamics for an open quantum
system** using the Lindblad master equation and non-Hermitian Schrödinger
equation (NHSE) solvers implemented in MQED.
The Lindblad solver is built on top of the
`QuTiP <https://qutip.org/>`_ package, while the NHSE solver uses a custom
implementation.
All input parameters are specified in a YAML configuration file, and the
output is written to an HDF5 file.

By the end you will know how to:

- run the ``mqed_lindblad`` and ``mqed_nhse`` commands with default and custom
  parameters,
- switch between the Lindblad and non-Hermitian solvers,
- locate and interpret the HDF5 output file.

.. seealso::

   See the companion theory documentation for the derivation of the Lindblad
   and non-Hermitian equations of motion (citation forthcoming).


Prerequisites
-------------

Make sure you have installed the package and activated the environment as
described in :doc:`/installation`.

This tutorial uses example data bundled under ``data/example/GF_data/``.
To generate your own Green's function cache, see
:ref:`tutorial-gf-sommerfeld`.

.. tip::

   The example data ships with the repository under ``data/example/``.
   No prior simulation run is required to follow this tutorial.


Quick start
-----------

Run the Lindblad solver from the repository root with all defaults:

.. code-block:: bash

   mqed_lindblad

Or use the non-Hermitian approach (recommended for large systems):

.. code-block:: bash

   mqed_nhse


Customising the run
-------------------

**Override individual parameters** on the command line using
`Hydra <https://hydra.cc/>`_ syntax:

.. code-block:: bash

   mqed_nhse simulation.Nmol=50 simulation.d_nm=4.0

**Use a different YAML** in the same config directory
(``configs/Lindblad/``):

.. code-block:: bash

   mqed_nhse --config-name=my_nhse_config

**Use a YAML from an arbitrary directory:**

.. code-block:: bash

   mqed_nhse --config-dir=/path/to/my/configs --config-name=my_nhse_config

The same flags work for ``mqed_lindblad``:

.. code-block:: bash

   mqed_lindblad --config-name=my_lindblad_config simulation.Nmol=20

.. tip::

   You can combine ``--config-name`` (or ``--config-dir``) with individual
   parameter overrides.  This is useful for running the same physical setup
   with different solver methods or lattice sizes.


Lindblad vs NHSE
-----------------

Both solvers propagate the excitonic dynamics of :math:`N` emitters
coupled through the dyadic Green's function.

- **Lindblad** — full density-matrix evolution via the Lindblad master
  equation, implemented with QuTiP's ``mesolve``.  Scales as
  :math:`O(N^2)` in memory.
- **NHSE** — state-vector evolution under a non-Hermitian effective
  Hamiltonian.  Much lighter for large :math:`N` and yields identical
  expectation values for single-excitation problems.

.. tip::

   For transport studies with :math:`N \gtrsim 50`, the NHSE solver is
   typically an order of magnitude faster.


Configuration reference
-----------------------

The NHSE configuration file
(``configs/Lindblad/quantum_dynamics_nhse.yaml``) is reproduced below.
The Lindblad configuration follows the same structure and is located at
``configs/Lindblad/quantum_dynamics.yaml``.

.. code-block:: yaml

   # ── Green's function input ─────────────────────────────
   greens:
     h5_path: ${oc.env:MQED_ROOT,${hydra:runtime.cwd}}/data/example/GF_data/BEM_GF_planar_Ag_665nm_height_2nm.hdf5

   height: 2                    # nm (for naming only; actual height is read from the GF file)

   # ── Material information ───────────────────────────────
   material:
     name: silver
     geometry: planar            # or 'sphere', 'nanorod', etc.

   # ── Observables ────────────────────────────────────────
   observables:
     - name: root_MSD
       kind: derived
       enabled: true

     - name: X_shift
       kind: operator

     - name: X_shift2
       kind: operator

     # Conditional versions (recommended for NHSE transport comparison)
     - name: X_shift_cond
       type: callable

     - name: X_shift2_cond
       type: callable

     - name: IPR_site
       kind: callable
       params:
         Nmol: ${simulation.Nmol}

     # Optionally: site populations
     # - name: pop_site
     #   kind: operator
     #   params:
     #     site: 1                   # 1-based site index

   # ── Simulation controls ────────────────────────────────
   simulation:
     # Time grid (picoseconds)
     t_ps:
       start: 0.0
       stop: 150.0
       output_step: 5e-3

     coupling_limit:
       enable: false
       V_hop_radius: 1
       keep_V_on_site: false
       Gamma_rule: leave           # "leave" | "same_as_V" | "diagonal_only" | "limit_by_hops"
       Gamma_hop_radius: null
       keep_Gamma_on_site: true

     Nmol: 100                     # number of emitter sites
     d_nm: 3                       # lattice spacing (nm)

     lambda_nm: 665                # emitter wavelength (nm)
     gf_method: BEM                # 'BEM' or 'Fresnel'

     # Dipole orientations
     mu_D_debye: 3.8
     mu_A_debye: 3.8
     theta_deg: 90.0
     phi_deg: 'magic'              # or a numeric value
     mode: stationary              # or 'disorder'
     # disorder_sigma_phi_deg: 15.0  # uncomment for disordered orientations

   # ── Initial condition ──────────────────────────────────
   initial_state:
     site_index: 1                 # start the exciton at site |1⟩

   # ── Solver ─────────────────────────────────────────────
   solver:
     method: NonHermitian           # 'Lindblad' or 'NonHermitian'

   # ── Output ─────────────────────────────────────────────
   output:
     filename: ${simulation.gf_method}_${material.name}_${material.geometry}_${simulation.lambda_nm}nm_N${simulation.Nmol}_height_${height}nm_inter_${simulation.d_nm}nm.hdf5

.. list-table:: Key parameters at a glance
   :header-rows: 1
   :widths: 30 50 20

   * - Parameter
     - Description
     - Default
   * - ``greens.h5_path``
     - Path to the cached dyadic Green's function HDF5 file.
     - See YAML above
   * - ``simulation.Nmol``
     - Number of emitter sites in the chain.
     - ``100``
   * - ``simulation.d_nm``
     - Intermolecular (lattice) spacing in nm.
     - ``3``
   * - ``simulation.mu_D_debye``
     - Donor transition dipole moment (Debye).  ``mu_A_debye`` is the same
       for the acceptor.
     - ``3.8``
   * - ``simulation.t_ps``
     - Time grid: ``start``, ``stop``, and ``output_step`` in picoseconds.
     - ``0 – 150 ps``, step ``5e-3``
   * - ``coupling_limit.enable``
     - Truncate dipole-dipole interaction to a finite hopping radius.
     - ``false``
   * - ``coupling_limit.V_hop_radius``
     - Number of nearest-neighbour hops retained (when ``enable: true``).
     - ``1``
   * - ``solver.method``
     - Solver backend: ``Lindblad`` (QuTiP) or ``NonHermitian``.
     - ``NonHermitian``


Expected output
---------------

After the simulation finishes, a success message is printed to the terminal
(and to the Hydra log file
``outputs/NonHermitian/.../NonHermitian.log``):

.. code-block:: text

   2025-10-24 11:40:02.818 | SUCCESS | mqed.Lindblad.run_quantum_dynamics:app_run:229
   - Simulation complete. Output saved to:
     /.../MacroscopicQED/outputs/NonHermitian/.../NAME.hdf5

The HDF5 file contains:

- **time grid** — the array of time points in picoseconds,
- **state vectors / density matrices** — the time-evolved quantum state,
- **observables** — MSD, IPR, site populations, and any other quantities
  listed in the ``observables`` section of the configuration.


What's next?
------------

- :ref:`tutorial-plotting` — visualise MSD, participation ratio, and other
  transport observables from the dynamics output.
