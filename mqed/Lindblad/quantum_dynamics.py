from __future__ import annotations
import numpy as np 
from qutip import *
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
from typing import Callable, Iterable, List, Optional, Sequence, Tuple, Dict, Any, Union
from loguru import logger


from mqed.utils.au_unit import au_to_eV, ps_to_au
from mqed.Lindblad.ddi_matrix import _phi_wrapped_normal_deg, build_ddi_matrix_from_Gslice
from mqed.utils.orientation import resolve_angle_deg, spherical_to_cartesian_dipole
from mqed.Lindblad.coupling_filter import enforce_coupling_range

ArrayLike = Union[np.ndarray, Sequence[float]]

@dataclass(frozen=True)
class CouplingLimitConfig:
    enable: bool = False
    V_hop_radius: Optional[int] = None     # e.g. 1 → NN only
    keep_V_on_site: bool = False
    Gamma_rule: str = "leave"              # "leave" | "same_as_V" | "diagonal_only" | "limit_by_hops"
    Gamma_hop_radius: Optional[int] = None
    keep_Gamma_on_site: bool = True

@dataclass(frozen=True)
class SimulationConfig:
    tlist: np.ndarray                        # time grid
    # Common physical params (extend as needed)
    emitter_frequency: float    # emitter frequency in eV
    Nmol: int                   #number of molecules
    Rx_nm: np.ndarray           # array of intermolecular distance, must contain 0, 1d, 2d, ...
    d_nm: float                 # intermolecular distance, unit is nm
    mu_D_debye: float           # donor dipole
    mu_A_debye: Union[None,float]            # acceptor dipole, default is none
    theta_deg: float            # polar angle of dipole
    phi_deg: Union[str,float]           # azimuthal angle of dipole, string used for magic angle
    disorder_sigma_phi_deg: Union[None,float] # standard deviation of azimuthal angle
    mode: str                   # The string for angle-ordered system or disorder system.
    coupling_limit: CouplingLimitConfig = field(default_factory=CouplingLimitConfig)

@dataclass
class SimulationResult:
    tlist: np.ndarray
    expectations: Dict[str, np.ndarray] = field(default_factory=dict)

class QuantumDynamics(ABC):
    """
    Abstract base for quantum dynamics solvers sharing parameters and Hamiltonian building.
    Subclasses implement the actual propogation method in the 'evolve'.
    """

    def __init__(self, config:SimulationConfig):
        self.cfg = config
        self.dim = self.cfg.Nmol + 1
        self.omega_M = self.cfg.emitter_frequency /au_to_eV
    
    def build_hamiltonian(self,Green, seed:Union[int,None] = None) -> Qobj:
        """
        Construct Hamiltonian of the system
        ..math::
            &\hat{H}_\mathrm{M} = \sum_{\alpha=1}^{N_\mathrm{M}} \hbar \omega_\mathrm{M} \hat{\sigma}^{(+)}_\alpha \hat{\sigma}^{(-)}_\alpha, \\
            &\hat{\mathcal{H}}_\mathrm{CP}^{\mathrm{Sc}} = \sum_{\alpha=1}^{N_\mathrm{M}} \Delta_{\alpha}^{\text{Sc}} \, \hat{\sigma}^{(+)}_\alpha \hat{\sigma}^{(-)}_\alpha , \\
            & \hat{\mathcal{H}}_\mathrm{RDDI} =  \sum_{\alpha,\beta (\alpha\neq \beta) }^{N_\mathrm{M}} V_{\alpha\beta} \, \hat{\sigma}^{(+)}_\alpha \hat{\sigma}^{(-)}_\beta. 
        """

        H_np = np.zeros((self.dim,self.dim), dtype=complex)
        
        if self.cfg.mode == 'stationary':
            logger.info("Building Hamiltonian for angle-ordered system.")
            self.phi_deg   = resolve_angle_deg(self.cfg.phi_deg)
            # Convert angles to Cartesian vectors
            p_donor = spherical_to_cartesian_dipole(self.cfg.theta_deg,
                                                    self.phi_deg)
            p_acceptor = spherical_to_cartesian_dipole(self.cfg.theta_deg,
                                                    self.phi_deg)
        
        # here V and Gamma have unit of eV. 
            self.V_ab, self.Gamma_ab = build_ddi_matrix_from_Gslice(
                G_slice = Green,
                Rx_nm = self.cfg.Rx_nm,
                energy_emitter = self.cfg.emitter_frequency,
                N_mol= self.cfg.Nmol,
                d_nm = self.cfg.d_nm,
                uD = p_donor,
                uA = p_acceptor,
                mu_D_debye = self.cfg.mu_D_debye,
                mu_A_debye = self.cfg.mu_A_debye,
                mode = self.cfg.mode,
                phi_deg = self.phi_deg,
                theta_deg= self.cfg.theta_deg,
                disorder_sigma_phi_deg= self.cfg.disorder_sigma_phi_deg,
                disorder_seed= None
            )
        else:  # disorder mode
            logger.info("Building Hamiltonian for orientation-disordered system.")
            self.V_ab, self.Gamma_ab = build_ddi_matrix_from_Gslice(
                G_slice = Green,
                Rx_nm = self.cfg.Rx_nm,
                energy_emitter = self.cfg.emitter_frequency,
                N_mol= self.cfg.Nmol,
                d_nm = self.cfg.d_nm,
                mu_D_debye = self.cfg.mu_D_debye,
                mu_A_debye = self.cfg.mu_A_debye,
                mode = self.cfg.mode,
                phi_deg = self.cfg.phi_deg,
                theta_deg= self.cfg.theta_deg,
                disorder_sigma_phi_deg= self.cfg.disorder_sigma_phi_deg,
                disorder_seed= seed
            )
        cl = getattr(self.cfg, "coupling_limit", None)
        if cl is not None and cl.enable:
            self.V_ab, self.Gamma_ab = enforce_coupling_range(
                self.V_ab, self.Gamma_ab,
                V_hop_radius    = cl.V_hop_radius,
                keep_V_on_site  = cl.keep_V_on_site,
                Gamma_rule      = cl.Gamma_rule,
                Gamma_hop_radius= cl.Gamma_hop_radius,
                keep_Gamma_on_site = cl.keep_Gamma_on_site,
            )
            if cl.V_hop_radius ==1:
                logger.info("Applied coupling limit: V to nearest-neighbours only.")
            elif cl.V_hop_radius ==2:
                logger.info("Applied coupling limit: V to next-nearest-neighbours only.")
            else:
                logger.info("Applied coupling limit to V with hop radius %s.", str(cl.V_hop_radius))

        H_np[1:, 1:] = self.V_ab / au_to_eV
        np.fill_diagonal(H_np[1:, 1:], self.omega_M)
        self.Gamma_ab = self.Gamma_ab / au_to_eV   # convert to a.u.
        # breakpoint()
        self.Hamiltonian = Qobj(H_np, dims=[[self.dim], [self.dim]])
        return self.Hamiltonian
    
    def build_collapse_ops(self) -> list[Qobj]:
        """
        Build the standard Lindblad collapse operator from Gamma. Folow the https://en.wikipedia.org/wiki/Lindbladian,
        the general Lindblad equation has the form:
        ..math::
            \dot{\rho} = -\frac{i}{\hbar}[H, \rho] + \sum_{n,m} h_{nm} \left( A_n \rho A_m^{\dagger} - \frac{1}{2} \{ A_m^{\dagger} A_n, \rho \} \right)
        This can be diagonalized through unitary transformation u:
        ..math::
            u^{\dagger} h u = 
            \begin{bmatrix} 
            \gamma_1 & 0 & \cdots & 0 \\
            0 & \gamma_2 & \cdots & 0 \\
            \vdots & \vdots & \ddots & \vdots \\
            0 & 0 & \cdots & \gamma_{N^2-1}
            \end{bmatrix}
        The eigenvalues are non-negative. If we define another orthonormal operator basis:
        ..math::
            L_i = \sum_{j} u_{ji} A_j
        We can get the standard Lindblad equation:
        ..math::
            \dot{\rho} = -\frac{i}{\hbar}[H, \rho] + \sum_{i} \gamma_i \left( L_i \rho L_i^{\dagger} - \frac{1}{2} \{ L_i^{\dagger} L_i, \rho \} \right).
        return:
            c_ops: list[Qobj]: a list of collapse operator.
        """

            # ---- basis: |0>, |1>, ..., |N_mol| ; σ_j^- = |0><j|
        sigma_minus = [projection(self.dim, 0, j+1) for j in range(self.cfg.Nmol)]

        # ---- make Γ exactly Hermitian (Lindblad requires h ≥ 0 Hermitian)
        G = 0.5 * (self.Gamma_ab + self.Gamma_ab.conj().T)

        # ---- diagonalize (Hermitian eigenproblem)
        evals, evecs = np.linalg.eigh(G)

        # ---- numerical tolerance
        # choose something relative to the norm of Γ
        gscale = max(np.linalg.norm(G, 2), 1.0)
        tol = 1e-12 * gscale

        # clip tiny negatives and check for real PSD
        evals = np.real_if_close(evals, tol=1000)
        evals[evals < 0 & (evals > -tol)] = 0.0
        if np.any(evals < -tol):
            raise ValueError(
                f"Gamma is not positive semidefinite (min eig={evals.min():.3e}). "
                "Dynamics would not be CP. Check units/derivation."
            )

        # ---- build L_k
        c_ops: List[Qobj] = []
        for k, gamma_k in enumerate(evals):
            if gamma_k <= 0.0:
                continue  # null rate → no jump
            v = evecs[:, k]  # complex components v_j^(k)
            Lk = qzero(self.dim)
            for j, vj in enumerate(v):
                if abs(vj) > 0:
                    Lk += vj * sigma_minus[j]
            c_ops.append(np.sqrt(gamma_k) * Lk)

        self.c_ops = c_ops
        return self.c_ops
    
    @abstractmethod
    def evolve(self,
                rho_or_psi: Qobj,
                e_ops: Optional[Dict[str, Qobj]] = None,
                options: Optional[Dict[str,Any]] = None,)->SimulationResult:
        """
        Run time evolution and return the physical quantities.
        Args:
            rho_or_psi: Qobject, the initial density matrix or wavefunction.
            e_ops: expectation values of the operators.
            options: optional restrict for differential equation solver. See ... 
        """
        raise NotImplementedError
    
    def _pack_result(self,
                    tlist:np.ndarray,
                    states: List[Qobj],
                    eops: Optional[Dict[str, Qobj]])-> SimulationResult:
        """
        Pack the simulation result. 

        """
        out = SimulationResult(tlist = tlist, states = states)
        if eops:
            out.expectations = {name: expect(op, states) for name, op in eops.items()}
        return out 


class LindbladDynamics(QuantumDynamics):
    """
    Standard Lindblad dynamics solver provide by qutip.mesolver(). 
    ..math::
        \dot{\rho} = -\frac{i}{\hbar}[H, \rho] + \sum_{i} \gamma_i \left( L_i \rho L_i^{\dagger} - \frac{1}{2} \{ L_i^{\dagger} L_i, \rho \} \right).
    """
    def __init__(self, config, GreensFunction, seed: Union[int, None]=None):
        super().__init__(config)
        self.cfg = config
        # self.dim = self.cfg.Nmol + 1
        # self.rho = fock_dm(self.dim, 1)
        self.t_evaluation = self.cfg.tlist * ps_to_au
        self.GF = GreensFunction
        self.seed = seed
    
    def evolve(self,
                rho_or_psi: Qobj,
                e_ops: Optional[Dict[str, Qobj]] = None,
                options: Optional[Dict[str,Any]] = None,)->SimulationResult:
        """
        Run time evolution and return the physical quantities.
        Args:
            rho_or_psi: Qobject, the initial density matrix or wavefunction.
            e_ops: expectation values of the operators.
            options: optional restrict for differential equation solver.  
        """

        Hamiltonian = self.build_hamiltonian(self.GF, seed = self.seed)
        c_ops = self.build_collapse_ops()

        result = mesolve(Hamiltonian, rho_or_psi, self.t_evaluation, c_ops = c_ops,
                        e_ops = list(e_ops.values()) if e_ops else None, options= options)
        
        # If e_ops provided, map back by name; otherwise compute later if needed
        if e_ops:
            named = {name: np.asarray(result.expect[n]) for n, name in enumerate(e_ops.keys())}
        else:
            named = {}
        return SimulationResult(tlist=self.cfg.tlist ,expectations=named)

class NonHermitianSchDynamics(QuantumDynamics):
    def __init__(self, config,GreensFunction, seed: Union[int, None]=None):
        super().__init__(config)
        self.cfg = config
        self.GF = GreensFunction
        self.t_evaluation = self.cfg.tlist * ps_to_au
        self.seed = seed
    def eff_Hamiltonian(self) -> Qobj:
        """
        Construct effective non-Hermitian Hamiltonian by:
        """
        H = self.build_hamiltonian(self.GF,seed = self.seed)

        H_full = H.full().copy()
        
        G = 0.5 * (self.Gamma_ab + self.Gamma_ab.conj().T)
        H_full[1:self.cfg.Nmol+1, 1:self.cfg.Nmol+1] += -0.5j*G

        self.Heff = Qobj(H_full, dims=H.dims)
        return self.Heff
    
    def evolve(self, 
            rho_or_psi, 
            e_ops = None, 
            options = None
            )-> SimulationResult:
        
        psi0 = rho_or_psi
        Heff = self.eff_Hamiltonian()

        result = sesolve(Heff, psi0, self.t_evaluation,
                        e_ops= list(e_ops.values()) if e_ops else None,
                        options = options)
        
        named = {name: np.asarray(result.expect[n]) for n, name in enumerate(e_ops.keys())} if e_ops else {}
        return SimulationResult(tlist=self.cfg.tlist, expectations=named)



