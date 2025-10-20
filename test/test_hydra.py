# tests/test_hydra_wiring.py
from pathlib import Path
from hydra import initialize_config_dir, compose
from mqed.Lindblad.run_quantum_dynamics import app_run

def _compose(cfg_dir: Path, config_name: str, overrides=None):
    overrides = overrides or []
    from omegaconf import OmegaConf
    with initialize_config_dir(config_dir=str(cfg_dir.resolve()), version_base=None):
        cfg = compose(config_name=config_name, overrides=overrides)
    return cfg

def test_config_composes():
    cfg_dir = Path(__file__).parent / "../configs/Lindblad"
    cfg = _compose(cfg_dir, "quantum_dynamics")
    assert "simulation" in cfg
    assert cfg.simulation.Nmol > 0
    assert cfg.solver.method in ("Lindblad", "NonHermitian")

def test_overrides_apply():
    cfg_dir = Path(__file__).parent / "../configs/Lindblad"
    cfg = _compose(
        cfg_dir, "quantum_dynamics",
        overrides=["solver.method=NonHermitian", "simulation.t_ps.stop=12.5", "output.filename=test.h5"]
    )
    assert cfg.solver.method == "NonHermitian"
    assert abs(cfg.simulation.t_ps.stop - 12.5) < 1e-12
    assert cfg.output.filename == "test.h5"

def test_smoke_writes_to_tmp(tmp_path):
    cfg_dir = Path(__file__).parent / "../configs/Lindblad"
    cfg = _compose(cfg_dir, "quantum_dynamics", overrides=[f"output.filename=test_out.h5",
                                                        f"greens.h5_path={(Path(__file__).resolve().parents[1] / 'outputs/Dyadic_GF_analytical/2025-10-13/11-29-07/result_Ag_2_nm.hdf5')}",
                                                        f"simulation.Nmol=5"])
    # IMPORTANT: app_run writes into output_dir we pass
    app_run(cfg, output_dir=tmp_path)
    assert (tmp_path / "test_out.h5").exists()
