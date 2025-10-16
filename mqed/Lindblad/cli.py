# mqed/cli.py
from hydra import initialize, compose
from mqed.Lindblad.run_quantum_dynamics import app_run

def mqed_lindblad() -> None:
    # Adjust config_path to where your QD_default.yaml lives relative to this module
    with initialize(config_path="../configs/Lindblad", version_base=None):
        cfg = compose(config_name="quantum_dynamics", overrides=["solver.method=Lindblad"])
    app_run(cfg)

def mqed_nhse() -> None:
    with initialize(config_path="../configs/Lindblad", version_base=None):
        cfg = compose(config_name="quantum_dynamics", overrides=["solver.method=NonHermitian"])
    app_run(cfg)
