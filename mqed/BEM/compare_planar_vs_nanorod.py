from __future__ import annotations
from pathlib import Path

import hydra
import numpy as np
from omegaconf import DictConfig
from matplotlib import pyplot as plt
from hydra.core.hydra_config import HydraConfig
from loguru import logger

from mqed.utils.orientation import spherical_to_cartesian_dipole
from mqed.utils.dgf_data import load_gf_h5
from mqed.utils.logging_utils import setup_loggers_hydra_aware


def _apply_rcparams(rcparams: dict | None) -> None:
    plt.rcParams.update(rcparams or {})


def _apply_axes_style(ax, cfg) -> None:
    ax.set_xlabel(cfg.xlabel)
    ax.set_ylabel(cfg.ylabel)

    if cfg.get("title"):
        ax.set_title(cfg.title)

    if cfg.get("xlim") is not None:
        ax.set_xlim(cfg.xlim[0], cfg.xlim[1])
    if cfg.get("ylim") is not None:
        ax.set_ylim(cfg.ylim[0], cfg.ylim[1])

    ax.tick_params(
        direction=cfg.ticks.direction,
        top=cfg.ticks.top,
        right=cfg.ticks.right,
        labelsize=cfg.ticks.labelsize,
        length=cfg.ticks.length,
        width=cfg.ticks.width,
    )

    if cfg.get("grid", False):
        ax.grid(True, alpha=cfg.get("grid_alpha", 0.25))

    for spine in ax.spines.values():
        spine.set_linewidth(cfg.get("spine_width", 1.2))


def _clip_xy(x: np.ndarray, y: np.ndarray, xlim) -> tuple[np.ndarray, np.ndarray]:
    """Clip x,y to xlim = [xmin, xmax] if provided."""
    if xlim is None:
        return x, y
    xmin, xmax = float(xlim[0]), float(xlim[1])
    m = (x >= xmin) & (x <= xmax)
    return x[m], y[m]


def _drop_nonfinite(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Remove NaN/Inf pairs to avoid matplotlib errors and ugly spikes."""
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def _plot_series(ax, x: np.ndarray, y: np.ndarray, s) -> None:
    ax.plot(
        x, y,
        label=s.label,
        linestyle=s.get("linestyle", "-"),
        linewidth=s.get("linewidth", 2.5),
        marker=s.get("marker", None),
        markersize=s.get("markersize", 7),
        markerfacecolor=s.get("markerfacecolor", "none"),
        markeredgewidth=s.get("markeredgewidth", 1.5),
        color=s.get("color", None),
    )


def _compute_enhancement(
    h5_path: Path,
    x_key: str,
    energy_index: int,
    p_donor: np.ndarray,
    p_acc: np.ndarray,
) -> dict[str, np.ndarray]:
    """
    Returns dict with:
      x_nm, enh_real, enh_imag
    """
    data = load_gf_h5(h5_path)
    Gtot = np.asarray(data["G_total"])  # (M,N,3,3)
    Gvac = np.asarray(data["G_vac"])    # (M,N,3,3)
    x_nm = np.asarray(data[x_key]).ravel()

    m = int(energy_index)

    g_vac = np.einsum("i,...ij,j->...", p_acc, Gvac[m], p_donor)
    g_tot = np.einsum("i,...ij,j->...", p_acc, Gtot[m], p_donor)

    enh_real = np.real(g_tot) / np.real(g_vac)
    enh_imag = np.imag(g_tot) / np.imag(g_vac)

    return {
        "x_nm": x_nm,
        "enh_real": enh_real,
        "enh_imag": enh_imag,
    }


@hydra.main(config_path="../../configs/BEM", config_name="compare_planar_vs_nanorod", version_base=None)
def main(cfg: DictConfig) -> None:
    output_dir = Path(HydraConfig.get().runtime.output_dir)
    setup_loggers_hydra_aware()
    _apply_rcparams(cfg.plot.get("rcParams", None))

    planar_path = Path(cfg.paths.planar_h5)
    nanorod_path = Path(cfg.paths.nanorod_h5)

    for p in (planar_path, nanorod_path):
        if not p.exists():
            raise FileNotFoundError(f"Missing input file: {p}")

    # dipoles
    p_donor = spherical_to_cartesian_dipole(cfg.dipoles.donor.theta_deg, cfg.dipoles.donor.phi_deg)
    p_acc   = spherical_to_cartesian_dipole(cfg.dipoles.acceptor.theta_deg, cfg.dipoles.acceptor.phi_deg)

    # compute enhancements
    planar = _compute_enhancement(planar_path, cfg.data.x_key, cfg.data.energy_index, p_donor, p_acc)
    rod    = _compute_enhancement(nanorod_path, cfg.data.x_key, cfg.data.energy_index, p_donor, p_acc)

    series_data = {
        "planar_real": (planar["x_nm"], planar["enh_real"]),
        "planar_imag": (planar["x_nm"], planar["enh_imag"]),
        "nanorod_real": (rod["x_nm"], rod["enh_real"]),
        "nanorod_imag": (rod["x_nm"], rod["enh_imag"]),
    }

    # plot
    fig, ax = plt.subplots(figsize=tuple(cfg.plot.figsize), dpi=cfg.plot.get("dpi", 120))
    _apply_axes_style(ax, cfg.plot.axes)

    # title formatting
    title = cfg.plot.axes.get("title")
    if title:
        ax.set_title(title.format(dipole_frequency=cfg.dipole_frequency))

    # y scale
    yscale = cfg.plot.axes.get("yscale", "log")
    if yscale == "symlog":
        ax.set_yscale("symlog", linthresh=cfg.plot.axes.get("linthresh", 1e-3))
    else:
        ax.set_yscale(yscale)

    # draw series in YAML order; clip to xlim per-series to avoid shape mismatch
    xlim = cfg.plot.axes.get("xlim", None)

    for key in cfg.plot.series_order:
        s = cfg.plot.series[key]
        x, y = series_data[s.source]

        x, y = _clip_xy(x, y, xlim)
        x, y = _drop_nonfinite(x, y)

        _plot_series(ax, x, y, s)
        logger.info(f"Plotted {s.label} (source={s.source}, n={len(x)})")

    # legend
    leg = ax.legend(
        loc=cfg.plot.legend.loc,
        bbox_to_anchor=tuple(cfg.plot.legend.bbox_to_anchor) if "bbox_to_anchor" in cfg.plot.legend else None,
        frameon=cfg.plot.legend.frameon,
        fancybox=cfg.plot.legend.fancybox,
        framealpha=cfg.plot.legend.framealpha,
        fontsize=cfg.plot.legend.fontsize,
    )
    leg.get_frame().set_edgecolor(cfg.plot.legend.get("edgecolor", "0.4"))

    fig.tight_layout()

    if cfg.plot.get("save", False):
        plot_filename = Path(cfg.plot.save_path)
        plot_filepath = output_dir / plot_filename
        fig.savefig(plot_filepath, bbox_inches="tight")
        logger.success(f"Plot saved to {plot_filepath}")

    plt.close(fig)
    logger.success(f"Planar vs Nanorod GF comparison complete. Logs saved to: {plot_filepath.absolute()}")


if __name__ == "__main__":
    main()
