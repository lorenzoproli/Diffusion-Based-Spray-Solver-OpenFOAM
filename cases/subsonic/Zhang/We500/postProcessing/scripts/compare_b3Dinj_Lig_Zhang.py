#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Three-curve comparison @ Δx = 60 mm, We500:
    Zhang (2025) Fig. 24 experimental reference  vs
    OpenFOAM b_3Dinj (Madabhushi standard, b = 3 × Dinj)  vs
    OpenFOAM We500_Lig (Madabhushi-Ligament, b = 3 × Dinj)

Reads CSVs already produced by postprocess_zhang_kernel_single.py from:
    Zhang_cases_cluster/We500/b_3Dinj/postProcessing_Zhang_OpenFOAM/profile_Zhang_b_3Dinj_x060mm.csv
    Zhang_cases_cluster/We500_Lig/postProcessing_Zhang_OpenFOAM/profile_Zhang_We500_Lig_x060mm.csv

Writes a single 2x2 figure + four single panels into
    Zhang_cases_cluster/postProcessing_KernelStudy/We500/Zhang_compare_b3Dinj_vs_Lig_x060mm*.png
"""
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT  = ROOT / "postProcessing_KernelStudy" / "We500"
OUT.mkdir(parents=True, exist_ok=True)

# Zhang (2025) Fig. 24 reference, We500 (Weg ≈ 508)
ZHANG_REF = {
    "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
    "count_droplets": np.array([1420, 1500, 1580, 1580, 1260, 680, 450, 420, 300, 180, 90, 55]),
    "Umean_m_s":      np.array([20.0, 23.0, 18.0, 13.0, 9.3, 8.2, 7.0, 7.0, 7.0, 7.5, 7.8, 8.5]),
    "D10_um":         np.array([85, 90, 98, 112, 130, 138, 135, 132, 140, 143, 142, 150]),
    "D32_um":         np.array([170, 165, 245, 285, 300, 335, 340, 335, 310, 275, 265, 200]),
    "SMD_global_um":  270.0,
}

CASES = [
    {
        "csv":   ROOT / "We500" / "b_3Dinj"
                      / "postProcessing_Zhang_OpenFOAM" / "profile_Zhang_b_3Dinj_x060mm.csv",
        "label": r"OpenFOAM Madabhushi, $b=3\,D_{inj}$",
        "color": "tab:blue",
        "marker": "o",
    },
    {
        "csv":   ROOT / "We500_Lig"
                      / "postProcessing_Zhang_OpenFOAM" / "profile_Zhang_We500_Lig_x060mm.csv",
        "label": r"OpenFOAM Madabhushi-Ligament, $b=3\,D_{inj}$",
        "color": "tab:red",
        "marker": "s",
    },
]

QUANTITIES = [
    ("count_droplets", "count [-]", "log"),
    ("Umean_m_s",      r"$\overline{U}_x$ [m/s]", "linear"),
    ("D10_um",         r"$D_{10}$ [µm]", "linear"),
    ("D32_um",         r"$D_{32}$ [µm]", "linear"),
]


def load(case):
    df = pd.read_csv(case["csv"])
    return df


def _plot_one(ax, quantity, ylabel, scale):
    ax.plot(ZHANG_REF["z_mm"], ZHANG_REF[quantity],
            "k^--", label="Zhang (2025) Fig. 24", markersize=6, linewidth=1.2)
    for c in CASES:
        df = load(c)
        ax.plot(df["z_center_mm"].to_numpy(), df[quantity].to_numpy(),
                color=c["color"], marker=c["marker"], linestyle="-",
                label=c["label"], markersize=5, linewidth=1.4)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel(ylabel)
    if scale == "log":
        ax.set_yscale("log")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=8)


def main():
    # Sanity check
    for c in CASES:
        if not c["csv"].is_file():
            sys.exit(f"missing CSV: {c['csv']}")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    titles = ["count", "Umean", "D10", "D32"]
    for ax, (q, ylab, sc), title in zip(axes.flat, QUANTITIES, titles):
        _plot_one(ax, q, ylab, sc)
        ax.set_title(f"We500, Δx = 60 mm — {title}")
    fig.suptitle("Zhang (2025) vs OpenFOAM Madabhushi vs Madabhushi-Ligament  "
                 r"($b = 3\,D_{inj}$)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_combined = OUT / "Zhang_compare_b3Dinj_vs_Lig_x060mm.png"
    fig.savefig(out_combined, dpi=150)
    plt.close(fig)
    print(f"[ok ] {out_combined.name}")

    for (q, ylab, sc), title in zip(QUANTITIES, titles):
        fig, ax = plt.subplots(figsize=(7, 5))
        _plot_one(ax, q, ylab, sc)
        ax.set_title(f"We500, Δx = 60 mm — {title}\n"
                     r"Zhang (2025) vs OpenFOAM $b=3D_{inj}$ vs Lig")
        fig.tight_layout()
        out_single = OUT / f"Zhang_compare_b3Dinj_vs_Lig_x060mm_{title}.png"
        fig.savefig(out_single, dpi=150)
        plt.close(fig)
        print(f"[ok ] {out_single.name}")

    print(f"\n[done] outputs in {OUT}")


if __name__ == "__main__":
    main()
