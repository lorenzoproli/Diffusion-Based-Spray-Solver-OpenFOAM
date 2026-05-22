#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kernel bandwidth sensitivity study for Zhang GLOP cases.

Reads CSVs produced by postprocess_zhang_kernel_single.py from each of:
    <root>/b_0_cellPoint/postProcessing_Zhang_OpenFOAM/profile_Zhang_b_0_cellPoint_x060mm.csv
    <root>/b_Dinj_over_16/...
    <root>/b_Dinj_over_8/...
    <root>/b_Dinj/...
    <root>/b_3Dinj/...

where <root> defaults to the parent of this script's directory (Zhang_cases/).

Produces:
    Zhang_kernel_comparison_{count,Umean,D10,D32}_x060mm.png  (all cases overlaid)
    Zhang_kernel_{case_name}_x060mm.png                        (one 2×2 panel per case)
    Zhang_kernel_comparison_SMD_global.png
    Zhang_kernel_comparison_summary.csv

Reference: Zhang (2025) Energy 335, 138192
    - Fig. 24: spatial profiles at Δx = 60 mm (overlaid for a selected Weg)
    - Fig. 23b: global SMD vs Weg (power-law fit + experimental points)

Usage:
    cd Zhang_cases
    python3 postProcessing_KernelStudy/compare_zhang_kernel.py
    python3 postProcessing_KernelStudy/compare_zhang_kernel.py --weg-ref We500
    python3 postProcessing_KernelStudy/compare_zhang_kernel.py --root /path/to/Zhang_cases
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Kernel cases
# =============================================================================
D_INJ_M = 1.6e-3   # jet diameter [m]

CASES = {
    "b_0_PCM": {
        "label":            "PCM (cellPoint), b = 0",
        "smoothBandwidth_m": 0.0,
        "color":            "tab:gray",
        "marker":           "D",
    },
    "b_Dinj": {
        "label":            r"DBM, $b = D_{inj}$",
        "smoothBandwidth_m": D_INJ_M,               # 1.6e-3 m
        "color":            "tab:blue",
        "marker":           "o",
    },
    "b_3Dinj": {
        "label":            r"DBM, $b = 3D_{inj}$",
        "smoothBandwidth_m": 3 * D_INJ_M,           # 4.8e-3 m
        "color":            "tab:orange",
        "marker":           "s",
    },
    "b_6Dinj": {
        "label":            r"DBM, $b = 6D_{inj}$",
        "smoothBandwidth_m": 6 * D_INJ_M,           # 9.6e-3 m
        "color":            "tab:red",
        "marker":           "^",
    },
}

# =============================================================================
# Zhang (2025) reference data
# =============================================================================
# Fig. 24: D32(z), Umean(z), D10(z), count(z) at Δx = 60 mm
ZHANG_REF_PROFILES = {
    "We100": {
        "Weg": 113,
        # No digitized profile data — only metadata for selection
    },
    "We500": {
        "Weg":            508,
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "count_droplets": np.array([1420, 1500, 1580, 1580, 1260, 680, 450, 420, 300, 180, 90, 55]),
        "Umean_m_s":      np.array([20.0, 23.0, 18.0, 13.0, 9.3, 8.2, 7.0, 7.0, 7.0, 7.5, 7.8, 8.5]),
        "D10_um":         np.array([85, 90, 98, 112, 130, 138, 135, 132, 140, 143, 142, 150]),
        "D32_um":         np.array([170, 165, 245, 285, 300, 335, 340, 335, 310, 275, 265, 200]),
    },
    "We3500": {
        "Weg":            3459,
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "count_droplets": np.array([14000, 17500, 20500, 18000, 14000, 12000, 8500, 6500, 4200, 3000, 2500, 1200]),
        "Umean_m_s":      np.array([62, 48, 37, 29, 28, 25, 25, 24, 24, 24, 24, 24]),
        "D10_um":         np.array([64, 68, 75, 79, 82, 84, 85, 87, 88, 89, 90, 94]),
        "D32_um":         np.array([73, 79, 89, 99, 101, 101, 101, 107, 108, 108, 113, 113]),
    },
    "We8500": {
        "Weg": 8624,
        # No digitized profile data published
    },
}

# Fig. 23b: experimental SMD points [Weg, SMD_um] + power-law fit
# SMD = 5825.53 * Weg^(-0.51)  (Zhang 2025 eq. in Fig. 23b)
EXP_SMD = np.array([
    [129,  430],   # Exp 1 (TL 10 %)
    [508,  270],   # Exp 2 (TL 30 %)
    [3136,  92],   # Exp 3 (TL 70 %)
    [6562,  75],   # Exp 4 (TL 100 %)
])
SIM_SMD = np.array([
    [113,  450],   # Sim 1
    [508,  260],   # Sim 2
    [3459,  88],   # Sim 3
    [8624,  72],   # Sim 4
])

def powerlaw_smd(weg):
    return 5825.53 * weg ** (-0.51)


# =============================================================================
# I/O helpers
# =============================================================================
def load_profile(root, case_name, weg_ref="We500", x_offset_mm=60):
    """Load profile CSV for a kernel case under <root>/<weg_ref>/<case_name>/...

    The nested layout Zhang_cases/<We*>/<b_*> is used (one Weg per parent
    directory, multiple kernel bandwidths inside it).
    """
    csv = (root / weg_ref / case_name / "postProcessing_Zhang_OpenFOAM"
           / f"profile_Zhang_{case_name}_x{x_offset_mm:03d}mm.csv")
    if not csv.is_file():
        return None, csv
    try:
        arr = np.genfromtxt(csv, delimiter=",", names=True)
    except Exception as e:
        print(f"[warn] cannot read {csv}: {e}")
        return None, csv
    return arr, csv


def global_smd_from_profile(arr):
    d32 = arr["D32_m"]
    nc  = arr["count_droplets"]
    ok  = np.isfinite(d32) & (d32 > 0) & (nc > 0)
    if not np.any(ok):
        return np.nan
    num = np.sum(nc[ok] * d32[ok] ** 3)
    den = np.sum(nc[ok] * d32[ok] ** 2)
    return (num / den * 1e6) if den > 0 else np.nan


# =============================================================================
# Plotting
# =============================================================================
def plot_profile_quantity(root, out_dir, key, ylabel, fname, ref, weg_ref, log=False):
    fig, ax = plt.subplots(figsize=(9, 5.5))

    if ref is not None and key in ref:
        ax.plot(ref["z_mm"], ref[key], "k-", lw=2.2, label="Zhang (2025) exp/sim",
                zorder=3)

    for case_name, meta in CASES.items():
        data, csv_path = load_profile(root, case_name, weg_ref=weg_ref)
        if data is None:
            print(f"[warn] missing {csv_path}")
            continue
        if key not in data.dtype.names:
            print(f"[warn] {key} not found in {csv_path}")
            continue
        z_mm = data["z_center_mm"]
        y    = data[key]
        ok   = np.isfinite(y) & (y > 0) if "count" not in key else np.isfinite(y)
        ax.plot(z_mm[ok], y[ok],
                color=meta["color"], marker=meta["marker"],
                lw=1.8, ms=5, ls="--",
                label=meta["label"])

    ax.set_xlabel("z [mm]")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, 65)
    if log:
        ax.set_yscale("log")
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(fontsize=8, loc="best")
    ax.set_title(f"Zhang kernel study — {weg_ref} @ $\\Delta x = 60$ mm")
    fig.tight_layout()
    out = out_dir / fname
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[ok ] {out.name}")


def plot_single_case_panel(root, out_dir, case_name, meta, ref, weg_ref):
    """2×2 figure (count, Umean, D10, D32) for one kernel case vs reference."""
    data, csv_path = load_profile(root, case_name, weg_ref=weg_ref)
    if data is None:
        print(f"[warn] missing {csv_path}, skipping single-case plot")
        return

    quantities = [
        ("count_droplets", "Equivalent droplet count [-]", True),
        ("Umean_m_s",      "Mean droplet velocity [m/s]",  False),
        ("D10_um",         r"$D_{10}$ [µm]",               False),
        ("D32_um",         r"$D_{32}$ [µm]",               False),
    ]

    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs = axs.flatten()

    for ax, (key, ylabel, log) in zip(axs, quantities):
        if ref is not None and key in ref:
            ax.plot(ref["z_mm"], ref[key], "k-", lw=2.2, label="Zhang (2025) exp/sim",
                    zorder=3)
        if key in data.dtype.names:
            z_mm = data["z_center_mm"]
            y    = data[key]
            ok   = np.isfinite(y) & (y > 0) if "count" not in key else np.isfinite(y)
            ax.plot(z_mm[ok], y[ok],
                    color=meta["color"], marker=meta["marker"],
                    lw=1.8, ms=5, ls="--", label=meta["label"])
        ax.set_xlabel("z [mm]")
        ax.set_ylabel(ylabel)
        ax.set_xlim(0, 65)
        if log:
            ax.set_yscale("log")
        ax.grid(True, ls=":", alpha=0.5)
        ax.legend(fontsize=8, loc="best")

    fig.suptitle(
        rf"Zhang kernel study — {weg_ref} — {meta['label']} @ $\Delta x = 60$ mm",
        fontsize=12,
    )
    fig.tight_layout()
    out = out_dir / f"Zhang_kernel_{weg_ref}_{case_name}_x060mm.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[ok ] {out.name}")


def plot_smd_global(root, out_dir, weg_ref):
    """Bar chart of global SMD for all kernel cases, with Fig. 23b reference."""
    ref = ZHANG_REF_PROFILES.get(weg_ref, {})
    weg_val = ref.get("Weg", None)

    fig, ax = plt.subplots(figsize=(10, 5))

    smd_values, smd_labels, colors = [], [], []
    for case_name, meta in CASES.items():
        data, _ = load_profile(root, case_name, weg_ref=weg_ref)
        smd = global_smd_from_profile(data) if data is not None else np.nan
        smd_values.append(smd)
        smd_labels.append(meta["label"])
        colors.append(meta["color"])

    x = np.arange(len(smd_labels))
    bars = ax.bar(x, smd_values, color=colors, alpha=0.8, edgecolor="k", linewidth=0.8)

    # Annotate values
    for bar, val in zip(bars, smd_values):
        if np.isfinite(val):
            ax.text(bar.get_x() + bar.get_width() / 2, val + 1,
                    f"{val:.1f}", ha="center", va="bottom", fontsize=8)

    # Reference lines
    if weg_val is not None:
        fit_val = powerlaw_smd(weg_val)
        ax.axhline(fit_val, color="tab:blue", ls="-", lw=1.8,
                   label=f"Power-law fit ({weg_ref}, Weg={weg_val}): {fit_val:.0f} µm")
    # Experimental points at this Weg
    if weg_val is not None:
        for row in EXP_SMD:
            if abs(row[0] - weg_val) / weg_val < 0.3:
                ax.axhline(row[1], color="tab:red", ls="--", lw=1.5,
                           label=f"EXP Zhang (2025) @ Weg≈{row[0]:.0f}: {row[1]:.0f} µm")
                break

    ax.set_xticks(x)
    ax.set_xticklabels(smd_labels, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel(r"Global $D_{32}$ [µm]  (nParticle-weighted)")
    ax.set_title(f"Zhang kernel study — global SMD @ Δx = 60 mm  (ref: {weg_ref})")
    ax.legend(fontsize=8)
    ax.grid(axis="y", ls=":", alpha=0.5)
    fig.tight_layout()
    out = out_dir / f"Zhang_kernel_{weg_ref}_SMD_global.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[ok ] {out.name}")


# =============================================================================
# Summary CSV
# =============================================================================
def write_summary(root, out_dir, weg_ref):
    ref = ZHANG_REF_PROFILES.get(weg_ref, {})
    ref_d32 = ref.get("D32_um")
    ref_z   = ref.get("z_mm")

    rows = []
    for case_name, meta in CASES.items():
        data, csv_path = load_profile(root, case_name, weg_ref=weg_ref)
        if data is None:
            print(f"[warn] missing {csv_path}, skipped in summary")
            continue
        smd_global = global_smd_from_profile(data)

        # RMSE vs reference profile
        rmse_d32 = np.nan
        if ref_d32 is not None:
            y_of = data["D32_um"]
            z_of = data["z_center_mm"]
            ok   = np.isfinite(y_of) & (y_of > 0)
            if np.any(ok):
                interp = np.interp(z_of[ok], ref_z, ref_d32,
                                   left=np.nan, right=np.nan)
                valid  = np.isfinite(interp)
                if np.any(valid):
                    rmse_d32 = float(np.sqrt(np.mean(
                        (y_of[ok][valid] - interp[valid]) ** 2)))

        rows.append({
            "case":               case_name,
            "label":              meta["label"],
            "smoothBandwidth_m":  meta["smoothBandwidth_m"],
            "smoothBandwidth_mm": meta["smoothBandwidth_m"] * 1e3,
            "b_over_Dinj":        meta["smoothBandwidth_m"] / D_INJ_M,
            "SMD_global_um":      round(smd_global, 2) if np.isfinite(smd_global) else "",
            "valid_bins":         int(np.sum(np.isfinite(data["D32_um"]) & (data["D32_um"] > 0))),
            f"RMSE_D32_vs_{weg_ref}_um": round(rmse_d32, 2) if np.isfinite(rmse_d32) else "",
        })

    out_csv = out_dir / f"Zhang_kernel_{weg_ref}_summary.csv"
    if rows:
        keys = list(rows[0].keys())
        with open(out_csv, "w") as f:
            f.write(",".join(keys) + "\n")
            for r in rows:
                f.write(",".join(str(r.get(k, "")) for k in keys) + "\n")
        print(f"[ok ] {out_csv.name}")
    else:
        print("[warn] no data for summary CSV")


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description="Zhang kernel bandwidth comparison")
    p.add_argument("--root",    default=None,
                   help="Zhang_cases root dir (default: parent of this script's dir)")
    p.add_argument("--weg-ref", default="We3500",
                   choices=list(ZHANG_REF_PROFILES.keys()),
                   help="reference Weg condition for profile overlay (default: We3500)")
    p.add_argument("--x-offset-mm", type=int, default=60,
                   help="downstream plane offset in mm (default: 60)")
    return p.parse_args()


def main():
    args    = parse_args()
    here    = Path(__file__).resolve().parent
    root    = Path(args.root).resolve() if args.root else here.parent
    out_dir = here / args.weg_ref
    out_dir.mkdir(parents=True, exist_ok=True)
    ref     = ZHANG_REF_PROFILES.get(args.weg_ref)

    print(f"[info] root      : {root}")
    print(f"[info] weg ref   : {args.weg_ref}")
    print(f"[info] plane     : Δx = {args.x_offset_mm} mm")
    print(f"[info] output    : {out_dir}\n")

    # Profile comparison plots
    for key, ylabel, fname, log in [
        ("count_droplets", "Equivalent droplet count [-]",
         f"Zhang_kernel_{args.weg_ref}_count_x060mm.png",  True),
        ("Umean_m_s",      "Mean droplet velocity [m/s]",
         f"Zhang_kernel_{args.weg_ref}_Umean_x060mm.png",  False),
        ("D10_um",         r"$D_{10}$ [µm]",
         f"Zhang_kernel_{args.weg_ref}_D10_x060mm.png",    False),
        ("D32_um",         r"$D_{32}$ [µm]",
         f"Zhang_kernel_{args.weg_ref}_D32_x060mm.png",    False),
    ]:
        plot_profile_quantity(root, out_dir, key, ylabel, fname, ref,
                              weg_ref=args.weg_ref, log=log)

    # Individual 2×2 panel per kernel case
    for case_name, meta in CASES.items():
        plot_single_case_panel(root, out_dir, case_name, meta, ref,
                               weg_ref=args.weg_ref)

    # Global SMD bar chart
    plot_smd_global(root, out_dir, args.weg_ref)

    # Summary CSV
    write_summary(root, out_dir, args.weg_ref)

    print(f"\n[done] outputs in {out_dir}")


if __name__ == "__main__":
    main()
