#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zhang-style post-processing for a single OpenFOAM kernel-bandwidth case.

Identical methodology to postprocess_zhang_single.py (WeComparison folder):
  - Snapshot-based Lagrangian statistics, nParticle-weighted, no flux weighting
  - Sectional statistics at Δx = 60 mm from jet orifice, binned in z
  - Min-droplets threshold per bin (default 25 equivalent droplets)

Usage:
    python3 postprocess_zhang_kernel_single.py \\
        --case ../b_Dinj \\
        --case-label b_Dinj \\
        --weg-ref We3500

Outputs inside <case>/postProcessing_Zhang_OpenFOAM/:
    profile_Zhang_b_Dinj_x060mm.csv
    profile_Zhang_b_Dinj_x060mm_{count,Umean,D10,D32}.png
    metadata_b_Dinj.json
"""

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Reference data — Zhang (2025) Energy 335, 138192 — Fig. 24, Δx = 60 mm
# =============================================================================
D_INJ_M = 1.6e-3   # jet diameter [m]

ZHANG_REF = {
    "We100": {
        # No digitized Fig. 24 profile data inserted yet; Fig. 23b metadata only
        "SMD_global_um": 430.0,
        "Weg": 113,
    },
    "We500": {
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "count_droplets": np.array([1420, 1500, 1580, 1580, 1260, 680, 450, 420, 300, 180, 90, 55]),
        "Umean_m_s":      np.array([20.0, 23.0, 18.0, 13.0, 9.3, 8.2, 7.0, 7.0, 7.0, 7.5, 7.8, 8.5]),
        "D10_um":         np.array([85, 90, 98, 112, 130, 138, 135, 132, 140, 143, 142, 150]),
        "D32_um":         np.array([170, 165, 245, 285, 300, 335, 340, 335, 310, 275, 265, 200]),
        # Fig. 23b: experimental SMD at 60 mm, Weg = 508
        "SMD_global_um":  270.0,
        "Weg":            508,
    },
    "We3500": {
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "count_droplets": np.array([14000, 17500, 20500, 18000, 14000, 12000, 8500, 6500, 4200, 3000, 2500, 1200]),
        "Umean_m_s":      np.array([62, 48, 37, 29, 28, 25, 25, 24, 24, 24, 24, 24]),
        "D10_um":         np.array([64, 68, 75, 79, 82, 84, 85, 87, 88, 89, 90, 94]),
        "D32_um":         np.array([73, 79, 89, 99, 101, 101, 101, 107, 108, 108, 113, 113]),
        # Fig. 23b: experimental SMD at 60 mm, Weg = 3136
        "SMD_global_um":  92.0,
        "Weg":            3459,
    },
    "We8500": {
        # No published digitized Fig. 24 data; Fig. 23b has one exp point
        "SMD_global_um":  75.0,
        "Weg":            8624,
    },
}

# =============================================================================
# Lagrangian readers (same as postprocess_zhang_single.py)
# =============================================================================
def _read_text(path):
    return Path(path).read_text(encoding="latin-1", errors="ignore")


def read_positions(path):
    txt = _read_text(path)
    triples = re.findall(
        r"\(\s*"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*\)",
        txt,
    )
    if not triples:
        return np.empty((0, 3), dtype=float)
    return np.array([[float(a), float(b), float(c)] for a, b, c in triples], dtype=float)


def read_scalar(path):
    txt = _read_text(path)
    m = re.search(r"\n\s*(\d+)\s*\n\s*\(\s*\n", txt)
    if not m:
        m = re.search(r"\n\s*(\d+)\s*\(\s*\n", txt)
    if not m:
        raise RuntimeError(f"Cannot find scalar list start in {path}")
    N = int(m.group(1))
    start = m.end()
    end = txt.find("\n)", start)
    if end == -1:
        end = txt.rfind(")")
    block = txt[start:end]
    vals = [float(l.strip()) for l in block.splitlines() if l.strip()]
    vals = np.array(vals, dtype=float)
    if len(vals) != N:
        raise RuntimeError(f"Scalar mismatch in {path}: expected {N}, got {len(vals)}")
    return vals


def read_vector(path):
    return read_positions(path)


# =============================================================================
# Time directory enumeration
# =============================================================================
def _fmt_time(t):
    s = f"{t:.10f}".rstrip("0").rstrip(".")
    return s if s else "0"


def build_time_list(t_start, t_end, t_step, case_dir):
    times = []
    n = int(round((t_end - t_start) / t_step)) + 1
    for k in range(n):
        t = t_start + k * t_step
        for nm in [_fmt_time(t), f"{t:.6g}"]:
            if (case_dir / nm).is_dir():
                times.append(nm)
                break
    return times


def load_ensemble(case_dir, cloud, time_names):
    all_pos, all_d, all_np, all_U, used_times = [], [], [], [], []
    for t in time_names:
        base = case_dir / t / "lagrangian" / cloud
        pth_pos = next((c for c in [base / "positions", base / "position0"] if c.is_file()), None)
        pth_d   = base / "d"
        pth_np  = base / "nParticle"
        pth_U   = base / "U"
        if pth_pos is None or not pth_d.is_file() or not pth_np.is_file():
            print(f"[skip] time {t}: missing positions/d/nParticle")
            continue
        pos   = read_positions(pth_pos)
        d     = read_scalar(pth_d)
        npart = read_scalar(pth_np)
        U     = read_vector(pth_U) if pth_U.is_file() else None
        if not (len(pos) == len(d) == len(npart)):
            print(f"[skip] time {t}: size mismatch")
            continue
        if U is not None and len(U) != len(d):
            U = None
        all_pos.append(pos); all_d.append(d); all_np.append(npart); all_U.append(U)
        used_times.append(t)
    if not used_times:
        raise RuntimeError("No valid time steps loaded.")
    pos   = np.vstack(all_pos)
    d     = np.hstack(all_d)
    npart = np.hstack(all_np)
    U     = np.vstack(all_U) if all(u is not None for u in all_U) else None
    return pos, d, npart, U, used_times


# =============================================================================
# Sectional statistics
# =============================================================================
def _wmean(v, w):
    den = np.sum(w)
    return np.sum(w * v) / den if den > 0 else np.nan


def _wD32(d, w):
    num = np.sum(w * d ** 3)
    den = np.sum(w * d ** 2)
    return num / den if den > 0 else np.nan


def compute_sectional(pos, d, npart, U, x_plane, slab_half, z_bins, min_droplets):
    mask  = np.abs(pos[:, 0] - x_plane) <= slab_half
    posp  = pos[mask]; dp = d[mask]; np_ = npart[mask]
    Up    = U[mask] if U is not None else None
    z     = posp[:, 2]
    centers = 0.5 * (z_bins[:-1] + z_bins[1:])
    out = {
        "z_center_m":     centers,
        "z_center_mm":    centers * 1e3,
        "count_parcels":  np.zeros(len(centers), dtype=int),
        "count_droplets": np.zeros(len(centers), dtype=float),
        "Umean_m_s":      np.full(len(centers), np.nan),
        "D10_m":          np.full(len(centers), np.nan),
        "D10_um":         np.full(len(centers), np.nan),
        "D32_m":          np.full(len(centers), np.nan),
        "D32_um":         np.full(len(centers), np.nan),
    }
    for i, (z0, z1) in enumerate(zip(z_bins[:-1], z_bins[1:])):
        m = (z >= z0) & (z < z1)
        di = dp[m]; ni = np_[m]
        out["count_parcels"][i]  = int(di.size)
        ndrop                    = float(np.sum(ni))
        out["count_droplets"][i] = ndrop
        if ndrop < min_droplets:
            continue
        out["D10_m"][i]  = D10 = _wmean(di, ni)
        out["D10_um"][i] = D10 * 1e6
        D32 = _wD32(di, ni)
        out["D32_m"][i]  = D32
        out["D32_um"][i] = D32 * 1e6
        if Up is not None:
            ui = np.linalg.norm(Up[m], axis=1)
            out["Umean_m_s"][i] = _wmean(ui, ni)
    return out


def global_smd(out):
    d_vals = out["D32_m"]
    n_vals = out["count_droplets"]
    mask   = np.isfinite(d_vals) & (d_vals > 0) & (n_vals > 0)
    if not np.any(mask):
        return np.nan
    # volume-weighted mean of D32 bins
    num = np.sum(n_vals[mask] * d_vals[mask] ** 3)
    den = np.sum(n_vals[mask] * d_vals[mask] ** 2)
    return (num / den * 1e6) if den > 0 else np.nan


def save_csv(stats, csv_path):
    keys = ["z_center_m", "z_center_mm", "count_parcels", "count_droplets",
            "Umean_m_s", "D10_m", "D10_um", "D32_m", "D32_um"]
    arr    = np.column_stack([stats[k] for k in keys])
    header = ",".join(keys)
    np.savetxt(csv_path, arr, delimiter=",", header=header, comments="")


def plot_quantities(stats, title, out_prefix, ref=None):
    z_mm = stats["z_center_mm"]
    quantities = [
        ("count_droplets", "Equivalent droplet count [-]",  "_count.png",
         ref["count_droplets"] if ref and "count_droplets" in ref else None),
        ("Umean_m_s",      r"Mean droplet velocity [m/s]",  "_Umean.png",
         ref["Umean_m_s"]      if ref and "Umean_m_s"      in ref else None),
        ("D10_um",         r"$D_{10}$ [µm]",                "_D10.png",
         ref["D10_um"]         if ref and "D10_um"         in ref else None),
        ("D32_um",         r"$D_{32}$ [µm]",                "_D32.png",
         ref["D32_um"]         if ref and "D32_um"         in ref else None),
    ]
    for key, ylabel, suffix, refvals in quantities:
        plt.figure(figsize=(8, 5))
        y  = stats[key]
        ok = np.isfinite(y) & (y > 0) if "count" not in key else np.ones(len(y), bool)
        plt.plot(z_mm[ok], y[ok], "ko-", lw=1.8, ms=5, label="OpenFOAM")
        if refvals is not None:
            plt.plot(ref["z_mm"], refvals, "rs--", lw=1.2, ms=5,
                     fillstyle="none", label="Zhang (2025) ref")
        plt.xlabel("z [mm]"); plt.ylabel(ylabel)
        plt.title(title); plt.grid(True); plt.legend()
        plt.tight_layout()
        plt.savefig(str(out_prefix) + suffix, dpi=200)
        plt.close()


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description="Zhang kernel-study single-case postprocessing")
    p.add_argument("--case",       required=True,
                   help="case directory (relative or absolute)")
    p.add_argument("--case-label", required=True,
                   help="label for output files, e.g. b_Dinj")
    p.add_argument("--weg-ref",    default="We3500",
                   choices=list(ZHANG_REF.keys()),
                   help="which Zhang (2025) reference to overlay (default: We3500)")
    p.add_argument("--x-offsets",  nargs="+", type=float, default=[0.060],
                   help="downstream offsets from x-jet [m] (default: 0.060)")
    p.add_argument("--times-start", type=float, default=0.02)
    p.add_argument("--times-end",   type=float, default=0.03)
    p.add_argument("--times-step",  type=float, default=0.0005)
    p.add_argument("--cloud",       default="sprayCloud")
    p.add_argument("--x-jet",       type=float, default=0.015,
                   help="x-coordinate of jet orifice exit [m] (default: 0.015)")
    p.add_argument("--slice-half-thickness", type=float, default=5e-4)
    p.add_argument("--min-droplets-per-bin", type=float, default=25.0)
    p.add_argument("--nz-stats",    type=int,   default=13)
    p.add_argument("--z-max",       type=float, default=0.065)
    p.add_argument("--z-min",       type=float, default=0.0)
    return p.parse_args()


def resolve_case(arg):
    p = Path(arg)
    if p.is_absolute() and p.is_dir():
        return p.resolve()
    if p.is_dir():
        return p.resolve()
    # Try relative to Zhang_cases (grandparent of this script)
    root = Path(__file__).resolve().parent.parent
    cand = root / arg
    if cand.is_dir():
        return cand.resolve()
    raise SystemExit(f"Case directory not found: {arg}")


def main():
    args     = parse_args()
    case_dir = resolve_case(args.case)
    out_root = case_dir / "postProcessing_Zhang_OpenFOAM"
    out_root.mkdir(parents=True, exist_ok=True)

    times = build_time_list(args.times_start, args.times_end, args.times_step, case_dir)
    if not times:
        raise SystemExit(f"No time directories found in {case_dir} "
                         f"between {args.times_start} and {args.times_end}")
    print(f"[info] case   : {case_dir.name}")
    print(f"[info] times  : {len(times)} ({times[0]} … {times[-1]})")

    pos, d, npart, U, used_times = load_ensemble(case_dir, args.cloud, times)
    print(f"[info] parcels: {len(d):,}")

    z_bins = np.linspace(args.z_min, args.z_max, args.nz_stats + 1)
    ref    = ZHANG_REF.get(args.weg_ref)

    summary = {
        "case_dir":               str(case_dir),
        "case_label":             args.case_label,
        "weg_ref":                args.weg_ref,
        "cloud":                  args.cloud,
        "x_jet_m":                args.x_jet,
        "x_offsets_m":            args.x_offsets,
        "slice_half_thickness_m": args.slice_half_thickness,
        "min_droplets_per_bin":   args.min_droplets_per_bin,
        "times_used":             used_times,
        "z_bins":                 z_bins.tolist(),
        "weighting":              "nParticle (snapshot, no flux)",
    }

    for dx in args.x_offsets:
        x_plane = args.x_jet + dx
        stats   = compute_sectional(
            pos, d, npart, U, x_plane,
            args.slice_half_thickness, z_bins, args.min_droplets_per_bin,
        )
        smd     = global_smd(stats)
        tag     = f"x{int(round(dx * 1e3)):03d}mm"
        lbl     = args.case_label
        csv_path = out_root / f"profile_Zhang_{lbl}_{tag}.csv"
        save_csv(stats, csv_path)

        ref_used = ref if abs(dx - 0.060) < 1e-9 else None
        plot_quantities(
            stats,
            title=f"{lbl}  Δx = {dx*1e3:.0f} mm  (ref: {args.weg_ref})",
            out_prefix=out_root / f"profile_Zhang_{lbl}_{tag}",
            ref=ref_used,
        )
        print(f"[ok ] {csv_path.name}  |  global SMD = {smd:.1f} µm" if np.isfinite(smd)
              else f"[ok ] {csv_path.name}")

        if np.isfinite(smd):
            summary[f"SMD_global_um_{tag}"] = round(smd, 2)

    (out_root / f"metadata_{args.case_label}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"[done] outputs in {out_root}")


if __name__ == "__main__":
    main()
