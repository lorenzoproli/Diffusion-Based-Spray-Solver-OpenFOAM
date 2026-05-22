#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Zhang-style post-processing for a single OpenFOAM Lagrangian case.

Methodology aligned with:
  Zhang et al. (2025), Energy 335, 138192     [throttling levels paper]
  Zhang et al. (2026), Aerospace Sci. & Tech. 168, 111179 [LMR paper]

Statistics at one or more downstream slabs |x - (x_jet + dx)| <= dx/2,
binned in z, weighted by nParticle (snapshot-based, NOT flux-weighted).

Outputs (per case, per x-offset):
  <case>/postProcessing_Zhang_OpenFOAM/profile_Zhang_<we_label>_xNNNmm.csv
  ... + per-quantity PNGs (count, Umean, D10, D32)
"""

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# Optional digitized references from Zhang (2025) Fig. 24 (Δx = 60 mm)
# =============================================================================
ZHANG_REF = {
    "We100": {
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40]),
        "count_droplets": np.array([200, 260, 380, 440, 400, 330, 280, 270, 220]),
        "U":              np.array([9.0, 11.5, 8.0, 6.0, 4.5, 3.5, 3.0, 2.7, 2.5]),
        "D10":            np.array([140, 148, 145, 155, 190, 215, 250, 265, 230]),
        "D32":            np.array([410, 445, 580, 555, 520, 535, 660, 640, 575]),
    },
    "We500": {
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "count_droplets": np.array([1420, 1500, 1580, 1580, 1260, 680, 450, 420, 300, 180, 90, 55]),
        "U":              np.array([20.0, 23.0, 18.0, 13.0, 9.3, 8.2, 7.0, 7.0, 7.0, 7.5, 7.8, 8.5]),
        "D10":            np.array([85, 90, 98, 112, 130, 138, 135, 132, 140, 143, 142, 150]),
        "D32":            np.array([170, 165, 245, 285, 300, 335, 340, 335, 310, 275, 265, 200]),
    },
    "We3500": {
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "count_droplets": np.array([14000, 17500, 20500, 18000, 14000, 12000, 8500, 6500, 4200, 3000, 2500, 1200]),
        "U":              np.array([62, 48, 37, 29, 28, 25, 25, 24, 24, 24, 24, 24]),
        "D10":            np.array([64, 68, 75, 79, 82, 84, 85, 87, 88, 89, 90, 94]),
        "D32":            np.array([73, 79, 89, 99, 101, 101, 101, 107, 108, 108, 113, 113]),
    },
    "We8500": {
        "z_mm":           np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "count_droplets": np.array([11500, 17500, 18500, 16000, 12500, 10000, 7500, 6200, 4500, 3500, 2700, 1800]),
        "U":              np.array([90, 95, 75, 55, 47, 45, 44, 44, 44, 44, 44, 44]),
        "D10":            np.array([50, 55, 60, 67, 70, 72, 74, 75, 76, 76, 76, 76]),
        "D32":            np.array([60, 60, 65, 80, 86, 88, 89, 89, 89, 89, 90, 88]),
    },
}


# =============================================================================
# OpenFOAM lagrangian readers (robust to cppcomment / scientific notation)
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


def read_vector(path):
    return read_positions(path)


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
    vals = []
    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            vals.append(float(line))
        except ValueError:
            pass
    vals = np.array(vals, dtype=float)
    if len(vals) != N:
        raise RuntimeError(f"Scalar length mismatch in {path}: expected {N}, got {len(vals)}")
    return vals


# =============================================================================
# Time list
# =============================================================================
def _fmt_time(t):
    """Match OpenFOAM directory naming: trailing zeros stripped."""
    s = f"{t:.10f}".rstrip("0").rstrip(".")
    return s if s else "0"


def build_time_list(t_start, t_end, t_step, case_dir):
    """Build the list of OpenFOAM time directory names that actually exist."""
    times = []
    n = int(round((t_end - t_start) / t_step)) + 1
    for k in range(n):
        t = t_start + k * t_step
        nm = _fmt_time(t)
        if (case_dir / nm).is_dir():
            times.append(nm)
        else:
            # try alternative formatting (e.g. 0.0250 vs 0.025)
            alt = f"{t:.6g}"
            if (case_dir / alt).is_dir():
                times.append(alt)
    return times


def load_ensemble(case_dir, cloud, time_names):
    all_pos, all_d, all_np, all_U, used_times = [], [], [], [], []
    for t in time_names:
        base = case_dir / t / "lagrangian" / cloud
        pos_candidates = [base / "positions", base / "position0"]
        pth_pos = next((c for c in pos_candidates if c.is_file()), None)
        pth_d  = base / "d"
        pth_np = base / "nParticle"
        pth_U  = base / "U"
        if pth_pos is None or not pth_d.is_file() or not pth_np.is_file():
            print(f"[skip] time {t}: missing positions/d/nParticle")
            continue
        pos   = read_positions(pth_pos)
        d     = read_scalar(pth_d)
        npart = read_scalar(pth_np)
        U     = read_vector(pth_U) if pth_U.is_file() else None
        if not (len(pos) == len(d) == len(npart)):
            print(f"[skip] time {t}: size mismatch pos={len(pos)} d={len(d)} nP={len(npart)}")
            continue
        if U is not None and len(U) != len(d):
            U = None
        all_pos.append(pos); all_d.append(d); all_np.append(npart); all_U.append(U)
        used_times.append(t)
    if not used_times:
        raise RuntimeError("No valid times loaded.")
    pos   = np.vstack(all_pos)
    d     = np.hstack(all_d)
    npart = np.hstack(all_np)
    U = np.vstack(all_U) if all(u is not None for u in all_U) else None
    return pos, d, npart, U, used_times


# =============================================================================
# Sectional statistics (snapshot, nParticle weighted, no flux)
# =============================================================================
def _wmean(v, w):
    den = np.sum(w)
    return np.sum(w * v) / den if den > 0 else np.nan


def _wD32(d, w):
    num = np.sum(w * d ** 3)
    den = np.sum(w * d ** 2)
    return num / den if den > 0 else np.nan


def compute_sectional(pos, d, npart, U, x_plane, slab_half, z_bins, min_droplets,
                      use_velocity_magnitude=True):
    mask = np.abs(pos[:, 0] - x_plane) <= slab_half
    posp = pos[mask]; dp = d[mask]; np_ = npart[mask]
    Up   = U[mask] if U is not None else None
    z    = posp[:, 2]
    centers = 0.5 * (z_bins[:-1] + z_bins[1:])
    out = {
        "z_center_m":     centers,
        "z_center_mm":    centers * 1e3,
        "count_parcels":  np.zeros_like(centers, dtype=int),
        "count_droplets": np.zeros_like(centers, dtype=float),
        "Umean_m_s":      np.full_like(centers, np.nan, dtype=float),
        "D10_m":          np.full_like(centers, np.nan, dtype=float),
        "D10_um":         np.full_like(centers, np.nan, dtype=float),
        "D32_m":          np.full_like(centers, np.nan, dtype=float),
        "D32_um":         np.full_like(centers, np.nan, dtype=float),
    }
    for i, (z0, z1) in enumerate(zip(z_bins[:-1], z_bins[1:])):
        m = (z >= z0) & (z < z1)
        di = dp[m]; ni = np_[m]
        out["count_parcels"][i]  = int(di.size)
        ndrop                    = float(np.sum(ni))
        out["count_droplets"][i] = ndrop
        if ndrop < min_droplets:
            continue
        D10 = _wmean(di, ni)
        D32 = _wD32(di, ni)
        out["D10_m"][i]  = D10; out["D10_um"][i] = D10 * 1e6
        out["D32_m"][i]  = D32; out["D32_um"][i] = D32 * 1e6
        if Up is not None:
            ui = np.linalg.norm(Up[m], axis=1) if use_velocity_magnitude else Up[m][:, 0]
            out["Umean_m_s"][i] = _wmean(ui, ni)
    return out


def save_csv(stats, csv_path):
    keys = [
        "z_center_m", "z_center_mm",
        "count_parcels", "count_droplets",
        "Umean_m_s",
        "D10_m", "D10_um",
        "D32_m", "D32_um",
    ]
    arr = np.column_stack([stats[k] for k in keys])
    header = ",".join(keys)
    np.savetxt(csv_path, arr, delimiter=",", header=header, comments="")


def plot_quantities(stats, title_prefix, out_prefix, ref=None):
    z_mm = stats["z_center_mm"]

    quantities = [
        ("count_droplets", "Equivalent droplet count [-]",
         "_count.png", ref["count_droplets"] if ref else None),
        ("Umean_m_s",      r"Mean droplet velocity [m/s]",
         "_Umean.png",      ref["U"] if ref else None),
        ("D10_um",         r"$D_{10}$ [µm]",
         "_D10.png",         ref["D10"] if ref else None),
        ("D32_um",         r"$D_{32}$ [µm]",
         "_D32.png",         ref["D32"] if ref else None),
    ]
    for key, ylabel, suffix, refvals in quantities:
        plt.figure(figsize=(8, 5))
        y = stats[key]
        m = np.isfinite(y) & (y != 0) if "count" not in key else np.ones_like(y, dtype=bool)
        plt.plot(z_mm[m], y[m], "ko-", lw=1.8, ms=5, label=f"OpenFOAM {title_prefix}")
        if refvals is not None and ref is not None:
            plt.plot(ref["z_mm"], refvals, "rs--", lw=1.2, ms=5,
                     fillstyle="none", label="Zhang (2025) ref")
        plt.xlabel("z [mm]")
        plt.ylabel(ylabel)
        plt.title(f"{title_prefix}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.savefig(str(out_prefix) + suffix, dpi=200)
        plt.close()


# =============================================================================
# CLI
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(description="Zhang-style post-processing (single OpenFOAM case)")
    p.add_argument("--case",   required=True,
                   help="case directory (relative or absolute)")
    p.add_argument("--we-label", required=True,
                   help="label used in output filenames, e.g. We500")
    p.add_argument("--x-offsets", nargs="+", type=float,
                   default=[0.060],
                   help="downstream offsets from x-jet (m), e.g. 0.030 0.045 0.060")
    p.add_argument("--times-start", type=float, default=0.02)
    p.add_argument("--times-end",   type=float, default=0.03)
    p.add_argument("--times-step",  type=float, default=0.0005)
    p.add_argument("--cloud", default="sprayCloud")
    p.add_argument("--x-jet", type=float, default=0.015)
    p.add_argument("--slice-half-thickness", type=float, default=5e-4)
    p.add_argument("--min-droplets-per-bin", type=float, default=25.0)
    p.add_argument("--nz-stats", type=int, default=13)
    p.add_argument("--z-max",  type=float, default=0.065)
    p.add_argument("--z-min",  type=float, default=0.0)
    return p.parse_args()


def resolve_case_dir(arg):
    p = Path(arg)
    if p.is_absolute() and p.is_dir():
        return p.resolve()
    # try relative to cwd
    if p.is_dir():
        return p.resolve()
    # try relative to Zhang_cases (parent of this script)
    here = Path(__file__).resolve().parent.parent
    cand = here / arg
    if cand.is_dir():
        return cand.resolve()
    raise SystemExit(f"Case directory not found: {arg}")


def main():
    args = parse_args()
    case_dir = resolve_case_dir(args.case)
    out_root = case_dir / "postProcessing_Zhang_OpenFOAM"
    out_root.mkdir(parents=True, exist_ok=True)

    # Build time list
    times = build_time_list(args.times_start, args.times_end, args.times_step, case_dir)
    if not times:
        raise SystemExit(f"No time directories found in {case_dir} between "
                         f"{args.times_start} and {args.times_end}")
    print(f"[info] case   : {case_dir}")
    print(f"[info] cloud  : {args.cloud}")
    print(f"[info] times  : {len(times)} (first={times[0]} last={times[-1]})")

    pos, d, npart, U, used_times = load_ensemble(case_dir, args.cloud, times)
    print(f"[info] merged parcels: {len(d)}")

    z_bins = np.linspace(args.z_min, args.z_max, args.nz_stats + 1)

    # Reference data only at Δx = 60 mm
    ref = ZHANG_REF.get(args.we_label)

    summary = {
        "case_dir": str(case_dir),
        "we_label": args.we_label,
        "cloud": args.cloud,
        "x_jet_m": args.x_jet,
        "x_offsets_m": list(args.x_offsets),
        "slice_half_thickness_m": args.slice_half_thickness,
        "min_droplets_per_bin": args.min_droplets_per_bin,
        "times_used": used_times,
        "z_bins": z_bins.tolist(),
        "weighting_method": "nParticle (snapshot, no flux)",
    }

    for dx in args.x_offsets:
        x_plane = args.x_jet + dx
        stats = compute_sectional(
            pos, d, npart, U, x_plane,
            args.slice_half_thickness, z_bins,
            args.min_droplets_per_bin,
        )
        tag = f"x{int(round(dx * 1e3)):03d}mm"
        csv_path = out_root / f"profile_Zhang_{args.we_label}_{tag}.csv"
        save_csv(stats, csv_path)
        plot_prefix = out_root / f"profile_Zhang_{args.we_label}_{tag}"
        # only attach reference at Δx = 60 mm to avoid mis-comparing planes
        ref_used = ref if abs(dx - 0.060) < 1e-9 else None
        plot_quantities(
            stats,
            title_prefix=f"{args.we_label} Δx={dx*1e3:.0f} mm",
            out_prefix=plot_prefix,
            ref=ref_used,
        )
        print(f"[ok ] wrote {csv_path.name}")

    (out_root / f"metadata_{args.we_label}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(f"[done] outputs in {out_root}")


if __name__ == "__main__":
    main()
