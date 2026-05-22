#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-processing for OpenFOAM lagrangian droplet data,
aligned with Zhang et al. (2025), Energy 335, 138192.

Methodology (from the paper, Section 3.4 and Appendices B/C):
  - Statistics collected at the y-z plane Δx = 60 mm downstream
    of the liquid jet outlet.
  - Droplet count, mean velocity, D10 (average diameter), and
    D32 (Sauter mean diameter) as functions of z-height.
  - Minimum 25 droplets per z-bin to ensure statistical reliability.
  - In the paper's VOF-LPT framework each Lagrangian particle is
    one physical droplet.  In OpenFOAM spray, each parcel represents
    nParticle droplets, so all statistics are weighted by nParticle
    (but NOT by velocity flux — the paper uses snapshot/presence-based
    counting, not flux-weighted MLDI-style sampling).

Reference data digitized from Fig. 24 of Zhang et al. (2025):
  - We_g = 500  (green curves, panels a₁/a₂, b₁/b₂, c₁, d₁)
  - We_g = 3500 (blue curves, panels a₁, b₁, c₂, d₂)
"""

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

try:
    from scipy.ndimage import gaussian_filter, label
    SCIPY_OK = True
except Exception:
    SCIPY_OK = False


# =============================================================================
# USER CONFIG
# =============================================================================

CASE_DIR = Path(".").resolve()
CLOUD = "sprayCloud"

# Geometry
x_jet = 0.015  # [m] liquid jet exit x-location
domain = {
    "x": (-0.004, 0.080),
    "y": (-0.016, 0.016),
    "z": (0.000, 0.065),
}

# Stable-stage time sampling
times = [
    "0.02",
    "0.0205",
    "0.021",
    "0.0215",
    "0.022",
    "0.0225",
    "0.023",
    "0.0235",
    "0.024",
    "0.0245",
    "0.025",
    "0.0255",
    "0.026",
    "0.0265",
    "0.027",
    "0.0275",
    "0.028",
    "0.0285",
    "0.029",
    "0.0295",
    "0.03",
]

# Statistics plane: ONLY Δx = 60 mm downstream of jet exit
# (per Zhang et al. 2025, Section 3.4 — the validated location)
x_offset = 0.060  # [m]
slice_half_thickness = 5e-4  # [m]

# z-binning for sectional statistics
nz_stats = 13
z_bins = np.linspace(domain["z"][0], domain["z"][1], nz_stats + 1)

# Minimum droplets per bin (paper: 25, per Section 3.4)
min_droplets_per_bin = 25

# 3D concentration grid
nx, ny, nz = 84, 32, 65

# Concentration-map processing
phi_thresh = 0.10
min_region_pixels = 10
gaussian_sigma = 3.0

# Projection method
projection_mode = "sum"  # "sum" or "max"

# Velocity statistic
use_velocity_magnitude = True

OUTDIR = CASE_DIR / "paper_style_postproc_t0200_t0300"
OUTDIR.mkdir(exist_ok=True)


# =============================================================================
# DIGITIZED REFERENCE DATA - Zhang et al. (2025), Energy 335, Fig. 24
# Δx = 60 mm, y-z plane
#
# We_g = 500  (red curves in Fig. 24)
# We_g = 3500 (green curves in Fig. 24)
#
# Values read from the figure panels at 300 DPI.  Where the right-hand
# zoom panels (a₂/b₂/c₂/d₂) give better resolution they were preferred.
# =============================================================================

ZHANG2025_REF = {
    # --- We_g = 500 (red curves in Fig. 24) -----------------------------------
    #
    # Panel a₂ (particle count, zoom 0–1500, red curve with square markers).
    # Markers at z ≈ 3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57 mm.
    "particle_count_We500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([1420, 1500, 1580, 1580, 1260, 680, 450, 420, 300, 180, 90, 55]),
    },
    # Panel b₂ (velocity, zoom 0–25 m/s, red curve).
    "U_We500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([20.0, 23.0, 18.0, 13.0, 9.3, 8.2, 7.0, 7.0, 7.0, 7.5, 7.8, 8.5]),
    },
    # Panel c₁ (D10, 0–300 µm, red curve).
    "D10_We500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([85, 90, 98, 112, 130, 138, 135, 132, 140, 143, 142, 150]),
    },
    # Panel d₁ (D32/SMD, 0–700 µm, red curve).
    "D32_We500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([170, 165, 245, 285, 300, 335, 340, 335, 310, 275, 265, 200]),
    },

    # --- We_g = 3500 (green curves in Fig. 24) --------------------------------
    #
    # Panel a₁ (particle count, full scale ×10⁴, green curve with triangles).
    # Data points at roughly 5 mm spacing.
    "particle_count_We3500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([14000, 17500, 20500, 18000, 14000, 12000, 8500, 6500, 4200, 3000, 2500, 1200]),
    },
    # Panel b₁ (velocity, 0–100 m/s, green curve).
    "U_We3500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([62, 48, 37, 29, 28, 25, 25, 24, 24, 24, 24, 24]),
    },
    # Panel c₂ (D10, zoom 40–100 µm, green curve).  Read from zoom panel.
    "D10_We3500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([64, 68, 75, 79, 82, 84, 85, 87, 88, 89, 90, 94]),
    },
    # Panel d₂ (D32/SMD, zoom 50–120 µm, green curve).  Read from zoom panel.
    "D32_We3500": {
        "z_mm": np.array([3, 7, 12, 17, 22, 27, 32, 35, 40, 45, 50, 57]),
        "value": np.array([73, 79, 89, 99, 101, 101, 101, 107, 108, 108, 113, 113]),
    },
}


# =============================================================================
# OPENFOAM READERS
# =============================================================================

def read_positions(path):
    """Robust reader for OpenFOAM vectorField files."""
    txt = Path(path).read_text(encoding="latin-1", errors="ignore")
    triples = re.findall(
        r"\(\s*"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*"
        r"\)",
        txt
    )
    if not triples:
        return np.empty((0, 3), dtype=float)
    return np.array([[float(a), float(b), float(c)] for a, b, c in triples], dtype=float)


def read_vector(path):
    return read_positions(path)


def read_scalar(path):
    """Robust reader for OpenFOAM scalarField files."""
    txt = Path(path).read_text(encoding="latin-1", errors="ignore")
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
    if end == -1:
        raise RuntimeError(f"Cannot find scalar list end in {path}")
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
# DATA LOADING
# =============================================================================

def load_ensemble(case_dir: Path, cloud: str, time_names):
    all_pos, all_d, all_np, all_U, used_times = [], [], [], [], []

    for t in time_names:
        base = case_dir / t / "lagrangian" / cloud
        pos_candidates = [base / "positions", base / "position0"]

        pth_pos = None
        for cand in pos_candidates:
            if cand.is_file():
                pth_pos = cand
                break

        pth_d = base / "d"
        pth_np = base / "nParticle"
        pth_U = base / "U"

        if pth_pos is None or (not pth_d.is_file()) or (not pth_np.is_file()):
            print(f"[skip] time {t}: missing position/d/nParticle")
            continue

        pos = read_positions(pth_pos)
        d = read_scalar(pth_d)
        npart = read_scalar(pth_np)

        if pth_U.is_file():
            U = read_vector(pth_U)
        else:
            U = None

        print(
            f"[debug] time {t} | posFile={pth_pos.name} | "
            f"pos={len(pos)} d={len(d)} nP={len(npart)} "
            f"U={len(U) if U is not None else 'None'}"
        )

        if not (len(pos) == len(d) == len(npart)):
            raise RuntimeError(
                f"Mismatch at time {t}: pos={len(pos)}, d={len(d)}, nP={len(npart)} "
                f"(position field used: {pth_pos.name})"
            )
        if U is not None and len(U) != len(d):
            raise RuntimeError(
                f"Mismatch at time {t}: U={len(U)}, d={len(d)} "
                f"(position field used: {pth_pos.name})"
            )

        all_pos.append(pos)
        all_d.append(d)
        all_np.append(npart)
        all_U.append(U)
        used_times.append(t)

    if not used_times:
        raise RuntimeError("No valid times found. Check time list and lagrangian paths.")

    pos = np.vstack(all_pos)
    d = np.hstack(all_d)
    npart = np.hstack(all_np)
    have_U = all(u is not None for u in all_U)
    U = np.vstack(all_U) if have_U else None

    return pos, d, npart, U, used_times


# =============================================================================
# SECTIONAL STATISTICS AT Δx = 60 mm
#
# Weighting: nParticle only (snapshot/presence-based).
# This matches the VOF-LPT paper methodology where each Lagrangian
# particle = 1 physical droplet.  In OpenFOAM spray, nParticle converts
# parcels → physical droplets.  NO flux (|u_x|) weighting is applied.
# =============================================================================

def weighted_mean(values, weights):
    den = np.sum(weights)
    if den <= 0:
        return np.nan
    return np.sum(weights * values) / den


def weighted_D10(d, weights):
    return weighted_mean(d, weights)


def weighted_D32(d, weights):
    num = np.sum(weights * d**3)
    den = np.sum(weights * d**2)
    if den <= 0:
        return np.nan
    return num / den


def compute_sectional_stats(pos, d, npart, U, x_plane, dx, z_bins, min_droplets):
    """
    Collect droplet statistics in z-bins within a slab |x - x_plane| <= dx.

    All averages are weighted by nParticle (snapshot weighting).
    Bins with fewer than min_droplets equivalent droplets are masked as NaN.
    """
    mask = np.abs(pos[:, 0] - x_plane) <= dx
    pos_p = pos[mask]
    d_p = d[mask]
    n_p = npart[mask]
    U_p = U[mask] if U is not None else None

    z = pos_p[:, 2]
    z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

    out = {
        "z_center_m": z_centers,
        "count_parcels": np.zeros_like(z_centers, dtype=int),
        "count_droplets": np.zeros_like(z_centers, dtype=float),
        "Umean_m_s": np.full_like(z_centers, np.nan, dtype=float),
        "D10_m": np.full_like(z_centers, np.nan, dtype=float),
        "D32_m": np.full_like(z_centers, np.nan, dtype=float),
    }

    for i, (z0, z1) in enumerate(zip(z_bins[:-1], z_bins[1:])):
        m = (z >= z0) & (z < z1)
        di = d_p[m]
        ni = n_p[m]

        out["count_parcels"][i] = di.size
        n_droplets = np.sum(ni)
        out["count_droplets"][i] = n_droplets

        # Apply paper threshold: min 25 equivalent droplets
        if n_droplets < min_droplets:
            continue

        # --- nParticle weighting (snapshot-based, NO flux) ---
        out["D10_m"][i] = weighted_D10(di, ni)
        out["D32_m"][i] = weighted_D32(di, ni)

        if U_p is not None:
            if use_velocity_magnitude:
                ui = np.linalg.norm(U_p[m], axis=1)
            else:
                ui = U_p[m][:, 0]
            out["Umean_m_s"][i] = weighted_mean(ui, ni)

    return out


def save_sectional_stats(stats, out_csv: Path):
    keys = ["z_center_m", "count_parcels", "count_droplets", "Umean_m_s", "D10_m", "D32_m"]
    arr = np.column_stack([stats[k] for k in keys])
    np.savetxt(out_csv, arr, delimiter=",", header=",".join(keys), comments="")


def plot_sectional_stats(stats, title_prefix: str, out_png_prefix: Path):
    z_mm = stats["z_center_m"] * 1e3

    plt.figure(figsize=(8, 5))
    plt.plot(z_mm, stats["count_parcels"], marker="o")
    plt.xlabel("z [mm]")
    plt.ylabel("Parcel count per bin [-]")
    plt.title(f"{title_prefix} - parcel count")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(str(out_png_prefix) + "_count.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(z_mm, stats["count_droplets"], marker="o")
    plt.xlabel("z [mm]")
    plt.ylabel("Equivalent droplet count per bin [-]")
    plt.title(f"{title_prefix} - droplet count (weighted by nParticle)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(str(out_png_prefix) + "_dropletCount.png", dpi=200)
    plt.close()

    if np.any(np.isfinite(stats["Umean_m_s"])):
        plt.figure(figsize=(8, 5))
        plt.plot(z_mm, stats["Umean_m_s"], marker="o")
        plt.xlabel("z [mm]")
        plt.ylabel("Mean droplet velocity [m/s]")
        plt.title(f"{title_prefix} - mean velocity")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(str(out_png_prefix) + "_Umean.png", dpi=200)
        plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(z_mm, stats["D10_m"] * 1e6, marker="o")
    plt.xlabel("z [mm]")
    plt.ylabel(r"$D_{10}$ [µm]")
    plt.title(f"{title_prefix} - D10")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(str(out_png_prefix) + "_D10.png", dpi=200)
    plt.close()

    plt.figure(figsize=(8, 5))
    plt.plot(z_mm, stats["D32_m"] * 1e6, marker="o")
    plt.xlabel("z [mm]")
    plt.ylabel(r"$D_{32}$ [µm]")
    plt.title(f"{title_prefix} - D32")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(str(out_png_prefix) + "_D32.png", dpi=200)
    plt.close()


# =============================================================================
# OVERLAY PLOTS WITH ZHANG ET AL. (2025) REFERENCE — Fig. 24 style
# =============================================================================

def plot_overlay_vs_zhang2025(outdir, of_stats, quantity_key, ylabel, out_png,
                              ref_data=ZHANG2025_REF):
    """
    Single-panel overlay: OpenFOAM line at Δx=60mm vs digitized Zhang (2025)
    reference curves for We=500 and We=3500.
    """
    fig, ax = plt.subplots(figsize=(8, 5.5))

    # --- OpenFOAM data ---
    z_mm = of_stats["z_center_m"] * 1e3
    y_of = of_stats[quantity_key]
    if quantity_key in ["D10_m", "D32_m"]:
        y_of = y_of * 1e6
    mask = np.isfinite(y_of)

    ax.plot(z_mm[mask], y_of[mask],
            color="black", lw=2.0, marker="o", ms=5, zorder=10,
            label=r"OpenFOAM $\Delta x = 60$ mm")

    # --- Zhang (2025) reference ---
    # Map quantity_key → reference dict keys
    ref_map = {
        "Umean_m_s":       ("U_We500",              "U_We3500"),
        "D10_m":           ("D10_We500",             "D10_We3500"),
        "D32_m":           ("D32_We500",             "D32_We3500"),
        "count_droplets":  ("particle_count_We500",  "particle_count_We3500"),
    }

    if quantity_key in ref_map:
        key500, key3500 = ref_map[quantity_key]

        if key500 in ref_data:
            rd = ref_data[key500]
            ax.plot(rd["z_mm"], rd["value"],
                    color="red", lw=1.5, ls="--", marker="s", ms=4,
                    fillstyle="none", zorder=5,
                    label=r"Zhang (2025) $We_g = 500$")

        if key3500 in ref_data:
            rd = ref_data[key3500]
            ax.plot(rd["z_mm"], rd["value"],
                    color="green", lw=1.5, ls="--", marker="^", ms=4,
                    fillstyle="none", zorder=5,
                    label=r"Zhang (2025) $We_g = 3500$")

    ax.set_xlabel("z (mm)", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.legend(frameon=True, fontsize=9)
    ax.set_xlim(0, 65)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(outdir / out_png, dpi=300)
    plt.close(fig)


# =============================================================================
# 3D CONCENTRATION FIELD + PROJECTIONS
# =============================================================================

def compute_3d_volume_fraction(pos, d, npart, domain, nx, ny, nz):
    x0, x1 = domain["x"]
    y0, y1 = domain["y"]
    z0, z1 = domain["z"]
    H, edges = np.histogramdd(
        sample=pos,
        bins=(nx, ny, nz),
        range=((x0, x1), (y0, y1), (z0, z1)),
        weights=npart * (np.pi * d**3 / 6.0)
    )
    dx = (x1 - x0) / nx
    dy = (y1 - y0) / ny
    dz = (z1 - z0) / nz
    phi = H / (dx * dy * dz)
    return phi, edges, (dx, dy, dz)


def log_normalize(phi):
    eps = np.finfo(float).tiny
    logphi = np.log10(phi + eps)
    vals = logphi[np.isfinite(logphi)]
    vmin = float(vals.min())
    vmax = float(vals.max())
    if vmax - vmin < 1e-16:
        phi_n = np.zeros_like(logphi)
    else:
        phi_n = (logphi - vmin) / (vmax - vmin)
    phi_n[~np.isfinite(phi_n)] = 0.0
    return phi_n, {"log10_min": vmin, "log10_max": vmax}


def project_field(phi_norm, mode="sum"):
    if mode == "sum":
        front = np.sum(phi_norm, axis=1)
        top   = np.sum(phi_norm, axis=2)
        left  = np.sum(phi_norm, axis=0)
    elif mode == "max":
        front = np.max(phi_norm, axis=1)
        top   = np.max(phi_norm, axis=2)
        left  = np.max(phi_norm, axis=0)
    else:
        raise ValueError("projection_mode must be 'sum' or 'max'")

    def _norm(a):
        amin = np.min(a)
        amax = np.max(a)
        if amax - amin < 1e-16:
            return np.zeros_like(a)
        return (a - amin) / (amax - amin)

    return _norm(front), _norm(top), _norm(left)


def remove_small_regions(mask, min_pixels):
    if not SCIPY_OK:
        return mask
    labels, nlab = label(mask.astype(np.int8))
    if nlab == 0:
        return mask
    out = np.zeros_like(mask, dtype=bool)
    for lab_id in range(1, nlab + 1):
        region = labels == lab_id
        if np.count_nonzero(region) >= min_pixels:
            out |= region
    return out


def process_projection_image(img_norm, phi_thresh, min_region_pixels, sigma):
    mask = img_norm > phi_thresh
    mask = remove_small_regions(mask, min_region_pixels)
    img_smooth = gaussian_filter(img_norm.astype(float), sigma=sigma) if SCIPY_OK else img_norm.astype(float)
    contour_band = (img_smooth >= 0.4) & (img_smooth <= 0.6)
    return {"mask": mask, "smooth": img_smooth, "contour_band": contour_band}


def save_projection_plot(img, title, out_png, extent, cmap="turbo"):
    plt.figure(figsize=(8, 6))
    plt.imshow(img.T, origin="lower", aspect="auto", extent=extent, cmap=cmap)
    plt.xlabel(title.split(" vs ")[0] + " [mm]")
    plt.ylabel(title.split(" vs ")[1] + " [mm]")
    plt.title(title)
    plt.colorbar(label="normalized projected concentration [-]")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def save_boundary_overlay(img, proc, title, out_png, extent):
    plt.figure(figsize=(8, 6))
    im = plt.imshow(img.T, origin="lower", aspect="auto", extent=extent, cmap="turbo")
    plt.contour(proc["mask"].T.astype(float), levels=[0.5], extent=extent, colors="w", linewidths=1.0)
    plt.contour(proc["contour_band"].T.astype(float), levels=[0.5], extent=extent, colors="k", linewidths=1.0)
    plt.xlabel(title.split(" vs ")[0] + " [mm]")
    plt.ylabel(title.split(" vs ")[1] + " [mm]")
    plt.title(title + " - boundary extraction")
    if np.nanmax(img) > np.nanmin(img):
        plt.colorbar(im, label="normalized projected concentration [-]")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


# =============================================================================
# MAIN
# =============================================================================

def main():
    pos, d, npart, U, used_times = load_ensemble(CASE_DIR, CLOUD, times)

    meta = {
        "case_dir": str(CASE_DIR),
        "cloud": CLOUD,
        "used_times": used_times,
        "x_jet_m": x_jet,
        "x_offset_m": x_offset,
        "slice_half_thickness_m": slice_half_thickness,
        "min_droplets_per_bin": min_droplets_per_bin,
        "domain": domain,
        "grid_shape": [nx, ny, nz],
        "projection_mode": projection_mode,
        "weighting_method": "nParticle (snapshot, no flux)",
        "reference": "Zhang et al. (2025), Energy 335, 138192, Fig. 24",
        "openfoam_adaptation": {"use_nParticle_weighting": True,
                                 "use_flux_weighting": False},
    }
    (OUTDIR / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    # ------------------------------------------------------------------
    # Sectional statistics at Δx = 60 mm
    # ------------------------------------------------------------------
    x_plane = x_jet + x_offset
    stats = compute_sectional_stats(
        pos, d, npart, U, x_plane, slice_half_thickness, z_bins,
        min_droplets_per_bin
    )
    tag = f"x_{int(round(x_offset*1e3)):03d}mm_downstream"
    save_sectional_stats(stats, OUTDIR / f"{tag}.csv")
    plot_sectional_stats(stats, tag, OUTDIR / tag)

    # ------------------------------------------------------------------
    # Overlay plots vs Zhang (2025) Fig. 24
    # ------------------------------------------------------------------
    plot_overlay_vs_zhang2025(
        OUTDIR, stats,
        quantity_key="Umean_m_s",
        ylabel=r"$U$ [m/s]",
        out_png="overlay_U_vs_zhang2025.png"
    )
    plot_overlay_vs_zhang2025(
        OUTDIR, stats,
        quantity_key="D10_m",
        ylabel=r"$D_{10}$ [µm]",
        out_png="overlay_D10_vs_zhang2025.png"
    )
    plot_overlay_vs_zhang2025(
        OUTDIR, stats,
        quantity_key="D32_m",
        ylabel=r"$D_{32}$ [µm]",
        out_png="overlay_D32_vs_zhang2025.png"
    )
    plot_overlay_vs_zhang2025(
        OUTDIR, stats,
        quantity_key="count_droplets",
        ylabel="Equivalent droplet count [-]",
        out_png="overlay_count_vs_zhang2025.png"
    )

    # ------------------------------------------------------------------
    # 3D concentration field & projections
    # ------------------------------------------------------------------
    phi, edges, spacing = compute_3d_volume_fraction(pos, d, npart, domain, nx, ny, nz)
    phi_norm, norm_info = log_normalize(phi)

    np.savez_compressed(
        OUTDIR / "concentration_3d.npz",
        phi=phi,
        phi_norm=phi_norm,
        x_edges=edges[0],
        y_edges=edges[1],
        z_edges=edges[2],
        spacing=np.array(spacing),
    )

    front, top, left = project_field(phi_norm, mode=projection_mode)

    x_mm = np.array(domain["x"]) * 1e3
    y_mm = np.array(domain["y"]) * 1e3
    z_mm = np.array(domain["z"]) * 1e3

    save_projection_plot(front, "x vs z", OUTDIR / "projection_front_xz.png",
                         [x_mm[0], x_mm[1], z_mm[0], z_mm[1]])
    save_projection_plot(top,   "x vs y", OUTDIR / "projection_top_xy.png",
                         [x_mm[0], x_mm[1], y_mm[0], y_mm[1]])
    save_projection_plot(left,  "y vs z", OUTDIR / "projection_left_yz.png",
                         [y_mm[0], y_mm[1], z_mm[0], z_mm[1]])

    front_proc = process_projection_image(front, phi_thresh, min_region_pixels, gaussian_sigma)
    top_proc   = process_projection_image(top,   phi_thresh, min_region_pixels, gaussian_sigma)
    left_proc  = process_projection_image(left,  phi_thresh, min_region_pixels, gaussian_sigma)

    save_boundary_overlay(front, front_proc, "x vs z",
                          OUTDIR / "boundary_front_xz.png",
                          [x_mm[0], x_mm[1], z_mm[0], z_mm[1]])
    save_boundary_overlay(top, top_proc, "x vs y",
                          OUTDIR / "boundary_top_xy.png",
                          [x_mm[0], x_mm[1], y_mm[0], y_mm[1]])
    save_boundary_overlay(left, left_proc, "y vs z",
                          OUTDIR / "boundary_left_yz.png",
                          [y_mm[0], y_mm[1], z_mm[0], z_mm[1]])

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\nUsed times: {used_times}")
    print(f"Total merged parcels: {len(d)}")
    print(f"Weighting: nParticle only (snapshot, no flux)")
    print(f"Min droplets/bin: {min_droplets_per_bin}")
    print(f"Statistics plane: Δx = {x_offset*1e3:.0f} mm")
    print(f"Outputs written to: {OUTDIR}")


if __name__ == "__main__":
    main()