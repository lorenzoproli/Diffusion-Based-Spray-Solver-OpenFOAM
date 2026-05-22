#!/usr/bin/env python3
"""
Lambert-style post-processing for OpenFOAM ParcelPlaneSampler data.

Purpose
-------
Reads crossing-event files written by the custom OpenFOAM ParcelPlaneSampler on
planes x/d = 30 and x/d = 60, computes nParticle-weighted profiles of:

    D32(y) = sum(n_i d_i^3) / sum(n_i d_i^2)
    Ux(y)  = sum(n_i Ux_i)  / sum(n_i)

and compares them against digitized Lambert/experiment reference curves.

Two metrics are reported per bin:

  * "total"     : all valid parcels in the bin
  * "filtered"  : parcels with ms != -1 only, i.e. excluding parcels that are
                  still in the first Madabhushi PE breakup (deforming blobs,
                  pre-tb). This matches the experimental observable of Sekar
                  et al. (PDPA / shadowgraphy), which only resolves stable
                  spherical droplets and not deforming blobs in transit.

Expected sampler columns, new format:
    time x y z Ux Uy Uz d nParticle user ms tc KHindex yState yDotState faceID

The script is backward-compatible with the older sampler format:
    time x y z Ux Uy Uz d nParticle faceID

Recommended execution from the case directory:
    python3 postProcessing/postprocess_lambert_openfoam.py --case .

If automatic file discovery fails, pass explicit directories/files:
    python3 postProcessing/postprocess_lambert_openfoam.py --case . \
        --plane30 postProcessing/.../parcelPlaneSampler_30d \
        --plane60 postProcessing/.../parcelPlaneSampler_60d
"""

from __future__ import annotations

from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# User configuration defaults
# -----------------------------------------------------------------------------

D_JET = 457e-6          # [m] jet diameter, Lambert J=10 reference
DZ_CENTER_DEFAULT = 0.001   # [m] centerline half-width in z, |z| <= 1 mm
N_BINS_Y_DEFAULT = 30
Y_MIN_DEFAULT = 0.0         # [m]
Y_MAX_DEFAULT = 0.015       # [m]
MIN_EVENTS_DEFAULT = 500    # conservative filtering for stable profile plots

PLANES = {
    "30d": {"x_nd": 30},
    "60d": {"x_nd": 60},
}

# -----------------------------------------------------------------------------
# Digitized reference data: [y_mm, value]
# -----------------------------------------------------------------------------

REF = {
    "30d": {
        "fluent_smd": np.array([
            [0.57, 15.37], [0.97, 16.04], [1.35, 17.09], [1.67, 17.92],
            [2.24, 19.20], [2.78, 20.42], [3.37, 21.63], [3.65, 22.40],
            [4.19, 24.28], [4.52, 25.69], [5.01, 27.54], [5.11, 27.76],
            [5.48, 28.47], [5.61, 28.63], [5.73, 28.69], [6.13, 28.47],
            [6.65, 27.38], [6.92, 26.74], [7.12, 26.10], [7.34, 25.59],
            [7.50, 25.21], [7.89, 24.92], [8.56, 24.86], [9.08, 24.89],
            [9.54, 25.91]
        ]),
        "exp_smd": np.array([
            [1.00, 17.54], [2.01, 19.52], [3.01, 22.01], [4.00, 25.50],
            [5.01, 25.34], [6.01, 22.68], [7.00, 21.76], [8.04, 22.84],
            [9.02, 24.47], [10.03, 24.92]
        ]),
        "fluent_ux": np.array([
            [0.57, 96.86], [0.98, 95.56], [1.36, 92.12], [1.68, 87.86],
            [2.26, 79.21], [2.80, 71.75], [3.37, 66.54], [3.66, 65.35],
            [4.19, 66.31], [4.54, 67.73], [5.00, 71.29], [5.10, 71.65],
            [5.46, 74.62], [5.57, 75.45], [5.71, 76.28], [6.16, 78.42],
            [6.64, 79.37], [6.90, 79.85], [7.10, 80.09], [7.31, 81.04],
            [7.50, 81.63], [7.86, 83.41], [8.53, 85.79], [9.07, 86.75],
            [9.53, 87.11]
        ]),
        "exp_ux": np.array([
            [1.01, 73.74], [2.03, 71.26], [3.01, 62.27], [4.02, 61.33],
            [5.00, 68.80], [6.00, 85.41], [7.00, 93.36], [8.03, 95.39],
            [9.02, 92.91], [9.99, 90.43]
        ]),
    },
    "60d": {
        "fluent_smd": np.array([
            [1.03, 16.67], [1.39, 17.33], [1.82, 18.06], [2.23, 18.98],
            [2.73, 19.65], [2.85, 19.81], [3.02, 19.97], [3.25, 20.29],
            [3.50, 20.79], [3.95, 21.49], [4.14, 21.81], [4.41, 22.35],
            [4.74, 22.95], [4.90, 23.37], [5.13, 24.00], [5.27, 24.29],
            [5.52, 25.08], [5.71, 25.52], [5.98, 26.25], [6.37, 27.27],
            [6.85, 27.94], [7.24, 28.22], [7.61, 28.38], [8.09, 28.38],
            [8.48, 27.84], [8.75, 27.08], [8.95, 26.44], [9.16, 25.81],
            [9.53, 25.21], [9.97, 25.08], [10.24, 24.98], [10.40, 25.05],
            [10.77, 25.43], [11.13, 25.43], [11.33, 25.78], [11.77, 27.65],
            [11.91, 28.10], [12.24, 26.25]
        ]),
        "exp_smd": np.array([
            [1.01, 17.17], [1.99, 18.89], [3.00, 20.29], [4.01, 22.67],
            [5.02, 24.35], [6.02, 23.65], [7.01, 23.33], [7.98, 23.11],
            [9.02, 23.49], [10.01, 24.86], [11.02, 25.65], [12.01, 26.41]
        ]),
        # Updated x/d=60 Ux Lambert/Fluent dataset provided by user.
        "fluent_ux": np.array([
            [1.0356164383561648, 99.2945205479452],
            [1.4000000000000021, 97.8690476190476],
            [1.8219178082191796, 95.61073059360729],
            [2.2246575342465764, 92.51891715590345],
            [2.7041095890410958, 88.47537508153945],
            [2.8383561643835620, 87.04794520547944],
            [3.0109589041095894, 85.26369863013697],
            [3.2602739726027394, 82.76581865622961],
            [3.5095890410958910, 80.62508153946509],
            [3.9506849315068493, 77.65264187866927],
            [4.1424657534246590, 76.70189171559032],
            [4.4301369863013700, 75.51386170906719],
            [4.7178082191780835, 74.80202217873449],
            [4.9095890410958920, 74.80365296803652],
            [5.1397260273972610, 74.80560991519894],
            [5.3698630136986300, 75.28375733855184],
            [5.5424657534246580, 75.99951076320939],
            [5.6958904109589050, 76.71510110893671],
            [5.9835616438356160, 78.02707110241354],
            [6.3671232876712320, 79.81604696673189],
            [6.8657534246575340, 81.9631441617743],
            [7.2493150684931500, 83.87116764514025],
            [7.5945205479452060, 85.1836268754077],
            [8.0931506849315070, 87.21167645140247],
            [8.4958904109589040, 89.83414872798433],
            [8.7643835616438360, 91.50309849967383],
            [8.9561643835616420, 92.57615786040442],
            [9.1671232876712310, 93.5303326810176],
            [9.4931506849315070, 94.6045335942596],
            [9.9534246575342460, 95.4417808219178],
            [10.164383561643834, 95.68166992824526],
            [10.394520547945204, 95.92172211350292],
            [10.797260273972600, 96.04419439008478],
            [11.104109589041094, 96.04680365296802],
            [11.334246575342465, 96.04876060013045],
            [11.717808219178082, 95.57583170254402],
            [11.909589041095890, 95.45841487279843],
            [12.235616438356164, 95.6992824527071],
        ]),
        # Updated x/d=60 Ux experiment dataset provided by user.
        "exp_ux": np.array([
            [0.9972602739726035, 82.38943248532289],
            [1.9753424657534246, 85.13584474885845],
            [3.0109589041095894, 81.81131767775602],
            [3.9890410958904106, 78.72439660795824],
            [5.0054794520547930, 79.20923026744944],
            [6.0027397260273965, 88.02723418134376],
            [7.0000000000000000, 95.65476190476188],
            [8.0164383561643820, 99.71102413568165],
            [9.0136986301369860, 99.48140900195693],
            [10.030136986301368, 99.7281474233529],
            [11.008219178082191, 98.78408349641225],
            [11.986301369863012, 97.0066862361383],
        ]),
    },
}

# -----------------------------------------------------------------------------
# Sampler formats
# -----------------------------------------------------------------------------

OLD_COLS = ["time", "x", "y", "z", "Ux", "Uy", "Uz", "d", "nParticle", "faceID"]
NEW_COLS = [
    "time", "x", "y", "z", "Ux", "Uy", "Uz", "d", "nParticle",
    "user", "ms", "tc", "KHindex", "yState", "yDotState", "faceID"
]

PLANE_TOKENS = {
    "30d": ["parcelplanesampler_30d", "sampler_30d", "plane_30d", "plane-30", "30d"],
    "60d": ["parcelplanesampler_60d", "sampler_60d", "plane_60d", "plane-60", "60d"],
}

VALID_SUFFIXES = {".dat", ".csv", ".smp", ".txt"}


def _time_sort_key(path: Path) -> tuple[int, float, str]:
    for part in reversed(path.parts):
        try:
            return (0, float(part), str(path))
        except ValueError:
            pass
    return (1, 0.0, str(path))


def _first_data_field_count(path: Path) -> int:
    try:
        with path.open("r", errors="ignore") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                return len(re.split(r"[\s,]+", s))
    except OSError:
        pass
    return 0


def find_plane_files(case_dir: Path, key: str, explicit: str | None = None, out_dir_name: str = "") -> list[Path]:
    """Return all files associated with one sampling plane.

    Option A assumes separate sampler objects/files for 30d and 60d. If an
    explicit path is a directory, all supported data files below it are read.
    """
    if explicit:
        p = Path(explicit)
        if not p.is_absolute():
            p = case_dir / p
        if p.is_file():
            return [p]
        if p.is_dir():
            return sorted(
                [q for q in p.rglob("*") if q.is_file() and q.suffix.lower() in VALID_SUFFIXES],
                key=_time_sort_key,
            )
        return []

    search_root = case_dir / "postProcessing" if (case_dir / "postProcessing").is_dir() else case_dir
    tokens = PLANE_TOKENS[key]
    matches: list[Path] = []
    for p in search_root.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in VALID_SUFFIXES:
            continue
        full = str(p).lower()
        if out_dir_name and out_dir_name.lower() in full:
            continue
        if any(tok in full for tok in tokens):
            matches.append(p)

    seen: set[Path] = set()
    out: list[Path] = []
    for p in sorted(matches, key=_time_sort_key):
        if p not in seen:
            out.append(p)
            seen.add(p)
    return out


def read_openfoam_sampler(path: Path) -> pd.DataFrame:
    """Read one ParcelPlaneSampler output file, new or old format."""
    if not path.exists():
        return pd.DataFrame()

    n_fields = _first_data_field_count(path)
    if n_fields >= len(NEW_COLS):
        names = NEW_COLS
    elif n_fields == len(OLD_COLS):
        names = OLD_COLS
    elif n_fields > 0:
        # Fallback: read first n known columns only.
        names = NEW_COLS[:n_fields]
    else:
        return pd.DataFrame()

    df = pd.read_csv(
        path,
        sep=r"\s+|,",
        comment="#",
        header=None,
        names=names,
        engine="python",
        on_bad_lines="skip",
    )
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    required = ["y", "Ux", "d", "nParticle"]
    return df.dropna(subset=[c for c in required if c in df.columns])


def read_plane_dataset(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for p in paths:
        df = read_openfoam_sampler(p)
        if not df.empty:
            df["sourceFile"] = str(p)
            frames.append(df)
    if not frames:
        return pd.DataFrame(columns=NEW_COLS)
    return pd.concat(frames, ignore_index=True, sort=False)


def filter_transits(df: pd.DataFrame, dz_center: float | None) -> pd.DataFrame:
    if df.empty:
        return df
    valid = (
        np.isfinite(df["y"]) & np.isfinite(df["Ux"]) &
        np.isfinite(df["d"]) & np.isfinite(df["nParticle"])
    )
    valid &= (df["d"] > 0.0) & (df["nParticle"] > 0.0)
    if "z" in df.columns and dz_center is not None:
        valid &= df["z"].abs() <= dz_center
    return df.loc[valid].reset_index(drop=True)


def _weighted_stats(d: np.ndarray, ux: np.ndarray, n: np.ndarray) -> tuple[float, float, float, float, float]:
    """Return D32_m, Ux_mean, n_phys_sum, frac_events_d_lt_2um, frac_mass_d_lt_2um."""
    if len(d) == 0:
        return (np.nan, np.nan, np.nan, np.nan, np.nan)
    snd2 = np.sum(n*d**2)
    snd3 = np.sum(n*d**3)
    sn = np.sum(n)
    small = d < 2.0e-6
    d32 = snd3/snd2 if snd2 > 0 else np.nan
    ux_mean = np.sum(n*ux)/sn if sn > 0 else np.nan
    frac_events_small = float(np.sum(small))/len(d)
    frac_mass_small = np.sum(n[small]*d[small]**3)/(snd3 + 1e-300)
    return d32, ux_mean, sn, frac_events_small, frac_mass_small


def bin_profile(df: pd.DataFrame, n_bins: int, y_min: float, y_max: float, min_events: int) -> pd.DataFrame:
    edges = np.linspace(y_min, y_max, n_bins + 1)
    centers = 0.5*(edges[:-1] + edges[1:])

    base_cols = {
        "y_center": centers,
        "y_center_mm": centers*1e3,
        "D32_m": np.nan,
        "D32_um": np.nan,
        "Ux_mean": np.nan,
        "n_events": 0,
        "n_phys_sum": np.nan,
        "frac_events_d_lt_2um": np.nan,
        "frac_mass_d_lt_2um": np.nan,
        # Filtered metric (ms != -1 only): excludes parcels still in the first
        # Madabhushi PE breakup (deforming blobs, pre-tb). Aligns the numerical
        # observable with shadowgraphy/PDPA which only resolves stable droplets.
        "D32_um_filtered": np.nan,
        "Ux_mean_filtered": np.nan,
        "n_events_stable": 0,
        "n_events_pre_pe": 0,
        "mass_fraction_pre_pe": np.nan,
    }

    # Diagnostic columns when new sampler states are available.
    for u in [0, 1, 2]:
        base_cols[f"n_events_user{u}"] = 0
        base_cols[f"mass_fraction_user{u}"] = np.nan
        base_cols[f"D32_um_user{u}"] = np.nan
        base_cols[f"Ux_mean_user{u}"] = np.nan
    for m_label, m_val in [("m10", -10), ("m20", -20), ("m2", 2)]:
        base_cols[f"n_events_{m_label}"] = 0
        base_cols[f"mass_fraction_{m_label}"] = np.nan

    out = pd.DataFrame(base_cols)
    if df.empty:
        return out

    y = df["y"].to_numpy(float)
    d = df["d"].to_numpy(float)
    ux = df["Ux"].to_numpy(float)
    n = df["nParticle"].to_numpy(float)
    bin_idx = np.digitize(y, edges) - 1
    valid = np.isfinite(y) & np.isfinite(d) & np.isfinite(ux) & np.isfinite(n) & (d > 0) & (n > 0)

    has_user = "user" in df.columns and df["user"].notna().any()
    has_ms = "ms" in df.columns and df["ms"].notna().any()
    user = np.rint(df["user"].to_numpy(float)).astype(int) if has_user else None
    ms = np.rint(df["ms"].to_numpy(float)).astype(int) if has_ms else None

    for k in range(n_bins):
        sel = (bin_idx == k) & valid
        n_events = int(sel.sum())
        if n_events < min_events:
            continue

        dk = d[sel]
        uk = ux[sel]
        nk = n[sel]
        d32, ux_mean, sn, fe_small, fm_small = _weighted_stats(dk, uk, nk)
        mass_total = np.sum(nk*dk**3)

        out.at[k, "n_events"] = n_events
        out.at[k, "n_phys_sum"] = sn
        out.at[k, "frac_events_d_lt_2um"] = fe_small
        out.at[k, "frac_mass_d_lt_2um"] = fm_small
        out.at[k, "D32_m"] = d32
        out.at[k, "D32_um"] = 1e6*d32 if np.isfinite(d32) else np.nan
        out.at[k, "Ux_mean"] = ux_mean

        # ----- Filtered D32: exclude parcels in pre-PE primary breakup ----- #
        # Sekar (PDPA / shadowgraphy) only resolves stable spherical droplets,
        # not deforming blobs in pre-breakup. Excluding ms == -1 aligns the
        # numerical metric with the experimental observable.
        if has_ms and ms is not None:
            mk_bin = ms[sel]
            sub_pre_pe = mk_bin == -1
            sub_stable = ~sub_pre_pe

            n_pre_pe = int(np.sum(sub_pre_pe))
            n_stable = int(np.sum(sub_stable))
            out.at[k, "n_events_pre_pe"] = n_pre_pe
            out.at[k, "n_events_stable"] = n_stable

            if n_pre_pe > 0:
                out.at[k, "mass_fraction_pre_pe"] = (
                    np.sum(nk[sub_pre_pe]*dk[sub_pre_pe]**3)/(mass_total + 1e-300)
                )

            if n_stable >= max(10, min_events // 10):
                d32f, uxf, _, _, _ = _weighted_stats(
                    dk[sub_stable], uk[sub_stable], nk[sub_stable]
                )
                out.at[k, "D32_um_filtered"] = (
                    1e6*d32f if np.isfinite(d32f) else np.nan
                )
                out.at[k, "Ux_mean_filtered"] = uxf

        if has_user and user is not None:
            uk_user = user[sel]
            for u in [0, 1, 2]:
                sub = uk_user == u
                out.at[k, f"n_events_user{u}"] = int(np.sum(sub))
                if np.any(sub):
                    d32u, uxu, _, _, _ = _weighted_stats(dk[sub], uk[sub], nk[sub])
                    out.at[k, f"D32_um_user{u}"] = 1e6*d32u if np.isfinite(d32u) else np.nan
                    out.at[k, f"Ux_mean_user{u}"] = uxu
                    out.at[k, f"mass_fraction_user{u}"] = np.sum(nk[sub]*dk[sub]**3)/(mass_total + 1e-300)

        if has_ms and ms is not None:
            mk = ms[sel]
            for m_label, m_val in [("m10", -10), ("m20", -20), ("m2", 2)]:
                sub = mk == m_val
                out.at[k, f"n_events_{m_label}"] = int(np.sum(sub))
                if np.any(sub):
                    out.at[k, f"mass_fraction_{m_label}"] = np.sum(nk[sub]*dk[sub]**3)/(mass_total + 1e-300)

    return out


def _plot_panel(ax, prof: pd.DataFrame, ref: dict, quantity: str, xlabel: str, title: str) -> None:
    if quantity == "D32_m":
        ref_f = ref.get("fluent_smd")
        ref_e = ref.get("exp_smd")
    else:
        ref_f = ref.get("fluent_ux")
        ref_e = ref.get("exp_ux")

    if ref_f is not None:
        ax.plot(ref_f[:, 1], ref_f[:, 0], "b-", lw=1.6, label="FLUENT (Lambert 2019)")
    if ref_e is not None:
        ax.plot(ref_e[:, 1], ref_e[:, 0], "ks", ms=5, label="EXP (Sekar et al. 2014)")

    if prof is not None and not prof.empty:
        y_mm = prof["y_center_mm"].to_numpy()

        # Total: all valid parcels in the bin (includes pre-PE blobs).
        x = prof["D32_um"].to_numpy() if quantity == "D32_m" else prof["Ux_mean"].to_numpy()
        ok = np.isfinite(x)
        ax.plot(x[ok], y_mm[ok], "ro--", lw=1.4, ms=4, label="OpenFOAM total")

        # Filtered: ms != -1 only (stable droplets, no pre-PE blobs).
        col_filtered = "D32_um_filtered" if quantity == "D32_m" else "Ux_mean_filtered"
        if col_filtered in prof.columns:
            xf = prof[col_filtered].to_numpy()
            okf = np.isfinite(xf)
            if np.any(okf):
                ax.plot(
                    xf[okf], y_mm[okf],
                    "g^-", lw=1.4, ms=4,
                    label="OpenFOAM filtered (ms ≠ -1)",
                )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("y [mm]")
    ax.set_title(title)
    ax.set_xlim(0, 40 if quantity == "D32_m" else 150)
    ax.set_ylim(0, 12 if "30" in title else 16)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")


def plot_all(profiles: dict[str, pd.DataFrame], out_png: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle("J=10, We=1500 — Lambert-style OpenFOAM plane statistics", fontsize=13)
    for key, row in [("30d", 0), ("60d", 1)]:
        prof = profiles.get(key, pd.DataFrame())
        x_nd = PLANES[key]["x_nd"]
        _plot_panel(axes[row, 0], prof, REF[key], "D32_m", "D32 [µm]", f"SMD — x/d = {x_nd}")
        _plot_panel(axes[row, 1], prof, REF[key], "Ux_mean", "Ux [m/s]", f"Mean Ux — x/d = {x_nd}")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def print_user_summary(df: pd.DataFrame, key: str) -> None:
    if df.empty or "user" not in df.columns:
        return
    user_counts = df["user"].round().astype("Int64").value_counts(dropna=True).sort_index()
    print(f" -> user counts ({key}, after filters):")
    for u, c in user_counts.items():
        print(f"      user={u}: {c}")
    if "ms" in df.columns:
        ms_counts = df["ms"].round().astype("Int64").value_counts(dropna=True).sort_index()
        print(f" -> ms counts ({key}, after filters, main states):")
        for m in [-20, -10, -1, 0, 2]:
            if m in ms_counts.index:
                print(f"      ms={m}: {ms_counts.loc[m]}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Lambert-style post-processing for OpenFOAM ParcelPlaneSampler files.")
    parser.add_argument("--case", type=Path, default=Path.cwd(), help="OpenFOAM case directory. Default: current directory")
    parser.add_argument("--plane30", type=str, default=None, help="File or directory for x/d=30 sampler data")
    parser.add_argument("--plane60", type=str, default=None, help="File or directory for x/d=60 sampler data")
    parser.add_argument("--dz", type=float, default=DZ_CENTER_DEFAULT, help="Centerline half-width in z [m]. Use negative value to disable filter.")
    parser.add_argument("--min-events", type=int, default=MIN_EVENTS_DEFAULT, help="Minimum events per y-bin")
    parser.add_argument("--bins", type=int, default=N_BINS_Y_DEFAULT, help="Number of y bins")
    parser.add_argument("--y-min", type=float, default=Y_MIN_DEFAULT, help="Minimum y [m]")
    parser.add_argument("--y-max", type=float, default=Y_MAX_DEFAULT, help="Maximum y [m]")
    args = parser.parse_args()

    case_dir = args.case.resolve()
    dz_center = None if args.dz < 0 else args.dz
    out_dir = case_dir / "postProcessing_Lambert_OpenFOAM"
    out_dir.mkdir(parents=True, exist_ok=True)

    explicit = {"30d": args.plane30, "60d": args.plane60}
    profiles: dict[str, pd.DataFrame] = {}

    print(f"Case directory: {case_dir}")
    if dz_center is None:
        print("Midline filter: disabled")
    else:
        print(f"Midline filter: |z| <= {dz_center*1e3:.2f} mm")
    print(f"MIN_EVENTS per bin: {args.min_events}")
    print("Mode: Option A, two separate samplers/files for x/d=30 and x/d=60")
    print()

    for key, info in PLANES.items():
        print(f"Processing plane {key} (x/d={info['x_nd']}) ...")
        files = find_plane_files(case_dir, key, explicit[key], out_dir.name)
        if not files:
            print(" -> WARNING: no file found; plotting references only.")
            profiles[key] = pd.DataFrame()
            continue

        print(f" -> files found: {len(files)}")
        for fp in files[:8]:
            try:
                rel = fp.relative_to(case_dir)
            except ValueError:
                rel = fp
            print(f"    - {rel}")
        if len(files) > 8:
            print(f"    ... {len(files)-8} more files")

        df_raw = read_plane_dataset(files)
        df_evt = filter_transits(df_raw, dz_center)
        prof = bin_profile(df_evt, args.bins, args.y_min, args.y_max, args.min_events)
        profiles[key] = prof

        csv_path = out_dir / f"profile_OpenFOAM_{key}.csv"
        raw_csv_path = out_dir / f"events_filtered_{key}.csv"
        prof.to_csv(csv_path, index=False)
        # Save filtered events for diagnostics; can be large but very useful.
        df_evt.to_csv(raw_csv_path, index=False)

        print(f" -> raw rows: {len(df_raw):8d} | filtered events: {len(df_evt):8d} | valid bins: {np.isfinite(prof['D32_m']).sum():2d}")
        print_user_summary(df_evt, key)
        print(f" -> saved profile: {csv_path}")
        print(f" -> saved filtered events: {raw_csv_path}")
        print()

    out_png = out_dir / "Lambert_OpenFOAM_Comparison.png"
    plot_all(profiles, out_png)
    print(f"Figure saved: {out_png}")


if __name__ == "__main__":
    main()
