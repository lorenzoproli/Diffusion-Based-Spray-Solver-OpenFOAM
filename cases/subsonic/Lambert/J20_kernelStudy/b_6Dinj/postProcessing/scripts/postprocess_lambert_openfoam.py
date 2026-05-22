#!/usr/bin/env python3
"""
Lambert-style post-processing for OpenFOAM ParcelPlaneSampler data — J=20 case.

Reads crossing-event files written by the custom OpenFOAM ParcelPlaneSampler on
planes x/d = 30 and x/d = 60, computes nParticle-weighted profiles of:
    D32(y) = sum(n_i d_i^3) / sum(n_i d_i^2)
    Ux(y)  = sum(n_i Ux_i)  / sum(n_i)
and compares them against Lambert / Sekar et al. reference data (J=20).

Usage:
    python3 postprocess_lambert_openfoam.py --case ../b_0_PCM_cell
"""
from __future__ import annotations

from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
D_JET = 4.57e-4        # [m] jet diameter J20
DZ_CENTER_DEFAULT = 0.001
N_BINS_Y_DEFAULT = 30
Y_MIN_DEFAULT = 0.0
Y_MAX_DEFAULT = 0.015
MIN_EVENTS_DEFAULT = 500

PLANES = {"30d": {"x_nd": 30}, "60d": {"x_nd": 60}}

# ---------------------------------------------------------------------------
# Reference data — Lambert J=20 — digitized from Lambert (ILASS 2019)
# Format: [y_mm, value]  (value = D32 [µm] or Ux [m/s])
# ---------------------------------------------------------------------------
REF = {
    "30d": {
        "fluent_smd": np.array([
            [0.597, 12.098], [0.995, 13.238], [1.368, 13.727], [1.666, 14.297],
            [2.288, 15.397], [2.810, 16.130], [3.382, 16.497], [3.680, 16.864],
            [4.178, 17.760], [4.526, 18.330], [5.023, 19.389], [5.098, 19.633],
            [5.446, 20.733], [5.570, 21.018], [5.694, 21.426], [6.142, 22.851],
            [6.639, 24.603], [6.888, 25.377], [7.137, 26.069], [7.311, 26.640],
            [7.510, 27.128], [7.883, 27.862], [8.554, 28.228], [9.052, 27.454],
            [9.524, 26.069], [9.897, 25.255], [10.071, 25.010], [10.345, 24.888],
            [10.593, 24.807], [10.817, 24.807], [11.364, 24.807], [11.861, 24.888],
            [12.309, 25.010], [12.657, 24.888], [12.931, 24.847], [13.229, 24.358],
        ]),
        "exp_smd": np.array([
            [2.014, 14.908], [3.009, 16.578], [4.004, 18.126], [5.023, 20.652],
            [6.018, 21.711], [7.037, 24.073], [8.007, 24.521], [9.002, 23.870],
            [9.996, 23.544], [11.016, 25.336], [12.011, 25.743],
        ]),
        "fluent_ux": np.array([
            [0.634, 100.759], [1.014, 100.169], [1.369,  98.841], [1.701,  96.481],
            [2.270,  90.875], [2.838,  83.646], [3.406,  76.122], [3.667,  72.582],
            [4.212,  67.861], [4.543,  65.943], [5.042,  65.943], [5.421,  66.828],
            [5.635,  66.976], [5.778,  67.566], [6.158,  69.336], [6.632,  71.254],
            [6.894,  72.139], [7.107,  73.172], [7.321,  74.204], [7.487,  75.237],
            [7.891,  77.155], [8.556,  80.105], [9.078,  80.991], [9.529,  81.581],
            [9.908,  82.171], [10.122, 82.613], [10.360, 83.203], [10.597, 83.793],
            [10.858, 84.679], [11.381, 85.859], [11.855, 86.449], [12.306, 86.596],
            [12.686, 87.334], [12.947, 87.187], [13.232, 87.334],
        ]),
        "exp_ux": np.array([
            [2.031,  79.958], [3.026,  74.204], [4.021,  65.353], [5.040,  57.387],
            [6.013,  57.829], [7.059,  68.894], [8.058,  79.958], [9.008,  87.482],
            [10.053, 92.940], [11.026, 91.465], [12.046, 88.957],
        ]),
    },
    "60d": {
        "fluent_smd": np.array([
            [1.023,  14.192], [1.406,  14.573], [1.760,  14.741], [2.203,  15.248],
            [2.675,  15.840], [2.852,  15.840], [3.029,  15.924], [3.206,  16.178],
            [3.501,  16.474], [3.914,  16.812], [4.179,  16.980], [4.415,  17.107],
            [4.652,  17.149], [4.858,  17.318], [5.124,  17.402], [5.448,  17.697],
            [5.684,  17.993], [5.979,  18.374], [6.392,  19.093], [6.834,  19.855],
            [7.276,  20.871], [7.600,  21.760], [8.130,  23.241], [8.484,  24.427],
            [8.749,  25.273], [8.955,  25.824], [9.191,  26.459], [9.485,  27.221],
            [9.899,  27.431], [10.194, 27.515], [10.430, 27.515], [10.755, 27.387],
            [11.139, 27.132], [11.316, 26.792], [11.700, 26.071], [11.878, 25.816],
            [12.114, 25.731], [12.291, 25.603], [12.468, 25.518], [12.616, 25.390],
            [13.000, 24.966], [13.236, 24.923], [13.354, 24.922], [13.531, 24.964],
            [13.945, 25.005], [14.564, 25.343], [14.978, 25.299], [15.273, 25.044],
        ]),
        "exp_smd": np.array([
            [1.997,  14.613], [2.970,  16.136], [4.001,  18.464], [5.004,  19.351],
            [5.977,  21.467], [6.979,  23.414], [8.012,  24.174], [8.985,  25.358],
            [9.989,  24.211], [10.992, 25.649], [11.966, 26.367], [12.999, 26.745],
            [13.942, 27.887], [14.976, 28.096],
        ]),
        "fluent_ux": np.array([
            [1.029,  102.461], [1.435,  102.019], [1.841,  100.988], [2.247,   99.516],
            [2.734,   97.308], [2.870,   96.719], [3.059,   95.247], [3.303,   93.922],
            [3.547,   92.303], [3.926,   89.211], [4.169,   87.150], [4.439,   84.206],
            [4.711,   81.851], [4.927,   80.526], [5.144,   79.054], [5.279,   78.170],
            [5.523,   76.698], [5.739,   75.668], [6.037,   74.784], [6.389,   74.343],
            [6.876,   74.343], [7.228,   75.373], [7.580,   76.257], [8.095,   78.170],
            [8.501,   79.642], [8.772,   80.673], [8.961,   81.262], [9.178,   82.292],
            [9.530,   84.059], [9.990,   85.089], [10.152,  85.678], [10.288,  86.414],
            [10.558,  87.150], [10.829,  88.475], [11.127,  89.211], [11.371,  90.389],
            [11.723,  92.156], [11.831,  92.450], [12.318,  93.333], [12.562,  94.217],
            [13.022,  95.247], [13.320,  95.247], [13.536,  95.394], [13.997,  95.836],
            [14.619,  95.689], [15.242,  96.130],
        ]),
        "exp_ux": np.array([
            [1.976,  91.420], [3.005,  90.978], [3.980,  86.267], [5.008,  82.881],
            [5.983,  78.023], [7.039,  80.820], [8.014,  83.617], [9.015,  88.623],
            [10.017, 95.100], [10.992, 97.897], [11.993, 98.927], [12.995, 98.044],
            [14.024, 97.014], [14.998, 95.983],
        ]),
    },
}

# ---------------------------------------------------------------------------
# Sampler column formats
# ---------------------------------------------------------------------------
OLD_COLS = ["time", "x", "y", "z", "Ux", "Uy", "Uz", "d", "nParticle", "faceID"]
NEW_COLS = [
    "time", "x", "y", "z", "Ux", "Uy", "Uz", "d", "nParticle",
    "user", "ms", "tc", "KHindex", "yState", "yDotState", "faceID",
]
PLANE_TOKENS = {
    "30d": ["parcelplanesampler_30d", "sampler_30d", "plane_30d", "30d"],
    "60d": ["parcelplanesampler_60d", "sampler_60d", "plane_60d", "60d"],
}
VALID_SUFFIXES = {".dat", ".csv", ".smp", ".txt"}


def _time_sort_key(path: Path) -> tuple:
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


def find_plane_files(case_dir: Path, key: str, explicit: str | None, out_dir_name: str) -> list[Path]:
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
    matches: list[Path] = []
    for p in search_root.rglob("*"):
        if not p.is_file() or p.suffix.lower() not in VALID_SUFFIXES:
            continue
        full = str(p).lower()
        if out_dir_name and out_dir_name.lower() in full:
            continue
        if any(tok in full for tok in PLANE_TOKENS[key]):
            matches.append(p)
    seen: set[Path] = set()
    out: list[Path] = []
    for p in sorted(matches, key=_time_sort_key):
        if p not in seen:
            out.append(p)
            seen.add(p)
    return out


def read_openfoam_sampler(path: Path) -> pd.DataFrame:
    n_fields = _first_data_field_count(path)
    if n_fields >= len(NEW_COLS):
        names = NEW_COLS
    elif n_fields == len(OLD_COLS):
        names = OLD_COLS
    elif n_fields > 0:
        names = NEW_COLS[:n_fields]
    else:
        return pd.DataFrame()
    df = pd.read_csv(
        path, sep=r"\s+|,", comment="#", header=None, names=names, engine="python", on_bad_lines="skip",
    )
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    required = ["y", "Ux", "d", "nParticle"]
    return df.dropna(subset=[c for c in required if c in df.columns])


def read_plane_dataset(paths: list[Path]) -> pd.DataFrame:
    frames = [read_openfoam_sampler(p) for p in paths]
    frames = [f for f in frames if not f.empty]
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


def _weighted_stats(d, ux, n):
    if len(d) == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    snd2 = np.sum(n * d**2)
    snd3 = np.sum(n * d**3)
    sn = np.sum(n)
    small = d < 2e-6
    d32 = snd3 / snd2 if snd2 > 0 else np.nan
    ux_mean = np.sum(n * ux) / sn if sn > 0 else np.nan
    return d32, ux_mean, sn, float(np.sum(small)) / len(d), np.sum(n[small] * d[small]**3) / (snd3 + 1e-300)


def bin_profile(df: pd.DataFrame, n_bins: int, y_min: float, y_max: float, min_events: int) -> pd.DataFrame:
    edges = np.linspace(y_min, y_max, n_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    base = {
        "y_center": centers, "y_center_mm": centers * 1e3,
        "D32_m": np.nan, "D32_um": np.nan, "Ux_mean": np.nan,
        "n_events": 0, "n_phys_sum": np.nan,
        "frac_events_d_lt_2um": np.nan, "frac_mass_d_lt_2um": np.nan,
    }
    for u in [0, 1, 2]:
        base[f"n_events_user{u}"] = 0
        base[f"mass_fraction_user{u}"] = np.nan
        base[f"D32_um_user{u}"] = np.nan
        base[f"Ux_mean_user{u}"] = np.nan
    for lbl, _ in [("m10", -10), ("m20", -20), ("m2", 2)]:
        base[f"n_events_{lbl}"] = 0
        base[f"mass_fraction_{lbl}"] = np.nan

    out = pd.DataFrame(base)
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
        if int(sel.sum()) < min_events:
            continue
        dk, uk, nk = d[sel], ux[sel], n[sel]
        d32, ux_mean, sn, fe_small, fm_small = _weighted_stats(dk, uk, nk)
        mass_total = np.sum(nk * dk**3)
        out.at[k, "n_events"] = int(sel.sum())
        out.at[k, "n_phys_sum"] = sn
        out.at[k, "frac_events_d_lt_2um"] = fe_small
        out.at[k, "frac_mass_d_lt_2um"] = fm_small
        out.at[k, "D32_m"] = d32
        out.at[k, "D32_um"] = 1e6 * d32 if np.isfinite(d32) else np.nan
        out.at[k, "Ux_mean"] = ux_mean

        if has_user and user is not None:
            uk_user = user[sel]
            for u in [0, 1, 2]:
                sub = uk_user == u
                out.at[k, f"n_events_user{u}"] = int(np.sum(sub))
                if np.any(sub):
                    d32u, uxu, *_ = _weighted_stats(dk[sub], uk[sub], nk[sub])
                    out.at[k, f"D32_um_user{u}"] = 1e6 * d32u if np.isfinite(d32u) else np.nan
                    out.at[k, f"Ux_mean_user{u}"] = uxu
                    out.at[k, f"mass_fraction_user{u}"] = np.sum(nk[sub] * dk[sub]**3) / (mass_total + 1e-300)

        if has_ms and ms is not None:
            mk = ms[sel]
            for lbl, mv in [("m10", -10), ("m20", -20), ("m2", 2)]:
                sub = mk == mv
                out.at[k, f"n_events_{lbl}"] = int(np.sum(sub))
                if np.any(sub):
                    out.at[k, f"mass_fraction_{lbl}"] = np.sum(nk[sub] * dk[sub]**3) / (mass_total + 1e-300)

    return out


def _plot_panel(ax, prof: pd.DataFrame, ref: dict, quantity: str, xlabel: str, title: str) -> None:
    ref_f = ref.get("fluent_smd" if quantity == "D32_m" else "fluent_ux")
    ref_e = ref.get("exp_smd" if quantity == "D32_m" else "exp_ux")
    if ref_f is not None:
        ax.plot(ref_f[:, 1], ref_f[:, 0], "b-", lw=1.6, label="FLUENT (Lambert 2019)")
    if ref_e is not None:
        ax.plot(ref_e[:, 1], ref_e[:, 0], "ks", ms=5, label="EXP (Sekar et al. 2014)")
    if prof is not None and not prof.empty:
        x = prof["D32_um" if quantity == "D32_m" else "Ux_mean"].to_numpy()
        y_mm = prof["y_center_mm"].to_numpy()
        ok = np.isfinite(x)
        ax.plot(x[ok], y_mm[ok], "ro--", lw=1.4, ms=4, label="OpenFOAM")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("y [mm]")
    ax.set_title(title)
    _cands = []
    if ref_f is not None and len(ref_f) > 0:
        _cands.append(float(np.nanmax(ref_f[:, 1])))
    if ref_e is not None and len(ref_e) > 0:
        _cands.append(float(np.nanmax(ref_e[:, 1])))
    if prof is not None and not prof.empty:
        _col = "D32_um" if quantity == "D32_m" else "Ux_mean"
        _v = pd.to_numeric(prof[_col], errors="coerce").max()
        if np.isfinite(_v):
            _cands.append(float(_v))
    _default_max = 40.0 if quantity == "D32_m" else 120.0
    _x_max = max([_default_max] + [c * 1.15 for c in _cands]) if _cands else _default_max
    ax.set_xlim(0, _x_max)
    ax.set_ylim(0, 15 if "30" in title else 16)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")


def plot_all(profiles: dict, case_label: str, out_png: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"J=20 — {case_label}", fontsize=13)
    for key, row in [("30d", 0), ("60d", 1)]:
        x_nd = PLANES[key]["x_nd"]
        prof = profiles.get(key, pd.DataFrame())
        _plot_panel(axes[row, 0], prof, REF[key], "D32_m", "D32 [µm]", f"SMD — x/d = {x_nd}")
        _plot_panel(axes[row, 1], prof, REF[key], "Ux_mean", "Ux [m/s]", f"Mean Ux — x/d = {x_nd}")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, default=Path.cwd())
    parser.add_argument("--plane30", type=str, default=None)
    parser.add_argument("--plane60", type=str, default=None)
    parser.add_argument("--dz", type=float, default=DZ_CENTER_DEFAULT)
    parser.add_argument("--min-events", type=int, default=MIN_EVENTS_DEFAULT)
    parser.add_argument("--bins", type=int, default=N_BINS_Y_DEFAULT)
    parser.add_argument("--y-min", type=float, default=Y_MIN_DEFAULT)
    parser.add_argument("--y-max", type=float, default=Y_MAX_DEFAULT)
    args = parser.parse_args()

    case_dir = args.case.resolve()
    dz_center = None if args.dz < 0 else args.dz
    out_dir = case_dir / "postProcessing_Lambert_OpenFOAM"
    out_dir.mkdir(parents=True, exist_ok=True)

    explicit = {"30d": args.plane30, "60d": args.plane60}
    profiles: dict[str, pd.DataFrame] = {}

    print(f"Case directory: {case_dir}")
    print(f"Midline filter: {'disabled' if dz_center is None else f'|z| <= {dz_center*1e3:.2f} mm'}")
    print(f"MIN_EVENTS per bin: {args.min_events}\n")

    for key, info in PLANES.items():
        print(f"Processing plane {key} (x/d={info['x_nd']}) ...")
        files = find_plane_files(case_dir, key, explicit[key], out_dir.name)
        if not files:
            print(" -> WARNING: no file found; plotting references only.")
            profiles[key] = pd.DataFrame()
            continue
        print(f" -> {len(files)} file(s) found")

        df_raw = read_plane_dataset(files)
        df_evt = filter_transits(df_raw, dz_center)
        prof = bin_profile(df_evt, args.bins, args.y_min, args.y_max, args.min_events)
        profiles[key] = prof
        prof.to_csv(out_dir / f"profile_OpenFOAM_{key}.csv", index=False)
        df_evt.to_csv(out_dir / f"events_filtered_{key}.csv", index=False)
        print(f" -> raw: {len(df_raw):8d} | filtered: {len(df_evt):8d} | valid bins: {np.isfinite(prof['D32_m']).sum():2d}")

    plot_all(profiles, case_dir.name, out_dir / "Lambert_OpenFOAM_Comparison.png")
    print(f"\nFigure and CSV saved in: {out_dir}")


if __name__ == "__main__":
    main()
