#!/usr/bin/env python3
"""Extract opto-aligned pivot-turn metrics for SMP544 inhibition and VES041 activation."""

from pathlib import Path
import re

import h5py
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "Behavior_analysis" / "smp544_pivot_results"
SMP_ROOT = Path("/Volumes/Neurobio/wilsonlab/Sophia/behavior/SMP544 unblinded data")
SMP_MAP = {
    1: ("2026_09_17", "fly01", "experimental"),
    2: ("2026_09_17", "fly02", "control"),
    3: ("2026_09_18", "fly01", "experimental"),
    5: ("2026_09_18", "fly03", "experimental"),
    8: ("2026_09_25", "fly01", "experimental"),
    10: ("2026_09_30", "fly01", "control"),
    12: ("2026_09_30", "fly03", "control"),
}
VES_ROOTS = {
    "ves041_activation": Path("/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers"),
}
VES_ACTIVATION_EXPERIMENTAL_IDS = {
    2, 7, 16, 19, 21, 26, 53, 58, 65, 70, 72, 74, 77, 79, 81, 85, 119, 130,
    14, 18, 24, 29, 40, 44, 66, 73, 75, 80, 90, 100, 102, 105, 108, 110,
}
VES_ACTIVATION_CONTROL_IDS = {
    3, 4, 5, 8, 12, 15, 25, 27, 69, 99, 101, 103, 111, 114, 121, 124,
    35, 52, 59, 63, 89, 93, 107, 113, 118, 122, 123, 125, 129, 132, 133, 137,
    17, 20, 23, 31, 34, 39, 46, 50, 54, 57, 62, 84, 87, 88, 91, 96,
}
WIN_SECONDS = 1.0


def vector(group, name):
    return np.asarray(group[name][...], dtype=float).reshape(-1)


def window_values(idx, heading_deg, x, y, angular_velocity, fwd, speed):
    heading = np.unwrap(np.deg2rad(heading_deg[idx]))
    heading_diff = np.rad2deg(np.diff(heading))
    signed_net = float(np.rad2deg(heading[-1] - heading[0]))
    total_heading = float(np.nansum(np.abs(heading_diff)))
    directionality = abs(signed_net) / total_heading if total_heading > 0 else np.nan
    xw, yw = x[idx], y[idx]
    steps = np.hypot(np.diff(xw), np.diff(yw))
    path_length = float(np.nansum(steps))
    displacement = float(np.hypot(xw[-1] - xw[0], yw[-1] - yw[0]))
    av = angular_velocity[idx]
    av = av[np.isfinite(av)]
    return {
        "signed_net_heading_deg": signed_net,
        "abs_net_heading_deg": abs(signed_net),
        "total_heading_change_deg": total_heading,
        "turn_directionality": directionality,
        "mean_abs_angular_velocity": float(np.mean(np.abs(av))) if av.size else np.nan,
        "p95_abs_angular_velocity": float(np.quantile(np.abs(av), .95)) if av.size else np.nan,
        "peak_abs_angular_velocity": float(np.max(np.abs(av))) if av.size else np.nan,
        "path_length": path_length,
        "displacement": displacement,
        "mean_forward_velocity": float(np.nanmean(fwd[idx])),
        "mean_total_speed": float(np.nanmean(speed[idx])),
    }


def process_file(path, dataset, fly, group):
    try:
        with h5py.File(path, "r") as handle:
            g = handle["exptData"]
            required = ["x", "y", "headingPosition", "optoStim_ds", "t_ds",
                        "totSpeed", "forwardVelocity", "smoothedangularVelocity"]
            if any(name not in g for name in required):
                return []
            x, y, heading, opto, time, speed, fwd, angular_velocity = [vector(g, k) for k in required]
    except OSError:
        # One legacy Chrimson control file is pre-v7.3 rather than HDF5.
        print(f"Skipping non-HDF5 MAT file: {path}")
        return []
    n = len(x)
    if not all(len(v) == n for v in [y, heading, opto, speed, fwd, angular_velocity]):
        return []
    rate = (len(time) - 1) / (time[-1] - time[0])
    win = max(3, int(round(WIN_SECONDS * rate)))
    half = max(3, int(round(.5 * rate)))
    active = opto != 0
    edges = np.diff(np.r_[False, active, False].astype(np.int8))
    onsets = np.flatnonzero(edges == 1)
    offsets = np.flatnonzero(edges == -1) - 1
    rows = []
    for pulse, (start, end) in enumerate(zip(onsets, offsets), start=1):
        if start - win < 0 or start + win > n or end - start + 1 < win:
            continue
        windows = {
            "pre_1s": np.arange(start - win, start),
            "opto_1s": np.arange(start, start + win),
            "pre_500ms": np.arange(start - half, start),
            "opto_500ms": np.arange(start, start + half),
        }
        row = {
            "dataset": dataset, "fly": fly, "group": group,
            "recording": path.stem.replace("_processed_final", ""),
            "pulse": pulse, "sample_rate_hz": rate,
        }
        for label, idx in windows.items():
            for metric, value in window_values(idx, heading, x, y, angular_velocity, fwd, speed).items():
                row[f"{metric}_{label}"] = value
        # Positive values indicate an opto-evoked increase.
        for metric in ["abs_net_heading_deg", "total_heading_change_deg",
                       "mean_abs_angular_velocity", "p95_abs_angular_velocity",
                       "peak_abs_angular_velocity", "path_length", "displacement",
                       "mean_forward_velocity", "mean_total_speed"]:
            row[f"delta_{metric}_1s"] = row[f"{metric}_opto_1s"] - row[f"{metric}_pre_1s"]
            row[f"delta_{metric}_500ms"] = row[f"{metric}_opto_500ms"] - row[f"{metric}_pre_500ms"]
        row["pivot_index_1s"] = row["abs_net_heading_deg_opto_1s"] / (row["displacement_opto_1s"] + .5)
        row["pivot_index_500ms"] = row["abs_net_heading_deg_opto_500ms"] / (row["displacement_opto_500ms"] + .25)
        rows.append(row)
    return rows


def process_smp_traces(path, fly, group):
    """Return -1 to +1.5 s traces for visual validation of pivot-like events."""
    with h5py.File(path, "r") as handle:
        g = handle["exptData"]
        x, y, heading, opto, time, speed, fwd, angular_velocity = [
            vector(g, k) for k in ["x", "y", "headingPosition", "optoStim_ds", "t_ds",
                                  "totSpeed", "forwardVelocity", "smoothedangularVelocity"]
        ]
    rate = (len(time) - 1) / (time[-1] - time[0])
    pre = int(round(rate))
    post = int(round(1.5 * rate))
    active = opto != 0
    edges = np.diff(np.r_[False, active, False].astype(np.int8))
    onsets = np.flatnonzero(edges == 1)
    rows = []
    for pulse, start in enumerate(onsets, start=1):
        if start - pre < 0 or start + post >= len(x):
            continue
        idx = np.arange(start - pre, start + post + 1)
        heading_unwrapped = np.rad2deg(np.unwrap(np.deg2rad(heading[idx])))
        for j, sample in enumerate(idx):
            rows.append({
                "fly": fly, "group": group, "recording": path.stem.replace("_processed_final", ""),
                "pulse": pulse, "time_from_onset_s": (sample - start) / rate,
                "x_relative": x[sample] - x[start], "y_relative": y[sample] - y[start],
                "heading_relative_deg": heading_unwrapped[j] - heading_unwrapped[pre],
                "angular_velocity": angular_velocity[sample],
                "forward_velocity": fwd[sample], "total_speed": speed[sample],
                "opto_on": bool(opto[sample] != 0),
            })
    return rows


def ves_activation_fly_id(path):
    for parent in path.parents:
        match = re.fullmatch(r"fly0*(\d+)", parent.name, re.I)
        if match:
            return int(match.group(1))
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    traces = []
    manifest = []
    for fly_id, (date, folder, group) in SMP_MAP.items():
        files = sorted((SMP_ROOT / date / folder).rglob("*_processed_final.mat"))
        for path in files:
            rows.extend(process_file(path, "smp544", f"smp544_{fly_id:02d}", group))
            traces.extend(process_smp_traces(path, f"smp544_{fly_id:02d}", group))
            manifest.append({"dataset": "smp544", "fly": f"smp544_{fly_id:02d}",
                             "group": group, "path": str(path)})
    for dataset, root in VES_ROOTS.items():
        for path in sorted(root.rglob("*_processed_final.mat")):
            fly_id = ves_activation_fly_id(path)
            if fly_id in VES_ACTIVATION_EXPERIMENTAL_IDS:
                group = "experimental"
            elif fly_id in VES_ACTIVATION_CONTROL_IDS:
                group = "control"
            else:
                continue
            fly = f"{dataset}:fly{fly_id:03d}"
            rows.extend(process_file(path, dataset, fly, group))
            manifest.append({"dataset": dataset, "fly": fly, "group": group, "path": str(path)})
    data = pd.DataFrame(rows)
    data.to_csv(OUT / "presentation_metrics.csv", index=False)
    pd.DataFrame(traces).to_csv(OUT / "smp544_opto_aligned_traces.csv", index=False)
    pd.DataFrame(manifest).to_csv(OUT / "recording_manifest.csv", index=False)
    print(data.groupby(["dataset", "group"]).agg(presentations=("pulse", "size"), flies=("fly", "nunique")))


if __name__ == "__main__":
    main()
