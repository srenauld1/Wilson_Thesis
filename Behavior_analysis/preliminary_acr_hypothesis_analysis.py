#!/usr/bin/env python3
"""Preliminary fly-level analysis of the ACR/rubyACR silencing hypothesis.

The script reads only behavior-rate datasets from *_processed_final.mat files,
computes matched 1 s pre/during-opto metrics for each pulse, averages pulses and
recordings within fly, and performs exploratory fly-level inference separately
for the unblinded and blinded-round-1 cohorts.
"""

from __future__ import annotations

import hashlib
import math
import re
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Behavior_analysis" / "preliminary_acr_results"
COHORTS = {
    "unblinded": ROOT / "unblinded RubyACR copy",
    "blinded_round1": ROOT / "round1 copy",
}
MIN_SPEED_RHO = 3.0
WALK_THRESHOLD = 1.0
WINDOW_SECONDS = 1.0
N_RESAMPLES = 20_000
SEED = 20260807


def vec(dataset: h5py.Dataset) -> np.ndarray:
    return np.asarray(dataset[...], dtype=float).reshape(-1)


def parse_fly(path: Path, cohort: str) -> dict[str, object]:
    folder = next(p for p in path.parts if re.match(r"^fly\d+_", p, re.I))
    lower = folder.lower()
    is_control = ("gfp" in lower) or ("empty" in lower)
    is_experimental = (not is_control) and ("acr" in lower)
    if "empty" in lower:
        genotype = "empty_control"
    elif "gfp" in lower:
        genotype = "gfp_control"
    elif "acr" in lower:
        genotype = "acr_experimental"
    else:
        genotype = "unclassified"
    return {
        "cohort": cohort,
        "fly": f"{cohort}:{folder}",
        "fly_folder": folder,
        "genotype_class": genotype,
        "experimental": bool(is_experimental),
        "control": bool(is_control),
    }


def window_metrics(idx: np.ndarray, heading_deg: np.ndarray, x: np.ndarray,
                   y: np.ndarray, speed: np.ndarray) -> tuple[float, ...]:
    heading = np.deg2rad(heading_deg[idx])
    moving_heading = heading[speed[idx] >= MIN_SPEED_RHO]
    moving_heading = moving_heading[np.isfinite(moving_heading)]
    rho = np.abs(np.mean(np.exp(1j * moving_heading))) if moving_heading.size >= 2 else np.nan

    finite_h = np.isfinite(heading)
    if finite_h.sum() >= 2:
        hchange = np.nansum(np.abs(np.rad2deg(np.diff(np.unwrap(heading[finite_h])))))
    else:
        hchange = np.nan

    xw, yw = x[idx], y[idx]
    good = np.isfinite(xw) & np.isfinite(yw)
    xw, yw = xw[good], yw[good]
    if xw.size >= 2:
        steps = np.hypot(np.diff(xw), np.diff(yw))
        path_length = np.sum(steps)
        displacement = np.hypot(xw[-1] - xw[0], yw[-1] - yw[0])
        straightness = displacement / path_length if path_length > 0 else np.nan
    else:
        path_length = displacement = straightness = np.nan
    return rho, hchange, straightness, path_length, displacement


def mean_or_nan(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(np.mean(x)) if x.size else np.nan


def process_file(path: Path, cohort: str) -> list[dict[str, object]]:
    meta = parse_fly(path, cohort)
    rows: list[dict[str, object]] = []
    with h5py.File(path, "r") as handle:
        group = handle["exptData"]
        needed = ["x", "y", "headingPosition", "optoStim_ds", "t_ds",
                  "totSpeed", "forwardVelocity"]
        if any(name not in group for name in needed):
            return rows
        x, y, heading, opto, time, speed, fwd = [vec(group[name]) for name in needed]

    n = len(x)
    if not all(len(a) == n for a in (y, heading, opto, speed, fwd)) or n < 3:
        return rows
    rate = (len(time) - 1) / (time[-1] - time[0])
    win = max(3, int(np.rint(WINDOW_SECONDS * rate)))
    on = opto != 0
    edges = np.diff(np.r_[False, on, False].astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1) - 1

    recording = path.stem.replace("_processed_final", "")
    for pulse, (start, end) in enumerate(zip(starts, ends), start=1):
        pre = np.arange(start - win, start)
        during = np.arange(start, start + win)
        if pre[0] < 0 or during[-1] >= n or (end - start + 1) < win:
            continue
        pre_metrics = window_metrics(pre, heading, x, y, speed)
        dur_metrics = window_metrics(during, heading, x, y, speed)
        onset_lookback = max(1, int(np.rint(0.25 * rate)))
        onset_pre = np.arange(max(0, start - onset_lookback), start)
        onset_speed = mean_or_nan(speed[onset_pre])
        onset_fwd = mean_or_nan(fwd[onset_pre])
        fwd_pre_walking = fwd[pre][fwd[pre] > WALK_THRESHOLD]
        fwd_dur_walking = fwd[during][fwd[during] > WALK_THRESHOLD]
        row = dict(meta)
        row.update({"recording": recording, "pulse": pulse, "sample_rate_hz": rate})
        for label, values in (("pre", pre_metrics), ("during", dur_metrics)):
            for metric, value in zip(
                ("rho", "hchange_deg", "straightness", "path_length", "displacement"), values
            ):
                row[f"{metric}_{label}"] = value
        row.update({
            "fwd_walking_pre": mean_or_nan(fwd_pre_walking),
            "fwd_walking_during": mean_or_nan(fwd_dur_walking),
            "fwd_all_pre": mean_or_nan(fwd[pre]),
            "fwd_all_during": mean_or_nan(fwd[during]),
            "tot_speed_pre": mean_or_nan(speed[pre]),
            "tot_speed_during": mean_or_nan(speed[during]),
            "onset_pre_250ms_tot_speed": onset_speed,
            "onset_pre_250ms_fwd_velocity": onset_fwd,
            "moving_at_onset": bool(np.isfinite(onset_speed) and onset_speed >= MIN_SPEED_RHO),
            "moving_pre_1s": bool(mean_or_nan(speed[pre]) >= MIN_SPEED_RHO),
        })
        rows.append(row)
    return rows


def bootstrap_mean_ci(x: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    x = x[np.isfinite(x)]
    if x.size < 2:
        return np.nan, np.nan
    draws = rng.choice(x, size=(N_RESAMPLES, x.size), replace=True).mean(axis=1)
    return tuple(np.quantile(draws, [0.025, 0.975]))


def sign_flip_p(x: np.ndarray, rng: np.random.Generator) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    observed = abs(x.mean())
    signs = rng.choice((-1.0, 1.0), size=(N_RESAMPLES, x.size))
    null = np.abs((signs * x).mean(axis=1))
    return float((1 + np.sum(null >= observed)) / (N_RESAMPLES + 1))


def permutation_group_p(exp: np.ndarray, ctrl: np.ndarray, rng: np.random.Generator) -> float:
    exp, ctrl = exp[np.isfinite(exp)], ctrl[np.isfinite(ctrl)]
    if exp.size == 0 or ctrl.size == 0:
        return np.nan
    observed = abs(exp.mean() - ctrl.mean())
    pooled = np.r_[exp, ctrl]
    null = np.empty(N_RESAMPLES)
    for i in range(N_RESAMPLES):
        shuffled = rng.permutation(pooled)
        null[i] = abs(shuffled[:exp.size].mean() - shuffled[exp.size:].mean())
    return float((1 + np.sum(null >= observed)) / (N_RESAMPLES + 1))


def spearman_permutation(x: np.ndarray, y: np.ndarray,
                         rng: np.random.Generator) -> tuple[float, float]:
    """Spearman rho with a two-sided Monte Carlo permutation p-value."""
    rx = pd.Series(x).rank(method="average").to_numpy(float)
    ry = pd.Series(y).rank(method="average").to_numpy(float)
    observed = float(np.corrcoef(rx, ry)[0, 1])
    null = np.empty(N_RESAMPLES)
    for i in range(N_RESAMPLES):
        null[i] = np.corrcoef(rx, rng.permutation(ry))[0, 1]
    pval = float((1 + np.sum(np.abs(null) >= abs(observed))) / (N_RESAMPLES + 1))
    return observed, pval


def hc3_interaction(model_df: pd.DataFrame, predictor: str,
                    outcome: str) -> tuple[float, float]:
    """OLS wiggliness-by-group interaction with an HC3 normal-reference p-value."""
    x = model_df[predictor].to_numpy(float)
    y = model_df[outcome].to_numpy(float)
    group = model_df.experimental.to_numpy(float)
    sd = np.std(x)
    z = (x - np.mean(x)) / sd if sd > 0 else np.zeros_like(x)
    design = np.column_stack([np.ones(len(x)), z, group, z * group])
    xtx_inv = np.linalg.pinv(design.T @ design)
    beta = xtx_inv @ design.T @ y
    residual = y - design @ beta
    leverage = np.einsum("ij,jk,ik->i", design, xtx_inv, design)
    scaled = residual**2 / np.maximum(1 - leverage, 1e-8) ** 2
    meat = design.T @ (design * scaled[:, None])
    covariance = xtx_inv @ meat @ xtx_inv
    se = math.sqrt(max(covariance[3, 3], 0))
    zstat = beta[3] / se if se > 0 else np.nan
    pval = math.erfc(abs(zstat) / math.sqrt(2)) if np.isfinite(zstat) else np.nan
    return float(beta[3]), float(pval)


def bh_adjust(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjustment, preserving missing values and row order."""
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().astype(float)
    if valid.empty:
        return out
    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out.loc[order] = np.minimum(adjusted, 1.0)
    return out


def fisher_two_sided(success_exp: int, total_exp: int,
                     success_ctrl: int, total_ctrl: int) -> float:
    """Two-sided Fisher exact p-value for a 2x2 joint-responder table."""
    total_success = success_exp + success_ctrl
    total = total_exp + total_ctrl
    lo = max(0, total_success - total_ctrl)
    hi = min(total_exp, total_success)
    denom = math.comb(total, total_success)
    probs = {
        k: math.comb(total_exp, k) * math.comb(total_ctrl, total_success - k) / denom
        for k in range(lo, hi + 1)
    }
    observed = probs[success_exp]
    return float(sum(prob for prob in probs.values() if prob <= observed + 1e-15))


def add_derived(fly: pd.DataFrame) -> pd.DataFrame:
    fly = fly.copy()
    for metric in ("rho", "hchange_deg", "straightness", "path_length", "displacement",
                   "fwd_walking", "fwd_all", "tot_speed"):
        fly[f"delta_{metric}"] = fly[f"{metric}_during"] - fly[f"{metric}_pre"]
    # All response columns point in the hypothesized direction: positive = support.
    fly["response_less_heading_change"] = -fly["delta_hchange_deg"]
    fly["response_more_straight"] = fly["delta_straightness"]
    fly["response_more_consistent"] = fly["delta_rho"]
    fly["response_faster_walking"] = fly["delta_fwd_walking"]
    fly["response_faster_all"] = fly["delta_fwd_all"]
    fly["wiggle_heading_change"] = fly["hchange_deg_pre"]
    fly["wiggle_low_straightness"] = 1 - fly["straightness_pre"]
    fly["wiggle_low_rho"] = 1 - fly["rho_pre"]
    return fly


def descriptive_tests(fly: pd.DataFrame) -> pd.DataFrame:
    outcomes = {
        "less heading change (deg)": "response_less_heading_change",
        "more path straightness": "response_more_straight",
        "more heading consistency (rho)": "response_more_consistent",
        "faster while walking": "response_faster_walking",
        "higher mean forward velocity": "response_faster_all",
    }
    records = []
    for cohort, cdf in fly.groupby("cohort", sort=False):
        for outcome_label, col in outcomes.items():
            for experimental, gdf in cdf.groupby("experimental", sort=False):
                values = gdf[col].to_numpy(float)
                values = values[np.isfinite(values)]
                salt = int(hashlib.sha1(f"{cohort}{col}{experimental}".encode()).hexdigest()[:8], 16)
                rng = np.random.default_rng(SEED + salt)
                lo, hi = bootstrap_mean_ci(values, rng)
                records.append({
                    "cohort": cohort,
                    "outcome": outcome_label,
                    "group": "experimental" if experimental else "control",
                    "n_flies": values.size,
                    "mean_hypothesis_oriented_change": np.mean(values) if values.size else np.nan,
                    "median_hypothesis_oriented_change": np.median(values) if values.size else np.nan,
                    "bootstrap_95ci_low": lo,
                    "bootstrap_95ci_high": hi,
                    "two_sided_sign_flip_p": sign_flip_p(values, rng),
                    "proportion_in_hypothesized_direction": np.mean(values > 0) if values.size else np.nan,
                })
            exp = cdf.loc[cdf.experimental, col].to_numpy(float)
            ctrl = cdf.loc[~cdf.experimental, col].to_numpy(float)
            exp_f, ctrl_f = exp[np.isfinite(exp)], ctrl[np.isfinite(ctrl)]
            salt = int(hashlib.sha1(f"group{cohort}{col}".encode()).hexdigest()[:8], 16)
            rng = np.random.default_rng(SEED + salt)
            records.append({
                "cohort": cohort,
                "outcome": outcome_label,
                "group": "experimental_minus_control",
                "n_flies": f"{exp_f.size} vs {ctrl_f.size}",
                "mean_hypothesis_oriented_change": exp_f.mean() - ctrl_f.mean(),
                "median_hypothesis_oriented_change": np.nan,
                "bootstrap_95ci_low": np.nan,
                "bootstrap_95ci_high": np.nan,
                "two_sided_sign_flip_p": permutation_group_p(exp_f, ctrl_f, rng),
                "proportion_in_hypothesized_direction": np.nan,
            })
    result = pd.DataFrame(records)
    result["bh_q_within_cohort_group"] = (
        result.groupby(["cohort", "group"], group_keys=False)["two_sided_sign_flip_p"]
        .apply(bh_adjust)
    )
    return result


def baseline_response_tests(fly: pd.DataFrame) -> pd.DataFrame:
    predictors = {
        "high pre heading change": "wiggle_heading_change",
        "low pre straightness": "wiggle_low_straightness",
        "low pre rho": "wiggle_low_rho",
    }
    outcomes = {
        "less heading change": "response_less_heading_change",
        "more straightness": "response_more_straight",
        "more heading consistency": "response_more_consistent",
        "faster while walking": "response_faster_walking",
        "higher mean forward velocity": "response_faster_all",
    }
    records = []
    for cohort, cdf in fly.groupby("cohort", sort=False):
        for pred_label, pred in predictors.items():
            for outcome_label, outcome in outcomes.items():
                for experimental, gdf in cdf.groupby("experimental", sort=False):
                    valid = gdf[[pred, outcome]].dropna()
                    if len(valid) >= 3 and valid[pred].nunique() > 1 and valid[outcome].nunique() > 1:
                        salt = int(hashlib.sha1(f"corr{cohort}{pred}{outcome}{experimental}".encode()).hexdigest()[:8], 16)
                        rho, pval = spearman_permutation(
                            valid[pred].to_numpy(float), valid[outcome].to_numpy(float),
                            np.random.default_rng(SEED + salt),
                        )
                    else:
                        rho = pval = np.nan
                    records.append({
                        "cohort": cohort, "predictor": pred_label, "outcome": outcome_label,
                        "group": "experimental" if experimental else "control",
                        "n_flies": len(valid), "spearman_rho": rho, "spearman_two_sided_p": pval,
                        "interaction_beta": np.nan, "interaction_hc3_p": np.nan,
                    })

                model_df = cdf[[pred, outcome, "experimental"]].dropna().copy()
                if len(model_df) >= 8 and model_df["experimental"].nunique() == 2:
                    beta, pval = hc3_interaction(model_df, pred, outcome)
                else:
                    beta = pval = np.nan
                records.append({
                    "cohort": cohort, "predictor": pred_label, "outcome": outcome_label,
                    "group": "experimental_vs_control_slope", "n_flies": len(model_df),
                    "spearman_rho": np.nan, "spearman_two_sided_p": np.nan,
                    "interaction_beta": beta, "interaction_hc3_p": pval,
                })
    result = pd.DataFrame(records)
    result["p_for_adjustment"] = np.where(
        result.group == "experimental_vs_control_slope",
        result.interaction_hc3_p, result.spearman_two_sided_p,
    )
    result["bh_q_within_cohort_group"] = (
        result.groupby(["cohort", "group"], group_keys=False)["p_for_adjustment"]
        .apply(bh_adjust)
    )
    return result.drop(columns="p_for_adjustment")


def response_space_counts(fly: pd.DataFrame) -> pd.DataFrame:
    d = fly.dropna(subset=["response_less_heading_change", "response_faster_walking"]).copy()
    d["faster_and_less_heading_change"] = (
        (d.response_less_heading_change > 0) & (d.response_faster_walking > 0)
    )
    result = (d.groupby(["cohort", "experimental"], as_index=False)
              .agg(n_flies=("fly", "size"), n_joint_responders=("faster_and_less_heading_change", "sum"),
                   proportion_joint_responders=("faster_and_less_heading_change", "mean")))
    result["experimental_vs_control_fisher_two_sided_p"] = np.nan
    for cohort, idx in result.groupby("cohort").groups.items():
        sub = result.loc[idx].set_index("experimental")
        pval = fisher_two_sided(
            int(sub.loc[True, "n_joint_responders"]), int(sub.loc[True, "n_flies"]),
            int(sub.loc[False, "n_joint_responders"]), int(sub.loc[False, "n_flies"]),
        )
        result.loc[idx, "experimental_vs_control_fisher_two_sided_p"] = pval
    return result


def make_plots(fly: pd.DataFrame) -> None:
    plot_script = Path(__file__).with_name("plot_preliminary_acr_results.R")
    if plot_script.exists():
        subprocess.run(["Rscript", str(plot_script), str(OUT)], check=True)


def write_report(fly: pd.DataFrame, summary: pd.DataFrame, baseline: pd.DataFrame,
                 joint: pd.DataFrame, pulse: pd.DataFrame) -> None:
    lines = [
        "# Preliminary ACR/rubyACR silencing analysis",
        "",
        "## Design and interpretation",
        "",
        "- Unit of inference: fly. Pulses and repeated recordings were averaged within fly.",
        "- Windows: 1 s immediately before opto onset and the first 1 s during opto.",
        "- Experimental: names containing ACR/rubyACR, except `empty`; controls: names containing GFP or `empty`.",
        "- Heading change is total absolute unwrapped angular displacement. Path straightness is displacement/path length. Heading consistency is circular mean resultant length (rho), excluding samples with total speed < 3.",
        "- Forward walking speed is averaged over samples with forward velocity > 1; mean forward velocity over all samples is included as a sensitivity measure.",
        "- All p-values are two-sided and exploratory. Benjamini-Hochberg q-values are provided within each cohort and test family.",
        "- Positive hypothesis-oriented changes mean less heading change, more straightness/consistency, or faster walking.",
        "",
        "## Sample sizes",
        "",
    ]
    counts = (fly.groupby(["cohort", "experimental"]).size().rename("n").reset_index())
    for _, r in counts.iterrows():
        lines.append(f"- {r.cohort}, {'experimental' if r.experimental else 'control'}: {r.n} flies")
    lines += ["", f"Usable opto pulses: {len(pulse)} across {len(fly)} flies.", "", "## Fly-level response summary", ""]
    chosen = summary[summary.group.isin(["experimental", "control", "experimental_minus_control"])].copy()
    for cohort in fly.cohort.unique():
        lines.append(f"### {cohort}")
        lines.append("")
        sub = chosen[chosen.cohort == cohort]
        for outcome in sub.outcome.unique():
            rows = sub[sub.outcome == outcome].set_index("group")
            e, c, d = rows.loc["experimental"], rows.loc["control"], rows.loc["experimental_minus_control"]
            lines.append(
                f"- {outcome}: experimental mean {e.mean_hypothesis_oriented_change:.3g} "
                f"(95% bootstrap CI {e.bootstrap_95ci_low:.3g} to {e.bootstrap_95ci_high:.3g}; "
                f"sign-flip p={e.two_sided_sign_flip_p:.3g}, BH q={e.bh_q_within_cohort_group:.3g}); control mean "
                f"{c.mean_hypothesis_oriented_change:.3g}; experimental−control "
                f"{d.mean_hypothesis_oriented_change:.3g} (permutation p={d.two_sided_sign_flip_p:.3g}, "
                f"BH q={d.bh_q_within_cohort_group:.3g})."
            )
        lines.append("")
    lines += ["## Baseline wiggliness relationships in experimental flies", ""]
    exp = baseline[baseline.group == "experimental"]
    for cohort in fly.cohort.unique():
        lines.append(f"### {cohort}")
        lines.append("")
        sub = exp[exp.cohort == cohort]
        for _, r in sub.iterrows():
            lines.append(f"- {r.predictor} → {r.outcome}: Spearman rho={r.spearman_rho:.3g}, p={r.spearman_two_sided_p:.3g}, BH q={r.bh_q_within_cohort_group:.3g}, n={int(r.n_flies)}.")
        lines.append("")
    lines += [
        "## Important caveats",
        "",
        "- This is exploratory and uses several correlated outcomes/predictors; isolated p<0.05 results should not be treated as confirmatory.",
        "- Relating a baseline metric to its own change score is vulnerable to regression to the mean. The baseline-response table therefore also includes an experimental-vs-control slope interaction (HC3 robust p-value) as a more informative check.",
        "- Averaging across pulses gives each fly equal weight but does not model pulse-to-pulse habituation or recording date/batch.",
        "- The walking-conditioned speed measure describes speed while moving forward; the all-sample forward-velocity sensitivity outcome also reflects stopping/turning.",
        "",
        "## Output files",
        "",
        "- `fly_level_metrics.csv`: analysis unit and derived responses.",
        "- `pulse_level_metrics.csv`: per-pulse audit table.",
        "- `fly_level_response_tests.csv`: directional changes and experimental-control permutation tests.",
        "- `baseline_wiggliness_response_tests.csv`: Spearman associations and group-slope interactions.",
        "- `joint_response_counts.csv`: flies that became both faster and less turning.",
        "- PNG figures show fly-level change distributions and the joint response space.",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    skipped: list[dict[str, str]] = []
    for cohort, folder in COHORTS.items():
        files = sorted(folder.rglob("*_processed_final.mat"))
        print(f"{cohort}: {len(files)} recordings", flush=True)
        for i, path in enumerate(files, start=1):
            try:
                file_rows = process_file(path, cohort)
                rows.extend(file_rows)
                if not file_rows:
                    skipped.append({"cohort": cohort, "path": str(path), "reason": "no usable pulses or missing fields"})
            except (OSError, KeyError, ValueError) as exc:
                skipped.append({"cohort": cohort, "path": str(path), "reason": repr(exc)})
                print(f"  SKIPPED {path}: {exc}", flush=True)
            if i % 5 == 0 or i == len(files):
                print(f"  processed {i}/{len(files)}", flush=True)
    pulse = pd.DataFrame(rows)
    if pulse.empty:
        raise RuntimeError("No usable opto pulses were found")
    if (~(pulse.experimental | pulse.control)).any():
        raise RuntimeError("Unclassified fly names are present")

    id_cols = ["cohort", "fly", "fly_folder", "genotype_class", "experimental", "control"]
    numeric_cols = [c for c in pulse.select_dtypes(include=[np.number]).columns if c != "pulse"]
    fly = pulse.groupby(id_cols, as_index=False)[numeric_cols].mean()
    counts = pulse.groupby("fly").agg(n_pulses=("pulse", "size"), n_recordings=("recording", "nunique"))
    fly = fly.merge(counts, left_on="fly", right_index=True, how="left")
    fly = add_derived(fly)

    summary = descriptive_tests(fly)
    baseline = baseline_response_tests(fly)
    joint = response_space_counts(fly)

    pulse.to_csv(OUT / "pulse_level_metrics.csv", index=False)
    pd.DataFrame(skipped, columns=["cohort", "path", "reason"]).to_csv(
        OUT / "skipped_recordings.csv", index=False
    )
    fly.to_csv(OUT / "fly_level_metrics.csv", index=False)
    summary.to_csv(OUT / "fly_level_response_tests.csv", index=False)
    baseline.to_csv(OUT / "baseline_wiggliness_response_tests.csv", index=False)
    joint.to_csv(OUT / "joint_response_counts.csv", index=False)
    make_plots(fly)
    write_report(fly, summary, baseline, joint, pulse)
    print(f"Wrote results to {OUT}")


if __name__ == "__main__":
    main()
