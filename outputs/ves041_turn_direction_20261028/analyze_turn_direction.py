"""VES041 turn-direction analysis from saved local pulse metrics.

Run: python3 analyze_turn_direction.py
Requirements: numpy, pandas. No acquisition server is needed.
Negative heading = left, positive = right, following findYawVelPeaksFT.m.
The convention is traced through position -> gradient -> angular velocity.
"""
from pathlib import Path
from math import comb
import hashlib
import json
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'source_data'
SEED = 20261028
BOOT = 100_000
PERM = 200_000
ORDER = ['ves041_1_Chrim', 'ves041_2_Chrim', 'ves041_1_gfp',
         'ves041_2_gfp', 'empty_Chrim']
LABELS = dict(zip(ORDER, ['VES041 line 1 + csChrimson',
    'VES041 line 2 + csChrimson', 'Line 1 + GFP', 'Line 2 + GFP',
    'Empty + csChrimson']))
mapping_text = (SRC / 'cschrim_behaviorcrunch.R').read_text()
mapping = {int(n): g for ids, g in re.findall(
    r'Fly %in% c\(([^)]+)\) ~ "([^"]+)"', mapping_text)
    for n in re.findall(r'\d+', ids)}
assert len(mapping) == 84

d = pd.read_csv(SRC / 'presentation_metrics.csv')
d = d.loc[d.dataset.eq('ves041_activation')].copy()
d['Fly'] = d.fly.str.extract(r'fly(\d+)$').astype(int)
d['genotype'] = d.Fly.map(mapping)
assert len(d) == 6480 and d.Fly.nunique() == 81
assert not d.duplicated(['Fly', 'recording', 'pulse']).any()
assert d.groupby('Fly').size().eq(80).all()
assert d.genotype.notna().all()

for window, label in [('pre_1s', 'pre'), ('opto_1s', 'light')]:
    signed = d['signed_net_heading_deg_' + window]
    absolute = d['abs_net_heading_deg_' + window]
    assert np.isfinite(signed).all()
    assert np.allclose(signed.abs(), absolute)
    d[label + '_signed_deg'] = signed
    d[label + '_right'] = (signed > 0).astype(int)
    d[label + '_left'] = (signed < 0).astype(int)
    d[label + '_zero'] = (signed == 0).astype(int)

f = d.groupby(['Fly', 'genotype']).agg(
    n_pulses=('pulse', 'size'),
    mean_signed_light_deg=('light_signed_deg', 'mean'),
    mean_signed_pre_deg=('pre_signed_deg', 'mean'),
    right_pulses=('light_right', 'sum'), left_pulses=('light_left', 'sum'),
    zero_pulses=('light_zero', 'sum'), pre_right_pulses=('pre_right', 'sum'),
    pre_left_pulses=('pre_left', 'sum'), pre_zero_pulses=('pre_zero', 'sum'),
).reset_index()
f['light_minus_pre_deg'] = f.mean_signed_light_deg - f.mean_signed_pre_deg
f['right_fraction'] = f.right_pulses / (f.right_pulses + f.left_pulses)
f['pre_right_fraction'] = f.pre_right_pulses / (f.pre_right_pulses + f.pre_left_pulses)
assert np.isfinite(f[['mean_signed_light_deg', 'mean_signed_pre_deg',
    'light_minus_pre_deg', 'right_fraction', 'pre_right_fraction']]).all().all()
assert f.n_pulses.eq(80).all()
assert f.groupby('genotype').size().reindex(ORDER).tolist() == [18, 16, 16, 16, 15]
f.to_csv(ROOT / 'fly_turn_direction.csv', index=False)

def bootstrap_mean(values, seed):
    """Resample independent flies; never treat pulses as independent flies."""
    x = np.asarray(values, dtype=float)
    rng = np.random.default_rng(seed)
    means = x[rng.integers(0, len(x), (BOOT, len(x)))].mean(axis=1)
    lo, hi = np.quantile(means, [.025, .975])
    return dict(mean=float(x.mean()), ci95_lo=float(lo), ci95_hi=float(hi))

def sign_flip(values, seed):
    """Two-sided Monte Carlo sign-flip test of mean=0; assumes symmetry.

    Labels are signs of independent fly estimates, not individual pulses.
    This is an exploratory analysis, with a plus-one Monte Carlo correction.
    """
    x = np.asarray(values, dtype=float)
    obs = abs(x.mean())
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(PERM // 1000):
        signs = rng.integers(0, 2, (1000, len(x))) * 2 - 1
        null = (signs * x).mean(axis=1)
        hits += np.count_nonzero(np.abs(null) >= obs - 1e-12)
    return (hits + 1) / (PERM + 1)

def binomial_two_sided(k, n):
    """Exact two-sided sign-count test under P(positive fly mean)=0.5."""
    return min(1., 2 * sum(comb(n, j) for j in range(min(k, n-k) + 1)) / 2**n)

def summarize(sub, name, seed):
    light = bootstrap_mean(sub.mean_signed_light_deg, seed)
    baseline = bootstrap_mean(sub.mean_signed_pre_deg, seed + 1)
    delta = bootstrap_mean(sub.light_minus_pre_deg, seed + 2)
    fraction = bootstrap_mean(sub.right_fraction * 100, seed + 3)
    right = int((sub.mean_signed_light_deg > 0).sum())
    left = int((sub.mean_signed_light_deg < 0).sum())
    return dict(genotype=name, label=LABELS.get(name, name), n_flies=len(sub),
        n_pulses=int(sub.n_pulses.sum()),
        light_mean_deg=light['mean'], light_ci95_lo=light['ci95_lo'],
        light_ci95_hi=light['ci95_hi'],
        baseline_mean_deg=baseline['mean'], baseline_ci95_lo=baseline['ci95_lo'],
        baseline_ci95_hi=baseline['ci95_hi'],
        light_minus_pre_mean_deg=delta['mean'],
        light_minus_pre_ci95_lo=delta['ci95_lo'], light_minus_pre_ci95_hi=delta['ci95_hi'],
        right_mean_percent=fraction['mean'], right_ci95_lo=fraction['ci95_lo'],
        right_ci95_hi=fraction['ci95_hi'],
        flies_mean_right=right, flies_mean_left=left,
        right_pulses=int(sub.right_pulses.sum()), left_pulses=int(sub.left_pulses.sum()),
        zero_pulses=int(sub.zero_pulses.sum()))

summary = []
for i, g in enumerate(ORDER):
    summary.append(summarize(f.loc[f.genotype.eq(g)], g, SEED + i * 10))
activation = f.loc[f.genotype.isin(ORDER[:2])]
assert len(activation) == 34 and activation.n_pulses.sum() == 2720
pooled = summarize(activation, 'VES041 activation pooled', SEED + 100)
summary.append(pooled)
pd.DataFrame(summary).to_csv(ROOT / 'group_turn_direction.csv', index=False)

# One primary exploratory population test. Line estimates are descriptive
# decompositions of the pooled VES041 cohort, not a new collection of tests.
primary_p = sign_flip(activation.mean_signed_light_deg, SEED + 200)
sign_p = binomial_two_sided(pooled['flies_mean_right'], len(activation))
delta_p = sign_flip(activation.light_minus_pre_deg, SEED + 201)
meta = dict(analysis_date='2026-10-28', seed=SEED, bootstrap_resamples=BOOT,
    sign_flip_samples=PERM, positive_direction='right', negative_direction='left',
    sign_convention_status='Lab analysis-code convention; no new hardware calibration.',
    primary_test=dict(estimand='Mean signed light-window heading across 34 VES041 flies',
        p=primary_p, method='two-sided fly-level sign-flip, plus-one correction',
        assumption='independent flies and symmetry of fly estimates under a zero-mean null'),
    sensitivity_tests=dict(exact_fly_sign_count_p=sign_p,
        light_minus_pre_sign_flip_p=delta_p, status='descriptive sensitivity; unadjusted'),
    group_order=ORDER, labels=LABELS, summary=summary,
    time_window='First 1 s of light; last minus first unwrapped heading sample, approximately 0.999 s apart.',
    mean_definition='Arithmetic mean of signed endpoint shifts over 80 pulses within each fly, then equal weight across flies.',
    direction_fraction='Per fly, number of positive endpoint pulses / number of nonzero endpoint pulses. No turn-amplitude cutoff.',
    uncertainty='Pointwise percentile 95% bootstrap intervals across flies; no within-fly error bars are drawn.',
    exclusions='Select dataset=ves041_activation only. No movement filter, outlier removal, imputation, or new fly exclusions.',
    limitations=['Net endpoint direction within a window does not count every individual turn bout.',
        'The analysis begins with saved exports; raw acquisition processing is not rerun.',
        'Nonsignificance does not establish equivalence or prove zero directional bias.',
        'The light-versus-baseline sensitivity uses paired fly means; zero initial heading is a coordinate reference.'],
    sources={p.name:dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size)
        for p in sorted(SRC.iterdir()) if p.is_file()})
(ROOT / 'analysis_results.json').write_text(json.dumps(meta, indent=2))
print(pd.DataFrame(summary).round(4).to_string(index=False))
print('\nPrimary mean-direction p:', primary_p)
print('Exact sign-count p:', sign_p)
print('Baseline-paired sensitivity p:', delta_p)
