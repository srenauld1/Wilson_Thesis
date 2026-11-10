"""Fly-level comparisons for native-rate imaging-defined inhibition events."""
from pathlib import Path
import json
import numpy as np
from analyze_ves041_inhibition_pivot_frequency import aggregate, inference

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/ves041_inhibition_imaging_rules_20261104'

def main():
    status=json.loads((OUT/'extraction_status.json').read_text())
    assert len(status)==28 and all(r['status']=='complete' for r in status)
    records=[];trials=[];events=[];examples=[]
    for path in sorted((OUT/'recordings').glob('*.json')):
        a=json.loads(path.read_text());r=a['record'];records.append({**r,**a['metadata']})
        for row in a['presentations']:
            trials.append({**r,**row})
            # Decompose the result to distinguish fewer slowdowns from a
            # smaller fraction of slowdowns accompanied by rotation.
            if row['definition']=='primary':
                trials.append({**r,**row,'definition':'all_slowdowns','events':row['all_slowdowns']})
                trials.append({**r,**row,'definition':'slowdowns_without_rotation','events':row['all_slowdowns']-row['events']})
        if isinstance(a['events'],dict):a['events']=[a['events']]
        for e in a['events']:events.append({**r,**e})
        if isinstance(a['examples'],dict):a['examples']=[a['examples']]
        for e in a['examples']:examples.append({**r,**e})
    assert len(records)==28 and sum(r['quality_subset'] for r in records)==17
    flies=aggregate(trials)
    assert all(r['n_presentations']==40 for r in flies if r['subset'] in ['all','quality_subset'])
    groups=inference(flies)
    rng=np.random.default_rng(20261104)
    for s in groups:
        selected=[r for r in flies if r['subset']==s['subset'] and r['definition']==s['definition']]
        ex=np.array([r['delta_rate_per_min'] for r in selected if r['group']=='experimental'])
        co=np.array([r['delta_rate_per_min'] for r in selected if r['group']=='control'])
        samples=ex[rng.integers(len(ex),size=(20000,len(ex)))].mean(axis=1)-co[rng.integers(len(co),size=(20000,len(co)))].mean(axis=1)
        s['difference_in_changes_bootstrap95']=np.quantile(samples,[.025,.975]).tolist()
    for name,data in [('manifest',records),('presentations',trials),('events',events),('examples',examples),('fly_summary',flies),('group_summary',groups)]:
        (OUT/f'{name}.json').write_text(json.dumps(data,indent=2,allow_nan=False))
    lines=[
      'VES041 INHIBITION: APPLYING THE IMAGING SLOWDOWN/REORIENTATION RULES',
      'November 4, 2026','',
      'COHORT: Same blinded round 1, 5-second recordings as the prior inhibition analysis. All 28 readable, genotype-assigned flies (11 experimental, 17 control); fixed good-walker subset 8 experimental/9 control. Forty presentations per fly. Unassigned genotypes and damaged ID15 remain excluded; no new behavioral-response exclusions.',
      'DETECTOR: Native forwardVelocity and angularVelocity, using the original MATLAB smoothdata movmedian over 0.25 s. Forward acceleration <=-2 mm/s^2 defines contiguous candidate episodes; align to the strongest deceleration in each episode. Median smoothed forward velocity -0.75:-0.25 s must be >=1 mm/s, and must fall >=1 mm/s relative to +0.25:+0.75 s. Accepted slowdowns spaced >=2 s chronologically. This spacing includes slowdowns without rotation, matching the original imaging pipeline.',
      'ROTATION: Integral of absolute 0.25-second median-smoothed yaw over 0:+0.5 s >=15 degrees. Not net heading, not a 30-degree rule, no displacement or directionality cutoff. Fast integral implementation verified against pizza_ves041_turn_angle on five actual events per recording.',
      'VALIDITY: Complete -2:+3 s event windows with finite kinematics and no time gaps. No fluorescence is measured in inhibition experiments, so the imaging fluorescence-validity requirement is inapplicable. Independent imaging-specific behaviorVolumeValid flags are also unavailable; native finite behavior/time checks are applied, not invented tracking-quality flags. No new smoothing of displayed event traces.',
      'ASSIGNMENT: Detect once over each continuous recording; assign the strongest-deceleration anchor to the full 5 s opto period or matched immediately preceding 5 s light-off period. The detection/display window can straddle light transitions, just as imaging windows can straddle heat transitions. Counting per phase does NOT require the whole -2:+3 s window to lie in that phase.',
      'BOUNDARY CHECK: within_phase_support restricts anchors using symmetric edge buffers of 0.875 s plus one native sample. This keeps the speed-comparison/turn windows and their median-smoothing support within the assigned phase; exposure is reduced equally in pre/opto periods. The full -2:+3 s validity window can still cross phase boundaries.',
      'RATES: Sum events and matched exposure within each fly, events per minute. Group means weight flies equally. No moving-at-onset restriction in the primary comparison; moving-onset subsets are descriptive.',
      'INFERENCE: Two-sided seeded permutation tests on fly-level rates and opto-minus-pre changes; within-group sign-flip tests on one change per fly. Exact enumeration where feasible, otherwise 49,999 resamples, matching the earlier inhibition analysis. BH adjustments separate between-group versus within-group tests, and all flies versus good walkers. Bootstrap intervals resample fly changes within genotype groups (20,000 resamples, seed 20261104); pointwise percentile intervals. Boundary sensitivity is descriptive. These are exploratory re-analyses of the same flies, not an independent replication.','',
      'RESULTS: IMAGING-DEFINED SLOWDOWNS WITH REORIENTATION'
    ]
    for subset in ['all','quality_subset']:
        s=next(s for s in groups if s['subset']==subset and s['definition']=='primary')
        lines+=['',subset]
        for group in ['experimental','control']:
            lines.append(f'{group}: n={s[f"n_{group}"]}; pre={s[f"{group}_pre_rate_per_min"]:.4f}; opto={s[f"{group}_inhibition_rate_per_min"]:.4f}; change={s[f"{group}_delta_rate_per_min"]:.4f} events/min; decreasing flies={s[f"{group}_flies_decrease"]}; increasing flies={s[f"{group}_flies_increase"]}')
        lines.append(f'Experimental-control difference in changes: {s["difference_in_changes"]:.4f} events/min; bootstrap 95% interval {s["difference_in_changes_bootstrap95"]}')
        for k,v in s.items():
            if k.endswith('_p') or k.endswith('_bh_q'):lines.append(f'{k}: {v:.6g}')
        support=next(x for x in groups if x['subset']==subset and x['definition']=='within_phase_support')
        lines.append(f'Within-phase-support check: experimental change={support["experimental_delta_rate_per_min"]:.4f}; control={support["control_delta_rate_per_min"]:.4f}; difference={support["difference_in_changes"]:.4f} events/min; bootstrap95={support["difference_in_changes_bootstrap95"]}')
        slow=next(x for x in groups if x['subset']==subset and x['definition']=='all_slowdowns')
        lines.append(f'All slowdowns, irrespective of turn: experimental pre/opto={slow["experimental_pre_rate_per_min"]:.4f}/{slow["experimental_inhibition_rate_per_min"]:.4f}; control={slow["control_pre_rate_per_min"]:.4f}/{slow["control_inhibition_rate_per_min"]:.4f} events/min')
    lines+=['','INDIVIDUAL FLIES: primary definition, all presentations']
    for r in sorted((x for x in flies if x['definition']=='primary' and x['subset']=='all'),key=lambda x:(x['group'],x['date'],x['daily_fly'])):
        lines.append(f'ID {",".join(r["cohort_ids"])}; {r["group"]}; good walker={next(x["quality_subset"] for x in records if x["fly_key"]==r["fly_key"])}; pre/opto counts={r["pre_events"]}/{r["inhibition_events"]}; rates={r["pre_rate_per_min"]:.4f}/{r["inhibition_rate_per_min"]:.4f}; change={r["delta_rate_per_min"]:.4f} events/min')
    (OUT/'results_and_methods.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines[10:35]),flush=True)

if __name__=='__main__':main()
