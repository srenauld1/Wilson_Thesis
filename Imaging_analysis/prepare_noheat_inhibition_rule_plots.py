from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]/'outputs/ves041_imaging_inhibition_rules_noheat_20261104'
a=json.loads((ROOT/'analysis.json').read_text())
time=np.asarray(a['time']);values=np.asarray(a['fly_means'],dtype=float)
mean=np.asarray(a['mean'],dtype=float);sem=np.asarray(a['sem'],dtype=float);counts=np.asarray(a['counts'])
assert values.shape==(101,6,3,2) and counts.shape==(6,2)
assert np.array_equal(counts.sum(axis=0),[259,238])
np.testing.assert_allclose(np.mean(values,axis=1),mean,atol=1e-10)
np.testing.assert_allclose(np.std(values,axis=1,ddof=1)/np.sqrt(6),sem,atol=1e-10)
events=a['events'];assert all(x['rotation_deg']>30 and x['displacement_mm']<=2 and x['directionality']>=.7 for x in events)
rows=[];summary=[];metrics=['dff','forward','rotation'];defs=['imaging','localized']
for j,definition in enumerate(defs):
 for v,metric in enumerate(metrics):
  for t in range(len(time)):
   summary.append(dict(definition=definition,metric=metric,time=float(time[t]),mean=float(mean[t,v,j]),sem=float(sem[t,v,j]),n=6,events=int(counts[:,j].sum())))
   for f,fly in enumerate(a['flies']):rows.append(dict(definition=definition,metric=metric,time=float(time[t]),fly=fly,value=float(values[t,f,v,j]),events=int(counts[f,j])))
(ROOT/'fly_traces.json').write_text(json.dumps(rows))
(ROOT/'mean_sem.json').write_text(json.dumps(summary))
report=['VES041 NO-HEAT IMAGING: APPLYING THE ORIGINAL INHIBITION LOCALIZED-TURN DEFINITION','',
 'No heat cohort was analyzed. Established pooled dark/bar/bar-jump recordings only.',
 '10 eligible recordings, 6 mapped flies; the February 17 animal identity is inferred, as in the previous pooled analysis.',
 'Two figure-only recordings (20260311-1 and 20260313-1) lack saved actual heading and remain unassessable for this strict rule. They are excluded from BOTH columns of the comparison.',
 'All six flies contribute to both definitions. Imaging slowdown+rotation: 259 events. Inhibition-style localized turns: 238 events.',
 'Detection uses recorded unwrapped ball heading and original physical ball XY at the saved imaging-aligned behavioral sampling rate; no heading reconstructed from yaw and no invented high-rate stream.',
 'Physical XY consistency: median displacement-derived speed / translational speed per recording: '+', '.join(f'{m["path_speed_ratio"]:.4f}' for m in a['metadata'])+'.',
 'New selection: strict >30-degree net heading change in a frame-rounded 1-second window; <=2 mm net XY displacement; directionality >=0.7 (net / cumulative absolute heading change). Sliding windows about 0.1 s apart; overlapping qualifying windows merged. No slowdown requirement or added 2-second event spacing.',
 'Actual window duration range: '+str([min(m['window_duration_s'] for m in a['metadata']),max(m['window_duration_s'] for m in a['metadata'])])+' s. Lower imaging-aligned sampling than inhibition can affect directionality estimates.',
 'Localized-turn alignment: CENTER of the representative qualifying window with maximum net rotation in the merged bout. It is NOT measured turn onset, strongest deceleration, or a fluorescence peak.',
 'Old imaging-rule alignment: strongest deceleration. Different anchors preclude direct cross-definition fluorescence-latency comparisons.',
 'Full valid -2:+3-second epochs. Original dF/F and kinematics retained; no fluorescence amplitude/peak gate, baseline subtraction, normalization, or display smoothing added. Linear time interpolation for display matches the original imaging plots.',
 'Events pooled across recordings within each fly, then equal-weight average over flies. Pointwise SEM: sample SD of the six fly means / sqrt(6). No hypothesis tests added.',
 'Old imaging events are the saved original group4 epochs, restricted to the same recordings and common full heading/XY/fluorescence validity.',
 '', 'PER-FLY COUNTS (imaging rule / localized-turn rule)']
report += [f'{fly}: {int(counts[i,0])} / {int(counts[i,1])}' for i,fly in enumerate(a['flies'])]
report+=['','EXTRACTION STATUS']+[f'{s["trial"]}: {s["status"]}; {s["message"]}' for s in a['status']]
(ROOT/'methods_and_counts.txt').write_text('\n'.join(report)+'\n')
print('\n'.join(report[:10]))
