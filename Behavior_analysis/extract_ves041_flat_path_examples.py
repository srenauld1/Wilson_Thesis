"""Actual paired flat paths for blinded good walkers, with explicit selection."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
import numpy as np
import h5py
from audit_ves041_long_stim_behavior import read_vector

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
KIN=ROOT/'outputs/ves041_blinded_good_walkers_speed_linked_turns_20261102/kinematic_cache'
OUT=SOURCE/'flat_path_examples'

def extract(rec):
    with np.load(KIN/f'{rec["recording_key"]}.npz') as z:
        speed=z['speed'];forward=z['forward'];light=z['light'];dt=float(z['dt'])
    with h5py.File(rec['source'],'r') as f:
        g=f['exptData'];n=g['totSpeed'].size
        tn=next(k for k in ['t_ds','tDS','t'] if k in g and g[k].size==n)
        t=g[tn];original_dt=(float(t[-1,-1])-float(t[0,0]))/(n-1)
        step=round(dt/original_dt)
        x=read_vector(g['x'],step);y=read_vector(g['y'],step)
    assert len(x)==len(light) and len(y)==len(light)
    edge=np.diff(np.r_[False,light,False].astype(np.int8));starts=np.flatnonzero(edge==1);ends=np.flatnonzero(edge==-1)
    rows=[]
    for pulse,(a,b) in enumerate(zip(starts,ends),1):
        p=a-(b-a)
        if p<0 or np.any(light[p:a]):continue
        moving_onset=bool(np.mean(speed[a-round(1/dt):a]>=3)>=.5)
        for phase,lo,hi in [('pre',p,a),('inhibition',a,b)]:
            idx=np.unique(np.r_[np.arange(lo,hi,max(1,round(.025/dt))),hi-1])
            rows.append(dict(pulse=pulse,phase=phase,moving_onset=moving_onset,
                mean_forward=float(np.mean(forward[lo:hi])),
                sample_indices=(idx-lo).tolist(),
                time_s=((idx-lo)*dt).tolist(),
                x=(x[idx]-x[lo]).tolist(),y=(y[idx]-y[lo]).tolist()))
    return dict(**rec,dt=dt,paths=rows)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((SOURCE/'manifest.json').read_text())
    trials=json.loads((SOURCE/'presentations.json').read_text())
    metrics={(t['recording_key'],t['pulse'],t['phase']):t for t in trials if t['definition']=='primary'}
    events=json.loads((SOURCE/'event_examples.json').read_text())
    by_event={}
    for e in events:by_event.setdefault((e['recording_key'],e['pulse'],e['phase']),[]).append(e)
    results=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for i,future in enumerate(as_completed([pool.submit(extract,r) for r in manifest]),1):
            r=future.result()
            for path in r['paths']:
                key=(r['recording_key'],path['pulse'],path['phase'])
                path['events']=metrics[key]['events']
                path['event_windows']=[dict(start=e['representative_start'],end=e['representative_end']) for e in by_event.get(key,[])]
            results.append(r);print(f'Extracted {i}/{len(manifest)}',flush=True)
    results.sort(key=lambda r:(r['group'],r['date'],r['daily_fly']))
    chosen=[]
    for fly_id in ['22C','8','13','30']:
        r=next(r for r in results if r['cohort_ids']==[fly_id])
        eligible=[p for p in r['paths'] if p['phase']=='pre' and p['events']>=1 and p['moving_onset']]
        median=float(np.median([p['mean_forward'] for p in eligible]))
        selected=min(eligible,key=lambda p:(abs(p['mean_forward']-median),p['pulse']))
        chosen.append(dict(fly_id=fly_id,recording_key=r['recording_key'],group=r['group'],pulse=selected['pulse'],
                           eligibility='At least one pre-opto reorientation and moving at light onset',
                           selection='Pre-opto mean forward velocity closest to the median among eligible trials; opto-on behavior not used'))
    (OUT/'all_flat_paths.json').write_text(json.dumps(results))
    (OUT/'selected_examples.json').write_text(json.dumps(chosen,indent=2))
    lines=['FLAT-PATH EXAMPLES: VES041 BLINDED GOOD WALKERS','',
      'These are measured paths reconstructed from stored x/y data, not drawings or averages. Every displayed segment is 5 s. Paths are translated independently to start at (0,0); orientations and relative geometry are preserved. No rotations or reflections are applied.',
      'The main figure includes experimental IDs 22 and 8 and control IDs 13 and 30. Fly 22 illustrates a strong experimental reduction; fly 8 a smaller reduction. Control 13 is included deliberately because it also shows a strong reduction; control 30 has the same mean pre-opto frequency as experimental fly 8. These examples do not establish a group-level difference.',
      'Within each fly, eligible examples have >=1 detected pre-opto reorientation and movement at light onset. Choose the trial with pre-opto mean forward velocity closest to the eligible-trial median; break ties by pulse number. The opto-on path and event count are not used to select the trial.',
      'Purple segments mark the representative 1 s windows classified as >30 degree turns, <=2 mm net displacement, and directionality >=0.7. This corresponds to the spatial definition in the current talk figure, not the later speed-linked definition.',
      'The atlas contains all 40 presentations from every included fly, including trials with no reorientations. Individual-fly mean rates remain the appropriate evidence for group comparisons.','']
    for c in chosen:
        r=next(r for r in results if r['recording_key']==c['recording_key'])
        p=next(p for p in r['paths'] if p['pulse']==c['pulse'] and p['phase']=='pre')
        o=next(p for p in r['paths'] if p['pulse']==c['pulse'] and p['phase']=='inhibition')
        lines.append(f'{c["group"]} ID {c["fly_id"]}, {r["date"]}, daily fly {r["daily_fly"]}, pulse {c["pulse"]}; pre/opto events {p["events"]}/{o["events"]}; mean forward {p["mean_forward"]:.2f}/{o["mean_forward"]:.2f} mm/s; source {r["source"]}')
    (OUT/'selection_and_methods.txt').write_text('\n'.join(lines)+'\n')
    print('Selected examples:',json.dumps(chosen),flush=True)

if __name__=='__main__':main()
