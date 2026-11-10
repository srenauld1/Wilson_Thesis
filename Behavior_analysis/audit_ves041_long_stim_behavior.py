"""Descriptive behavior screening of the verified long-stimulus inventory.

Selection uses pre-light movement, without inspecting genotype-specific responses.
No hypothesis tests or automatic experiment exclusions are produced.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import numpy as np
import h5py

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/ves041_inhibition_duration_audit_20261102'
THRESHOLD = 3.0  # Existing lab minVel used in master_behavior_analysis.m.

def read_vector(ds, step=1):
    if ds.ndim == 1:
        return np.asarray(ds[::step], dtype=float).reshape(-1)
    if ds.shape[0] == 1:
        return np.asarray(ds[0, ::step], dtype=float).reshape(-1)
    return np.asarray(ds[::step, 0], dtype=float).reshape(-1)

def mean_or_none(a):
    a=np.asarray(a);a=a[np.isfinite(a)]
    return float(np.mean(a)) if a.size else None

def write_report(results):
    lines=['VES041 LONG-STIMULUS RECORDINGS: DESCRIPTIVE WALKING SCREEN',
           'Audit date: November 2, 2026', '',
           'What good behavior means here: enough spontaneous walking and enough presentations with a moving pre-light baseline to support kinematic analysis. Turning style and the size of the inhibition response are not selection criteria.',
           'Movement threshold: total speed >=3 mm/s, following the existing lab minVel setting. Baseline: the 1 s immediately before each light onset. A moving-baseline presentation has >=50% moving samples and >=95% finite speed, forward velocity, and heading.',
           'Prioritize is a descriptive convenience flag: >=50% moving baseline samples, >=50% moving during light-off periods, and >=10 moving-baseline presentations. These coverage cutoffs are screening choices, not validated biological thresholds or exclusion rules. Recordings below them can still contribute their moving subset.',
           'The 240 s protocol has only two blocks, so its assessment uses light-off movement and is listed separately. One-second baselines before just two long blocks are a limited sample.',
           'All 130 readable recordings have finite speed, forward velocity, and heading at the sampled frames. This does not establish raw-video or FicTrac tracking quality. The January 21, 2026 cohort-ID 15 long-stimulus processed files are damaged; the light command is readable in raw data, but behavior was not reconstructed.',
           'Repeated cells/trials belong to the same fly. Selection should refer to the listed recording: a fly can walk well in one trial and poorly in another. No recordings were excluded, and no hypothesis tests were run.', '',
           'Columns: cohort ID; date/daily fly/cell/trial; group; pre-light moving %; light-off moving %; moving-baseline presentations; mean pre-light forward velocity (mm/s); review flag.']
    for cohort,duration in sorted({(r['cohort'],r['duration_s']) for r in results}):
        lines += ['',f'{cohort.upper()} / {duration} SECOND LIGHT BLOCKS']
        subset=sorted([r for r in results if r['cohort']==cohort and r['duration_s']==duration],
                      key=lambda r:r.get('baseline_moving_fraction',-1),reverse=True)
        for r in subset:
            name=f'ID {",".join(r["cohort_ids"]) or "unassigned"}; {r["date"]} fly {r["daily_fly"]} cell {r["cell"]} trial {r["trial"]}; {r["group"]}'
            if r['status']!='ok':lines += [name+'; CANNOT ASSESS: processed behavior unreadable'];continue
            if duration==240:flag='long-block review; only 2 onsets'
            elif r['baseline_moving_fraction']>=.5 and r['off_moving_fraction']>=.5 and r['n_moving_baseline_pulses']>=10:flag='prioritize'
            elif r['n_moving_baseline_pulses']>=10:flag='usable moving subset'
            else:flag='sparse moving baselines'
            q='/'.join(f'{100*x:.0f}' if x is not None else 'NA' for x in r['off_movement_by_quarter'])
            lines += [name+f'; {100*r["baseline_moving_fraction"]:.1f}%; {100*r["off_moving_fraction"]:.1f}%; {r["n_moving_baseline_pulses"]}/{r["n_pulses"]}; {r["baseline_mean_forward"]:.2f}; {flag}',
                      '  Light-off moving % across successive quarters: '+q,
                      '  File: '+r['behavior_file']]
    (OUT/'behavior_quality_report.txt').write_text('\n'.join(lines)+'\n')

def audit(rec):
    result={**rec}
    candidates=sorted(rec['paths'],key=lambda p:(not p.endswith('_processed_final.mat'),len(p)))
    errors=[]
    for path in candidates:
        if path.endswith('_raw.mat'):continue
        try:
            with h5py.File(path,'r') as f:
                g=f['exptData']; n=g['totSpeed'].size
                tn=next(k for k in ['t_ds','tDS','t'] if k in g and g[k].size==n)
                td=g[tn];start=float(td[0,0]);end=float(td[-1,-1]);dt=(end-start)/(n-1)
                # Approximately 100 Hz is sufficient for movement coverage and QC.
                step=max(1,int(round(.01/dt))); dt*=step
                speed=read_vector(g['totSpeed'],step)
                forward=read_vector(g['forwardVelocity'],step)
                heading=read_vector(g['headingPosition'],step)
                on=next(k for k in ['optoStim_ds','optoStim'] if k in g)
                scale=g[on].size/n
                if abs(scale-round(scale))>1e-5:raise ValueError('Unaligned stimulus and behavior lengths')
                opto=read_vector(g[on],step*round(scale))
            n=min(map(len,[speed,forward,heading,opto]));speed=speed[:n];forward=forward[:n];heading=heading[:n];opto=opto[:n]
            finite=np.isfinite(speed)&np.isfinite(forward)&np.isfinite(heading)
            light=np.isfinite(opto)&(opto!=0)
            edges=np.diff(np.r_[False,light,False].astype(np.int8));onsets=np.flatnonzero(edges==1)
            moving=finite&(speed>=THRESHOLD)
            prefrac=[];prevalid=[];prespeed=[];preforward=[];preindices=[]
            for onset in onsets:
                i=np.arange(max(0,onset-round(1/dt)),onset)
                if len(i)<round(.95/dt):continue
                preindices.extend(i)
                prevalid.append(float(np.mean(finite[i])))
                prefrac.append(float(np.mean(moving[i])))
                prespeed.append(mean_or_none(speed[i]))
                preforward.append(mean_or_none(forward[i]))
            pi=np.asarray(preindices,dtype=int)
            valid_windows=np.asarray(prevalid)>=.95
            active_windows=valid_windows&(np.asarray(prefrac)>=.5)
            off=~light
            quarters=np.minimum(3,np.arange(n)*4//n)
            qmovement=[float(np.mean(moving[off&(quarters==q)])) if np.any(off&(quarters==q)) else None for q in range(4)]
            baseline_fraction=float(np.mean(moving[pi])) if pi.size else None
            coverage=('frequent baseline movement' if baseline_fraction>=.5 else
                      'intermittent baseline movement' if baseline_fraction>=.2 else 'little baseline movement')
            finite_fraction=float(np.mean(finite))
            if finite_fraction<.99:coverage+='; review missing values'
            per_pulse=[{'pulse':i+1,'pre_moving_fraction':p,'pre_valid_fraction':v,'pre_mean_speed':s,'pre_mean_forward':fw} for i,(p,v,s,fw) in enumerate(zip(prefrac,prevalid,prespeed,preforward))]
            bins=[]
            for b in range(int(np.ceil(n*dt/10))):
                i=np.arange(round(b*10/dt),min(n,round((b+1)*10/dt)))
                if len(i):bins.append({'time_s':b*10+5,'moving_fraction':float(np.mean(moving[i])),'light_fraction':float(np.mean(light[i])),'mean_forward':mean_or_none(forward[i])})
            result.update(status='ok',behavior_file=path,duration_recording_s=n*dt,
                finite_fraction=finite_fraction,off_moving_fraction=float(np.mean(moving[off])),
                off_mean_speed=mean_or_none(speed[off]),off_median_speed=float(np.nanmedian(speed[off])),
                off_mean_forward=mean_or_none(forward[off]),
                baseline_moving_fraction=baseline_fraction,baseline_mean_speed=mean_or_none(speed[pi]),
                baseline_mean_forward=mean_or_none(forward[pi]),
                n_pulses=len(prefrac),n_moving_baseline_pulses=int(np.sum(active_windows)),
                fraction_moving_baseline_pulses=float(np.mean(active_windows)),
                n_valid_baseline_pulses=int(np.sum(valid_windows)),
                off_movement_by_quarter=qmovement,
                speed_p999=float(np.nanquantile(speed,.999)),
                classification=coverage,pulses_detail=per_pulse,time_bins=bins,
                notes='Only two long blocks; baseline presentation count is inherently limited.' if rec['duration_s']==240 else '')
            return result
        except Exception as e:errors.append(f'{path}: {e}')
    return {**result,'status':'unreadable behavior','classification':'cannot assess behavior','errors':errors}

def main():
    records=json.loads((OUT/'behavior_quality_records.json').read_text())
    results=[]
    with ProcessPoolExecutor(max_workers=3) as pool, (OUT/'behavior_quality_results.jsonl').open('w') as handle:
        for i,future in enumerate(as_completed([pool.submit(audit,r) for r in records]),1):
            row=future.result();results.append(row);handle.write(json.dumps(row,allow_nan=False)+'\n');handle.flush()
            if i%10==0:print(f'Completed {i}/{len(records)}',flush=True)
    results.sort(key=lambda r:(r['date'],r['daily_fly'],r['cell'],r['trial']))
    (OUT/'behavior_quality_results.json').write_text(json.dumps(results,indent=2,allow_nan=False))
    summary=[{k:v for k,v in r.items() if k not in ['pulses_detail','time_bins','paths','errors']} for r in results]
    (OUT/'behavior_quality_summary.json').write_text(json.dumps(summary,indent=2))
    write_report(results)
    print('Finished',len(results),'recordings;',sum(r['status']=='ok' for r in results),'readable',flush=True)

if __name__=='__main__':main()
