"""At least 30-degree turns flanked by forward-velocity changes.

Restricted to the previously selected blinded good walkers. Speed-linked turns
replace the spatial pivot definition; no low-displacement criterion is applied.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import numpy as np
import h5py
from scipy.ndimage import gaussian_filter1d
from audit_ves041_long_stim_behavior import read_vector
from analyze_ves041_inhibition_pivot_frequency import aggregate,inference

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
OUT=ROOT/'outputs/ves041_blinded_good_walkers_speed_linked_turns_20261102'
DEFINITIONS=[
    ('primary',.5,20,'flanking',False,None),
    ('window025',.25,20,'flanking',False,None),
    ('window1',1.,20,'flanking',False,None),
    ('yaw10',.5,10,'flanking',False,None),
    ('yaw30',.5,30,'flanking',False,None),
    ('turn_vs_flank',.5,20,'turn_vs_flank',False,None),
    ('same_phase_flanks',.5,20,'flanking',True,None),
    ('max1s_turn',.5,20,'flanking',False,1.),
]

def write_report(flies,groups,events,counts):
    lines=['VES041 BLINDED GOOD WALKERS: SPEED-LINKED TURNS',
      'Analysis date: November 2, 2026','',
      'SCOPE',
      'The same 8 experimental and 9 control good walkers from blinded round 1; 5 s inhibition. Good-walker selection is unchanged and independent of turn detection. All 40 presentations per selected fly enter the main comparison, including trials starting while stationary. No unblinded or SMP544 flies are included.', '',
      'DEFINITION',
      'A turn qualifies when its absolute net heading change is AT LEAST 30 degrees and EITHER the preceding forward slowdown is AT LEAST 2 mm/s OR the immediately following forward speed-up is AT LEAST 2 mm/s. Both conditions can pass, but a turn is counted only once. The prior displacement and directionality restrictions are removed.',
      'Forward velocity is the signed forwardVelocity signal, not total speed or absolute forward velocity. Using turn onset as 0 and turn end as T:',
      '  Slowdown = mean forward velocity [-1,-0.5) s minus mean [-0.5,0) s >=2 mm/s.',
      '  Speed-up = mean forward velocity [T+0.5,T+1) s minus mean [T,T+0.5) s >=2 mm/s.',
      'This is a literal speed change before or after a turn. A sensitivity analysis instead compares the immediate pre-turn 0.5 s with the first up-to-0.5 s of the turn, or the final up-to-0.5 s of the turn with the immediate post-turn 0.5 s.',
      'Turn boundaries: unwrapped heading is Gaussian-smoothed with sigma 30 ms, then differentiated. A directional turn bout is a contiguous interval with angular velocity >=20 deg/s or <=-20 deg/s. Same-direction bouts separated by <=0.1 s are merged. Opposite-direction bouts are separate. Net heading change is measured across the resulting bout. Primary analysis has no maximum bout duration; <=1 s turns are a sensitivity check. Angular boundary thresholds 10 and 30 deg/s are also checked.', '',
      'PHASE ASSIGNMENT AND EXPOSURE',
      'Count turns fully contained in each 5 s inhibition period or the equal-duration 5 s immediately preceding it. Primary speed flanks may extend outside the phase but must be observed in the recording. This permits full-period assessment, including light-onset/offset responses.',
      'The same_phase_flanks check requires the entire turn and speed-change sequence to stay within its assigned phase. It uses symmetric 1 s edge buffers, leaving 3 s of eligible interior exposure per phase. This guards against a light-induced speed change qualifying a preceding turn. Primary and conservative checks are reported separately.',
      'Rate = qualified turns / eligible elapsed minute. Flies are weighted equally; repeated presentations are aggregated within each fly. Between-group comparisons examine inhibition rate and the difference in inhibition-minus-pre changes. Paired within-group tests use one change per fly. Two-sided exact permutations when feasible, otherwise 49,999 seeded permutations. BH adjustment is separate for the two between-group and two within-group primary comparisons. Secondary definitions are descriptive.', '',
      'RESULTS']
    for r in groups:
        if r['definition']!='primary':continue
        lines+=['',r['subset']]
        for group in ['experimental','control']:
            lines.append(f'{group}: n={r[f"n_{group}"]}; pre={r[f"{group}_pre_rate_per_min"]:.3f}; inhibition={r[f"{group}_inhibition_rate_per_min"]:.3f}; change={r[f"{group}_delta_rate_per_min"]:.3f} turns/min; decreasing flies={r[f"{group}_flies_decrease"]}; increasing flies={r[f"{group}_flies_increase"]}')
        lines.append(f'Experimental minus control inhibition rate: {r["opto_difference"]:.3f}; difference in changes: {r["difference_in_changes"]:.3f} turns/min')
        for k,v in r.items():
            if v is not None and (k.endswith('_p') or k.endswith('_bh_q')):lines.append(f'{k}: {v:.6g}')
    lines+=['','SENSITIVITY (all selected flies)']
    for r in groups:
        if r['subset']=='all':lines.append(f'{r["definition"]}: experimental pre/on={r["experimental_pre_rate_per_min"]:.3f}/{r["experimental_inhibition_rate_per_min"]:.3f}; control pre/on={r["control_pre_rate_per_min"]:.3f}/{r["control_inhibition_rate_per_min"]:.3f}; difference in changes={r["difference_in_changes"]:.3f} turns/min')
    lines+=['','INDIVIDUAL FLIES (primary definition)']
    for r in sorted([r for r in flies if r['definition']=='primary' and r['subset']=='all'],key=lambda r:(r['group'],r['date'],r['daily_fly'])):
        lines.append(f'ID {",".join(r["cohort_ids"])}; {r["group"]}; {r["date"]} daily fly {r["daily_fly"]}; pre/on events={r["pre_events"]}/{r["inhibition_events"]}; pre/on rates={r["pre_rate_per_min"]:.3f}/{r["inhibition_rate_per_min"]:.3f}; change={r["delta_rate_per_min"]:.3f}')
    lines+=['','INTERPRETATION',
      'The experimental mean rate decreases more than the control mean rate under the selected 0.5 s definition, but the between-group difference in changes is uncertain. Restricting to moving-onset trials preserves the descriptive direction. Requiring all speed flanks within phase preserves the direction; widening the speed windows to 1 s greatly attenuates the group difference. The result is sensitive to timing choices and is not firm evidence for selective suppression.',
      'These are speed-linked turning bouts, not necessarily in-place pivots. Event identity depends on angular bout boundaries and averaging windows. Synthetic checks and example kinematic traces support implementation, but raw tracking video has not been manually scored.']
    (OUT/'results_and_methods.txt').write_text('\n'.join(lines)+'\n')

def turn_candidates(heading,dt,threshold):
    smoothed=gaussian_filter1d(heading,max(.03/dt,.5))
    angular=np.gradient(smoothed,dt)
    sign=np.where(angular>=threshold,1,np.where(angular<=-threshold,-1,0))
    edges=np.flatnonzero(np.r_[True,sign[1:]!=sign[:-1],True])
    segments=[]
    for a,b in zip(edges[:-1],edges[1:]):
        direction=int(sign[a])
        if not direction:continue
        if segments and direction==segments[-1][2] and (a-segments[-1][1])*dt<=.1:
            segments[-1]=(segments[-1][0],int(b),direction)
        else:segments.append((int(a),int(b),direction))
    turns=[]
    for a,b,direction in segments:
        end=min(b,len(heading)-1)
        rotation=abs(float(smoothed[end]-smoothed[a]))
        if rotation>=30:
            turns.append(dict(start=a,end=end,anchor=(a+end)//2,
                              rotation_deg=rotation,duration_s=(end-a)*dt,
                              peak_angular_velocity=float(np.max(np.abs(angular[a:end+1])))))
    return turns

def speed_gate(event,forward,dt,window,mode):
    w=round(window/dt);a=event['start'];b=event['end'];n=len(forward)
    if mode=='flanking':
        if a<2*w or b+2*w>n:return None
        early_pre=float(np.mean(forward[a-2*w:a-w]));late_pre=float(np.mean(forward[a-w:a]))
        early_post=float(np.mean(forward[b:b+w]));late_post=float(np.mean(forward[b+w:b+2*w]))
        slow=early_pre-late_pre;fast=late_post-early_post
        support_start=a-2*w;support_end=b+2*w
    else:
        if a<w or b+w>n:return None
        early_pre=float(np.mean(forward[a-w:a]));late_pre=float(np.mean(forward[a:min(a+w,b+1)]))
        early_post=float(np.mean(forward[max(a,b-w):b+1]));late_post=float(np.mean(forward[b:b+w]))
        slow=early_pre-late_pre;fast=late_post-early_post
        support_start=a-w;support_end=b+w
    return dict(slowdown_mm_s=slow,speedup_mm_s=fast,slowdown_pass=slow>=2,
                speedup_pass=fast>=2,qualifies=slow>=2 or fast>=2,
                support_start=support_start,support_end=support_end,
                early_pre_mean=early_pre,late_pre_mean=late_pre,
                early_post_mean=early_post,late_post_mean=late_post)

def process(rec):
    cache=OUT/'kinematic_cache'/f'{rec["recording_key"]}.npz'
    if cache.exists():
        with np.load(cache) as f:
            h=f['heading'];forward=f['forward'];speed=f['speed'];light=f['light'];dt=float(f['dt'])
    else:
        with h5py.File(rec['source'],'r') as f:
            g=f['exptData'];n=g['totSpeed'].size
            tn=next(k for k in ['t_ds','tDS','t'] if k in g and g[k].size==n)
            td=g[tn];dt=(float(td[-1,-1])-float(td[0,0]))/(n-1)
            step=max(1,round(.01/dt));dt*=step
            speed=read_vector(g['totSpeed'],step);forward=read_vector(g['forwardVelocity'],step)
            h=read_vector(g['headingPosition'],step)
            on=next(k for k in ['optoStim_ds','optoStim'] if k in g)
            scale=g[on].size/n
            if abs(scale-round(scale))>1e-5:raise ValueError('Unaligned stimulus and behavior')
            opto=read_vector(g[on],step*round(scale))
        n=min(map(len,[speed,forward,h,opto]));speed=speed[:n];forward=forward[:n]
        h=np.rad2deg(np.unwrap(np.deg2rad(h[:n])));light=opto[:n]!=0
        np.savez_compressed(cache,heading=h,forward=forward,speed=speed,light=light,dt=dt)
    base={k:rec[k] for k in ['cohort','date','daily_fly','cell','trial','cohort_ids','duration_s','group','fly_key','recording_key','quality_subset']}
    edges=np.diff(np.r_[False,light,False].astype(np.int8));starts=np.flatnonzero(edges==1);ends=np.flatnonzero(edges==-1)
    candidates={threshold:turn_candidates(h,dt,threshold) for threshold in [10,20,30]}
    trials=[];eventrows=[];candidate_counts=[]
    for definition,window,threshold,mode,within_phase,max_duration in DEFINITIONS:
        gate_events=[]
        for e in candidates[threshold]:
            if max_duration is not None and e['duration_s']>max_duration:continue
            gate=speed_gate(e,forward,dt,window,mode)
            if gate is not None:gate_events.append({**e,**gate})
        for pulse,(a,b) in enumerate(zip(starts,ends),1):
            length=b-a;p=a-length
            if p<0 or np.any(light[p:a]):continue
            one=round(1/dt);moving_onset=np.mean(speed[a-one:a]>=3)>=.5
            for phase,lo,hi in [('pre',p,a),('inhibition',a,b)]:
                # Same-phase flanks exclude a symmetric edge buffer from each phase.
                pad=round((2*window if mode=='flanking' else window)/dt) if within_phase else 0
                eligible_start=lo+pad;eligible_end=hi-pad
                exposure=(eligible_end-eligible_start)*dt
                if exposure<=0:continue
                available=[e for e in gate_events if e['start']>=eligible_start and e['end']<=eligible_end]
                accepted=[e for e in available if e['qualifies']]
                trials.append({**base,'pulse':pulse,'phase':phase,'definition':definition,
                  'events':len(accepted),'seconds':exposure,'moving_onset':bool(moving_onset),
                  'moving_fraction':float(np.mean(speed[eligible_start:eligible_end]>=3)),
                  'mean_forward':float(np.mean(forward[eligible_start:eligible_end])),
                  'onset_1s_pivot':any(e['anchor']<lo+one for e in accepted)})
                candidate_counts.append({**base,'pulse':pulse,'phase':phase,'definition':definition,
                  'turns_ge30':len(available),'slowdown_only':sum(e['slowdown_pass'] and not e['speedup_pass'] for e in accepted),
                  'speedup_only':sum(e['speedup_pass'] and not e['slowdown_pass'] for e in accepted),
                  'both':sum(e['speedup_pass'] and e['slowdown_pass'] for e in accepted)})
                if definition=='primary':
                    for e in available:
                        aa=e['start'];bb=e['end'];idx=np.arange(max(0,aa-round(1.1/dt)),min(len(h),bb+round(1.1/dt)),max(1,round(.025/dt)))
                        eventrows.append({**base,'pulse':pulse,'phase':phase,**e,
                          'relative_time_to_light_s':(aa-a)*dt,
                          'trace_time_s':[(int(i)-aa)*dt for i in idx],
                          'trace_forward':[float(v) for v in forward[idx]],
                          'trace_heading':[float(v) for v in h[idx]-h[aa]]})
    return trials,eventrows,candidate_counts

def main():
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'kinematic_cache').mkdir(exist_ok=True)
    records=json.loads((SOURCE/'manifest.json').read_text())
    assert len(records)==17 and all(r['cohort']=='blinded round 1' and r['quality_subset'] for r in records)
    trials=[];events=[];counts=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(process,r) for r in records]
        for i,future in enumerate(as_completed(futures),1):
            t,e,c=future.result();trials.extend(t);events.extend(e);counts.extend(c)
            print(f'Extracted {i}/{len(records)}',flush=True)
    for name,obj in [('presentations',trials),('events',events),('candidate_counts',counts),('manifest',records)]:
        (OUT/f'{name}.json').write_text(json.dumps(obj,allow_nan=False))
    flies=[f for f in aggregate(trials) if f['subset']!='quality_subset']
    (OUT/'fly_summary.json').write_text(json.dumps(flies,indent=2))
    groups=inference(flies);(OUT/'group_summary.json').write_text(json.dumps(groups,indent=2))
    write_report(flies,groups,events,counts)
    (OUT/'parameters.json').write_text(json.dumps(dict(angle_min_deg=30,speed_change_mm_s=2,
       speed_signal='forwardVelocity',primary_window_s=.5,angular_boundary_deg_s=20,
       heading_smoothing_sigma_s=.03,max_merge_gap_s=.1,definitions=DEFINITIONS),indent=2))
    print('Completed:',len(records),'recordings;',sum(e['qualifies'] for e in events),'primary speed-linked turns',flush=True)

if __name__=='__main__':main()
