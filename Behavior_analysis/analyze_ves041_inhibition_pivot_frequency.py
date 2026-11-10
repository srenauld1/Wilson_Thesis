"""Pivot-like bout frequency in verified VES041 long-stimulus recordings.

Reads source data only. Hypothesis tests are performed on fly aggregates,
never on pooled presentations from different flies.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from collections import defaultdict
import argparse
import json
import numpy as np
import h5py
from scipy.stats import permutation_test
from audit_ves041_long_stim_behavior import read_vector

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'outputs/ves041_inhibition_duration_audit_20261102'
OUT=ROOT/'outputs/ves041_inhibition_pivot_frequency_20261102'
DEFINITIONS=[('primary',45,2,'displacement'),('turn30',30,2,'displacement'),
             ('turn60',60,2,'displacement'),('disp1',45,1,'displacement'),
             ('disp3',45,3,'displacement'),('path2',45,2,'path')]

def window_features(h,x,y,dt):
    w=round(1/dt)
    starts=np.arange(0,len(h)-w,max(1,round(.1/dt)))
    ends=starts+w
    total=np.r_[0,np.cumsum(np.abs(np.diff(h)))]
    path=np.r_[0,np.cumsum(np.hypot(np.diff(x),np.diff(y)))]
    net=np.abs(h[ends]-h[starts])
    turn=total[ends]-total[starts]
    directionality=np.divide(net,turn,out=np.zeros_like(net),where=turn>0)
    displacement=np.hypot(x[ends]-x[starts],y[ends]-y[starts])
    length=path[ends]-path[starts]
    return starts,ends,net,directionality,displacement,length

def bouts(features,name,angle,mm,kind,dt,strict=False):
    starts,ends,net,di,displacement,length=features
    translation=displacement if kind=='displacement' else length
    rotation_ok=net>angle if strict else net>=angle
    selected=np.flatnonzero(rotation_ok&(di>=.7)&(translation<=mm))
    events=[]
    for i in selected:
        if events and starts[i]<=events[-1]['last_end']:
            e=events[-1];e['last_end']=int(ends[i])
            if net[i]>e['rotation_deg']:
                e.update(representative_start=int(starts[i]),representative_end=int(ends[i]),rotation_deg=float(net[i]),displacement_mm=float(displacement[i]),path_mm=float(length[i]),directionality=float(di[i]))
        else:
            events.append(dict(first_start=int(starts[i]),last_end=int(ends[i]),representative_start=int(starts[i]),representative_end=int(ends[i]),rotation_deg=float(net[i]),displacement_mm=float(displacement[i]),path_mm=float(length[i]),directionality=float(di[i])))
    return events

def process(rec,definitions=None,strict=False):
    definitions=definitions or DEFINITIONS
    primary_angle=next(d[1] for d in definitions if d[0]=='primary')
    path=rec['behavior_file']
    with h5py.File(path,'r') as f:
        g=f['exptData'];n=g['totSpeed'].size
        tn=next(k for k in ['t_ds','tDS','t'] if k in g and g[k].size==n)
        td=g[tn];dt=(float(td[-1,-1])-float(td[0,0]))/(n-1)
        step=max(1,round(.01/dt));dt*=step
        speed=read_vector(g['totSpeed'],step);forward=read_vector(g['forwardVelocity'],step)
        h=read_vector(g['headingPosition'],step)
        x=read_vector(g['x'],step);y=read_vector(g['y'],step)
        on=next(k for k in ['optoStim_ds','optoStim'] if k in g)
        scale=g[on].size/n
        if abs(scale-round(scale))>1e-5:raise ValueError('Unaligned stimulus')
        opto=read_vector(g[on],step*round(scale))
    n=min(map(len,[speed,forward,h,x,y,opto]));speed=speed[:n];forward=forward[:n];h=h[:n];x=x[:n];y=y[:n];opto=opto[:n]
    h=np.rad2deg(np.unwrap(np.deg2rad(h)))
    light=opto!=0;ed=np.diff(np.r_[False,light,False].astype(np.int8))
    starts=np.flatnonzero(ed==1);ends=np.flatnonzero(ed==-1)
    base={k:rec[k] for k in ['cohort','date','daily_fly','cell','trial','cohort_ids','duration_s','group']}
    base['fly_key']=f'{rec["cohort"]}:{rec["date"]}:fly{rec["daily_fly"]:02d}'
    base['recording_key']=f'{rec["date"]}_fly{rec["daily_fly"]:02d}_cell{rec["cell"]:02d}_trial{rec["trial"]}'
    base['quality_subset']=rec['baseline_moving_fraction']>=.5 and rec['off_moving_fraction']>=.5 and rec['n_moving_baseline_pulses']>=10
    trials=[];event_examples=[];dropped=0
    for pulse,(a,b) in enumerate(zip(starts,ends),1):
        # Actual on duration determines equal-length pre and inhibition exposures.
        length=b-a
        available_off=a-(ends[pulse-2] if pulse>1 else 0)
        # At 50% duty cycle the sampled on/off lengths can differ by one frame.
        if available_off<length and (length-available_off)*dt<=.03:
            length=available_off;b=a+length
        p=a-length
        if p<0 or np.any(light[p:a]):dropped+=1;continue
        valid=np.isfinite(h[p:b])&np.isfinite(x[p:b])&np.isfinite(y[p:b])&np.isfinite(speed[p:b])
        if np.mean(valid)<.99:dropped+=1;continue
        one=round(1/dt)
        moving_onset=float(np.mean(speed[max(0,a-one):a]>=3))>=.5
        for phase,lo,hi in [('pre',p,a),('inhibition',a,b)]:
            features=window_features(h[lo:hi],x[lo:hi],y[lo:hi],dt)
            firstlo=a-one if phase=='pre' else a
            firsthi=a if phase=='pre' else a+one
            hh=h[firstlo:firsthi];xx=x[firstlo:firsthi];yy=y[firstlo:firsthi]
            net=abs(float(hh[-1]-hh[0]));total=float(np.sum(np.abs(np.diff(hh))))
            onset_flag=(net>primary_angle if strict else net>=primary_angle) and np.hypot(xx[-1]-xx[0],yy[-1]-yy[0])<=2 and total>0 and net/total>=.7
            for name,angle,mm,kind in definitions:
                events=bouts(features,name,angle,mm,kind,dt,strict)
                trials.append({**base,'pulse':pulse,'phase':phase,'definition':name,
                  'events':len(events),'seconds':length*dt,'moving_onset':moving_onset,
                  'moving_fraction':float(np.mean(speed[lo:hi]>=3)),
                  'mean_forward':float(np.mean(forward[lo:hi])),
                  'onset_1s_pivot':bool(onset_flag)})
                if name=='primary':
                    for event in events:
                        s=event['representative_start']+lo;t=event['representative_end']+lo
                        indices=np.arange(s,t+1,max(1,round(.05/dt)))
                        event_examples.append({**base,'pulse':pulse,'phase':phase,**event,
                          'relative_time_s':(s-a)*dt,
                          'x':[float(v) for v in x[indices]-x[s]],
                          'y':[float(v) for v in y[indices]-y[s]],
                          'heading':[float(v) for v in h[indices]-h[s]]})
    return {'record':{**base,'source':path,'dropped_presentations':dropped},'trials':trials,'examples':event_examples}

def avg(x):return float(np.mean(x))

def bh(p):
    p=np.array(p);order=np.argsort(p);q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty_like(q);out[order]=np.minimum(q,1);return out.tolist()

def between(x,y):
    if min(len(x),len(y))<2:return None
    if np.all(np.asarray(x)==0) and np.all(np.asarray(y)==0):return 1.0
    return float(permutation_test((np.array(x),np.array(y)),lambda a,b:np.mean(a)-np.mean(b),
       permutation_type='independent',n_resamples=49999,alternative='two-sided',rng=np.random.default_rng(544),vectorized=False).pvalue)

def within(x):
    if len(x)<3:return None
    if np.all(np.asarray(x)==0):return 1.0
    return float(permutation_test((np.array(x),),lambda a:np.mean(a),permutation_type='samples',
       n_resamples=49999,alternative='two-sided',rng=np.random.default_rng(41),vectorized=False).pvalue)

def aggregate(trials):
    flyrows=[]
    for state in ['all','moving_onset','stationary_onset','quality_subset']:
        groups=defaultdict(list)
        for t in trials:
            if state=='moving_onset' and not t['moving_onset']:continue
            if state=='stationary_onset' and t['moving_onset']:continue
            if state=='quality_subset' and not t['quality_subset']:continue
            key=(t['cohort'],t['duration_s'],t['fly_key'],t['definition'])
            groups[key].append(t)
        for key,ts in groups.items():
            r={k:ts[0][k] for k in ['cohort','duration_s','fly_key','date','daily_fly','group','cohort_ids']}
            r.update(subset=state,definition=key[-1])
            for phase in ['pre','inhibition']:
                tt=[t for t in ts if t['phase']==phase]
                seconds=sum(t['seconds'] for t in tt);events=sum(t['events'] for t in tt)
                r[f'{phase}_events']=events;r[f'{phase}_seconds']=seconds
                r[f'{phase}_rate_per_min']=events/seconds*60
                r[f'{phase}_incidence']=avg([t['events']>0 for t in tt])
                r[f'{phase}_moving_fraction']=sum(t['seconds']*t['moving_fraction'] for t in tt)/seconds
                r[f'{phase}_forward']=sum(t['seconds']*t['mean_forward'] for t in tt)/seconds
                r[f'{phase}_onset1s_incidence']=avg([t['onset_1s_pivot'] for t in tt])
            r['n_presentations']=len(ts)//2
            r['delta_rate_per_min']=r['inhibition_rate_per_min']-r['pre_rate_per_min']
            r['delta_incidence']=r['inhibition_incidence']-r['pre_incidence']
            r['delta_onset1s_incidence']=r['inhibition_onset1s_incidence']-r['pre_onset1s_incidence']
            flyrows.append(r)
    return flyrows

def inference(flyrows):
    results=[]
    conditions=sorted({(r['cohort'],r['duration_s'],r['definition'],r['subset']) for r in flyrows})
    for cohort,dur,definition,subset in conditions:
        data=[r for r in flyrows if (r['cohort'],r['duration_s'],r['definition'],r['subset'])==(cohort,dur,definition,subset)]
        ex=[r for r in data if r['group'].startswith('experimental')];co=[r for r in data if r['group'].startswith('control')]
        r=dict(cohort=cohort,duration_s=dur,definition=definition,subset=subset,n_experimental=len(ex),n_control=len(co))
        for label,group in [('experimental',ex),('control',co)]:
            if not group:continue
            for m in ['pre_rate_per_min','inhibition_rate_per_min','delta_rate_per_min','pre_incidence','inhibition_incidence','pre_moving_fraction','inhibition_moving_fraction','pre_forward','inhibition_forward','pre_onset1s_incidence','inhibition_onset1s_incidence','delta_onset1s_incidence']:
                r[f'{label}_{m}']=avg([x[m] for x in group])
            r[f'{label}_flies_decrease']=int(sum(x['delta_rate_per_min']< -1e-8 for x in group))
            r[f'{label}_flies_increase']=int(sum(x['delta_rate_per_min']>1e-8 for x in group))
        if ex and co:
            r['opto_difference']=r['experimental_inhibition_rate_per_min']-r['control_inhibition_rate_per_min']
            r['difference_in_changes']=r['experimental_delta_rate_per_min']-r['control_delta_rate_per_min']
        # Fly-level comparisons for all flies and the previously defined activity screen.
        if definition=='primary' and subset in ['all','quality_subset'] and dur in [5,10] and ex and co:
            for metric in ['inhibition_rate_per_min','delta_rate_per_min']:
                r[f'{metric}_between_p']=between([x[metric] for x in ex],[x[metric] for x in co])
            r['experimental_within_p']=within([x['delta_rate_per_min'] for x in ex])
            r['control_within_p']=within([x['delta_rate_per_min'] for x in co])
        results.append(r)
    for subset in ['all','quality_subset']:
        for keys in [['inhibition_rate_per_min_between_p','delta_rate_per_min_between_p'],['experimental_within_p','control_within_p']]:
            targets=[(r,k) for r in results for k in keys if r['subset']==subset and r.get(k) is not None]
            for (r,k),q in zip(targets,bh([r[k] for r,k in targets])):r[k[:-2]+'_bh_q']=q
    return results

def write_report(flyrows,summaries,manifest):
    lines=['VES041 INHIBITION: PIVOT-LIKE BOUT FREQUENCY',
      'Analysis date: November 2, 2026','',
      'INTERPRETATION',
      'There is no consistent experimental-control difference in pivot-like frequency during inhibition. In blinded round 1 the mean rates during light are nearly identical. Experimental flies show a descriptively larger pre-to-light decrease in round 1, whereas controls show the larger decrease in the unblinded 10 s cohort.',
      'The earlier activity-screened subset has lower rates during light in experimental flies, but the difference in within-fly changes is uncertain. For the 10 s active subset, the experimental and control reductions are nearly the same. A stricter in-place definition using total path length instead of net displacement does not support selective suppression.', '',
      'METHODS',
      'Primary definition: at least 45 degrees net heading change in a sliding 1 s window, net displacement <=2 mm, and absolute net heading change / total heading change >=0.7. Candidate windows are evaluated every ~0.1 s on behavior sampled at ~100 Hz. Overlapping qualifying windows are merged into one bout, so sustained turning is not counted repeatedly at each frame.',
      'Counted throughout the full 5, 10, or 240 s inhibition period and compared with an immediately preceding light-off window of equal duration. At exact 50% duty cycle, matched windows are shortened by at most one sampled frame to avoid including preceding light in the baseline. Rate = bouts / elapsed minute. These are pivot-like bouts defined from kinematics, not manually scored turns.',
      'Fly is the independent unit. Repeated cells/trials are aggregated within date/fly/cohort/protocol. Across-fly comparisons examine inhibition rate and the experimental-control difference in within-fly inhibition-minus-pre changes. Group means weight flies equally. Within-group paired tests use one rate change per fly; no per-fly significance claims are made from treating repeated presentations as independent animals.',
      'Tests: two-sided label permutations between groups and sign-flip permutations of within-fly changes. Exact enumeration is used when feasible, otherwise 49,999 seeded permutations. Benjamini-Hochberg adjustment is separate for between-group and within-group tests, and separate for all-recording versus activity-screened analyses.',
      'Activity-screened sensitivity uses the previously computed walking screen: >=50% of 1 s pre-onset samples moving, >=50% of light-off samples moving, and >=10 moving-baseline presentations. The moving-at-onset sensitivity uses >=50% samples at total speed >=3 mm/s in the preceding 1 s. No pivot measurements were used to select activity-screened recordings.',
      'Sensitivity definitions: 30 or 60 degree rotation, 1 or 3 mm net displacement, and a stricter <=2 mm total path length. The other thresholds stay fixed. The stricter path rule rejects walking loops that return near their starting point.',
      'Only genotype-assigned, readable recordings enter group comparisons. Round 1 ID 15 was excluded because processed behavior in the 5 s trial is damaged. The previously reported list contained 12 experimental IDs, but only 11 have readable long-stimulus behavior. Genotype-unassigned summer and other recordings are excluded. Unblinded 5 s data have two experimental flies and no assigned controls, and are descriptive only. The 240 s protocol has two blocks per fly and is descriptive only.',
      f'Processed recordings: {len(manifest)}. No SMP544 data were analyzed.', '',
      'GROUP SUMMARY',
      'Rates below are mean bouts per minute across flies. Delta is inhibition minus matched pre-light. Negative experimental-control difference in changes indicates a larger reduction in experimental flies.']
    for r in summaries:
        if r['definition']!='primary' or r['subset'] not in ['all','quality_subset','moving_onset']:continue
        lines+=['',f'{r["cohort"]}, {r["duration_s"]} s, subset={r["subset"]}']
        for group in ['experimental','control']:
            if not r[f'n_{group}']:continue
            lines.append(f'{group}: n={r[f"n_{group}"]}; pre={r[f"{group}_pre_rate_per_min"]:.3f}; inhibition={r[f"{group}_inhibition_rate_per_min"]:.3f}; delta={r[f"{group}_delta_rate_per_min"]:.3f}; flies decreasing={r[f"{group}_flies_decrease"]}; flies increasing={r[f"{group}_flies_increase"]}')
        if 'opto_difference' in r:
            lines.append(f'Experimental minus control inhibition rate: {r["opto_difference"]:.3f}; difference in changes: {r["difference_in_changes"]:.3f}')
        for key,value in r.items():
            if value is not None and (key.endswith('_p') or key.endswith('_bh_q')):lines.append(f'  {key}: {value:.5g}')
    lines+=['','STRICT IN-PLACE SENSITIVITY (<=2 MM TOTAL PATH LENGTH)']
    for r in summaries:
        if r['definition']=='path2' and r['subset']=='all' and 'difference_in_changes' in r:
            lines.append(f'{r["cohort"]}, {r["duration_s"]} s: experimental pre/on={r["experimental_pre_rate_per_min"]:.3f}/{r["experimental_inhibition_rate_per_min"]:.3f}; control pre/on={r["control_pre_rate_per_min"]:.3f}/{r["control_inhibition_rate_per_min"]:.3f}; difference in changes={r["difference_in_changes"]:.3f}')
    lines+=['','INDIVIDUAL FLIES (PRIMARY DEFINITION, ALL READABLE RECORDINGS)',
            'Each entry gives counts and rates from the same fly before and during inhibition.']
    for r in sorted([r for r in flyrows if r['definition']=='primary' and r['subset']=='all'],key=lambda r:(r['cohort'],r['duration_s'],r['group'],r['date'],r['daily_fly'])):
        lines.append(f'{r["cohort"]} {r["duration_s"]} s; ID {",".join(r["cohort_ids"]) or "unassigned"}; {r["date"]} daily fly {r["daily_fly"]}; {r["group"]}; {r["n_presentations"]} presentations; pre/on counts {r["pre_events"]}/{r["inhibition_events"]}; pre/on rates {r["pre_rate_per_min"]:.3f}/{r["inhibition_rate_per_min"]:.3f}; delta {r["delta_rate_per_min"]:.3f}')
    lines+=['','LIMITATIONS',
      'A low net displacement can include tight walking loops. The path-length sensitivity addresses this, but raw tracking videos have not been reviewed and event examples are illustrative rather than a blinded manual validation set. Faster forward walking can reduce classification as a low-translation pivot; both experimental and control groups accelerate during light. Cohorts differ in protocol and are not pooled. Activity-screened and threshold sensitivity analyses remain exploratory.']
    (OUT/'results_and_methods.txt').write_text('\n'.join(lines)+'\n')

def write_filtered_report(fly,groups,manifest,angle,strict,good_walkers=False):
    comparison='>' if strict else '>='
    lines=['VES041 BLINDED INHIBITION: PIVOT-LIKE BOUT FREQUENCY',
      f'Primary rotation criterion: {comparison}{angle:g} degrees within 1 s.',
      'Only blinded round 1 recordings with 5 s inhibition are included.',
      'Other primary criteria: net displacement <=2 mm and directionality >=0.7. Overlapping sliding-window detections are merged into distinct bouts. Matched pre-light and inhibition windows have equal exposure.',
      'Fly is the independent unit. All readable, genotype-assigned flies are included; activity-screened and moving-at-onset subsets are sensitivity analyses. ID 15 has damaged processed behavior and is omitted. No unblinded or SMP544 data are included.',
      'The activity screen is the previously defined >=50% pre-onset and light-off movement with >=10 moving-baseline presentations. Movement is total speed >=3 mm/s.',
      'Between-group two-sided permutation tests compare inhibition rates and within-fly inhibition-minus-pre changes. Within-group sign-flip tests use one change per fly. Up to 49,999 seeded permutations; exact enumeration when feasible. BH corrections are separate for between-group and within-group tests and for all-flies versus activity-screened analyses.',
      'The <=2 mm total-path-length sensitivity excludes tight walking loops. Raw tracking video has not been manually scored.', '', 'GROUP COMPARISONS (fly-weighted bouts/minute)']
    for r in groups:
        if r['definition']!='primary':continue
        lines+=['',r['subset']]
        for group in ['experimental','control']:
            lines.append(f'{group}: n={r[f"n_{group}"]}; pre={r[f"{group}_pre_rate_per_min"]:.3f}; inhibition={r[f"{group}_inhibition_rate_per_min"]:.3f}; change={r[f"{group}_delta_rate_per_min"]:.3f}; decreasing flies={r[f"{group}_flies_decrease"]}; increasing flies={r[f"{group}_flies_increase"]}')
        lines.append(f'Experimental minus control: inhibition rate={r["opto_difference"]:.3f}; difference in within-fly changes={r["difference_in_changes"]:.3f}')
        for k,v in r.items():
            if v is not None and (k.endswith('_p') or k.endswith('_bh_q')):lines.append(f'{k}={v:.6g}')
    lines+=['','DEFINITION SENSITIVITY (all flies)']
    for r in groups:
        if r['subset']=='all':lines.append(f'{r["definition"]}: inhibition experimental/control={r["experimental_inhibition_rate_per_min"]:.3f}/{r["control_inhibition_rate_per_min"]:.3f}; difference in changes={r["difference_in_changes"]:.3f}')
    lines+=['','INDIVIDUAL FLIES (primary definition, all presentations)']
    for r in sorted([r for r in fly if r['definition']=='primary' and r['subset']=='all'],key=lambda r:(r['group'],r['date'],r['daily_fly'])):
        lines.append(f'ID {",".join(r["cohort_ids"])}; {r["group"]}; {r["date"]} daily fly {r["daily_fly"]}; pre/on events={r["pre_events"]}/{r["inhibition_events"]}; pre/on bouts/min={r["pre_rate_per_min"]:.3f}/{r["inhibition_rate_per_min"]:.3f}; change={r["delta_rate_per_min"]:.3f}')
    if good_walkers:
        lines=[line.replace('All readable, genotype-assigned flies are included; activity-screened and moving-at-onset subsets are sensitivity analyses.',
            'Only the previously identified good walkers are included: 8 experimental and 9 control flies. Moving-at-onset and stationary-at-onset subsets are sensitivity analyses.')
            .replace('DEFINITION SENSITIVITY (all flies)','DEFINITION SENSITIVITY (selected good walkers)') for line in lines]
        lines.insert(2,'Good walkers only. Selection uses the preceding behavior screen and does not use pivot measurements or inhibition responses.')
    (OUT/'results_and_methods.txt').write_text('\n'.join(lines)+'\n')

def main():
    global OUT
    parser=argparse.ArgumentParser()
    parser.add_argument('--refresh-cohort')
    parser.add_argument('--summarize',action='store_true')
    parser.add_argument('--cohort')
    parser.add_argument('--angle',type=float,default=45)
    parser.add_argument('--strict-angle',action='store_true')
    parser.add_argument('--output',type=Path)
    parser.add_argument('--good-walkers',action='store_true')
    parser.add_argument('--from-cache',type=Path)
    args=parser.parse_args()
    if args.output:OUT=args.output.resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    definitions=DEFINITIONS if args.angle==45 else [
      ('primary',args.angle,2,'displacement'),('turn45',45,2,'displacement'),
      ('turn60',60,2,'displacement'),('disp1',args.angle,1,'displacement'),
      ('disp3',args.angle,3,'displacement'),('path2',args.angle,2,'path')]
    (OUT/'parameters.json').write_text(json.dumps({'cohort':args.cohort,'primary_angle_deg':args.angle,'strict_angle':args.strict_angle,'good_walkers':args.good_walkers,'definitions':definitions},indent=2))
    quality=json.loads((SOURCE/'behavior_quality_results.json').read_text())
    records=[r for r in quality if r['status']=='ok' and (r['group'].startswith('experimental') or r['group'].startswith('control'))]
    if args.cohort:records=[r for r in records if r['cohort']==args.cohort]
    if args.good_walkers:records=[r for r in records if r['baseline_moving_fraction']>=.5 and r['off_moving_fraction']>=.5 and r['n_moving_baseline_pulses']>=10]
    trials=[];examples=[];manifest=[];errors=[]
    if args.from_cache:
        cached_parameters=json.loads((args.from_cache/'parameters.json').read_text())
        if cached_parameters['primary_angle_deg']!=args.angle or cached_parameters['strict_angle']!=args.strict_angle:
            raise ValueError('Cached pivot definition does not match requested definition')
        def eligible(r):return (not args.cohort or r['cohort']==args.cohort) and (not args.good_walkers or r['quality_subset'])
        trials=[r for r in json.loads((args.from_cache/'presentations.json').read_text()) if eligible(r)]
        examples=[r for r in json.loads((args.from_cache/'event_examples.json').read_text()) if eligible(r)]
        manifest=[r for r in json.loads((args.from_cache/'manifest.json').read_text()) if eligible(r)]
        records=[]
    if args.refresh_cohort or args.summarize:
        trials=json.loads((OUT/'presentations.json').read_text())
        examples=json.loads((OUT/'event_examples.json').read_text())
        manifest=json.loads((OUT/'manifest.json').read_text())
        if args.summarize:records=[]
        else:
            records=[r for r in records if r['cohort']==args.refresh_cohort]
            trials=[r for r in trials if r['cohort']!=args.refresh_cohort]
            examples=[r for r in examples if r['cohort']!=args.refresh_cohort]
            manifest=[r for r in manifest if r['cohort']!=args.refresh_cohort]
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(process,r,definitions,args.strict_angle):r for r in records}
        for i,future in enumerate(as_completed(futures),1):
            try:
                result=future.result();trials.extend(result['trials']);examples.extend(result['examples']);manifest.append(result['record'])
            except Exception as e:errors.append({'path':futures[future]['behavior_file'],'error':str(e)})
            if i%10==0:print(f'Extracted {i}/{len(records)}',flush=True)
    for name,obj in [('presentations',trials),('event_examples',examples),('manifest',manifest),('errors',errors)]:
        (OUT/f'{name}.json').write_text(json.dumps(obj,allow_nan=False))
    fly=aggregate(trials)
    if args.good_walkers:fly=[r for r in fly if r['subset']!='quality_subset']
    (OUT/'fly_summary.json').write_text(json.dumps(fly,indent=2))
    print('Starting fly-level comparisons',flush=True)
    groups=inference(fly)
    (OUT/'group_summary.json').write_text(json.dumps(groups,indent=2))
    if args.cohort:write_filtered_report(fly,groups,manifest,args.angle,args.strict_angle,args.good_walkers)
    else:write_report(fly,groups,manifest)
    print('Finished',len(records),'recordings;',len(errors),'errors;',len(examples),'primary event examples',flush=True)

if __name__=='__main__':main()
