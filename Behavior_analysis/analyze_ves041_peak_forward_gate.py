"""Good-walker turns gated by forward velocity at peak absolute yaw."""
from pathlib import Path
from itertools import combinations, product
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
import csv
import json
import time
import numpy as np
import h5py

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
CACHE=ROOT/'outputs/ves041_blinded_good_walkers_speed_linked_turns_20261102/kinematic_cache'
OUT=ROOT/'outputs/ves041_peak_forward_gate_20261104'
MODES=['raw_signed','raw_absolute','processed_signed','angle_direction_only']


def read(ds,start=None,stop=None,step=1):
    sl=slice(start,stop,step)
    if ds.ndim==1:return np.asarray(ds[sl],float)
    return np.asarray(ds[0,sl] if ds.shape[0]==1 else ds[sl,0],float)


def candidates(h,dt):
    w=round(1/dt);stride=max(1,round(.1/dt))
    starts=np.arange(0,len(h)-w,stride);ends=starts+w
    net=abs(h[ends]-h[starts]);cum=np.r_[0,np.cumsum(abs(np.diff(h)))];total=cum[ends]-cum[starts]
    direction=np.divide(net,total,out=np.zeros_like(net),where=total>0)
    selected=np.flatnonzero((net>30)&(direction>=.7));events=[]
    for i in selected:
        props=dict(representative_start=int(starts[i]),representative_end=int(ends[i]),rotation_deg=float(net[i]),directionality=float(direction[i]))
        if events and starts[i]<=events[-1]['last_end']:
            events[-1]['last_end']=int(ends[i])
            if net[i]>events[-1]['rotation_deg']:events[-1].update(props)
        else:events.append(dict(first_start=int(starts[i]),last_end=int(ends[i]),**props))
    return events


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def between(a,b):
    vals=np.r_[a,b];n=len(a);obs=float(np.mean(a)-np.mean(b));all_sum=float(vals.sum());extreme=0;total=0
    for ix in combinations(range(len(vals)),n):
        sm=float(vals[list(ix)].sum());stat=sm/n-(all_sum-sm)/(len(vals)-n)
        extreme+=abs(stat)>=abs(obs)-1e-10;total+=1
    return dict(difference=obs,p=extreme/total,n_permutations=total)


def process(rec):
    with np.load(CACHE/(rec['recording_key']+'.npz')) as z:
        h=np.rad2deg(np.unwrap(np.deg2rad(z['heading'])));fw=z['forward'];speed=z['speed'];light=z['light'];dt=float(z['dt'])
    with h5py.File(rec['source'],'r') as f:
        g=f['exptData'];n=g['forwardVelocity'].size
        td=next(g[k] for k in ['t_ds','tDS'] if k in g and g[k].size==n)
        original_dt=float(read(td,1,2)[0]-read(td,0,1)[0]);t0=float(read(td,0,1)[0]);step=round(dt/original_dt)
        angular=read(g['angularVelocity'],step=step)
        assert len(angular)==len(h)
        raw_t0=float(read(g['t'],0,1)[0]);raw_t_dt=float(read(g['t'],1,2)[0]-raw_t0)
        ratio=g['t'].size/g['forwardVelocity_raw'].size
        assert np.isclose(ratio,round(ratio)) and g['angularVelocity_raw'].size==g['forwardVelocity_raw'].size
        raw_dt=raw_t_dt*ratio
        ed=np.diff(np.r_[False,light,False].astype(int));ons=np.flatnonzero(ed==1);offs=np.flatnonzero(ed==-1)
        assert len(ons)==40
        trials=[];events=[]
        for pulse,(a,b) in enumerate(zip(ons,offs),1):
            length=b-a;available=a-(offs[pulse-2] if pulse>1 else 0)
            if available<length and (length-available)*dt<=.03:length=available;b=a+length
            pre=a-length
            assert pre>=0 and not np.any(light[pre:a])
            moving=bool(np.mean(speed[max(0,a-round(1/dt)):a]>=3)>=.5)
            for phase,lo,hi in [('pre',pre,a),('light_on',a,b)]:
                es=candidates(h[lo:hi],dt)
                accepted={mode:0 for mode in MODES}
                for e in es:
                    s=lo+e['first_start'];z=lo+e['last_end']
                    raw_a=int(round((t0+s*dt-raw_t0)/raw_dt));raw_b=int(round((t0+z*dt-raw_t0)/raw_dt))+1
                    rv=read(g['angularVelocity_raw'],raw_a,raw_b)
                    assert len(rv)>0 and np.isfinite(rv).all()
                    pi=int(np.argmax(abs(rv)));raw_index=raw_a+pi
                    raw_fwd=float(read(g['forwardVelocity_raw'],raw_index,raw_index+1)[0])
                    processed_index=s+int(np.argmax(abs(angular[s:z+1])))
                    processed_fwd=float(fw[processed_index])
                    flags=dict(raw_signed=raw_fwd<=1,raw_absolute=abs(raw_fwd)<=1,
                               processed_signed=processed_fwd<=1,angle_direction_only=True)
                    for mode,flag in flags.items():accepted[mode]+=int(flag)
                    events.append(dict(fly_id=','.join(rec['cohort_ids']),group=rec['group'],recording_key=rec['recording_key'],
                        pulse=pulse,phase=phase,start_s=t0+s*dt,end_s=t0+z*dt,
                        representative_start_s=t0+(lo+e['representative_start'])*dt,
                        representative_end_s=t0+(lo+e['representative_end'])*dt,
                        rotation_deg=e['rotation_deg'],directionality=e['directionality'],
                        peak_yaw_time_s=raw_t0+raw_index*raw_dt,peak_yaw_deg_s=float(rv[pi]),forward_at_peak_mm_s=raw_fwd,
                        processed_peak_time_s=t0+processed_index*dt,processed_forward_at_peak_mm_s=processed_fwd,
                        **flags))
                for mode,count in accepted.items():
                    trials.append(dict(fly_id=','.join(rec['cohort_ids']),group=rec['group'],recording_key=rec['recording_key'],
                        pulse=pulse,phase=phase,definition=mode,events=count,seconds=(hi-lo)*dt,moving_onset=moving))
    return trials,events


def main():
    OUT.mkdir(exist_ok=True)
    records=json.loads((SRC/'manifest.json').read_text());assert len(records)==17
    assert sum(r['group']=='experimental' for r in records)==8
    old=json.loads((SRC/'fly_summary.json').read_text())
    old={','.join(r['cohort_ids']):r for r in old if r['subset']=='all' and r['definition']=='primary'}
    trials=[];events=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        pending={pool.submit(process,r):r for r in records}
        for i,future in enumerate(as_completed(pending),1):
            r=pending[future];tr,ev=future.result();trials.extend(tr);events.extend(ev)
            print(f'{i}/17 ID {r["cohort_ids"]}: {len(ev)} candidate bouts, {sum(e["raw_signed"] for e in ev)} pass raw signed cutoff',flush=True)
            (OUT/'progress.json').write_text(json.dumps(dict(completed=i,flies=17,last_id=r['cohort_ids'])))
    flies=[]
    for mode in MODES:
        for rec in records:
            id=','.join(rec['cohort_ids']);tt=[t for t in trials if t['fly_id']==id and t['definition']==mode]
            row=dict(fly_id=id,group=rec['group'],definition=mode,presentations=len(tt)//2)
            for phase in ['pre','light_on']:
                ts=[t for t in tt if t['phase']==phase];count=sum(t['events'] for t in ts);sec=sum(t['seconds'] for t in ts)
                row[phase+'_events']=count;row[phase+'_seconds']=sec;row[phase+'_rate_per_min']=60*count/sec
            row['change_per_min']=row['light_on_rate_per_min']-row['pre_rate_per_min']
            row['old_pre_events']=old[id]['pre_events'];row['old_light_on_events']=old[id]['inhibition_events']
            row['old_pre_rate_per_min']=old[id]['pre_rate_per_min'];row['old_light_on_rate_per_min']=old[id]['inhibition_rate_per_min']
            row['old_change_per_min']=old[id]['delta_rate_per_min']
            flies.append(row)
    groups=[]
    for mode in MODES:
        for group in ['experimental','control']:
            ff=[r for r in flies if r['definition']==mode and r['group']==group]
            ge=[e for e in events if e['group']==group and e[mode]]
            groups.append(dict(group=group,definition=mode,n_flies=len(ff),
                pre_events=sum(r['pre_events'] for r in ff),light_on_events=sum(r['light_on_events'] for r in ff),
                pre_mean_rate=float(np.mean([r['pre_rate_per_min'] for r in ff])),
                light_on_mean_rate=float(np.mean([r['light_on_rate_per_min'] for r in ff])),
                mean_change=float(np.mean([r['change_per_min'] for r in ff])),
                flies_decreased=sum(r['change_per_min']<-1e-8 for r in ff),
                flies_increased=sum(r['change_per_min']>1e-8 for r in ff),
                flies_unchanged=sum(abs(r['change_per_min'])<=1e-8 for r in ff),
                accepted_peaks_forward_below_minus1=sum(e['forward_at_peak_mm_s']<-1 for e in ge)))
    ex=[r for r in flies if r['definition']=='raw_signed' and r['group']=='experimental']
    co=[r for r in flies if r['definition']=='raw_signed' and r['group']=='control']
    stats={metric:between([r[metric] for r in ex],[r[metric] for r in co]) for metric in ['light_on_rate_per_min','change_per_min']}
    pv=np.array([v['p'] for v in stats.values()]);order=np.argsort(pv);q=np.minimum.accumulate((pv[order]*len(pv)/np.arange(1,len(pv)+1))[::-1])[::-1]
    for i,ix in enumerate(order):list(stats.values())[ix]['bh_q']=float(min(1,q[i]))
    for name,rows in [('fly_summary',flies),('group_summary',groups),('presentations',trials),('candidate_events',events)]:
        write_csv(OUT/(name+'.csv'),rows);(OUT/(name+'.json')).write_text(json.dumps(rows,indent=2,default=lambda x:x.item()))
    (OUT/'inference.json').write_text(json.dumps(stats,indent=2,default=lambda x:x.item()))
    (OUT/'parameters.json').write_text(json.dumps(dict(cohort='same 17 blinded round-1 good walkers; 5 s protocol',angle='>30 degrees in nominal 1 s windows',directionality_min=.7,displacement_requirement='removed',candidate_window_step_s=.1,merge='merge overlapping angle/directionality windows before applying forward-velocity gate once per bout',primary='signed raw forward velocity <=1 mm/s at largest absolute raw yaw peak within merged bout',raw_sampling='native stored raw velocity sampling; no additional smoothing',exposure='matched 5 s pre-light and light-on phases; each candidate window fully contained in phase',sensitivity=['absolute raw forward velocity <=1 mm/s','signed processed forward <=1 mm/s at processed yaw peak','angle/directionality only']),indent=2))
    print(json.dumps(groups,indent=2,default=lambda x:x.item()),flush=True);print(json.dumps(stats,indent=2,default=lambda x:x.item()),flush=True)

if __name__=='__main__':main()
