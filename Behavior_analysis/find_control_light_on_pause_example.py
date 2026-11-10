"""Review original detected control light-on reorientations for a brief pause."""
from pathlib import Path
import csv,json
import numpy as np
import h5py
from analyze_ves041_peak_forward_gate import read

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
CACHE=ROOT/'outputs/ves041_blinded_good_walkers_speed_linked_turns_20261102/kinematic_cache'
OUT=ROOT/'outputs/ves041_original_control_pause_example_20261104'

def main():
    OUT.mkdir(exist_ok=True)
    manifest=json.loads((SRC/'manifest.json').read_text())
    events=json.loads((SRC/'event_examples.json').read_text())
    controls=[r for r in manifest if r['group']=='control']
    rows=[]
    for i,rec in enumerate(controls,1):
        with np.load(CACHE/(rec['recording_key']+'.npz')) as z:
            dt=float(z['dt']);light=z['light']
        edge=np.diff(np.r_[False,light,False].astype(int));ons=np.flatnonzero(edge==1);offs=np.flatnonzero(edge==-1)
        with h5py.File(rec['source'],'r') as f:
            g=f['exptData'];t=g['t'];r0=float(read(t,0,1)[0]);rdt=float(read(t,1,2)[0]-r0)
            ratio=t.size/g['forwardVelocity_raw'].size;rdt*=ratio
            step=max(1,round(.001/rdt));sample_dt=step*rdt
            for j,e in enumerate([v for v in events if v['recording_key']==rec['recording_key'] and v['phase']=='inhibition'],1):
                a=int(ons[e['pulse']-1]);b=int(offs[e['pulse']-1]);start=(a+e['representative_start'])*dt;end=(a+e['representative_end'])*dt
                lower=start-1.5;upper=end+1.5;aa=round((lower-r0)/rdt);bb=round((upper-r0)/rdt)+1
                forward=read(g['forwardVelocity_raw'],aa,bb,step);yaw=read(g['angularVelocity_raw'],aa,bb,step)
                tt=r0+(aa+np.arange(len(forward))*step)*rdt
                inside=(tt>=start)&(tt<=end);ix=np.flatnonzero(inside)[np.argmax(abs(yaw[inside]))];peak_t=float(tt[ix])
                def avg(lo,hi,absolute=False):
                    v=forward[(tt>=peak_t+lo)&(tt<peak_t+hi)]
                    return float(np.mean(abs(v) if absolute else v))
                before=avg(-1,-.4);after=avg(.4,1)
                row=dict(fly_id=','.join(rec['cohort_ids']),recording_key=rec['recording_key'],pulse=e['pulse'],within_fly_light_event=j,
                    start_s=start,end_s=end,light_start_s=a*dt,light_end_s=b*dt,
                    rotation_deg=e['rotation_deg'],displacement_mm=e['displacement_mm'],directionality=e['directionality'],
                    peak_yaw_time_s=peak_t,peak_yaw_deg_s=float(yaw[ix]),forward_at_yaw_peak=float(forward[ix]),
                    mean_abs_forward_near_peak=avg(-.1,.1,True),mean_forward_near_peak=avg(-.1,.1),
                    before_forward=before,after_forward=after,
                    pre_clearance_s=peak_t-a*dt,post_clearance_s=b*dt-peak_t,
                    source=rec['source'])
                # A pause assessment is for illustration selection only; detector criteria are unchanged.
                smooth=np.convolve(abs(forward),np.ones(100)/100,mode='same')
                ixlo=ix
                while ixlo>0 and smooth[ixlo]<1.5:ixlo-=1
                ixhi=ix
                while ixhi+1<len(smooth) and smooth[ixhi]<1.5:ixhi+=1
                row['pause_duration_s']=(ixhi-ixlo)*sample_dt if smooth[ix]<1.5 else 0
                rows.append(row)
        print(f'{i}/9 control ID {rec["cohort_ids"]}: inspected original light-on events',flush=True)
    with (OUT/'all_control_light_on_candidates.csv').open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (OUT/'all_control_light_on_candidates.json').write_text(json.dumps(rows,indent=2))
    eligible=[r for r in rows if r['mean_abs_forward_near_peak']<=2 and r['before_forward']>=3 and r['after_forward']>=3 and r['pre_clearance_s']>=1 and r['post_clearance_s']>=1 and 40<=r['rotation_deg']<=120]
    median_angle=float(np.median([v['rotation_deg'] for v in eligible])) if eligible else None
    eligible.sort(key=lambda r:(abs(r['rotation_deg']-median_angle),r['mean_abs_forward_near_peak']))
    print('Eligible:',len(eligible));print(json.dumps(eligible[:10],indent=2))
    print('Best pause candidates:')
    shortlist=sorted([r for r in rows if 35<=r['rotation_deg']<=150],key=lambda r:-(min(r['before_forward'],r['after_forward'])-r['mean_abs_forward_near_peak']))
    print(json.dumps(shortlist[:8],indent=2))

if __name__=='__main__':main()
