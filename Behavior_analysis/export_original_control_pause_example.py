"""Package an illustrative control event from the original spatial detector."""
from pathlib import Path
import json
import h5py
import numpy as np
from analyze_ves041_peak_forward_gate import read

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
OUT=ROOT/'outputs/ves041_original_control_pause_example_20261104'

def main():
    reviewed=json.loads((OUT/'all_control_light_on_candidates.json').read_text())
    chosen=next(e for e in reviewed if e['fly_id']=='30' and e['pulse']==31)
    es=json.loads((SRC/'event_examples.json').read_text())
    stored=next(e for e in es if e['cohort_ids']==['30'] and e['phase']=='inhibition' and e['pulse']==31)
    meta=dict(**chosen,selection='Illustrative original control light-on event with a brief forward pause, forward locomotion before and after, and a moderate turn angle. Not selected to estimate a group effect.',
       source_analysis=str(SRC),original_criteria=dict(net_turn_deg='>30 within nominal 1 s',net_displacement_mm='<=2',directionality='>=0.7'),
       velocity='Stored forwardVelocity_raw and angularVelocity_raw sampled at 1 kHz; no added smoothing',
       angle_percentile_among_control_light_on_events=100*sum(e['rotation_deg']<=stored['rotation_deg'] for e in es if e['group']=='control' and e['phase']=='inhibition')/66)
    with h5py.File(chosen['source'],'r') as f:
        g=f['exptData'];td=next(g[k] for k in ['t_ds','tDS'] if k in g and g[k].size==g['headingPosition'].size)
        odt=float(read(td,1,2)[0]-read(td,0,1)[0]);step=round(.0105/odt)
        h=np.rad2deg(np.unwrap(np.deg2rad(read(g['headingPosition'],step=step))));x=read(g['x'],step=step);y=read(g['y'],step=step)
        dt=odt*step;s=round(chosen['start_s']/dt);z=round(chosen['end_s']/dt)
        net=abs(h[z]-h[s]);disp=float(np.hypot(x[z]-x[s],y[z]-y[s]));di=float(net/sum(abs(np.diff(h[s:z+1]))))
        assert np.isclose(net,stored['rotation_deg']) and np.isclose(disp,stored['displacement_mm']) and np.isclose(di,stored['directionality'])
        assert net>30 and disp<=2 and di>=.7
        rt=g['t'];rdt=float(read(rt,1,2)[0]-read(rt,0,1)[0]);r0=float(read(rt,0,1)[0]);rawstep=round(.001/rdt)
        lo=chosen['light_start_s']-7.5;hi=lo+25;aa=round((lo-r0)/rdt);bb=round((hi-r0)/rdt)+1
        t=read(rt,aa,bb,rawstep)-r0;fw=read(g['forwardVelocity_raw'],aa,bb,rawstep);av=read(g['angularVelocity_raw'],aa,bb,rawstep)
        np.savetxt(OUT/'Fly30_pulse31_raw_velocity.csv',np.column_stack([t,fw,av]),delimiter=',',header='time_s,forward_mm_s,yaw_deg_s',comments='',fmt='%.12g')
    meta['window_25s']=[lo,hi];meta['window_zoom']=[chosen['light_start_s']-2.5,chosen['light_start_s']+7.5]
    meta['verification']=dict(recomputed_rotation_deg=float(net),recomputed_displacement_mm=disp,recomputed_directionality=di,passed=True)
    (OUT/'selected_example.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
