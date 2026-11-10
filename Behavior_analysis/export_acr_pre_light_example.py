"""Package a detected ACR pre-light turn followed by a zero-event light pulse."""
import json
from pathlib import Path
import h5py
import numpy as np
from analyze_ves041_peak_forward_gate import read
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
OUT=ROOT/'outputs/ves041_original_acr_pre_light_example_20261104'

def main():
 candidates=json.loads((OUT/'candidate_review.json').read_text());r=next(e for e in candidates if e['fly_id']=='22C' and e['pulse']==31)
 es=json.loads((SRC/'event_examples.json').read_text());original=next(e for e in es if e['cohort_ids']==['22C'] and e['phase']=='pre' and e['pulse']==31 and np.isclose(e['rotation_deg'],r['rotation_deg']))
 trials=json.loads((SRC/'presentations.json').read_text());on=next(e for e in trials if e['cohort_ids']==['22C'] and e['phase']=='inhibition' and e['pulse']==31 and e['definition']=='primary')
 assert on['events']==0
 with h5py.File(r['source'],'r') as f:
  g=f['exptData'];td=next(g[k] for k in ['t_ds','tDS'] if k in g and g[k].size==g['headingPosition'].size);dt0=float(read(td,1,2)[0]-read(td,0,1)[0]);step=round(.0105/dt0);dt=step*dt0
  h=np.rad2deg(np.unwrap(np.deg2rad(read(g['headingPosition'],step=step))));x=read(g['x'],step=step);y=read(g['y'],step=step);a=round(r['start_s']/dt);b=round(r['end_s']/dt)
  net=float(abs(h[b]-h[a]));disp=float(np.hypot(x[b]-x[a],y[b]-y[a]));direction=float(net/sum(abs(np.diff(h[a:b+1]))))
  assert np.isclose(net,original['rotation_deg']) and np.isclose(disp,original['displacement_mm']) and np.isclose(direction,original['directionality'])
  assert net>30 and disp<=2 and direction>=.7
  rt=g['t'];rdt=float(read(rt,1,2)[0]-read(rt,0,1)[0]);rs=round(.001/rdt)
  lo=r['light_start_s']-7.5;hi=lo+25;aa=round(lo/rdt);bb=round(hi/rdt)+1
  time=read(rt,aa,bb,rs);fw=read(g['forwardVelocity_raw'],aa,bb,rs);yaw=read(g['angularVelocity_raw'],aa,bb,rs)
 np.savetxt(OUT/'ACR_Fly22C_pulse31_raw_velocity.csv',np.column_stack([time,fw,yaw]),delimiter=',',header='time_s,forward_mm_s,yaw_deg_s',comments='',fmt='%.12g')
 metadata=dict(**r,window_zoom=[r['light_start_s']-2.5,r['light_start_s']+7.5],window_25s=[lo,hi],
  light_on_detected_events=0,source_analysis=str(SRC),criteria='Original >30 deg net heading change in nominal 1 s, <=2 mm net displacement, directionality >=0.7',
  verification=dict(net_turn_deg=net,net_displacement_mm=disp,directionality=direction,passed=True),
  selection='Illustrative original-detected pre-light event, zero events during the following pulse, high forward walking during light, all visible on the same 10 s scale as the control. Selected to illustrate this pattern, not to estimate group effects.',
  gaussian_sigma_s=.010,raw_velocity_sampling_hz=1/(rdt*rs))
 (OUT/'selected_example.json').write_text(json.dumps(metadata,indent=2));print(json.dumps(metadata,indent=2))

if __name__=='__main__':main()
