"""Find an original ACR pre-light event followed by high-forward light-on walking."""
from pathlib import Path
import csv,json
import numpy as np
import h5py
from analyze_ves041_peak_forward_gate import read
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
CACHE=ROOT/'outputs/ves041_blinded_good_walkers_speed_linked_turns_20261102/kinematic_cache'
OUT=ROOT/'outputs/ves041_original_acr_pre_light_example_20261104'

def main():
 OUT.mkdir(exist_ok=True)
 records=json.loads((SRC/'manifest.json').read_text());ev=json.loads((SRC/'event_examples.json').read_text());tr=json.loads((SRC/'presentations.json').read_text())
 oncounts={(r['recording_key'],r['pulse']):r['events'] for r in tr if r['definition']=='primary' and r['phase']=='inhibition'}
 rows=[]
 for rec in records:
  if rec['group']!='experimental':continue
  z=np.load(CACHE/(rec['recording_key']+'.npz'));fw=z['forward'];l=z['light'];dt=float(z['dt']);ed=np.diff(np.r_[False,l,False].astype(int));ons=np.flatnonzero(ed==1);offs=np.flatnonzero(ed==-1)
  for e in ev:
   if e['recording_key']!=rec['recording_key'] or e['phase']!='pre' or oncounts[e['recording_key'],e['pulse']]!=0:continue
   a=ons[e['pulse']-1];b=offs[e['pulse']-1];lo=a-(b-a)
   start=(lo+e['representative_start'])*dt;end=(lo+e['representative_end'])*dt
   if start<a*dt-2.5 or end>a*dt:continue
   v=fw[a:b]
   rows.append(dict(fly_id=','.join(e['cohort_ids']),recording_key=rec['recording_key'],pulse=e['pulse'],start_s=start,end_s=end,light_start_s=a*dt,light_end_s=b*dt,
      rotation_deg=e['rotation_deg'],displacement_mm=e['displacement_mm'],directionality=e['directionality'],
      on_mean_forward=float(np.mean(v)),on_fraction_forward_gt3=float(np.mean(v>3)),on_median_forward=float(np.median(v)),on_min_forward=float(min(v)),source=rec['source']))
 rows.sort(key=lambda e:-(e['on_fraction_forward_gt3']*5+min(e['on_mean_forward'],15)))
 for r in rows[:12]:
  with h5py.File(r['source']) as f:
   g=f['exptData'];rdt=float(read(g['t'],1,2)[0]-read(g['t'],0,1)[0]);step=round(.001/rdt)
   a=round((r['light_start_s']-2.5)/rdt);b=round((r['light_start_s']+7.5)/rdt)+1
   t=read(g['t'],a,b,step);v=read(g['forwardVelocity_raw'],a,b,step);yaw=read(g['angularVelocity_raw'],a,b,step)
   kernel=np.exp(-.5*(np.arange(-40,41)/10)**2);kernel/=kernel.sum();smooth=np.convolve(v,kernel,'same');ys=np.convolve(yaw,kernel,'same');on=(t>=r['light_start_s'])&(t<r['light_end_s'])
   r.update(raw_on_mean=float(np.mean(v[on])),smoothed_on_fraction_forward_gt3=float(np.mean(smooth[on]>3)),smoothed_forward_max=float(max(smooth)),smoothed_forward_min=float(min(smooth)),smoothed_yaw_abs_max=float(max(abs(ys))),raw_on_median=float(np.median(v[on])))
  print(json.dumps(r),flush=True)
 (OUT/'candidate_review.json').write_text(json.dumps(rows,indent=2))
 print('Candidates:',len(rows),flush=True)

if __name__=='__main__':main()
