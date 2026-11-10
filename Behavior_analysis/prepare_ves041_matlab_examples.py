"""Portable inputs for MATLAB versions of the selected 25-second examples."""
import csv
import json
from pathlib import Path
import numpy as np
import h5py
import render_ves041_reorientations_clean as p

OUT=p.ROOT/'outputs/ves041_matlab_examples_20261104'

def main():
    (OUT/'source_data').mkdir(exist_ok=True,parents=True)
    verification=json.loads((p.SRC/'verification.json').read_text())
    configs=[]
    for fly,event_id,pulse in [(8,1,None),(8,6,None),(45,31,None),(45,None,32)]:
        r=p.load(fly)
        chosen=next((e for e in r['e'] if int(e['event_id'])==event_id),None)
        pulse=pulse or int(chosen['pulse'])
        on=next(a for a in r['light'] if int(a['pulse'])==pulse)
        lo=on['start_s']-7.5;hi=lo+25
        name=f'Fly{fly:02d}_'+(f'event{event_id:03d}' if event_id else f'light_bout{pulse:02d}')+'_25s_square'
        source=next(a['local_source'] for a in verification if int(a['fly_id'])==fly)
        with h5py.File(source,'r') as f:
            g=f['exptData'];t=g['t'];dt=float(t[1,0]-t[0,0]);t0=float(t[0,0])
            step=round(.001/dt);a=int(np.ceil(lo/dt));b=int(np.floor(hi/dt))+1
            tt=np.asarray(t[a:b:step,0])-t0
            raw=np.column_stack([tt,g['forwardVelocity_raw'][a:b:step,0],g['angularVelocity_raw'][a:b:step,0]])
        np.savetxt(OUT/'source_data'/f'{name}_velocity.csv',raw,delimiter=',',header='time_s,forward_mm_s,yaw_deg_s',comments='',fmt='%.12g')
        d=r['d'];take=(d['time_s']>=lo)&(d['time_s']<=hi)
        np.savetxt(OUT/'source_data'/f'{name}_heading.csv',np.column_stack([d['time_s'][take],d['heading_change_deg'][take]]),delimiter=',',header='time_s,heading_change_deg',comments='',fmt='%.12g')
        events=[e for e in r['e'] if e['start_s']<=hi and e['end_s']>=lo]
        light=[a for a in r['light'] if a['start_s']<hi and a['end_s']>lo]
        configs.append(dict(name=name,fly=fly,genotype=r['label'],event_id=event_id or 0,pulse=pulse,
          start_s=lo,end_s=hi,light=light,events=events,source=source,
          velocity_csv=f'{name}_velocity.csv',heading_csv=f'{name}_heading.csv'))
        print(name,lo,hi,'events',[e['event_id'] for e in events],flush=True)
    (OUT/'source_data/config.json').write_text(json.dumps(configs,indent=2))

if __name__=='__main__':main()
