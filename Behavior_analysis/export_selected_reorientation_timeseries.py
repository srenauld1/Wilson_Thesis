"""Selected high-resolution event traces, including the following light pulse."""
import numpy as np
import h5py
import json
import render_ves041_reorientations_clean as p


def export(fly_id,event_id,raw_velocity=True,window_s=25):
    r=p.load(fly_id)
    e=next(e for e in r['e'] if int(e['event_id'])==event_id)
    current=next(a for a in r['light'] if int(a['pulse'])==int(e['pulse']))
    following=next(a for a in r['light'] if a['start_s']>e['end_s'])
    lo=current['start_s']-7.5
    hi=lo+window_s
    assert lo<=e['start_s']<e['end_s']<=hi
    assert lo<=current['start_s']<current['end_s']<=hi
    velocity_record=r
    if raw_velocity:
        verification=json.loads((p.SRC/'verification.json').read_text())
        source=next(a['local_source'] for a in verification if int(a['fly_id'])==fly_id)
        with h5py.File(source,'r') as f:
            g=f['exptData'];t=g['t']
            dt=float(t[1,0]-t[0,0]);t0=float(t[0,0])
            step=max(1,round(.001/dt))
            a=max(0,int(np.floor((lo-t0)/dt)));b=min(t.size,int(np.ceil((hi-t0)/dt))+1)
            data=np.zeros(len(range(a,b,step)),dtype=[('time_s',float),('forward_mm_s',float),('yaw_deg_s',float)])
            data['time_s']=t[a:b:step,0]-t0
            data['forward_mm_s']=g['forwardVelocity_raw'][a:b:step,0]
            data['yaw_deg_s']=g['angularVelocity_raw'][a:b:step,0]
        assert np.isfinite(data['forward_mm_s']).all() and np.isfinite(data['yaw_deg_s']).all()
        assert data['forward_mm_s'].min()>=-20 and data['forward_mm_s'].max()<=50
        assert abs(data['yaw_deg_s']).max()<=1000
        velocity_record={**r,'d':data}
    s=p.SVG(1050,1050)
    p.title(s,r,f'Event #{event_id}  ·  {e["marker_s"]:.2f} s  ·  25-second window · full light-on pulse')
    p.legend(s,108)
    if raw_velocity:s.text(80,141,'Raw velocity fields · sampled at 1 kHz · no additional smoothing',16,p.GRAY)
    p.trace(s,velocity_record,lo,hi,'forward_mm_s',(110,200,870,155),(-20,50),[-20,0,20,40],'Forward velocity (mm/s)')
    p.trace(s,velocity_record,lo,hi,'yaw_deg_s',(110,450,870,155),(-1000,1000),[-1000,-500,0,500,1000],'Yaw velocity (deg/s)',event_labels=True)
    d=r['d'];h=d['heading_change_deg'][(d['time_s']>=lo)&(d['time_s']<=hi)]
    lower=np.floor(min(h)/30)*30;upper=np.ceil(max(h)/30)*30
    p.trace(s,r,lo,hi,'heading_change_deg',(110,700,870,155),(lower-10,upper+10),np.linspace(lower,upper,4),'Heading change (deg; original processed heading)',xticks=True)
    s.text(980,925,'Time from recording start (s)',17,p.GRAY,'end')
    s.text(80,1010,f'Net turn {e["rotation_deg"]:.1f}°  ·  Net displacement {e["displacement_mm"]:.2f} mm  ·  Directionality {e["directionality"]:.2f}',18)
    suffix='_less_smoothed_25s_square' if raw_velocity else '_25s_square'
    target=p.OUT/'events'/f'Fly{fly_id:02d}_event{event_id:03d}_timeseries_with_following_opto{suffix}.svg'
    s.save(target)
    if raw_velocity:
        np.savetxt(target.with_suffix('.csv'),np.column_stack([data[k] for k in data.dtype.names]),delimiter=',',header=','.join(data.dtype.names),comments='')
        target.with_suffix('.json').write_text(json.dumps(dict(source=source,velocity_fields=['forwardVelocity_raw','angularVelocity_raw'],sampling_hz=1/(dt*step),window_s=window_s,canvas_px=[1050,1050],light_onset_s=current['start_s'],additional_smoothing='none',heading='original processed heading',detections='existing primary >30-degree events, unchanged',start_s=lo,end_s=hi),indent=2))
    print(f'{target}: {lo:.3f}–{hi:.3f} s; following pulse {following["start_s"]:.3f}–{following["end_s"]:.3f} s')


if __name__=='__main__':
    export(8,1)
    export(8,6)
    export(45,31)
