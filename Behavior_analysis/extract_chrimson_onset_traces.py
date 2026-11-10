"""Raw-velocity onset traces: average presentations within fly, then flies."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import csv
import json
import re
import traceback
import h5py
import numpy as np

REPO=Path(__file__).resolve().parents[1]
SRC=REPO/'.build_ves041_activation/reproducibility/source_data'
OUT=REPO/'outputs/chrimson_onset_traces_20261103'
CACHE=OUT/'fly_cache'
SERVER=Path('/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers')

def read(ds,start=None,stop=None,step=1):
    sl=slice(start,stop,step)
    if ds.ndim==1:return np.asarray(ds[sl],dtype=float).ravel()
    return np.asarray(ds[0,sl] if ds.shape[0]==1 else ds[sl,0],dtype=float).ravel()

def prepare():
    text=(SRC/'cschrim_behaviorcrunch.R').read_text()
    mapping={int(n):g for ids,g in re.findall(r'Fly %in% c\(([^)]+)\) ~ "([^"]+)"',text) for n in re.findall(r'\d+',ids)}
    manifest={}
    with (SRC/'recording_manifest.csv').open() as f:
        for r in csv.DictReader(f):
            if r['dataset']=='ves041_activation':manifest[int(r['fly'].split('fly')[-1])]=Path(r['path'])
    manifest[1]=SERVER/'fly01/cell01_processed/trial01Astim2/2026_03_04_fly01_cell01_trial01A_processed_final.mat'
    manifest[6]=SERVER/'fly06/cell01_processed/trial01Astim2/2026_03_06_fly02_cell01_trial01A_processed_final.mat'
    records=[]
    for id,genotype in sorted(mapping.items()):
        final=manifest[id]
        prefix=final.name.replace('_processed_final.mat','')
        base=final.parent.parent/(prefix+' pattern20 stim2 function1_pro.mat')
        records.append(dict(fly=id,genotype=genotype,source=str(base),alternate=str(final)))
    assert len(records)==84
    return records

def extract(rec):
    id=rec['fly'];npz=CACHE/f'fly{id:03d}.npz';meta=CACHE/f'fly{id:03d}.json'
    if npz.exists() and meta.exists():return json.loads(meta.read_text())
    source=Path(rec['source'])
    if not source.is_file():source=Path(rec['alternate'])
    with h5py.File(source,'r') as f:
        g=f['exptData']
        fw=g['forwardVelocity_raw'];av=g['angularVelocity_raw'];op=g['optoStim'];t=g['t']
        n=fw.size;m=op.size
        assert n==av.size and m==t.size
        light_dt=float(np.median(np.diff(read(t,0,100))))
        light_t0=float(read(t,0,1)[0]);light_last=float(read(t,m-1,m)[0])
        assert np.isclose(light_t0+(m-1)*light_dt,light_last,atol=.001)
        ratio=m/n;assert np.isclose(ratio,round(ratio))
        ratio=int(round(ratio));signal_dt=light_dt*ratio
        # Coarse light-edge search only. Every onset/offset is refined in the
        # original acquisition-rate light command before segment extraction.
        stride=max(1,int(round(m/800000)))
        coarse=read(op,step=stride)>0
        edges=np.diff(np.r_[False,coarse,False].astype(np.int8))
        starts=np.flatnonzero(edges==1);ends=np.flatnonzero(edges==-1)
        bins=300;fw_trials=[];rot_trials=[];pulses=[];durations=[];skipped=[];missing=0
        for pulse,(a,b) in enumerate(zip(starts,ends),1):
            lo=max(0,a*stride-stride);hi=min(m,a*stride+1)
            onset=lo+int(np.flatnonzero(read(op,lo,hi)>0)[0])
            lo=max(onset,b*stride-stride);hi=min(m,b*stride+1)
            off=read(op,lo,hi)<=0
            offset=lo+int(np.flatnonzero(off)[0]) if np.any(off) else m
            duration=(offset-onset)*light_dt
            if not .95<=duration<=1.05:
                raise ValueError(f'Fly {id} pulse {pulse}: unexpected light duration {duration}')
            onset_signal=onset/ratio
            start=int(np.ceil(onset_signal-1/signal_dt));stop=int(np.ceil(onset_signal+2/signal_dt))
            if start<0 or stop>n:
                skipped.append(dict(pulse=pulse,reason='Insufficient full -1 to +2 s window'));continue
            rel=(np.arange(start,stop)-onset_signal)*signal_dt
            bi=np.floor((rel+1)/.01+1e-8).astype(int)
            valid=(bi>=0)&(bi<bins)
            forward=read(fw,start,stop);rotation=np.abs(read(av,start,stop))
            missing+=int(np.sum(~np.isfinite(forward))+np.sum(~np.isfinite(rotation)))
            def bin_mean(values):
                keep=valid&np.isfinite(values)
                sums=np.bincount(bi[keep],weights=values[keep],minlength=bins)
                counts=np.bincount(bi[keep],minlength=bins)
                return np.divide(sums,counts,out=np.full(bins,np.nan),where=counts>0)
            fw_trials.append(bin_mean(forward));rot_trials.append(bin_mean(rotation))
            pulses.append(pulse);durations.append(duration)
    forward=np.asarray(fw_trials);rotation=np.asarray(rot_trials)
    assert forward.shape==rotation.shape and forward.shape[0]>0 and forward.shape[1]==300
    time=np.arange(300)*.01-.995
    assert np.isfinite(forward).all() and np.isfinite(rotation).all(), 'Missing entire time bins require review'
    np.savez_compressed(npz,time=time,forward_trials=forward,rotation_trials=rotation,pulses=pulses,
                        forward_mean=np.mean(forward,axis=0),rotation_mean=np.mean(rotation,axis=0))
    result={**rec,'source_used':str(source),'n_pulses_detected':len(starts),'n_pulses_included':len(pulses),
            'skipped_pulses':skipped,'signal_rate_hz':1/signal_dt,'light_duration_median_s':float(np.median(durations)),
            'nonfinite_samples':missing,'cache':str(npz)}
    meta.write_text(json.dumps(result,indent=2))
    return result

def main():
    CACHE.mkdir(parents=True,exist_ok=True)
    records=prepare();(OUT/'source_manifest.json').write_text(json.dumps(records,indent=2))
    results=[];errors=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        pending={pool.submit(extract,r):r for r in records}
        for future in as_completed(pending):
            rec=pending[future]
            try:
                r=future.result();results.append(r)
                print(f'EXTRACTED {len(results)}/84 fly {r["fly"]}: {r["n_pulses_included"]} presentations',flush=True)
            except Exception:
                errors.append({**rec,'error':traceback.format_exc()});print(f'ERROR fly {rec["fly"]}: {errors[-1]["error"]}',flush=True)
            (OUT/'extraction_results.json').write_text(json.dumps(sorted(results,key=lambda x:x['fly']),indent=2))
            (OUT/'errors.json').write_text(json.dumps(errors,indent=2))
    rows=[];summary=[]
    for r in sorted(results,key=lambda x:x['fly']):
        with np.load(r['cache']) as a:
            for metric in ['forward','rotation']:
                values=a[metric+'_mean']
                for t,y in zip(a['time'],values):rows.append(dict(fly=r['fly'],genotype=r['genotype'],metric=metric,time=float(t),value=float(y),presentations=r['n_pulses_included']))
    for genotype in sorted({r['genotype'] for r in results}):
        group=[r for r in results if r['genotype']==genotype]
        for metric in ['forward','rotation']:
            matrix=np.asarray([np.load(r['cache'])[metric+'_mean'] for r in group])
            means=np.mean(matrix,axis=0);sem=np.std(matrix,axis=0,ddof=1)/np.sqrt(len(group))
            for j in range(300):summary.append(dict(genotype=genotype,metric=metric,time=float(time_grid()[j]),mean=float(means[j]),sem=float(sem[j]),n=len(group)))
    (OUT/'fly_traces.json').write_text(json.dumps(rows))
    (OUT/'genotype_mean_sem.json').write_text(json.dumps(summary))
    methods=[
        'CSCHRIMSON ACTIVATION: ONSET-ALIGNED RAW VELOCITY TRACES',
        'Source cohort: all 84 genotype-assigned flies in the archived activation table. One existing keeper recording per fly, matching the established manifest; IDs 1 and 6 added from their intact keeper recordings.',
        'Existing original function1_pro.mat files preferred. No master rerun or modification of server sources. This also recovers ID39 from its intact original processed file.',
        'Fields: forwardVelocity_raw (signed forward velocity, mm/s); abs(angularVelocity_raw) (rotational speed, deg/s). Absolute value taken BEFORE trial/bin averaging.',
        'Onsets and offsets detected from optoStim; coarse scan refined at original acquisition rate. Window -1 to +2 s relative to each measured light onset. Only full windows are eligible.',
        'Consecutive 10-ms bins; mean of original samples in each bin. No additional smoothing, baseline subtraction, normalization, moving-only filter or outlier rejection.',
        'Average all eligible presentation traces within fly. Then arithmetic mean and sample-SD/sqrt(n) across fly averages within genotype at each bin. These are pointwise SEM bands, not confidence intervals or presentation-level uncertainty.',
        f'Extracted {len(results)}/84 flies; {sum(r["n_pulses_included"] for r in results)} presentations; errors: {len(errors)}.',
        'Retain source manifest, per-fly NPZ presentation traces, extraction metadata and plotted values for audit. Shaded interval 0-1 s is nominal stimulation; actual median durations retained per recording.',
        'No statistical tests added.'
    ]
    (OUT/'methods.txt').write_text('\n'.join(methods)+'\n')
    print(methods[-3],flush=True)
    if errors:raise RuntimeError('Some flies failed; review errors before labeling this an all-fly plot.')

def time_grid():return np.arange(300)*.01-.995

if __name__=='__main__':main()
