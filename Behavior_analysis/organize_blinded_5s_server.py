"""Copy-only server organization. Never modifies an original recording."""
import argparse
import hashlib
import json
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SERVER = Path('/Volumes/Neurobio/wilsonlab/Sophia/behavior/RubyACR Data')
TARGET = SERVER / '5s stim'
LOCAL = REPO / 'outputs/blinded_5s_server_organization_20261102'

CRITERIA = '''GOOD WALKING: blinded round 1, 5-second optogenetic stimulation

This is the pre-existing descriptive activity screen, applied per recording.
Movement means total speed >=3 mm/s (not forward velocity alone).
A recording must meet ALL THREE criteria:
1. >=50% of pooled samples in the 1-second pre-light baseline windows moving.
2. >=50% of samples during all light-off periods moving.
3. >=10 presentations with a moving baseline: >=50% moving samples in the
   preceding 1 second, with >=95% finite speed, forward velocity and heading.

These are screening choices, not validated biological thresholds. They do not
use reorientation frequency, genotype, or the size/direction of the optogenetic
response. Classification is recording-specific, not a permanent fly phenotype.
Genotype-unassigned flies are included and explicitly identified in the manifest.
Passing this screen does not establish raw-video tracking quality.

The 'not good walking' folder contains recordings failing this screen. Its
'unassessable' subfolder contains unreadable files, NOT confirmed poor walkers.
All originals remain untouched. Only existing processed files were copied;
no raw acquisition files were copied or reconstructed. Only blinded round 1
trial01B/stim5 recordings are included, not trial01A/stim2 or unblinded data.
'''

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def exclusive_text(path, text):
    if path.exists():
        if path.read_text() != text:
            raise RuntimeError(f'Refusing to overwrite different existing file: {path}')
    else:
        path.write_text(text)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--execute',action='store_true'); args=ap.parse_args()
    audit=json.loads((REPO/'outputs/ves041_inhibition_duration_audit_20261102/behavior_quality_results.json').read_text())
    records=[]
    for r in audit:
        if r['cohort']!='blinded round 1' or r['duration_s']!=5:
            continue
        candidates=[Path(p) for p in r['paths'] if '/RubyACR Data/round1/' in p and p.endswith('_pro.mat')]
        if not candidates: # damaged ID15 was represented by raw path in old audit
            candidates=[SERVER/'round1/fly015/cell01_processed/2026_01_21_fly05_cell01_trial01B pattern20 stim5 function1_pro.mat']
            assert r['cohort_ids']==['15']
        assert len(candidates)==1 and candidates[0].is_file()
        source=candidates[0]
        good=r['status']=='ok' and r['baseline_moving_fraction']>=.5 and r['off_moving_fraction']>=.5 and r['n_moving_baseline_pulses']>=10
        classification='good walking' if good else 'not good walking'
        parent=TARGET/classification
        if r['status']!='ok': parent=parent/'unassessable'
        dest=parent/('fly'+r['cohort_ids'][0])/source.parent.name/source.name
        record={k:r.get(k) for k in ['cohort_ids','date','daily_fly','cell','trial','group','labels','status','baseline_moving_fraction','off_moving_fraction','n_moving_baseline_pulses']}
        record.update(source=str(source),destination=str(dest),classification=classification,source_bytes=source.stat().st_size,analysis_status='pending' if r['status']=='ok' else 'unassessable: truncated processed source')
        records.append(record)
    assert len(records)==42 and len({r['destination'] for r in records})==42
    LOCAL.mkdir(parents=True,exist_ok=True)
    (LOCAL/'plan.json').write_text(json.dumps(records,indent=2))
    print(f'Plan: {len(records)} files, {sum(r["source_bytes"] for r in records)/1e9:.2f} GB',flush=True)
    if not args.execute:return
    for name in ['good walking','not good walking']:(TARGET/name).mkdir(parents=True,exist_ok=True)
    exclusive_text(TARGET/'good walking/GOOD_WALKING_CRITERIA.txt',CRITERIA)
    exclusive_text(TARGET/'README.txt',CRITERIA+'\nSee recording_manifest.json and master_status.json for copy integrity and analysis status.\n')
    # Copy all source files unchanged, with checksums and no overwrite.
    for i,r in enumerate(records,1):
        source=Path(r['source']);dest=Path(r['destination']);dest.parent.mkdir(parents=True,exist_ok=True)
        digest=sha(source)
        if dest.exists():
            if sha(dest)!=digest:raise RuntimeError(f'Destination mismatch; will not overwrite {dest}')
        else:
            shutil.copy2(source,dest)
            if sha(dest)!=digest:raise RuntimeError(f'Copy integrity failure: {dest}')
        r['sha256']=digest
        print(f'COPY VERIFIED {i}/42 ID {r["cohort_ids"][0]} -> {r["classification"]}',flush=True)
        exclusive_text(dest.parent/'recording_info.json',json.dumps(r,indent=2)+'\n')
        (LOCAL/'copy_manifest.json').write_text(json.dumps(records,indent=2))
    exclusive_text(TARGET/'recording_manifest.json',json.dumps(records,indent=2)+'\n')
    (LOCAL/'copy_manifest.json').write_text(json.dumps(records,indent=2))
    # Save the exact scripts used, without changing any existing originals.
    provenance=TARGET/'analysis provenance';provenance.mkdir(exist_ok=True)
    for name in ['master_behavior_analysis.m','run_blinded_5s_server_master.m','process_fictrac_panels.m','plotExpt_sliced_opto_velocity.m','plot_flatpath_opto.m','plot_flatpath_aligned_opto.m','compute_absolute_circular_diff.m','detect_local_peaks.m','save_plot_with_title_as_filename.m','organize_blinded_5s_server.py']:
        p=REPO/'Behavior_analysis'/name;d=provenance/name
        if d.exists():
            if sha(d)!=sha(p):raise RuntimeError(f'Provenance mismatch: {d}')
        else:shutil.copy2(p,d)
    print('COPY PHASE COMPLETE',flush=True)

if __name__=='__main__':main()
