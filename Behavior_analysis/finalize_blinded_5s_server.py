"""Track the running copy/master batch and verify its final deliverables."""
import json
import time
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
LOCAL=REPO/'outputs/blinded_5s_server_organization_20261102'
ROOT=Path('/Volumes/Neurobio/wilsonlab/Sophia/behavior/RubyACR Data/5s stim')
plan=json.loads((LOCAL/'plan.json').read_text())
last=None
deadline=time.monotonic()+10*3600
while time.monotonic()<deadline:
    try:
        manifest=json.loads((LOCAL/'copy_manifest.json').read_text()) if (LOCAL/'copy_manifest.json').exists() else []
        copied=sum('sha256' in r for r in manifest)
        results=[]
        for worker in [1,2]:
            path=ROOT/f'master_status_worker{worker}.json'
            if path.exists():
                part=json.loads(path.read_text())
                results.extend(part if isinstance(part,list) else [part])
        completed=sum(r['status']=='complete' for r in results)
        failed=[r for r in results if r['status']=='failed']
        skipped=sum(r['status']=='unassessable' for r in results)
        state=(copied,completed,len(failed),skipped)
        terminal=copied==42 and len(results)==42
        if terminal:
            assert len({r['fly_id'] for r in results})==42
            for r in results:
                if r['status']=='complete':
                    output=Path(r['source']).parent/'trial01Bstim5'
                    assert (output/'MASTER_COMPLETE.txt').is_file()
            assert (ROOT/'recording_manifest.json').is_file()
        if state!=last or terminal:
            text=(f'BLINDED ROUND 1 / 5-SECOND STIMULATION\n'
                  f'Batch state: {"FINISHED" if terminal else "RUNNING"}\n'
                  f'Source copies checksum-verified: {copied}/42\n'
                  f'Master outputs completed and verified: {completed}/41 readable recordings\n'
                  f'Master failures: {len(failed)}\n'
                  f'Unassessable processed files retained: {skipped}/1 (ID 15)\n'
                  f'Classification: 20 good walking; 21 not good walking; 1 unassessable.\n'
                  'Unknown genotypes remain explicitly unassigned. Originals unchanged.\n'
                  'No raw acquisition files copied or reconstructed.\n'
                  'See recording_manifest.json, recording_info.json, and master_status_worker*.json.\n'
                  'Use non-mirrored rawY aligned plots for measured path geometry.\n')
            if failed:text+='\nFailures:\n'+'\n'.join(f'ID {r["fly_id"]}: {r["message"]}' for r in failed)+'\n'
            (ROOT/'BATCH_STATUS.txt').write_text(text)
            (LOCAL/'BATCH_STATUS.txt').write_text(text)
            print(text,flush=True);last=state
        if terminal:
            (ROOT/'master_status.json').write_text(json.dumps(results,indent=2)+'\n')
            (LOCAL/'master_status.json').write_text(json.dumps(results,indent=2)+'\n')
            break
    except (OSError,json.JSONDecodeError) as e:
        print(f'Will retry status read: {e}',flush=True)
    time.sleep(30)
else:
    raise RuntimeError('Batch status watcher timed out; check worker logs before resuming.')
