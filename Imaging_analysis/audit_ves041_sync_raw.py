"""Read-only raw CSV/DAQ audit; MATLAB renders the resulting small MAT files.

Recompute exact coded-marker matching, then compare every saved row, edge and
update timestamp. Raw recordings and alignment sidecars are never modified.
"""
import csv
import json
import sys
from pathlib import Path
import traceback
import numpy as np
import h5py
from scipy.io import loadmat, savemat

sys.path.insert(0, '/Users/sophiarenauld/Documents/GitHub/scopa')
from bergii_sync import csv_native, match_markers


def matlab_value(file, node):
    # Same small-struct decoder as bergii_import, without its TIFF dependency.
    if isinstance(node, h5py.Group):
        return {key: matlab_value(file, value) for key, value in node.items()}
    if node.attrs.get('MATLAB_empty', 0): return []
    value = node[()]
    if h5py.check_dtype(ref=value.dtype):
        return [matlab_value(file, file[ref]) for ref in value.ravel()]
    if node.attrs.get('MATLAB_class') == b'char':
        return ''.join(chr(int(c)) for c in value.ravel()).rstrip('\x00')
    return value.item() if value.size == 1 else value.ravel().tolist()


def audit(item, out):
    rec = item['recording']
    folder = Path(item['input']).parent
    receipt_path = folder/(rec+'_bergii_behavior_import.json')
    if not receipt_path.exists():
        receipt_path = folder/(rec+'_bergii_import.json')
    receipt = json.loads(receipt_path.read_text())
    run = receipt['runId']
    daq = folder/(run+'_imaging_daq.mat')
    source_csv = folder/(run+'_fictrac_frames.csv')
    sidecar = folder/(rec+'_bergii_.mat')
    protected = [daq, source_csv, sidecar, Path(item['input'])]
    before = {str(p): (p.stat().st_size, p.stat().st_mtime_ns) for p in protected}
    saved = loadmat(sidecar, simplify_cells=True)['bergii']
    rows, emitted, cf, cs, native = csv_native(source_csv)
    with h5py.File(daq, 'r') as h:
        meta = matlab_value(h, h['metadata'])
        names = meta['channels']['name']
        roles = ['fictrac_update_toggle_v', 'fictrac_sync_code_v',
                 'fictrac_forward_position_wrapped_v']
        assert all(names.count(s) == 1 for s in roles), 'Missing/ambiguous sync channel'
        inds = [names.index(s) for s in roles]
        ds = h['samples']; count = ds.shape[1]
        assert count == int(meta['actualScans'])
        assert meta['sharedAcquisition'] and not meta.get('simulated')
        for i in inds:
            assert int(meta['channels']['matrixIndex'][i])-int(meta['channels']['indexBase'][i]) == i
        rate = float(meta.get('actualRateHz', meta['rateHz']))
        fast = np.empty(count, dtype=bool); slow = np.empty(count, dtype=bool)
        forward = np.empty(count, dtype=np.float32)
        # Read each compressed chunk once, rather than separately per channel.
        for start in range(0, count, 1000000):
            stop = min(count, start+1000000)
            block = ds[:, start:stop]
            assert np.all(np.isfinite(block[inds])), 'Nonfinite raw sync/forward sample'
            fast[start:stop] = block[inds[0]] > 2.5
            slow[start:stop] = block[inds[1]] > 2.5
            forward[start:stop] = block[inds[2]]
            assert np.allclose(h['timeSeconds'][0, start:stop], np.arange(start, stop)/rate, rtol=0, atol=1e-8)
    edges_all = np.flatnonzero(np.diff(fast.astype(np.int8)) != 0)+1
    off, trim, n = match_markers(fast[edges_all], slow[np.minimum(edges_all+3, count-1)], cf, cs)
    edges = edges_all[:n]; ix = emitted[off:off+n]
    fdaq = fast[edges].astype(float); sdaq = slow[np.minimum(edges+3, count-1)].astype(float)
    mismatch_fast = int(np.count_nonzero(fdaq != cf[off:off+n]))
    mismatch_slow = int(np.count_nonzero(sdaq != cs[off:off+n]))
    first = int(saved['originalSampleStart'])-1
    times = (edges-first)/rate
    stored_edges = np.asarray(saved['markerEdgesOriginalSamplesOneBased']).reshape(-1)
    stored_rows = np.asarray(saved['csvRowIndicesOneBased']).reshape(-1)
    stored_times = np.asarray(saved['updateTimes']).reshape(-1)
    assert len(stored_edges) == n == len(stored_rows) == len(stored_times), 'Saved match length differs'
    edge_errors = int(np.count_nonzero(edges+1 != stored_edges))
    row_errors = int(np.count_nonzero(ix+1 != stored_rows))
    timing_error = np.max(np.abs(times-stored_times))
    assert edge_errors == 0 and row_errors == 0 and timing_error < 1e-9, 'Saved mapping differs from raw reconstruction'
    # Cross-check the exact traces actually used by the imaging analysis.
    with h5py.File(item['input'], 'r') as h:
        imported = matlab_value(h, h['dq'])
        if isinstance(imported, list): imported = imported[0]
        analysis_sync = imported['bergii']
    for key in ['csvRowIndicesOneBased', 'markerEdgesOriginalSamplesOneBased', 'updateTimes',
                'behaviorVolumeValid', 'trackingGapVolumesOneBased', 'incompleteVolumesOneBased']:
        assert np.array_equal(np.asarray(analysis_sync.get(key, [])).reshape(-1), np.asarray(saved.get(key, [])).reshape(-1)), 'A2P input has stale mapping or coverage: '+key
    mid = (edges[:-1]+edges[1:])//2
    forward_error = np.angle(np.exp(1j*(forward[mid]*2*np.pi/10-native['bf'][off:off+n-1])))
    p95 = float(np.percentile(np.abs(forward_error), 95))
    assert p95 <= .05, 'Forward-voltage cross-check exceeds importer QC limit'
    gaps = np.diff(times); gap_ix = np.flatnonzero(gaps > .15)
    invalid = np.asarray(saved['incompleteVolumesOneBased']).reshape(-1)
    gap_volumes = np.asarray(saved.get('trackingGapVolumesOneBased', [])).reshape(-1)
    zoom = {}
    for label, center in [('start', max(0., times[0])+.1), ('end', times[-1]-1.1)]:
        take = (times >= center) & (times <= center+1)
        lo = max(0, int((center+first/rate)*rate)); hi = min(count, lo+int(rate))
        z = np.arange(lo, hi, max(1, int(rate/2000)))
        zoom[label+'RawTime'] = (z-first)/rate
        zoom[label+'RawFast'] = fast[z].astype(float)
        zoom[label+'RawSlow'] = slow[z].astype(float)
        zoom[label+'CsvTime'] = times[take]
        zoom[label+'CsvFast'] = cf[off:off+n][take]
        zoom[label+'CsvSlow'] = cs[off:off+n][take]
    # Local offset error curve supplements the exhaustive exact-match check.
    shifts = np.arange(-20, 21); errors = np.full(len(shifts), np.nan)
    for j, shift in enumerate(shifts):
        start = off+shift
        if 0 <= start <= len(cf)-n:
            errors[j] = np.mean((fdaq != cf[start:start+n]) | (sdaq != cs[start:start+n]))
    flags = []
    if not meta.get('complete'): flags.append('incomplete DAQ acquisition')
    if receipt['acquisitionStatus'] != 'complete': flags.append('acquisition '+receipt['acquisitionStatus'])
    if len(gap_ix): flags.append(str(len(gap_ix))+' tracking gaps >150 ms')
    if len(invalid): flags.append(str(len(invalid))+' volumes lack complete behavior coverage')
    if trim: flags.append('one verified terminal cleanup edge excluded')
    status = 'PASS with coverage notes' if flags else 'PASS'
    result = dict(recording=rec, status=status, matchedUpdates=n, fastMismatches=mismatch_fast,
        slowMismatches=mismatch_slow, savedRowMismatches=row_errors, savedEdgeMismatches=edge_errors,
        maxSavedTimeErrorSec=float(timing_error), forwardErrorP95Rad=p95,
        csvOffsetEmittedRows=off, excludedCleanupEdges=trim, trackingGaps=len(gap_ix),
        maximumUpdateGapSec=float(np.max(gaps)), incompleteBehaviorVolumes=len(invalid),
        trackingGapVolumes=len(gap_volumes), sourceDaqComplete=bool(meta.get('complete')),
        acquisitionStatus=receipt['acquisitionStatus'], notes='; '.join(flags), error='')
    savemat(out/(rec+'_sync_audit.mat'), dict(**zoom, recording=rec, updateTimes=times,
        intervalTimes=times[1:], updateIntervals=gaps, timingError=times-stored_times,
        forwardErrorTimes=times[:-1], forwardErrorRad=forward_error,
        offsetShifts=shifts, offsetMismatchFraction=errors, notes=result['notes'],
        matchedUpdates=n, fastMismatches=mismatch_fast, slowMismatches=mismatch_slow,
        incompleteVolumes=invalid, trackingGapVolumes=gap_volumes), do_compression=True)
    for p in protected:
        assert before[str(p)] == (p.stat().st_size, p.stat().st_mtime_ns), 'Source changed during audit'
    return result


def main():
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    manifest = Path('/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis/hand_drawn_15s_ave_f_20261026/run_report.csv')
    items = list(csv.DictReader(manifest.open())); assert len(items) == 16
    prior_path = out/'sync_audit_report.json'
    prior = ({r['recording']: r for r in json.loads(prior_path.read_text())}
             if prior_path.exists() and '--retry-failed' in sys.argv else {})
    results = []
    for i, item in enumerate(items):
        print('AUDIT', i+1, '/16', item['recording'], flush=True)
        try:
            cached = prior.get(item['recording'], {})
            if cached.get('status', '').startswith('PASS'):
                result = cached
            else:
                result = audit(item, out)
        except Exception:
            result = dict(recording=item['recording'], status='FAIL', error=traceback.format_exc())
        results.append(result)
        (out/'sync_audit_report.json').write_text(json.dumps(results, indent=2)+'\n')
        keys = list(dict.fromkeys(k for r in results for k in r))
        with (out/'sync_audit_report.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=keys); writer.writeheader(); writer.writerows(results)
        print(result, flush=True)
    if any(r['status'] == 'FAIL' for r in results): raise SystemExit(1)


if __name__ == '__main__': main()
