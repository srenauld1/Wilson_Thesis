"""Map existing >30-degree detections to the two example recordings."""
from pathlib import Path
import csv
import json
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / 'outputs/ves041_blinded_good_walkers_pivots_gt30_20261102'
LOCAL = ROOT / 'outputs/ves041_fly8_fly45_master_rerun_20261102'
OUT = ROOT / 'outputs/ves041_example_reorientation_timeseries_20261104'


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    manifest = json.loads((ANALYSIS / 'manifest.json').read_text())
    saved_events = json.loads((ANALYSIS / 'event_examples.json').read_text())
    presentations = json.loads((ANALYSIS / 'presentations.json').read_text())
    reports = []
    for fly_id in ['8', '45']:
        rec = next(r for r in manifest if r['cohort_ids'] == [fly_id])
        key = rec['recording_key']
        source = next((LOCAL / key).rglob('*processed_final.mat'))
        with h5py.File(source, 'r') as f:
            g = f['exptData']
            raw_t = np.asarray(g['t_ds'][:, 0])
            original_dt = (raw_t[-1] - raw_t[0]) / (len(raw_t) - 1)
            step = max(1, round(.01 / original_dt))
            dt = step * original_dt
            arrays = {k: np.asarray(g[k][::step, 0]) for k in
                      ['x', 'y', 'headingPosition', 'forwardVelocity',
                       'angularVelocity', 'totSpeed', 'optoStim_ds']}
        t = raw_t[::step] - raw_t[0]
        light = arrays['optoStim_ds'] != 0
        heading = np.rad2deg(np.unwrap(np.deg2rad(arrays['headingPosition'])))
        edges = np.diff(np.r_[False, light, False].astype(int))
        starts = np.flatnonzero(edges == 1)
        ends = np.flatnonzero(edges == -1)
        assert len(starts) == len(ends) == 40
        cache_path = ROOT / 'outputs/ves041_blinded_good_walkers_speed_linked_turns_20261102/kinematic_cache' / (key + '.npz')
        with np.load(cache_path) as cached:
            assert np.array_equal(light, cached['light'])
            assert np.allclose(arrays['forwardVelocity'], cached['forward'])
            assert np.allclose(arrays['headingPosition'], cached['heading'])
            assert np.isclose(dt, cached['dt'])
        windows, phase_starts = [], {}
        for pulse, (a, b) in enumerate(zip(starts, ends), 1):
            length = b - a
            available_off = a - (ends[pulse - 2] if pulse > 1 else 0)
            if available_off < length and (length - available_off) * dt <= .03:
                length = available_off
                b = a + length
            for phase, lo, hi in [('pre', a-length, a), ('inhibition', a, b)]:
                phase_starts[pulse, phase] = lo
                windows.append(dict(pulse=pulse, phase=phase,
                                    start_s=float(t[lo]), end_s=float(t[hi-1]+dt)))
        events = []
        for e in saved_events:
            if e['recording_key'] != key:
                continue
            lo = phase_starts[e['pulse'], e['phase']]
            s, z = lo + e['representative_start'], lo + e['representative_end']
            net = abs(heading[z] - heading[s])
            distance = np.hypot(arrays['x'][z]-arrays['x'][s], arrays['y'][z]-arrays['y'][s])
            turn_sum = np.sum(np.abs(np.diff(heading[s:z+1])))
            directionality = net / turn_sum
            assert np.isclose(net, e['rotation_deg'], atol=1e-7)
            assert np.isclose(distance, e['displacement_mm'], atol=1e-7)
            assert np.isclose(directionality, e['directionality'], atol=1e-7)
            assert net > 30 and distance <= 2 and directionality >= .7
            # The stored example path must match the native local coordinates.
            idx = np.arange(s, z+1, max(1, round(.05/dt)))
            assert np.allclose(arrays['x'][idx]-arrays['x'][s], e['x'])
            assert np.allclose(arrays['y'][idx]-arrays['y'][s], e['y'])
            mid = (s + z) // 2
            events.append(dict(pulse=e['pulse'], phase=e['phase'],
                start_s=float(t[s]), end_s=float(t[z]), marker_s=float(t[mid]),
                bout_start_s=float(t[lo+e['first_start']]),
                bout_end_s=float(t[lo+e['last_end']]),
                x_mm=float(arrays['x'][mid]), y_mm=float(arrays['y'][mid]),
                rotation_deg=net, displacement_mm=distance,
                directionality=directionality))
        events.sort(key=lambda e: e['start_s'])
        for i, e in enumerate(events, 1):
            e['event_id'] = i
        expected = sum(r['events'] for r in presentations
                       if r['recording_key'] == key and r['definition'] == 'primary')
        assert len(events) == expected
        write_csv(OUT / f'fly{fly_id}_events.csv', events)
        write_csv(OUT / f'fly{fly_id}_scored_windows.csv', windows)
        write_csv(OUT / f'fly{fly_id}_light_intervals.csv', [
            dict(pulse=i, start_s=float(t[a]), end_s=float(t[b-1]+dt))
            for i, (a,b) in enumerate(zip(starts,ends),1)])
        columns = [t, arrays['x'], arrays['y'], arrays['forwardVelocity'],
                   arrays['angularVelocity'], heading-heading[0], arrays['totSpeed'], light]
        np.savetxt(OUT / f'fly{fly_id}_timeseries.csv', np.column_stack(columns),
                   delimiter=',', header='time_s,x_mm,y_mm,forward_mm_s,yaw_deg_s,heading_change_deg,total_speed_mm_s,light_on',
                   comments='', fmt='%.9g')
        counts = {p: sum(e['phase'] == p for e in events) for p in ['pre','inhibition']}
        reports.append(dict(fly_id=fly_id, label='VES041_1 x '+('ACR' if fly_id=='8' else 'GFP'),
                            events=len(events), counts=counts, dt_s=dt,
                            local_source=str(source), original_source=rec['source']))
        print(f'Fly {fly_id}: {len(events)} verified events; {counts}', flush=True)
    (OUT / 'verification.json').write_text(json.dumps(reports, indent=2))
    (OUT / 'README.txt').write_text(
        'Fly 8 (ACR) and Fly 45 (GFP): full recordings and existing reorientation detections.\n'
        'Detection source: ves041_blinded_good_walkers_pivots_gt30_20261102, primary definition.\n'
        'Rule: >30 degrees net heading change in a nominal 1 s window, <=2 mm net displacement, directionality >=0.7.\n'
        'Directionality is absolute net heading change divided by summed absolute heading changes.\n'
        'Windows advance by approximately 0.1 s; overlapping detections merge into one bout.\n'
        'Orange segments show the stored representative window (largest net turn in the bout).\n'
        'Numbers mark the midpoint of that window, with the same ID on path, traces, and event CSV.\n'
        'Purple shading marks every light-on interval. Gray bars show matched pre-light scoring windows.\n'
        'Only the matched 5 s pre-light and 5 s light-on windows were scored. Other times are unscored, not zero-event intervals.\n'
        'The full 20-minute recording is plotted, including unscored periods.\n'
        'Native x/y geometry is preserved, with no rotation, reflection, or translation.\n'
        'Stored velocity traces are plotted at the existing analysis sampling (~95 Hz), without additional smoothing.\n'
        'Heading is unwrapped and shown relative to the initial heading.\n'
        'Every saved event was checked against local x/y and heading values and the saved counts.\n'
        'The time series match the existing analysis cache. No new events were detected or selected.\n'
        'Detail panels cover consecutive 200 s intervals; event numbering is retained across panels.\n'
        'The events folder contains one SVG per event, with 3 s context before and after its representative window.\n')


if __name__ == '__main__':
    main()
