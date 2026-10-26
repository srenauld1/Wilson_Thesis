"""Compute full-recording per-pixel maps from saved MATLAB v7.3 stacks.

HDF5 stores MATLAB [y,x,z,t] in reversed [t,z,x,y] order. Process
one x band at a time, retaining every time sample and every z plane.
"""
from pathlib import Path
import argparse
import json
import time
import h5py
import numpy as np


def compute(source, output, band_width=32):
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(source, 'r') as src:
        stack = src['stack']
        nt, nz, nx, ny = stack.shape
        assert stack.dtype == np.uint16, 'Expected finite unsigned registered images'
        maps = {name: np.empty((nz, nx, ny), np.float64)
                for name in ('mean_fluorescence', 'max_fluorescence')}
        cache = output / 'pixel_fluorescence_maps.h5'
        cached = False
        if cache.exists():
            with h5py.File(cache, 'r') as previous:
                cached = (previous.attrs.get('source_file') == str(source.resolve())
                          and previous.attrs.get('source_mtime_ns') == source.stat().st_mtime_ns
                          and all(name in previous and previous[name].shape == (ny, nx, nz) for name in maps))
                if cached:
                    for name in maps:
                        maps[name] = previous[name][:].transpose(2, 1, 0)
            if cached:
                print(f'{source.parent.name}: using complete mean/max cache', flush=True)
        for x0 in ([] if cached else range(0, nx, band_width)):
            x1 = min(nx, x0 + band_width)
            # Contiguous time dimension permits efficient temporal reductions.
            data = np.ascontiguousarray(np.moveaxis(stack[:, :, x0:x1, :], 0, -1))
            maps['mean_fluorescence'][:, x0:x1] = np.mean(data, axis=-1, dtype=np.float64)
            maps['max_fluorescence'][:, x0:x1] = np.max(data, axis=-1)
            del data
            print(f'{source.parent.name}: x {x1}/{nx}; {time.monotonic()-start:.1f} s', flush=True)
        maps['f0'] = maps['mean_fluorescence'].copy()
        valid = maps['f0'] > 0
        dff = np.full_like(maps['f0'], np.nan)
        np.divide(maps['max_fluorescence'] - maps['f0'], maps['f0'],
                  out=dff, where=valid)
        assert np.all(dff[valid] >= 0)
        maps['max_dff'] = dff
        # Transpose exactly once into explicit y,x,z arrays in the map file.
        with h5py.File(output / 'pixel_fluorescence_maps.h5', 'w') as dst:
            dst.attrs['source_file'] = str(source.resolve())
            dst.attrs['source_mtime_ns'] = source.stat().st_mtime_ns
            for name, array in maps.items():
                ds = dst.create_dataset(name, data=array.transpose(2, 1, 0), compression='gzip')
                ds.attrs['axis_order'] = 'y,x,z'
            dst.create_dataset('valid_baseline', data=valid.transpose(2, 1, 0))
        mean_upper = float(np.percentile(maps['mean_fluorescence'], 99.5, method='hazen'))
        dff_upper = float(np.percentile(dff[valid], 99.5, method='hazen'))
        meta = dict(source=str(source), source_shape_tzxy=list(stack.shape),
                    output_axis_order='y,x,z', volumes=nt, slices=nz,
                    baseline='Each pixel temporal mean fluorescence over every volume. Fifth-percentile baselines were checked in trial 11 and were zero at all pixels, so the mean definition is used consistently for all three recordings.',
                    mean='Mean fluorescence over all volumes, no additional smoothing',
                    max_dff='(max_t(F) - F0) / F0; exact temporal maximum, no additional smoothing',
                    invalid_baseline='F0 <= 0 is NaN in max_dff; no epsilon division',
                    valid_baseline_pixels=int(valid.sum()), total_pixels=int(valid.size),
                    display='One shared scale over all four slices per recording and map; upper limit is 99.5th percentile, numerical maps retain all values',
                    mean_display_limits=[0, mean_upper], max_dff_display_limits=[0, dff_upper],
                    max_dff_full_range=[float(np.nanmin(dff)), float(np.nanmax(dff))])
        (output / 'pixel_fluorescence_map_settings.json').write_text(json.dumps(meta, indent=2)+'\n')
        print(f'COMPLETE {source.parent.name}: {nt} volumes, {nz} slices, {time.monotonic()-start:.1f} s', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path('/Users/sophiarenauld/stacks/2026-10-20'))
    args = parser.parse_args()
    files = sorted(p for p in args.root.glob('*fly1*/*_s_.mat')
                   if any(s in p.parent.name for s in ('heatbarrier', 'stripes')))
    assert len(files) == 3, f'Expected three selected heat recordings, found {len(files)}'
    for source in files:
        compute(source, source.parent / 'smr_analysis' / 'pixel_fluorescence_maps')


if __name__ == '__main__':
    main()
