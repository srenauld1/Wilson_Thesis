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


def compute(source, output, band_width=32, bin_factor=1):
    start = time.monotonic()
    output.mkdir(parents=True, exist_ok=True)
    with h5py.File(source, 'r') as src:
        stack = src['stack']
        nt, nz, nx, ny = stack.shape
        assert nx % bin_factor == 0 and ny % bin_factor == 0 and band_width % bin_factor == 0
        bx, by = nx // bin_factor, ny // bin_factor
        assert stack.dtype == np.uint16, 'Expected finite unsigned registered images'
        maps = {name: np.empty((nz, bx, by), np.float64)
                for name in ('mean_fluorescence', 'max_fluorescence', 'f0')}
        cache = output / 'pixel_fluorescence_maps.h5'
        cached = False
        if cache.exists():
            with h5py.File(cache, 'r') as previous:
                cached = (previous.attrs.get('spatial_bin', 1) == bin_factor
                          and previous.attrs.get('baseline_percentile') == 5
                          and previous.attrs.get('source_file') == str(source.resolve())
                          and previous.attrs.get('source_mtime_ns') == source.stat().st_mtime_ns
                          and all(name in previous and previous[name].shape == (by, bx, nz) for name in maps))
                if cached:
                    for name in maps:
                        maps[name] = previous[name][:].transpose(2, 1, 0)
            if cached:
                print(f'{source.parent.name}: using complete fifth-percentile map cache', flush=True)
        for x0 in ([] if cached else range(0, nx, band_width)):
            x1 = min(nx, x0 + band_width)
            # Contiguous time dimension permits efficient temporal reductions.
            raw = stack[:, :, x0:x1, :]
            if bin_factor > 1:
                raw = raw.reshape(nt, nz, (x1-x0)//bin_factor, bin_factor, by, bin_factor).mean(axis=(3, 5), dtype=np.float64)
            data = np.ascontiguousarray(np.moveaxis(raw, 0, -1))
            del raw
            b0, b1 = x0 // bin_factor, x1 // bin_factor
            maps['mean_fluorescence'][:, b0:b1] = np.mean(data, axis=-1, dtype=np.float64)
            maps['max_fluorescence'][:, b0:b1] = np.max(data, axis=-1)
            maps['f0'][:, b0:b1] = np.percentile(
                data, 5, axis=-1, method='hazen', overwrite_input=True)
            del data
            print(f'{source.parent.name}: x {x1}/{nx}; {time.monotonic()-start:.1f} s', flush=True)
        valid = maps['f0'] > 0
        dff = np.full_like(maps['f0'], np.nan)
        np.divide(maps['max_fluorescence'] - maps['f0'], maps['f0'],
                  out=dff, where=valid)
        assert np.all(dff[valid] >= 0)
        maps['max_dff'] = dff
        # Transpose exactly once into explicit y,x,z arrays in the map file.
        with h5py.File(output / 'pixel_fluorescence_maps.h5', 'w') as dst:
            dst.attrs['spatial_bin'] = bin_factor
            dst.attrs['baseline_percentile'] = 5
            dst.attrs['source_file'] = str(source.resolve())
            dst.attrs['source_mtime_ns'] = source.stat().st_mtime_ns
            for name, array in maps.items():
                ds = dst.create_dataset(name, data=array.transpose(2, 1, 0), compression='gzip')
                ds.attrs['axis_order'] = 'y,x,z'
            dst.create_dataset('valid_baseline', data=valid.transpose(2, 1, 0))
        mean_upper = float(np.percentile(maps['mean_fluorescence'], 99.5, method='hazen'))
        dff_upper = float(np.percentile(dff[valid], 99.5, method='hazen')) if valid.any() else 1.0
        meta = dict(source=str(source), source_shape_tzxy=list(stack.shape),
                    output_axis_order='y,x,z', output_shape_yxz=[by, bx, nz],
                    spatial_bin=bin_factor, binning='Nonoverlapping block mean fluorescence at each time point, before temporal percentile or maximum; no temporal smoothing',
                    volumes=nt, slices=nz,
                    baseline='Each pixel temporal 5th percentile over every volume, including zeros. Hazen interpolation matches MATLAB prctile. No alternative baseline or epsilon substitution.',
                    mean='Mean fluorescence over all volumes, no additional smoothing',
                    max_dff='(max_t(F) - F0) / F0; exact temporal maximum, no additional smoothing',
                    invalid_baseline='F0 <= 0 is NaN in max_dff; no epsilon division',
                    valid_baseline_pixels=int(valid.sum()), total_pixels=int(valid.size),
                    display='One shared scale over all four slices per recording and map; upper limit is 99.5th percentile, numerical maps retain all values',
                    mean_display_limits=[0, mean_upper], max_dff_display_limits=[0, dff_upper],
                    max_dff_full_range=[float(np.nanmin(dff)), float(np.nanmax(dff))] if valid.any() else None)
        (output / 'pixel_fluorescence_map_settings.json').write_text(json.dumps(meta, indent=2)+'\n')
        print(f'COMPLETE {source.parent.name}: {nt} volumes, {nz} slices, {time.monotonic()-start:.1f} s', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path('/Users/sophiarenauld/stacks/2026-10-20'))
    parser.add_argument('--bin', type=int, default=1, choices=(1, 4), dest='bin_factor')
    args = parser.parse_args()
    files = sorted(p for p in args.root.glob('*fly1*/*_s_.mat')
                   if any(s in p.parent.name for s in ('heatbarrier', 'stripes')))
    assert len(files) == 3, f'Expected three selected heat recordings, found {len(files)}'
    for source in files:
        folder = 'pixel_fluorescence_maps' if args.bin_factor == 1 else 'pixel_fluorescence_maps_binned4x4'
        compute(source, source.parent / 'smr_analysis' / folder, bin_factor=args.bin_factor)


if __name__ == '__main__':
    main()
