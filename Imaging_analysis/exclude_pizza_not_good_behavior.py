"""Apply the user's four recording exclusions, preserving removed outputs."""
import csv
import shutil
from pathlib import Path

root = Path('/Users/sophiarenauld/stacks/pizza_talk_analysis')
base = root / 'ves041_noheat_imaging'
archive = root / 'excluded_not_good_behavior'
archive.mkdir(exist_ok=True)


def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f))


def write(p, rows):
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


rows = read(base / 'manifest.csv')
removed = []
for r in rows:
    trial = r['trial']
    if not (trial.startswith(('20251001', '20260213')) or trial == '20261020_3_1'):
        continue
    dest = archive / trial
    dest.mkdir(exist_ok=True)
    for stem_key in ('eventFigureStem', 'pathFigureStem'):
        old = Path(r[stem_key])
        new = dest / old.name
        for ext in ('.png', '.fig'):
            src = Path(str(old) + ext)
            target = Path(str(new) + ext)
            if src.exists() and src != target:
                assert not target.exists()
                shutil.move(str(src), str(target))
        r[stem_key] = str(new)
    for ext in ('.csv', '.mat'):
        src = base / 'analysis_data' / f'{trial}_events{ext}'
        if src.exists():
            shutil.move(str(src), str(dest / src.name))
    r['status'] = 'excluded: user marked not good behavior'
    r['category'] = 'excluded'
    removed.append(trial)
assert len(removed) == 4
write(base / 'manifest.csv', rows)

mapping = [r for r in read(base / 'bar_fly_mapping.csv') if r['trial'] not in removed]
for r in mapping:
    for date in ('20261019', '20261024'):
        if r['trial'].startswith(date):
            r['flyId'] = f'{date}_fly1'
            r['identityBasis'] = f'User confirmed: all {date} VES041 recordings are one fly'
assert len(mapping) == 9 and len({r['flyId'] for r in mapping}) == 7
write(base / 'bar_fly_mapping.csv', mapping)

prior = archive / 'previous_pooled_with_excluded_trials'
prior.mkdir(exist_ok=True)
for p in (base / 'bar').glob('ves041_bar_pooled_equal_fly_peaks_slowdowns*'):
    assert not (prior / p.name).exists()
    shutil.move(str(p), str(prior / p.name))
(archive / 'exclusions.txt').write_text(
    'User marked these recordings not good behavior on 2026-10-28.\n'
    + '\n'.join(removed)
    + '\nIndividual outputs and superseded pooled outputs are preserved here, '
      'outside the included-analysis folders.\n')
print('Excluded 4 recordings; updated bar cohort: 9 recordings, 7 flies.')
