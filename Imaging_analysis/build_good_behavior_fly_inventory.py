"""Prepare source-backed fly/condition records for the pizza-talk workbook."""
import csv
import json
import re
from pathlib import Path

ROOT = Path('/Users/sophiarenauld/stacks')
REPO = Path(__file__).resolve().parents[1]
DEST = REPO / 'outputs/pizza_talk_20261028_fly_inventory'
DEST.mkdir(parents=True, exist_ok=True)


def read(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def identity(trial, folder, cell_type):
    date = re.match(r'(\d{8})', trial).group(1)
    if cell_type == 'VES041' and date in ('20261019', '20261024'):
        num, basis = '1', 'User confirmed: all recordings on this date are one fly'
    elif cell_type == 'VES041' and date == '20261020':
        number = int(trial.split('_')[1])
        num = '2' if 3 <= number <= 6 else '1'
        basis = 'User confirmed: trials 3–6 fly 2; trials 7–13 fly 1'
    else:
        match = re.search(r'fly(\d+)', folder, re.I)
        if match:
            num, basis = match.group(1), 'Explicit fly number in recording folder'
        else:
            num = '1'
            basis = ('Inferred first fly; later same-date folders explicitly label fly 2'
                     if date in ('20260213', '20260217') else
                     'Inferred single unnumbered preparation on this date; confirm identity')
    return f'{cell_type}_{date}_F{num}', date, basis


def condition(trial, folder, context=''):
    text = f'{trial} {folder} {context}'.lower()
    if 'odor' in text:
        return 'Odor'
    if 'dark' in text:
        return 'Dark'
    if 'pyramid' in text:
        return 'Heat pyramid'
    if 'stripe' in text:
        return 'Heat stripe'
    if 'olheat' in text or 'heatol' in text or 'heat-ol' in text:
        return 'Heat OL'
    if 'heat' in text or 'barrier' in text:
        return 'Heat barrier'
    if 'grating' in text:
        return 'Gratings'
    if 'jump' in text or 'jummp' in text:
        return 'Bar jumps'
    if 'bar' in text:
        return 'CL bar'
    raise ValueError(f'Unknown condition: {trial}')


records = []
ves_source = ROOT / 'imaging_analysis_artifacts/reports/ves041_trial_inventory_20261028.csv'
for r in read(ves_source):
    fly, date, basis = identity(r['trial'], r['folder'], 'VES041')
    core = r['coreGoodBehavior'] == 'true'
    records.append(dict(trial=r['trial'], cellType='VES041', date=date, flyId=fly,
                        good='Y' if core else 'N',
                        condition=condition(r['trial'], r['folder'], r['context']),
                        review=r['status'], identityBasis=basis,
                        source=str(ves_source), folder=r['folder']))

smp_source = REPO / 'Imaging_analysis/compare_smp544_ves041_reorientation.m'
text = smp_source.read_text().split('smpTrials=[', 1)[1].split('];', 1)[0]
smp_trials = re.findall(r'"([^"\n]+smp544[^"\n]+)"', text)
assert len(smp_trials) == 6
for trial in smp_trials:
    folder = str(ROOT / trial)
    assert Path(folder).is_dir()
    fly, date, basis = identity(trial, folder, 'SMP544')
    records.append(dict(trial=trial, cellType='SMP544', date=date, flyId=fly,
                        good='Y', condition=condition(trial, folder),
                        review='Previously behavior-screened: sustained, varied walking',
                        identityBasis=basis, source=str(smp_source), folder=folder))

# Preserve the recent SMP544 records as unreviewed, not as failed behavior.
smp_recent = ROOT / 'smp544_20261026_analysis_20261028/run_report.csv'
for r in read(smp_recent):
    folder = str(Path(r['output']).parent)
    trial = Path(folder).name
    explicit = re.search(r'fly(\d+)', trial, re.I)
    fly = f'SMP544_20261026_F{explicit.group(1)}' if explicit else ''
    records.append(dict(trial=trial, cellType='SMP544', date='20261026', flyId=fly,
                        good='Unreviewed', condition=condition(trial, folder, r['context']),
                        review='Processed only; no recorded good-behavior decision; registration QC held',
                        identityBasis='Explicit fly number' if explicit else 'Unassigned; fly identity absent',
                        source=str(smp_recent), folder=folder))

records.sort(key=lambda x: (x['cellType'], x['date'], x['trial']))
accepted = [r for r in records if r['good'] == 'Y']
assert len(accepted) == 40
flies = []
# Retain explicitly rejected flies so their good-behavior flag is visibly N.
reviewed = [r for r in records if r['good'] != 'Unreviewed']
for fly in sorted({r['flyId'] for r in reviewed}):
    rows = [r for r in reviewed if r['flyId'] == fly]
    assert len({r['date'] for r in rows}) == 1
    bases = sorted({r['identityBasis'] for r in rows})
    flies.append(dict(flyId=fly, date=rows[0]['date'],
                      identity='Inferred' if any('Inferred' in b for b in bases) else
                      'User confirmed' if any('User confirmed' in b for b in bases) else 'Explicit',
                      identityBasis='; '.join(bases)))

payload = dict(records=records, flies=flies,
               conditionHeaders=['Dark', 'CL bar', 'Bar jumps', 'Gratings', 'Heat OL',
                                 'Heat barrier', 'Heat pyramid', 'Heat stripe'])
(DEST / 'inventory_source.json').write_text(json.dumps(payload, indent=2))
print(f'{len(flies)} flies; {len(accepted)} accepted trials; {len(records)} source records')
print('Accepted by cell type:', {c: sum(r['cellType'] == c for r in accepted) for c in ['VES041', 'SMP544']})
