from pathlib import Path
import csv,json
root=Path('/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis')
rows=list(csv.DictReader(open(root/'manifest.csv')))
assert len(rows)==24
assert sum(r['status']=='awaiting ROI extraction' for r in rows)==1
rows=[r for r in rows if r['status']=='complete']
assert len(rows)==23
report=[]
for r in rows:
 assert r['status']=='complete',r
 p=root/r['trial'];e=list(csv.DictReader(open(p/'events.csv')))
 assert sum(x['heatOn']=='1' for x in e)==int(r['heatOnGroupEntries'])
 assert sum(x['heatOn']=='0' for x in e)==int(r['heatOffGroupEntries'])
 assert all((float(x['heatAtEvent'])>float(r['heatThreshold']))==(x['heatOn']=='1') for x in e)
 for cond in ['heat_on','heat_off']:
  for ext in ['png','fig']:assert (p/f"{r['trial']}_{cond}_peaks_slowdowns.{ext}").stat().st_size>1000
 assert (p/'analysis.mat').stat().st_size>1000
 report.append({'trial':r['trial'],'groupEntries':len(e),'transitionGroupEntries':sum(x['windowCrossesHeatState']=='1' for x in e)})
assert len(list(root.rglob('*.png')))==46
assert len(list(root.rglob('*.fig')))==46
Path('outputs/pizza_heat_20261028/verification.json').write_text(json.dumps(report,indent=2))
print('Verified 23 completed trials; 46 PNGs; 46 FIGs; event-state assignments and table totals match.')
