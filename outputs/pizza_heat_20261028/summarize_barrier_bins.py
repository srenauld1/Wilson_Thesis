from pathlib import Path
import csv,json,collections
root=Path('/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier')
rows=list(csv.DictReader(open(root/'control_heat_bins_manifest.csv')))
assert len(rows)==14 and sum(r['status']=='complete' for r in rows)==13
assert sum(r['status']=='awaiting ROI extraction' for r in rows)==1
summary=[];pngs=figs=pathpanels=0
for r in rows:
 if r['status']!='complete':continue
 p=Path(r['output']);E=list(csv.DictReader(open(p/'events.csv')))
 assert len(E)==int(r['includedGroupEntries'])+int(r['excludedGroupEntries'])
 for e in E:
  b=int(e['heatBin']);h=float(e['medianHeatV']);control=float(e['controlTrial'])
  if b==1:assert control==1
  if b>1:assert control==0
  if b==2:assert -.05<=h<1
  if b==3:assert 1<=h<2.5
  if b==4:assert h>=2.5
 for group in range(1,9):
  for b,name in enumerate(['control_barrier','heat_0_to_lt1V','heat_1_to_lt2p5V','heat_ge2p5V'],1):
   selected=[e for e in E if int(e['eventGroup'])==group and int(e['heatBin'])==b]
   summary.append(dict(recording=r['trial'],eventGroup=group,condition=name,events=len(selected)))
 T=list(csv.DictReader(open(p/'flat_path_trials.csv')));assert len(T)==int(r['flatPathTrials']);assert len({x['barrierTrial'] for x in T})==len(T);pathpanels+=len(T)
 n=(len(T)+11)//12;assert len(list(p.glob('flat_paths_dff_page_*.png')))==n
 assert len(list(p.glob('*.png')))==4+n and len(list(p.glob('*.fig')))==4+n
 pngs+=4+n;figs+=4+n
with (root/'control_heat_bin_counts.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=summary[0].keys());w.writeheader();w.writerows(summary)
report=dict(completeRecordings=13,awaitingROI=1,conditionFigures=52,flatPathPanels=pathpanels,flatPathPages=pngs-52,pngs=pngs,matlabFigures=figs)
(root/'control_heat_bins_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
