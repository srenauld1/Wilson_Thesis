from pathlib import Path
import csv,hashlib,shutil,json,collections
root=Path('/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis')
rows=list(csv.DictReader(open(root/'cohort.csv')));locations=[];counts=collections.Counter()
for r in rows:
 text=' '.join(r.get(k,'') for k in ['trial','folder','context']).lower()
 geometry=next((g for token,g in [('pyramid','pyramid'),('stripe','stripes'),('gradient','gradient'),('multiple','multiple_barriers')] if token in text),'barrier')
 src=root/r['trial'];dst=root/geometry/r['trial'];assert src.is_dir() and not dst.exists(),str(src)
 before={str(p.relative_to(src)):hashlib.sha256(p.read_bytes()).hexdigest() for p in src.rglob('*') if p.is_file()}
 dst.parent.mkdir(exist_ok=True);src.rename(dst)
 after={str(p.relative_to(dst)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dst.rglob('*') if p.is_file()}
 assert before==after,r['trial']
 locations.append(dict(trial=r['trial'],geometry=geometry,analysisFolder=str(dst)));counts[geometry]+=1
with open(root/'trial_locations.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=['trial','geometry','analysisFolder']);w.writeheader();w.writerows(locations)
shutil.copyfile('Imaging_analysis/plot_pizza_ves041_heat_trials.m',root/'plot_pizza_ves041_heat_trials.m')
p=root/'README.txt';s=p.read_text();s+='\nFolder organization: barrier/, multiple_barriers/, gradient/, stripes/, and pyramid/. Individual trial folders retain both heat-on and heat-off figures. Geometry follows recording names and existing context labels; legacy heat recordings retain their existing barrier classification. trial_locations.csv lists the new paths. April 24 trial 2 remains awaiting ROI extraction under barrier/.\n';p.write_text(s)
p=root/'analysis_settings.json';s=json.loads(p.read_text());s['folderLayout']='geometry/trial; geometry folders: barrier, multiple_barriers, gradient, stripes, pyramid';p.write_text(json.dumps(s,indent=2))
assert len(list(root.rglob('*.png')))==46 and len(list(root.rglob('*.fig')))==46
print(dict(counts));print('Verified all moved file hashes and all 46 PNG/46 FIG files.')
