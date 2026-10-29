import csv,pathlib
base=pathlib.Path('/Users/sophiarenauld/stacks/imaging_analysis_artifacts/reports')
out=base/'ves041_behavior_review_20261028'
read=lambda p:list(csv.DictReader(open(p,newline='')))
metrics=read(out/'behavior_metrics.csv');inv=read(base/'ves041_trial_inventory_20261028.csv')
qc={r['trial']:r for r in read(out/'source_qc.csv')}
notes={
'20261019_3_2':'Repeated walking and bidirectional turning throughout. Acquisition marked failed, but raw sync audit found complete DAQ and matching behavior; retain acquisition-status warning.',
'20261019_4_1':'Repeated forward bouts and turns across the recording.',
'20261019_5_1':'Intermittent but repeated locomotion across nearly all 30-second bins; 43-second longest quiet interval.',
'20261019_6_1':'Behavior is usable in the captured segment, but acquisition/DAQ is incomplete. Retain only verified covered samples; do not count as a complete stripe experiment.',
'20261019_7_1':'Repeated locomotion and turns in the shorter, approximately 180-second recording.',
'20261019_8_2':'Repeated walking and turns throughout; exclude the recorded brief tracking gap.',
'20261020_3_1':'Sustained locomotion and frequent turns. Acquisition marked failed, but raw sync audit found complete DAQ and matching behavior; retain acquisition-status warning.',
'20261020_4_2':'Sustained locomotion and turns in dark; some directional bias, with turns in both directions.',
'20261020_5_1':'Repeated walking and turns throughout; directional bias is present, not unidirectional rotation. Exclude the brief tracking gap.',
'20261020_6_1':'Strong sustained walking with frequent turns. Negative yaw predominates but positive turns are also present; retain direction-bias note.',
'20261020_7_1':'Repeated strong forward bouts and turns in both directions throughout.',
'20261020_8_1':'Intermittent locomotion with an approximately 85-second quiet interval and activity concentrated in bouts. Suitable for selected active epochs; not a core whole-recording good-behavior addition.',
'20261020_9_1':'Repeated forward bouts and bidirectional turns across the recording.',
'20261020_11_2':'Repeated walking and turns before and after a roughly 57-second quiet interval. Retain, but note intermittency and exclude tracking gaps.',
'20261020_12_1':'Lower activity but repeated forward bouts and turns distributed across the recording. Retain as lower-activity behavior, with limited heat exposure; do not infer heat-response absence.',
'20261020_13_1':'Repeated strong forward bouts and bidirectional turns; exclude brief tracking gap.',
'20261024_2_3_ves041_barjumps':'Weakest walking trial: 10.2% above 1 mm/s, six walking bouts lasting at least one second, and only 44% active complete 30-second bins. Turning is present; defer inclusion in core locomotion cohort. This is not a finding of tracking failure or absence of visual responses.',
'20261024_2_4_ves041_barjumps':'Repeated forward bouts and turns throughout; clearly stronger locomotor sampling than bar-jump 2_3.',
'20261024_3_2_ves041_barrier_bar':'Repeated low-to-moderate forward locomotion and bidirectional turns; no prolonged quiet interval. Limited barrier coverage is a protocol-sampling issue.',
'20261024_4_1_ves041_barrier_bar_multiple':'Repeated forward bouts and bidirectional turns across the full ten minutes. Geometry conditions remain separate for interpretation.',
'20261024_5_2_ves041_pyramid_bar':'Repeated forward locomotion and bidirectional turns. Exclude the roughly 2.4-second tracking gap and masked boundary samples; limited pyramid coverage does not invalidate behavior.',
'20261024_6_3_ves041_stripes_ball':'Repeated forward bouts and turns across the full ten minutes; strong candidate for heat/no-heat behavioral comparisons.',
'20261024_7_1_ves041_stripes_bar':'Repeated forward bouts and bidirectional turns throughout. Sparse heat slowdowns limit that comparison, not overall behavioral quality.'}
restricted={'20261019_6_1':'restricted: verified partial recording','20261020_8_1':'restricted: active epochs only'}
hold={'20261024_2_3_ves041_barjumps'}
decisions=[]
for m in metrics:
 if m['recent']!='1':continue
 trial=m['trial'];dec=restricted.get(trial,'hold: low locomotor sampling' if trial in hold else 'add: good behavior')
 q=qc.get(trial) or qc.get(str(pathlib.Path(m['sourceFile']).parent.parent),{})
 row=dict(trial=trial,decision=dec,includeCoreGoodBehavior=str(dec.startswith('add')).lower(),reason=notes[trial])
 row.update({k:m[k] for k in ['durationSec','validPct','walkingPct','turningPct','walkingBoutsAtLeast1Sec','active30SecBinsPct','longestStillSec','longestInvalidSec','trackingMaskAvailable','sourceFile']})
 row.update(acquisitionStatus=q.get('acquisitionStatus',''),sourceQcNotes=q.get('notes',''))
 decisions.append(row)
by={r['trial']:r for r in decisions}
for r in inv:
 r['behaviorReviewDecision']='previously tagged; not reclassified'
 r['behaviorReviewReason']='Existing Finder good behavior tag retained.'
 r['coreGoodBehavior']='true'
 if r['trial'] in by:
  d=by[r['trial']];r['behaviorReviewDecision']=d['decision'];r['behaviorReviewReason']=d['reason'];r['coreGoodBehavior']=d['includeCoreGoodBehavior'];r['status']='reviewed: '+d['decision'];r['evidence']=str(out/'behavior_decisions.csv')
def write(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write(out/'behavior_decisions.csv',decisions);write(base/'ves041_trial_inventory_20261028.csv',inv);write(base/'ves041_good_behavior_trials_20261028.csv',[r for r in inv if r['coreGoodBehavior']=='true'])
lines=['VES041 behavior review — October 28, 2026','','Reviewed 23 processed recordings from October 19, 20, and 24. Recommendation: 20 core additions, two restricted recordings, and one held out. Combined with the 18 existing Finder-tagged trials, the core list contains 38 recordings. Finder tags were not edited.','','Review basis','','Decisions use full-recording forward velocity, yaw velocity, trajectories, available tracking-validity masks, and acquisition/synchronization logs. No calcium trace, peak amplitude, or desired heat/visual response was used for selection. No raw camera video was reviewed. These decisions concern behavioral suitability, not imaging quality or proof of correct ball tracking.','','Thirteen of the 18 previously tagged trials had saved fly_data.mat files available at the standard location and supplied a descriptive reference. Their walking fractions span 12.3–78.1%; this argues against a single high walking-fraction cutoff. The five missing standard analysis files do not alter their existing tag classification. The recent decisions are an exploratory manual review, not a validated automatic classifier.','','Metrics','','Forward and yaw signals are aligned to imaging volumes and smoothed with a 0.25-second moving median. Walking means forward velocity >=1 mm/s. Turning means absolute yaw >=30 degrees/s. Quiet means abs(forward)<0.5 mm/s and abs(yaw)<10 degrees/s. Fractions use finite samples that pass the saved behavior-validity mask where available. Valid percentage is coverage of the aligned behavior, not an independent camera-tracking accuracy score. An active 30-second bin has >=10% walking among its valid samples; trailing partial bins are omitted. Walking bouts require >=1 second continuously over threshold. XY uses saved protocol coordinates when available, otherwise forward/yaw-derived coordinates, so old and new path frames should not be compared as a tracking-accuracy test.','','Restrictions','','October 19 6_1 has usable behavior in an approximately 166-second segment but incomplete acquisition/DAQ. Keep it outside the complete-recording core list; use only verified behavior-covered samples. October 20 8_1 has repeated locomotion but an approximately 85-second quiet interval and extended sparse periods; retain as an active-epoch candidate. October 24 2_3 has only 10.2% walking and six sustained bouts and is held out of the core locomotion cohort; turns are still present, so it is not a universally unusable trial.','','October 19 3_2 and October 20 3_1 retain acquisition-failed warnings. Their raw synchronization audit showed complete DAQ and matching behavior timestamps, so the behavioral review supports inclusion without clearing the acquisition-status flag. Brief tracking gaps in other accepted trials must remain masked. After a gap, missing displacement cannot be recovered from the saved coordinate assumption. Geometry-dependent heat comparisons need their own coverage checks. Low heat occupancy, not reaching a barrier, and lack of a desired stimulus response do not by themselves cause behavioral exclusion.','','Per-recording decisions','']
for r in decisions:
 lines.append(f"{r['trial']}: {r['decision']}. Walking {float(r['walkingPct']):.1f}%; turning {float(r['turningPct']):.1f}%; {r['walkingBoutsAtLeast1Sec']} sustained walking bouts. {r['reason']}")
(out/'review_notes.txt').write_text('\n'.join(lines)+'\n')
assert len(decisions)==23 and sum(d['includeCoreGoodBehavior']=='true' for d in decisions)==20
assert len([r for r in inv if r['coreGoodBehavior']=='true'])==38
print('23 decisions; 20 added; 2 restricted; 1 held; 38 core trials total.')
