function apply_inhibition_rules_to_noheat_imaging()
% No heat cohort. Strict localized-turn selection uses recorded heading/XY.
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';
out=fullfile(pwd,'outputs','ves041_imaging_inhibition_rules_noheat_20261104');
if ~isfolder(out),mkdir(out);end
F=readtable(fullfile(root,'dark_bar_barjump_fly_mapping.csv'),'Delimiter',',','TextType','string');
X={};oldX={};usedFly=strings(0,1);metadata=struct([]);eventRows=struct([]);status=struct([]);rel=(-2:.05:3)';
for k=1:height(F)
 name=char(F.trial(k));A=load(fullfile(root,'analysis_data',[name '_events.mat']),'d','meta','X','groups','rel');
 row=struct('trial',name,'flyId',char(F.flyId(k)),'status','','message','');
 if ~endsWith(A.meta.source,'fly_data.mat')||~isfile(A.meta.source)
  row.status='unassessable';row.message='Only archived plot traces available; no saved actual ball heading. No heading is inferred from yaw.';
  status=[status;row];continue %#ok<AGROW>
 end
 S=load(A.meta.source,'fly_data');q=S.fly_data.dq(1);d=A.d;n=numel(d.t);
 if ~all(isfield(q,{'bh','pxb','pyb'}))||any([numel(q.bh),numel(q.pxb),numel(q.pyb)]~=n)
  row.status='unassessable';row.message='Recorded heading or ball XY is missing or not aligned.';status=[status;row];continue %#ok<AGROW>
 end
 assert(max(abs(double(q.t(:))-d.t))<1e-8,'Heading/position times do not match saved fluorescence.');
 assert(max(abs(double(q.pxb(:))-d.x))<1e-8 && max(abs(double(q.pyb(:))-d.y))<1e-8,'Saved flat path differs from source ball XY.');
 valid=d.valid&isfinite(q.bh(:))&isfinite(d.x)&isfinite(d.y);
 B=blocks(d.t,valid);heading=nan(n,1);for j=1:size(B,1),ix=B(j,1):B(j,2);heading(ix)=rad2deg(unwrap(double(q.bh(ix))));end
 dt=median(diff(d.t));window=max(1,round(1/dt));stride=max(1,round(.1/dt));
 detected=struct('first_start',{},'last_end',{},'representative_start',{},'representative_end',{},'rotation_deg',{},'displacement_mm',{},'directionality',{});
 for j=1:size(B,1)
  lo=B(j,1);hi=B(j,2);h=heading(lo:hi);x=d.x(lo:hi);y=d.y(lo:hi);
  starts=(1:stride:numel(h)-window)';ends=starts+window;if isempty(starts),continue,end
  total=[0;cumsum(abs(diff(h)))];net=abs(h(ends)-h(starts));amount=total(ends)-total(starts);
  direction=zeros(size(net));take=amount>0;direction(take)=net(take)./amount(take);
  dispmm=hypot(x(ends)-x(starts),y(ends)-y(starts));passing=find(net>30&direction>=.7&dispmm<=2);
  for v=passing'
   a=starts(v)+lo-1;b=ends(v)+lo-1;
   if ~isempty(detected)&&a<=detected(end).last_end
    detected(end).last_end=b;
    if net(v)>detected(end).rotation_deg
     detected(end).representative_start=a;detected(end).representative_end=b;
     detected(end).rotation_deg=net(v);detected(end).displacement_mm=dispmm(v);detected(end).directionality=direction(v);
    end
   else
    detected(end+1)=struct('first_start',a,'last_end',b,'representative_start',a,'representative_end',b, ...
     'rotation_deg',net(v),'displacement_mm',dispmm(v),'directionality',direction(v)); %#ok<AGROW>
   end
  end
 end
 Y=nan(numel(rel),0,3);keep=false(numel(detected),1);centers=zeros(numel(detected),1);
 for j=1:numel(detected)
  e=detected(j);center=(d.t(e.representative_start)+d.t(e.representative_end))/2;centers(j)=center;
  keep(j)=complete(d.t,valid,center);
  r=e;r.trial=name;r.flyId=char(F.flyId(k));r.anchor_time=center;r.bout_start_time=d.t(e.first_start);r.bout_end_time=d.t(e.last_end);r.full_epoch_valid=keep(j);
  eventRows=[eventRows;r]; %#ok<AGROW>
  if keep(j),Y(:,end+1,:)=epoch(d,center,rel);end %#ok<AGROW>
 end
 old=A.X{4};oldkeep=arrayfun(@(i)complete(d.t,valid,d.t(i)),A.groups{4});old=old(:,oldkeep,:);
 assert(isequal(A.rel(:),rel)&&all(isfinite(Y(:)))&&all(isfinite(old(:))));
 % Physical consistency check: translation only; the historical totalspeed
 % field includes yaw and is not an XY calibration reference.
 pathSpeed=hypot(diff(d.x),diff(d.y))./diff(d.t);v=hypot(double(q.bvf(:)),double(q.bvs(:)));v=(v(1:end-1)+v(2:end))/2;
 ok=isfinite(pathSpeed)&isfinite(v)&v>1;calibration=median(pathSpeed(ok)./v(ok));
 m=struct('trial',name,'flyId',char(F.flyId(k)),'source',A.meta.source,'heading_source','dq(1).bh radians, unwrapped within valid blocks', ...
  'xy_source','dq(1).pxb/pyb, original ball flat-path coordinates','path_speed_ratio',calibration,'sample_rate_hz',1/dt, ...
  'window_samples',window,'window_duration_s',window*dt,'stride_samples',stride,'stride_s',stride*dt, ...
  'detected_bouts',numel(detected),'valid_epochs',nnz(keep),'old_imaging_events',size(old,2),'old_events_rejected_common_validity',nnz(~oldkeep));
 metadata=[metadata;m];usedFly(end+1,1)=F.flyId(k);X{end+1}=Y;oldX{end+1}=old; %#ok<AGROW>
 row.status='complete';row.message='Recorded heading and physical ball XY; source dF/F unchanged.';status=[status;row]; %#ok<AGROW>
 fprintf('%s: %d old-rule events, %d localized-turn events, path calibration ratio %.3f\n',name,size(old,2),size(Y,2),calibration);
end
flies=sort(unique(usedFly));nf=numel(flies);counts=zeros(nf,2);flyMeans=nan(numel(rel),nf,3,2);pooled=cell(nf,2);
for f=1:nf
 for k=find(usedFly==flies(f))'
  for definition=1:2
   if definition==1,a=oldX{k};else,a=X{k};end
   if isempty(pooled{f,definition}),pooled{f,definition}=a;else,pooled{f,definition}=cat(2,pooled{f,definition},a);end
  end
 end
 for definition=1:2
  a=pooled{f,definition};counts(f,definition)=size(a,2);
  if counts(f,definition)>0,flyMeans(:,f,:,definition)=mean(a,2);end
 end
end
meanTrace=nan(numel(rel),3,2);sem=meanTrace;contributing=sum(counts>0,1);
for definition=1:2
 for variable=1:3
  a=flyMeans(:,:,variable,definition);nvalid=sum(isfinite(a),2);
  meanTrace(:,variable,definition)=mean(a,2,'omitnan');sem(:,variable,definition)=std(a,0,2,'omitnan')./sqrt(nvalid);
  sem(nvalid<2,variable,definition)=nan;
 end
end
methods=struct('cohort','Only the established pooled dark/bar/bar-jump VES041 no-heat cohort. Heat cohort not loaded or analyzed.', ...
 'selection','>30 deg net unwrapped recorded heading change in ~1 s; <=2 mm net XY displacement; net/accumulated absolute heading change >=0.7. Sliding starts every ~0.1 s, rounded to available imaging-aligned frames. Overlapping/touching qualifying windows merged. No slowdown or 2-s refractory requirement.', ...
 'alignment','Center of the representative 1-second window with greatest net rotation within each merged bout. This is not a neural peak, slowdown time, or experimentally measured turn onset.', ...
 'comparator','Saved original imaging-rule slowdown+rotation events on exactly the same available recordings, restricted only to common complete heading/XY/fluorescence/behavior validity.', ...
 'comparator_alignment','Strongest deceleration. Latencies cannot be directly compared across definitions because their alignment anchors differ.', ...
 'fluorescence','Original stored dF/F; no peak/amplitude selection, added baseline subtraction, z-score, or extra smoothing. Linear display interpolation onto -2:0.05:+3 s, matching original figures.', ...
 'validity','Full -2:+3 s event epochs must be finite and behavior-valid without time gaps. No heat state, no opto state, no added visual/bar-jump exclusion.', ...
 'sampling','Use recorded heading and XY at the saved imaging-aligned behavioral rate, without inventing a 100-Hz stream. Lower sampling can alter directionality estimates compared with the inhibition data.', ...
 'pooling','Pool events over recordings within each mapped fly; equal weight per contributing fly. SEM is pointwise sample SD across fly means divided by sqrt(n). Zero-event flies are excluded for that definition, not replaced by zeros.', ...
 'missing_data','Two older figure-only recordings lack recoverable saved actual heading. They remain explicitly unassessable, not substituted by integration of yaw.', ...
 'identities','Existing fly mapping retained, including inferred identities. See original dark_bar_barjump_fly_mapping.csv and source metadata.', ...
 'statistics','Descriptive plots only. No new hypothesis tests.');
save(fullfile(out,'analysis.mat'),'X','oldX','pooled','flyMeans','meanTrace','sem','counts','flies','rel','metadata','status','eventRows','methods','contributing','-v7.3');
payload=struct('flies',flies,'counts',counts,'time',rel,'fly_means',flyMeans,'mean',meanTrace,'sem',sem,'metadata',metadata,'status',status,'events',eventRows,'methods',methods,'contributing',contributing);
fid=fopen(fullfile(out,'analysis.json'),'w');fprintf(fid,'%s',jsonencode(payload));fclose(fid);
fprintf('COMPLETE: %d recordings, %d estimated flies; old/new events %d/%d\n',numel(metadata),nf,sum(counts,1));
end

function z=epoch(d,t0,rel)
z=nan(numel(rel),1,3);vars={d.dff,d.fwd,abs(d.yaw)};
for k=1:3,z(:,1,k)=interp1(d.t,vars{k},t0+rel);end
end
function tf=complete(t,valid,t0)
fs=1/median(diff(t));q=find(t<=t0-2,1,'last');r=find(t>=t0+3,1,'first');
if isempty(q)||isempty(r),tf=false;else,tf=all(valid(q:r))&&all(diff(t(q:r))<=2/fs);end
end
function B=blocks(t,valid)
cut=[true;diff(t)>2*median(diff(t))];s=find(valid&(cut|[true;~valid(1:end-1)]));e=find(valid&([~valid(2:end);true]|[cut(2:end);true]));B=[s e];
end
