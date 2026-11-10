function result=crunch_pizza_four_barrier_original_turns()
% Requested four recordings, old +/-0.5-second >=30-degree turn selection.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier';
ids=["20260427-2_ves041_fly2_heat","20260427-7_ves041_fly3_heat_20mm","20261019_4_1","20261020_11_2"];
flies=["20260427_fly2","20260427_fly3","20261019_fly1","20261020_fly1"];
figures=["20260427-2_ves041_fly2_heat_heat_off_peaks_slowdowns.png","20260427-7_ves041_fly3_heat_20mm_heat_off_peaks_slowdowns.png","20261019_4_1_heat_off_peaks_slowdowns.png","20261020_11_2_heat_off_peaks_slowdowns.png"];
out=fullfile(root,'pooled_four_requested_trials_oldcriteria_20261102');if ~isfolder(out),mkdir(out);end
rel=(-2:.05:3)';conditions=["control","heat"];X=cell(4,2);counts=zeros(4,2);flyMeans=nan(numel(rel),4,3,2);allEvents=table;sourceRows=cell(0,7);
for k=1:4
 base=load(fullfile(root,ids(k),'analysis.mat'),'d','meta');
 if endsWith(string(base.meta.source),'fly_data.mat')
  a=load(base.meta.source,'fly_data');a=a.fly_data;q=a.dq(1);d=struct('t',double(q.t(:)),'dff',double(a.roi.ts(1,:))','fwd',double(q.bvf(:)));
  if isfield(q,'bvy_deg'),d.yaw=double(q.bvy_deg(:));else,d.yaw=rad2deg(double(q.bvy(:)));end
  d.heat=double(q.heat(:));d.valid=true(size(d.t));if isfield(q,'bergii')&&isfield(q.bergii,'behaviorVolumeValid'),d.valid=logical(q.bergii.behaviorVolumeValid(:));end
  P=struct('x',a.protocol.pos_x_mm(:),'control',a.protocol.control_trial(:),'trial',a.protocol.trial_number(:));
  contextSource=string(base.meta.source)+" protocol";
 else
  d=base.d;folder=fileparts(base.meta.source);v=dir(fullfile(folder,'*virmenData.mat'));assert(numel(v)==1);V=load(fullfile(v.folder,v.name),'vmData');V=V.vmData;tk=V.Seconds{1};if isduration(tk),tk=seconds(tk);end;tk=double(tk(:));
  type=col(V,'currTt');ctl=col(V,'control_trial');assert(all(abs(type-round(type))<1e-6|~isfinite(type)));
  trial=cumsum([true;diff(type)~=0]);barrier=ismember(type,[1 2]);ctl(~barrier)=nan;trial=double(trial);trial(~barrier)=nan;
  P=struct('x',interp1(tk,col(V,'pos_x'),d.t,'linear',nan),'control',interp1(tk,ctl,d.t,'previous',nan),'trial',interp1(tk,trial,d.t,'previous',nan));
  d.heat=interp1(tk,col(V,'heatDAQ'),d.t,'linear',nan);contextSource=string(fullfile(v.folder,v.name))+" currTt/control_trial/heatDAQ";
 end
 [E,arrays]=select(d,P,ids(k),flies(k),rel);allEvents=[allEvents;E];
 % Reproduce the last requested single-recording result exactly.
 if ids(k)=="20261020_11_2"
  old=load(fullfile(root,ids(k),'presentation_control_vs_heat_20261102','original_criteria_turns','comparison_data.mat'),'X','counts','selected');
  assert(isequal(cellfun(@(z)size(z,2),arrays),old.counts));
  for s=1:2
   take=E.selected&E.context==conditions(s);previous=old.selected.eventTime(old.selected.context==conditions(s));assert(max(abs(E.eventTime(take)-previous))<1e-6);
   delta=abs(arrays{s}-old.X{s});assert(all(delta(:)<1e-3),'Original turn arrays differ beyond CSV timestamp rounding');arrays{s}=old.X{s};
  end
 end
 for s=1:2,X{k,s}=arrays{s};counts(k,s)=size(arrays{s},2);if counts(k,s)>0,flyMeans(:,k,:,s)=mean(arrays{s},2);end,end
 writetable(E,fullfile(out,ids(k)+"_events.csv"));
 sourceRows(end+1,:)={ids(k),flies(k),figures(k),string(base.meta.source),contextSource,counts(k,1),counts(k,2)};
 fprintf('%s: %d control turns; %d heat turns\n',ids(k),counts(k,1),counts(k,2));
end
assert(numel(unique(flies))==4);meanTrace=nan(numel(rel),3,2);sem=meanTrace;contributing=sum(counts>0,1);
for s=1:2
 for row=1:3
  q=flyMeans(:,:,row,s);n=sum(isfinite(q),2);meanTrace(:,row,s)=mean(q,2,'omitnan');sem(:,row,s)=std(q,0,2,'omitnan')./sqrt(n);sem(n<2,row,s)=nan;
  % Verify equal-fly, not event-weighted, averages and between-fly SEM.
  expected=zeros(numel(rel),1);for k=find(counts(:,s)>0)',expected=expected+squeeze(mean(X{k,s}(:,:,row),2))/contributing(s);end
  assert(max(abs(expected-meanTrace(:,row,s)))<1e-10);
 end
end
fig=figure('Visible','off','Color','w','Position',[20 20 1120 1050]);tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};light={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};titles={'Control turns','Heat turns'};
for row=1:3
 z=reshape(flyMeans(:,:,row,:),[],1);band=[reshape(meanTrace(:,row,:)-sem(:,row,:),[],1);reshape(meanTrace(:,row,:)+sem(:,row,:),[],1)];z=[z;band];z=z(isfinite(z));pad=max(.01,.06*range(z));lim=[min(z)-pad max(z)+pad];
 for s=1:2
  ax=nexttile(tl);hold(ax,'on');
  for k=find(counts(:,s)>0)',plot(ax,rel,flyMeans(:,k,row,s),'Color',light{row},'LineWidth',1.1,'Tag','flyMean','UserData',flies(k));end
  mu=meanTrace(:,row,s);se=sem(:,row,s);fill(ax,[rel;flipud(rel)],[mu-se;flipud(mu+se)],colors{row},'FaceAlpha',.2,'EdgeColor','none','Tag','flySem');plot(ax,rel,mu,'Color',colors{row},'LineWidth',3,'Tag','equalFlyMean');
  xline(ax,0,':','Color',[.45 .45 .45],'LineWidth',1);xlim(ax,[-2 3]);ylim(ax,lim);ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=-2:1:3;
  if s==1,ylabel(ax,labels{row},'FontSize',21);end
  if row==3,xlabel(ax,'Time from strongest deceleration (s)','FontSize',17);else,ax.XTickLabel=[];end
  if row==1,title(ax,titles{s},'FontSize',23,'FontWeight','bold');subtitle(ax,sprintf('%d flies | %d events',contributing(s),sum(counts(:,s))),'FontSize',15);end
 end
end
subtitle(tl,{'4 requested recordings | Original turn criteria','Thin: fly means | Thick: equal-fly mean | Shading: SEM across flies'},'FontSize',15);
assert(numel(findall(fig,'Type','axes'))==6&&numel(findall(fig,'Tag','flyMean'))==3*sum(contributing));stem=fullfile(out,'ves041_four_trials_original_criteria_control_heat_turns');exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);
F=cell2table(sourceRows,'VariableNames',{'recording','flyId','requestedFigure','traceSource','contextSource','controlTurns','heatTurns'});writetable(F,fullfile(out,'cohort_and_counts.csv'));writetable(allEvents,fullfile(out,'all_slowdown_events.csv'));
methods=struct('selection','Deceleration<=-2 mm/s^2 after .25s median smoothing; forward-before>=1 and speed-drop>=1 mm/s; >=2 s accepted slowdown spacing; accumulated absolute smoothed yaw>=30 deg in +/-0.5 s', ...
 'speedWindows','Median forward -0.75:-0.25 versus +0.25:+0.75 s','eventWindow','-2:+3 s finite behavior, fluorescence, heat monitor and logged position; no time gaps', ...
 'control','Core +/-0.5 s entirely control_trial=1 and one barrier trial; no extra full-window same-control or placed-barrier requirement', ...
 'heat','Core +/-0.5 s entirely measured command>0.05 V, one barrier trial and not control; no-heat and transition events excluded', ...
 'legacyProtocol','Legacy currTt1=heat, currTt2=control; derive consecutive trial IDs from type changes; accommodation currTt3 explicitly excluded, including when saved control flag is stale', ...
 'averaging','Pool events within each fly and condition, then equal weight across contributing fly means; one recording per fly here', ...
 'uncertainty','Pointwise sample SD of fly means /sqrt(number of contributing flies); missing conditions remain NaN, not zero', ...
 'normalization','Original source dF/F and physical kinematic units; no added normalization, baseline subtraction, z-score, smoothing of display, or amplitude selection', ...
 'cohort','Only the four requested recordings; heat-off PNG names identify recordings, not the event subset used', ...
 'oct20Validation','All selected event arrays identical to previous original-criteria 15-control/18-heat comparison');
save(fullfile(out,'crunch_data.mat'),'X','counts','flyMeans','meanTrace','sem','contributing','rel','ids','flies','F','allEvents','methods','-v7.3');fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('folder',out,'recordings',4,'flies',4,'controlEvents',sum(counts(:,1)),'heatEvents',sum(counts(:,2)),'contributingFlies',contributing);disp(result);
end
function x=col(V,key)
x=V.(key);if iscell(x),x=x{1};end;x=double(x(:));
end
function [E,arrays]=select(d,P,id,fly,rel)
t=d.t(:);fs=1/median(diff(t));n=max(3,round(.25*fs));f=smoothdata(d.fwd,'movmedian',n);r=smoothdata(d.yaw,'movmedian',n);acc=gradient(f,t);h=d.heat(:);
valid=isfinite(f)&isfinite(r)&isfinite(h)&isfinite(d.dff)&isfinite(P.x)&d.valid;
candidate=valid&acc<=-2;beg=find(diff([false;candidate])==1);last=find(diff([candidate;false])==-1);rows=cell(0,11);lastAccepted=-Inf;
for k=1:numel(beg)
 [~,j]=min(acc(beg(k):last(k)));q=beg(k)+j-1;et=t(q);if et-lastAccepted<2||et<t(1)+2||et>t(end)-3,continue,end
 win=t>=et-2&t<=et+3;pre=t>=et-.75&t<et-.25;post=t>=et+.25&t<=et+.75;core=t>=et-.5&t<=et+.5;
 if ~all(valid(win))||any(diff(t(win))>2/fs),continue,end
 before=median(f(pre));after=median(f(post));drop=before-after;if before<1||drop<1,continue,end
 angle=trapz(t(core),abs(r(core)));on=h(core)>.05;
 if all(on),ctx="heat";elseif all(~on),ctx="no heat";else,ctx="transition";end
 if any(~isfinite(P.control(core)))||any(~isfinite(P.trial(core)))||numel(unique(P.trial(core)))~=1,ctx="transition";elseif all(P.control(core)==1),ctx="control";end
 rows(end+1,:)={id,fly,q,et,ctx,angle>=30,before,after,drop,angle,acc(q)};lastAccepted=et;
end
E=cell2table(rows,'VariableNames',{'recording','flyId','sampleIndex','eventTime','context','reorientation','speedBefore','speedAfter','speedDrop','absoluteTurnDeg','accelerationMmS2'});
E.context=string(E.context);E.selected=E.reorientation==1&ismember(E.context,["control","heat"]);arrays=cell(1,2);vars={d.dff,d.fwd,abs(d.yaw)};
for s=1:2
 context="control";if s==2,context="heat";end;events=E.eventTime(E.selected&E.context==context);arrays{s}=nan(numel(rel),numel(events),3);
 for j=1:numel(events),for v=1:3,arrays{s}(:,j,v)=interp1(t,vars{v},events(j)+rel,'linear',NaN);end,end
 assert(all(isfinite(arrays{s}(:))));
end
end
