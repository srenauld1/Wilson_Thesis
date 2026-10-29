function manifest=plot_pizza_ves041_noheat_imaging(out,onlyTrials)
% Per-recording 4x4 peak/slowdown comparisons and whole-recording flat paths.
% Selection is the reviewed core list; entire heat experiments are excluded.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
base=fileparts(mfilename('fullpath'));
input=fullfile(imaging_analysis_artifact_root(),'reports','ves041_good_behavior_trials_20261028.csv');
I=readtable(input,'Delimiter',',','TextType','string');
prior=table();
if nargin>=2&&~isempty(onlyTrials)
 prior=readtable(fullfile(out,'manifest.csv'),'Delimiter',',','TextType','string');
 I=I(ismember(I.trial,string(onlyTrials)),:);
end
if ~isfolder(out),mkdir(out);end
dataOut=fullfile(out,'analysis_data');if ~isfolder(dataOut),mkdir(dataOut);end
settings=struct('inventory',input,'windowSec',[-2 3],'displayDff','original source dF/F; no event baseline subtraction', ...
 'peakDetector','Standard plotting_ves041.m findpeaks on original dF/F', ...
 'noiseEstimator','1.4826 * MAD(dF/F - 30-second moving-median baseline)', ...
 'peakProminenceSigma',3,'peakDistanceSec',2,'peakWidthSec',.25, ...
 'slowdownSmoothSec',.25,'maximumAccelerationMmS2',-2, ...
 'minimumSpeedBeforeMmS',1,'minimumSpeedDropMmS',1,'slowdownSeparationSec',2, ...
 'speedBeforeWindowSec',[-.75 -.25],'speedAfterWindowSec',[.25 .75], ...
 'turnWindowSec',[0 .5],'turnAbsoluteYawDeg',15, ...
 'rotationDisplay','absolute yaw speed in deg/s; signed yaw retained in MAT', ...
 'eventOrder','chronological within each column','eventGridSec',.05, ...
 'heatExclusion','entire configured heat recordings and recordings with measured heat >0.05 V excluded', ...
 'stimulusEpochs','visual and odor epochs retained','colorLimits','full finite dF/F range; shared between event columns', ...
 'eventWindowValidity','full window must have finite dF/F, forward, yaw and saved valid behavior, with no time gaps');
settings.flatPath=struct('rendering','scatter','markerSize',14,'colorPercentiles',[5 95],'drawOrder','ascending original dF/F; higher values last','samplesRemovedForColor',false);
settings.folderLayout='bar and gratings in subfolders; dark and odor in parent';
fid=fopen(fullfile(out,'analysis_settings.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(settings,PrettyPrint=true));fclose(fid);
rows=cell(0,13);
for k=1:height(I)
 name=char(I.trial(k));folder=char(I.folder(k));analysis=I.analysis_folder(k);
 if ismissing(analysis)||strlength(analysis)==0,analysis=fullfile(folder,'smr_analysis');end
 descriptor=lower(string(name)+" "+I.context(k)+" "+folder);if ismissing(descriptor),descriptor=lower(string(name)+" "+folder);end
 if contains(descriptor,'heat')||contains(descriptor,'barrier')||contains(descriptor,'pyramid')||contains(descriptor,'stripes')
  rows(end+1,:)={name,'excluded: heat experiment',folder,'',0,0,0,0,0,0,'','',''};continue
 end
 fprintf('PIZZA %s\n',name);
 try
  [d,meta]=load_trial(folder,char(analysis));
  if meta.isHeat
   rows(end+1,:)={name,'excluded: heat experiment',folder,meta.source,0,0,0,0,0,0,meta.pathSource,'',''};continue
  end
  validate(d);
  [peakIdx,detector]=peaks(d);
  usablePeakIdx=peakIdx(arrayfun(@(i)complete_window(d,i),peakIdx));
  [slowIdx,turnDeg]=slowdowns(d);withTurn=turnDeg>=15;
  groups={usablePeakIdx,slowIdx,slowIdx(~withTurn),slowIdx(withTurn)};
  assert(numel(groups{2})==numel(groups{3})+numel(groups{4}));
  assert(isequal(sort([groups{3};groups{4}]),groups{2}),'Slowdown partition mismatch');
  rel=(-2:.05:3)';X=cell(1,4);
  for c=1:4
   X{c}=nan(numel(rel),numel(groups{c}),3);
   vars={d.dff,d.fwd,abs(d.yaw)};
   for j=1:numel(groups{c})
    tq=d.t(groups{c}(j))+rel;
    for v=1:3,X{c}(:,j,v)=interp1(d.t,vars{v},tq,'linear',NaN);end
   end
   assert(all(isfinite(X{c}(:))),'Event interpolation produced nonfinite samples');
  end
  category=pizza_ves041_trial_category(name,folder);figureOut=string(out);
  if any(category==["bar","gratings"]),figureOut=fullfile(out,category);if ~isfolder(figureOut),mkdir(figureOut);end,end
  stem=fullfile(figureOut,[name '_peaks_slowdowns']);event_figure(X,rel,name,stem);
  pathStem=fullfile(figureOut,[name '_flat_path_dff']);path_figure(d,name,pathStem,meta.pathSource);
  eventType=[repmat("peak",numel(peakIdx),1);repmat("slowdown",numel(slowIdx),1)];
  sampleIndex=[peakIdx;slowIdx];eventTime=d.t(sampleIndex);
  included=[ismember(peakIdx,usablePeakIdx);true(numel(slowIdx),1)];
  absoluteTurnDeg=[nan(numel(peakIdx),1);turnDeg];hasTurn=[nan(numel(peakIdx),1);double(withTurn)];
  E=table(eventType,sampleIndex,eventTime,included,absoluteTurnDeg,hasTurn);
  writetable(E,fullfile(dataOut,[name '_events.csv']));
  save(fullfile(dataOut,[name '_events.mat']),'d','meta','detector','groups','rel','X','E','settings','-v7.3');
  rows(end+1,:)={name,'complete',folder,meta.source,numel(peakIdx),numel(usablePeakIdx),numel(slowIdx),nnz(~withTurn),nnz(withTurn),100*mean(d.valid),meta.pathSource,stem,pathStem};
  fprintf('  %d/%d peaks, %d slowdowns = %d no turn + %d turn\n',numel(usablePeakIdx),numel(peakIdx),numel(slowIdx),nnz(~withTurn),nnz(withTurn));
 catch ME
  rows(end+1,:)={name,['FAILED: ' ME.message],folder,'',0,0,0,0,0,0,'','',''};
  warning('pizza:TrialFailed','%s: %s',name,ME.message);
 end
 manifest=cell2table(rows,'VariableNames',{'trial','status','trialFolder','source','detectedPeaks','plottedPeaks','slowdowns','slowdownsWithoutTurn','slowdownsWithTurn','validBehaviorPct','pathSource','eventFigureStem','pathFigureStem'});
if ~isempty(prior),manifest=[prior(~ismember(prior.trial,string(manifest.trial)),manifest.Properties.VariableNames);convertvars(manifest,@iscell,'string')];end
 writetable(manifest,fullfile(out,'manifest.csv'));
end
manifest=cell2table(rows,'VariableNames',{'trial','status','trialFolder','source','detectedPeaks','plottedPeaks','slowdowns','slowdownsWithoutTurn','slowdownsWithTurn','validBehaviorPct','pathSource','eventFigureStem','pathFigureStem'});
if ~isempty(prior),manifest=[prior(~ismember(prior.trial,string(manifest.trial)),manifest.Properties.VariableNames);convertvars(manifest,@iscell,'string')];end
writetable(manifest,fullfile(out,'manifest.csv'));
assert(~any(startsWith(string(manifest.status),'FAILED')),'Some requested trials failed; see manifest.csv');
end

function [d,m]=load_trial(folder,analysis)
f=fullfile(analysis,'fly_data.mat');m=struct('source',f,'isHeat',false,'pathSource','');
if ~isfile(f)
 [d,m]=from_saved_figures(analysis);return
end
S=load(f,'fly_data');a=S.fly_data;q=a.dq(1);
d.t=double(q.t(:));d.dff=double(a.roi(1).ts(1,:))';d.fwd=double(q.bvf(:));
odor=contains(lower(folder),'odor');
if odor
 assert(isfield(q,'vvy'),'Missing odor-specific yaw channel');d.yaw=rad2deg(double(q.vvy(:)));
else
 if isfield(q,'bvy_deg'),d.yaw=double(q.bvy_deg(:));else,d.yaw=rad2deg(double(q.bvy(:)));end
end
d.valid=isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
if isfield(q,'bergii')&&isfield(q.bergii,'behaviorVolumeValid'),d.valid=d.valid&logical(q.bergii.behaviorVolumeValid(:));end
if isfield(a,'protocol')&&isfield(a.protocol,'kind'),m.isHeat=startsWith(a.protocol.kind,'heat');end
if isfield(q,'heat')&&~isempty(q.heat),m.isHeat=m.isHeat||any(q.heat(:)>.05);end
if odor,fields={'px','py'};else,fields={'pxb','pyb'};end
if all(cellfun(@(key)isfield(q,key)&&numel(q.(key))==numel(d.t),fields))
 d.x=double(q.(fields{1})(:));d.y=double(q.(fields{2})(:));m.pathSource=['dq(1).' fields{1} '/' fields{2}];
else
 [d.x,d.y]=path_from_figure(analysis,numel(d.t));m.pathSource=fullfile(analysis,'x_color_linear_dff_vs_y_color_linear_dff.fig');
end
if odor,m.yawSource='dq(1).vvy radians/s converted to deg/s (odor channel correction)';else,m.yawSource='dq(1).bvy_deg or converted bvy';end
end

function [d,m]=from_saved_figures(analysis)
f=fullfile(analysis,'fwd_rot_vs_dff_separate.fig');assert(isfile(f),'No saved data or recoverable trace figure');
fig=openfig(f,'invisible');clean=onCleanup(@()close(fig));axs=findall(fig,'Type','axes');assert(numel(axs)==3,'Expected exactly three trace axes');
pos=cell2mat(get(axs,'Position'));[~,order]=sort(pos(:,2),'descend');vars=cell(1,3);times=cell(1,3);
for j=1:3
 lines=findall(axs(order(j)),'Type','line');[~,ix]=max(arrayfun(@(h)numel(h.XData),lines));vars{j}=double(lines(ix).YData(:));times{j}=double(lines(ix).XData(:));
end
assert(all(cellfun(@(t)all(isfinite(t))&&all(diff(t)>0),times)),'Invalid recovered time axes');
d=struct('t',times{1},'dff',vars{1},'fwd',recover_interp(times{2},vars{2},times{1}),'yaw',recover_interp(times{3},vars{3},times{1}));d.valid=isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
[d.x,d.y]=path_from_figure(analysis,numel(d.t));
m=struct('source',f,'isHeat',false,'pathSource',fullfile(analysis,'x_color_linear_dff_vs_y_color_linear_dff.fig'),'yawSource','Original saved signed rotational velocity in deg/s','recoveredFromSavedFigures',true,'behaviorResampledToDffTime',~isequal(times{1},times{2}));
end

function y=recover_interp(t,x,tq)
if isequal(t,tq),y=x;return,end
y=interp1(t,x,tq,'linear',NaN);g=find(diff(t)>2*median(diff(t)));
for k=g',y(tq>t(k)&tq<t(k+1))=nan;end
end

function [x,y]=path_from_figure(analysis,n)
f=fullfile(analysis,'x_color_linear_dff_vs_y_color_linear_dff.fig');assert(isfile(f),'No saved flat path');
fig=openfig(f,'invisible');clean=onCleanup(@()close(fig));objs=findall(fig,'Type','scatter');assert(~isempty(objs),'Expected numeric scatter path in saved figure');
[~,j]=max(arrayfun(@(h)numel(h.XData),objs));x=double(objs(j).XData(:));y=double(objs(j).YData(:));
assert(numel(x)==n&&numel(y)==n,'Saved path length differs from trace length');
end

function validate(d)
assert(all(isfinite(d.t))&&all(diff(d.t)>0),'Time must be finite and increasing');
for key={'dff','fwd','yaw','x','y','valid'},assert(numel(d.(key{1}))==numel(d.t),'Source lengths differ');end
assert(nnz(d.valid)>20,'Insufficient valid behavior');
end

function B=blocks(t,valid)
breaks=[true;diff(t)>2*median(diff(t))];
s=find(valid&(breaks|[true;~valid(1:end-1)]));e=find(valid&([~valid(2:end);true]|[breaks(2:end);true]));B=[s e];
end

function [idx,info]=peaks(d)
fs=1/median(diff(d.t));base=movmedian(d.dff,max(3,round(30*fs)),'omitmissing');res=d.dff-base;res=res(isfinite(res));sigma=1.4826*mad(res,1);prom=3*sigma;
info=struct('noiseSigma',sigma,'prominence',prom,'minimumDistanceSamples',max(1,round(2*fs)),'minimumWidthSamples',max(1,round(.25*fs)));
idx=zeros(0,1);B=blocks(d.t,isfinite(d.dff));
for j=1:size(B,1)
 ix=(B(j,1):B(j,2))';if numel(ix)<3,continue,end
 [~,loc]=findpeaks(d.dff(ix),'MinPeakProminence',prom,'MinPeakDistance',info.minimumDistanceSamples,'MinPeakWidth',info.minimumWidthSamples);idx=[idx;ix(loc)]; %#ok<AGROW>
end
% Enforce the same separation if a finite segment boundary is shorter than 2s.
if any(diff(d.t(idx))<info.minimumDistanceSamples/fs)
 [~,rank]=sort(d.dff(idx),'descend');keep=false(size(idx));
 for j=rank',if ~any(abs(d.t(idx(j))-d.t(idx(keep)))<info.minimumDistanceSamples/fs),keep(j)=true;end,end
 idx=idx(keep);
end
end

function tf=complete_window(d,i)
t0=d.t(i);fs=1/median(diff(d.t));ix=find(d.t>=t0-2&d.t<=t0+3);
tf=t0-2>=d.t(1)&&t0+3<=d.t(end)&&all(d.valid(ix))&&all(diff(d.t(ix))<=2/fs);
% Also require the interpolation bracketing samples, never crossing a masked gap.
q=find(d.t<=t0-2,1,'last');r=find(d.t>=t0+3,1,'first');
if isempty(q)||isempty(r),tf=false;else,tf=tf&&all(d.valid(q:r));end
end

function [idx,angle]=slowdowns(d)
fs=1/median(diff(d.t));n=max(3,round(.25*fs));f=nan(size(d.fwd));yaw=nan(size(d.yaw));acc=f;B=blocks(d.t,d.valid);
for j=1:size(B,1)
 ix=B(j,1):B(j,2);f(ix)=smoothdata(d.fwd(ix),'movmedian',n);yaw(ix)=smoothdata(d.yaw(ix),'movmedian',n);
 if numel(ix)>2,acc(ix)=gradient(f(ix),d.t(ix));end
end
candidate=d.valid&acc<=-2;change=diff([false;candidate;false]);s=find(change==1);e=find(change==-1)-1;
idx=zeros(0,1);angle=zeros(0,1);last=-inf;
for j=1:numel(s)
 [~,loc]=min(acc(s(j):e(j)));i=s(j)+loc-1;t0=d.t(i);
 if t0-last<2||~complete_window(d,i),continue,end
 pre=d.t>=t0-.75&d.t<t0-.25;post=d.t>=t0+.25&d.t<=t0+.75;core=d.t>=t0&d.t<=t0+.5;
 before=median(f(pre));drop=before-median(f(post));if before<1||drop<1,continue,end
 idx(end+1,1)=i;angle(end+1,1)=pizza_ves041_turn_angle(d,i,[0 .5]);last=t0; %#ok<AGROW>
end
end

function event_figure(X,rel,name,stem)
fig=figure('Visible','off','Color','w','Position',[10 10 1700 1250]);clean=onCleanup(@()close(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
labels={'VES041 peaks','All slowdowns','Slowdowns without turns','Slowdowns with turns'};
lims=cell(1,3);
for v=1:3,z=[];for c=1:2,q=X{c}(:,:,v);z=[z;q(:)];end;lims{v}=limits(z);end
for row=1:4
 for c=1:4
  ax=nexttile(tl);n=size(X{c},2);hold(ax,'on');
  if row==1
   if n>0,imagesc(ax,rel,1:n,X{c}(:,:,1)');set(ax,'YDir','normal');ylim(ax,[.5 n+.5]);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');ylim(ax,[.5 1.5]);end
   caxis(ax,lims{1});colormap(ax,parula(256));cb=colorbar(ax);cb.Label.String='Original \DeltaF/F';ylabel(ax,'Event (chronological)');title(ax,sprintf('%s | n=%d',labels{c},n));
  else
   v=row-1;if n>0,plot(ax,rel,X{c}(:,:,v),'Color',[.8 .8 .8],'LineWidth',.4);plot(ax,rel,mean(X{c}(:,:,v),2),'k','LineWidth',2);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
   ylim(ax,lims{v});names={'Original \DeltaF/F','Forward velocity (mm/s)','Rotational speed |yaw| (deg/s)'};ylabel(ax,names{v});grid(ax,'on');
  end
  xline(ax,0,':','Color',[.8 .2 .2]);xlim(ax,[-2 3]);
  if row==4,xlabel(ax,'Time from event (s)');end
 end
end
title(tl,strrep(name,'_',' '),'FontWeight','bold');
subtitle(tl,{'Peaks: standard pipeline findpeaks | Slowdowns: maximum deceleration | Turn: accumulated |yaw| >=15 deg over 0:+0.5 s after slowdown', ...
 'Original dF/F; thin lines: events; black: mean | Columns 3 + 4 partition column 2 | Visual and odor epochs retained'},'FontSize',10);
exportgraphics(fig,[stem '.png'],'Resolution',170);savefig_visible_on_open(fig,[stem '.fig']);
end

function path_figure(d,name,stem,pathSource)
plot_pizza_flat_path_scatter(d,name,stem,pathSource);
end

function lim=limits(z)
z=z(isfinite(z));if isempty(z),lim=[0 1];else,lim=[min(z) max(z)];if lim(1)==lim(2),lim=lim+[-.5 .5];else,pad=.025*diff(lim);lim=lim+[-pad pad];end,end
end
