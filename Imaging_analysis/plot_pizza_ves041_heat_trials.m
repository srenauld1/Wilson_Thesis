function M=plot_pizza_ves041_heat_trials(out)
% Match the noheat eight-column figures; classify by measured heat at t=0.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis';end
I=readtable(fullfile(imaging_analysis_artifact_root(),'reports','pizza_ves041_heat_cohort.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');
if ~isfolder(out),mkdir(out);end
rows=cell(0,7);rel=(-2:.05:3)';prior=table();if isfile(fullfile(out,'manifest.csv')),prior=readtable(fullfile(out,'manifest.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');end
for k=1:height(I)
 name=I.trial(k);if ~isempty(prior),old=prior(prior.trial==name & prior.status=="complete",:);if ~isempty(old),rows(end+1,:)=table2cell(old);continue,end,end
 name=I.trial(k);fprintf('HEAT %s\n',name);dest=fullfile(out,heat_geometry(I(k,:)),name);if ~isfolder(dest),mkdir(dest);end
 try
  [d,meta]=load_heat(I(k,:));validate_heat(d);
  [p,peakDetector]=peaks(d);p=p(arrayfun(@(i)complete_window(d,i),p));
  [s,turn]=slowdowns(d);[r,drop,~,~,rotationDetector]=pizza_heat_rotation_peaks(d);
  [u,angle,~,~,speedupDetector]=pizza_ves041_speedups(d);
  groups={p,s,s(turn<15),s(turn>=15),r(drop<1),r(drop>=1),u(angle<15),u(angle>=15)};
  X=cell(1,8);vars={d.dff,d.fwd,abs(d.yaw)};
  for c=1:8
   ix=groups{c};ix=ix(isfinite(d.heat(ix)));groups{c}=ix;X{c}=nan(numel(rel),numel(ix),3);
   for j=1:numel(ix),for v=1:3,X{c}(:,j,v)=interp1(d.t,vars{v},d.t(ix(j))+rel);end,end
   assert(all(isfinite(X{c}(:))));
  end
  % Shared limits across heat-on/off preserve visual comparability.
  lims=cell(1,3);for v=1:3,z=[];for c=1:8,q=X{c}(:,:,v);z=[z;q(:)];end;z=z(isfinite(z));if isempty(z),lims{v}=[0 1];else,lo=min(z);hi=max(z);pad=max(.01,.025*(hi-lo));lims{v}=[lo-pad hi+pad];end,end
  eventRows=cell(0,7);counts=zeros(2,8);
  for state=0:1
   Y=cell(1,8);
   for c=1:8
    ix=groups{c};sel=(d.heat(ix)>meta.threshold)==logical(state);Y{c}=X{c}(:,sel,:);counts(state+1,c)=nnz(sel);
    for j=find(sel)'
     win=d.t>=d.t(ix(j))-2&d.t<=d.t(ix(j))+3;h=d.heat(win);transition=any(~isfinite(h))||any((h>meta.threshold)~=logical(state));
     eventRows(end+1,:)={c,ix(j),d.t(ix(j)),state,d.heat(ix(j)),transition,mean(h>meta.threshold)};
    end
   end
   label='heat_off';if state,label='heat_on';end
   stem=fullfile(dest,name+"_"+label+"_peaks_slowdowns");
   pizza_heat_render(Y,rel,strrep(char(name),'_',' '), ...
    {sprintf('%s at event time | %s > %g | Shared heat-on/off axis limits',strrep(label,'_',' '),meta.heatSource,meta.threshold), ...
    'Original dF/F | Gray: events; black: mean | Window -2:+3 s; transitions retained and flagged | Turn: slowdown 0:+0.5 s; speed-up -0.5:0 s'},char(stem),false,strings(0,1),cellfun(@(x)size(x,2),Y),lims);
  end
  E=cell2table(eventRows,'VariableNames',{'group','sampleIndex','eventTime','heatOn','heatAtEvent','windowCrossesHeatState','heatOnFraction'});writetable(E,fullfile(dest,'events.csv'));
  assert(all(sum(counts,1)==cellfun(@numel,groups)));
  save(fullfile(dest,'analysis.mat'),'d','meta','groups','X','rel','E','counts','peakDetector','rotationDetector','speedupDetector','-v7.3');
  rows(end+1,:)={name,'complete',meta.source,meta.heatSource,meta.threshold,sum(counts(2,:)),sum(counts(1,:))};
 catch ME
  status=['FAILED: ' ME.message];if contains(ME.message,'Awaiting ROI extraction'),status='awaiting ROI extraction';end
  rows(end+1,:)={name,status,'','',nan,0,0};fprintf('%s\n',getReport(ME));
 end
 M=cell2table(rows,'VariableNames',{'trial','status','source','heatSource','heatThreshold','heatOnGroupEntries','heatOffGroupEntries'});writetable(M,fullfile(out,'manifest.csv'));
end
M=cell2table(rows,'VariableNames',{'trial','status','source','heatSource','heatThreshold','heatOnGroupEntries','heatOffGroupEntries'});writetable(M,fullfile(out,'manifest.csv'));
end

function [d,m]=load_heat(R)
folder=char(R.folder);a=R.analysis_folder;if ismissing(a)||strlength(a)==0,a=fullfile(folder,'smr_analysis');end
f=fullfile(a,'fly_data.mat');native=false;
if ~isfile(f)&&isempty(dir(fullfile(folder,'*virmenData.mat'))),z=dir(fullfile(folder,'*_or_z0z1_s_.mat'));if numel(z)==1,f=fullfile(folder,z.name);native=true;end,end
m=struct('source',char(f),'heatSource','DAQ heat command','threshold',.05);
if isfile(f)
 if native,S=load(f,'dq','roi');a=S;else,S=load(f,'fly_data');a=S.fly_data;end;q=a.dq(1);d.t=double(q.t(:));d.dff=double(a.roi(1).ts(1,:))';d.fwd=double(q.bvf(:));
 if isfield(q,'bvy_deg'),d.yaw=double(q.bvy_deg(:));else,d.yaw=rad2deg(double(q.bvy(:)));end
 d.heat=double(q.heat(:));d.valid=isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
 if isfield(q,'bergii')&&isfield(q.bergii,'behaviorVolumeValid'),d.valid=d.valid&logical(q.bergii.behaviorVolumeValid(:));end
else
 ff=dir(fullfile(folder,'*virmenData.mat'));ii=dir(fullfile(folder,'*_or_z0z1_s_.mat'));assert(~isempty(ii),'Awaiting ROI extraction: no processed fluorescence trace');assert(numel(ff)==1&&numel(ii)==1,'Missing or ambiguous legacy source files');
 S=load(fullfile(folder,ff.name));V=S.vmData;A=load(fullfile(folder,ii.name),'roi');
 t=S.daqVolumeTime.daqVolumeTime{1,1};t=double(t(:));d.t=t;d.dff=double(A.roi.ts(1,:))';
 tk=V.Seconds{1,1};if isduration(tk),tk=seconds(tk);end;tk=double(tk(:));[tk,ix]=unique(tk,'stable');
 f=V.vel_for{1,1};y=rad2deg(V.vel_yaw{1,1});h=V.heat{1,1};
 d.fwd=interp1(tk,double(f(ix)),t,'linear',nan);d.yaw=interp1(tk,double(y(ix)),t,'linear',nan);d.heat=interp1(tk,double(h(ix)),t,'previous',nan);
 d.valid=isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);m.source=fullfile(folder,ii.name);m.heatSource='Virmen heat command';m.threshold=.0001;
end
end
function validate_heat(d)
assert(all(isfinite(d.t))&&all(diff(d.t)>0),'Invalid time');for key={'dff','fwd','yaw','heat','valid'},assert(numel(d.(key{1}))==numel(d.t),'Length mismatch');end
assert(any(isfinite(d.heat)),'Missing heat channel');
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
 ix=(B(j,1):B(j,2))';if numel(ix)<3||d.t(ix(end))-d.t(ix(1))<5,continue,end
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


function pizza_heat_render(X,rel,titleText,subtitleText,stem,isNormalized,flies,eventCounts,sharedLimits)
% X contains original event traces or fly-mean traces, one cell per column.
if nargin<7,flies=strings(0,1);end
if nargin<8,eventCounts=cellfun(@(x)size(x,2),X);end
pooled=~isempty(flies);assert(numel(X)==8);
fig=figure('Visible','off','Color','w','Position',[10 10 3300 1400]);cleanup=onCleanup(@()close(fig));tl=tiledlayout(fig,4,8,'TileSpacing','compact','Padding','loose');
labels={{'VES041 peaks'},{'All slowdowns'},{'Slowdowns without turns'},{'Slowdowns with turns'}, ...
 {'Rotation peaks','Forward drop <1 mm/s'},{'Rotation peaks','Forward drop >=1 mm/s'}, ...
 {'Speed-ups','Without rotation (<15 deg)'},{'Speed-ups','With rotation (>=15 deg)'}};
lims=cell(1,3);allZ=[];
for v=1:3
 z=[];for c=1:8,q=X{c}(:,:,v);z=[z;q(:)];end;z=z(isfinite(z)); %#ok<AGROW>
 if isNormalized,allZ=[allZ;z];elseif isempty(z),lims{v}=[0 1];else,lims{v}=[min(z) max(z)];if diff(lims{v})==0,lims{v}=lims{v}+[-.5 .5];else,lims{v}=lims{v}+[-1 1]*.025*diff(lims{v});end,end %#ok<AGROW>
end
if nargin>=9,lims=sharedLimits;end
if isNormalized,bound=max(1,ceil(2*max(abs(allZ)))/2);for v=1:3,lims{v}=[-bound bound];end,end
flyLabel=replace(flies,'_fly',' F');flyLabel=replace(flyLabel,'_unlabeled_inferred',' *');flyLabel=replace(flyLabel,'_inferred','*');
for row=1:4,for c=1:8
 ax=nexttile(tl);hold(ax,'on');n=size(X{c},2);
 if row==1
  if n>0,im=imagesc(ax,rel,1:n,X{c}(:,:,1)');im.AlphaData=isfinite(im.CData);set(ax,'YDir','normal','Color',[.9 .9 .9]);ylim(ax,[.5 n+.5]);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');ylim(ax,[.5 1.5]);end
  caxis(ax,lims{1});cb=colorbar(ax);if isNormalized,cb.Label.String='Fly mean \DeltaF/F (z-score)';else,cb.Label.String='Original \DeltaF/F';end
  if pooled,yticks(ax,1:n);yticklabels(ax,flyLabel);ax.YAxis.FontSize=8;ylabel(ax,'Fly (* inferred identity)');count=sprintf('%d flies | %d events',nnz(any(isfinite(X{c}(:,:,1)),1)),eventCounts(c));else,ylabel(ax,'Event (chronological)');count=sprintf('n=%d',eventCounts(c));end
  title(ax,[labels{c} {count}]);
 else
  v=row-1;if n>0,plot(ax,rel,X{c}(:,:,v),'Color',[.75 .75 .75],'LineWidth',.5);plot(ax,rel,mean(X{c}(:,:,v),2,'omitnan'),'k','LineWidth',2.2);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
  if isNormalized,names={'\DeltaF/F (z-score)','Forward velocity (z-score)','|Yaw| speed (z-score)'};yline(ax,0,':','Color',[.45 .45 .45]);else,names={'Original \DeltaF/F','Forward velocity (mm/s)','|Yaw| speed (deg/s)'};end
  ylabel(ax,names{v});ylim(ax,lims{v});grid(ax,'on');
 end
 xline(ax,0,':','Color',[.8 .2 .2]);xlim(ax,[-2 3]);if row==4,xlabel(ax,'Time from event (s)');end
end,end
if isNormalized,cmap=[linspace(.15,1,128)' linspace(.3,1,128)' ones(128,1);ones(128,1) linspace(1,.2,128)' linspace(1,.15,128)'];colormap(fig,cmap);else,colormap(fig,parula(256));end
title(tl,titleText,'FontWeight','bold');subtitle(tl,subtitleText,'FontSize',10);assert(numel(findall(fig,'Type','axes'))==32);
exportgraphics(fig,[stem '.png'],'Resolution',170);savefig_visible_on_open(fig,[stem '.fig']);
end

function [idx,drop,before,after,info]=pizza_heat_rotation_peaks(d)
t=d.t;fs=1/median(diff(t));valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);cut=[true;diff(t)>2/fs];
s=find(valid&(cut|[true;~valid(1:end-1)]));e=find(valid&([~valid(2:end);true]|[cut(2:end);true]));
r=nan(size(t));f=r;for k=1:numel(s),ix=s(k):e(k);r(ix)=abs(smoothdata(d.yaw(ix),'gaussian',max(3,round(.12*fs))));f(ix)=smoothdata(d.fwd(ix),'movmedian',max(3,round(.25*fs)));end
baseline=movmedian(r,max(3,round(5*fs)),'omitmissing');res=r-baseline;res=res(isfinite(res));prom=max(20,2*1.4826*mad(res,1));
info=struct('prominence',prom,'minimumHeight',20,'minimumDistanceSec',1.5);candidates=zeros(0,1);
for k=1:numel(s)
 ix=(s(k):e(k))';if numel(ix)<3||t(ix(end))-t(ix(1))<5,continue,end
 [~,loc]=findpeaks(r(ix),'MinPeakHeight',20,'MinPeakProminence',prom,'MinPeakDistance',max(1,round(1.5*fs)));candidates=[candidates;ix(loc)]; %#ok<AGROW>
end
idx=zeros(0,1);drop=idx;before=idx;after=idx;
for i=candidates'
 lo=find(t<=t(i)-2,1,'last');hi=find(t>=t(i)+3,1,'first');
 if isempty(lo)||isempty(hi)||~all(valid(lo:hi))||any(diff(t(lo:hi))>2/fs),continue,end
 pre=t>=t(i)-.75&t<t(i)-.25;post=t>=t(i)+.25&t<=t(i)+.75;
 a=median(f(pre));b=median(f(post));idx(end+1,1)=i;before(end+1,1)=a;after(end+1,1)=b;drop(end+1,1)=a-b; %#ok<AGROW>
end
end

function geometry=heat_geometry(R)
text=lower(R.trial+" "+R.folder+" "+fillmissing(R.context,'constant',""));
if contains(text,"pyramid"),geometry="pyramid";
elseif contains(text,"stripe"),geometry="stripes";
elseif contains(text,"gradient"),geometry="barrier";
elseif contains(text,"multiple"),geometry="multiple_barriers";
else,geometry="barrier";end
end
