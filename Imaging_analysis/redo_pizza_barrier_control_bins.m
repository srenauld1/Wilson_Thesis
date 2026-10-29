function R=redo_pizza_barrier_control_bins(root)
% Retain existing event detectors; use explicit control labels and voltage bins.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis';end
L=readtable(fullfile(root,'trial_locations.csv'),'Delimiter',',','TextType','string');L=L(L.geometry=="barrier",:);rows=cell(0,6);
for k=1:height(L)
 id=L.trial(k);folder=char(L.analysisFolder(k));dest=fullfile(folder,'control_and_heat_bins');
 if ~isfile(fullfile(folder,'analysis.mat')),rows(end+1,:)={id,'awaiting ROI extraction',0,0,0,dest};continue,end
 fprintf('BARRIER BINS %s\n',id);
 try
  if isfile(fullfile(dest,'complete.json')),j=jsondecode(fileread(fullfile(dest,'complete.json')));rows(end+1,:)={id,'complete',j.includedEntries,j.excludedEntries,j.pathTrials,dest};continue,end
  if ~isfolder(dest),mkdir(dest);end
  S=load(fullfile(folder,'analysis.mat'));d=S.d;P=protocol(S);rel=S.rel;
  names=["control_barrier","heat_0_to_lt1V","heat_1_to_lt2p5V","heat_ge2p5V"];
  C=cell(4,8);eventRows=cell(0,9);
  for c=1:8
   ix=S.groups{c};assign=zeros(numel(ix),1);
   for j=1:numel(ix)
    i=ix(j);lo=find(d.t<=d.t(i)-2,1,'last');hi=find(d.t>=d.t(i)+3,1,'first');win=lo:hi;core=d.t>=d.t(i)-.5&d.t<=d.t(i)+.5;
    voltage=median(P.heat(core));reason="";b=0;
    if isempty(lo)||isempty(hi)||~all(P.active(win))||any(~isfinite(P.control(win)))||any(~isfinite(P.trial(win)))
     reason="unplaced barrier or missing protocol";
    elseif numel(unique(P.trial(win)))~=1||numel(unique(P.control(win)))~=1||numel(unique(P.geometry(win)))~=1
     reason="window crosses trial, control, or geometry boundary";
    elseif P.control(i)==1,b=1;
    elseif P.control(i)~=0,reason="unknown control flag";
    elseif ~all(isfinite(P.heat(core))),reason="missing measured heat command";
    elseif voltage<-.05,reason="unexpected negative measured command";
    elseif max(0,voltage)<1,b=2;
    elseif voltage<2.5,b=3;
    else,b=4;end
    assign(j)=b;eventRows(end+1,:)={c,i,d.t(i),P.trial(i),P.control(i),P.geometry(i),voltage,b,reason};
   end
   for b=1:4,C{b,c}=S.X{c}(:,assign==b,:);end
  end
  E=cell2table(eventRows,'VariableNames',{'eventGroup','sampleIndex','eventTime','barrierTrial','controlTrial','geometry','medianHeatV','heatBin','exclusionReason'});
  writetable(E,fullfile(dest,'events.csv'));lims=cell(1,3);
  for v=1:3,z=[];for b=1:4,for c=1:8,q=C{b,c}(:,:,v);z=[z;q(:)];end,end;z=z(isfinite(z));if isempty(z),lims{v}=[0 1];else,pad=max(.01,.025*(max(z)-min(z)));lims{v}=[min(z)-pad max(z)+pad];end,end
  for b=1:4
   X=C(b,:);stem=fullfile(dest,char(id+"_"+names(b)+"_peaks_slowdowns"));
   pizza_barrier_render(X,rel,strrep(char(id),'_',' '), ...
    {strrep(char(names(b)),'_',' '),'Only placed barriers; full -2:+3 s window within one trial and geometry | Heat: median measured command +/-0.5 s; bins 0/1/2.5 V', ...
    'Original dF/F | Gray: events; black: mean | Shared scales across all four conditions | Zero-voltage heat trials remain in low-heat bin'},stem,false,strings(0,1),cellfun(@(x)size(x,2),X),lims);
  end
  [pathTrials,pathPages]=paths(d,P,id,dest);
  settings=struct('controlDefinition','control_trial==1 during placed barrier only; accommodation and lead-in excluded', ...
   'heatBinEdgesV',[0 1 2.5 Inf],'heatBinWindowSec',[-.5 .5],'minimumAllowedMedianHeatV',-.05,'eventWindowSec',[-2 3], ...
   'fullWindowSameTrialGeometryAndControl',true,'sourceProtocol',P.source,'activeDefinition',P.activeDefinition,'pathPages',pathPages, ...
   'pathDisplay','Separate active barrier episodes, shared 5th/95th dF/F color limits within recording; all finite samples retained', ...
   'eventDetectors','Existing analysis.mat groups and arrays reused unchanged; turn threshold 15 deg and original time windows retained');
  save(fullfile(dest,'binned_analysis.mat'),'C','E','rel','P','settings','-v7.3');
  report=struct('includedEntries',nnz(E.heatBin>0),'excludedEntries',nnz(E.heatBin==0),'pathTrials',pathTrials);
  fid=fopen(fullfile(dest,'complete.json'),'w');fprintf(fid,'%s',jsonencode(report));fclose(fid);
  % Previous broad voltage split is retained away from the revised figures.
  old=fullfile(folder,'previous_voltage_on_off');if ~isfolder(old),mkdir(old);end
  files=[dir(fullfile(folder,'*_heat_on_peaks_slowdowns.*'));dir(fullfile(folder,'*_heat_off_peaks_slowdowns.*'))];
  for j=1:numel(files),movefile(fullfile(files(j).folder,files(j).name),fullfile(old,files(j).name));end
  rows(end+1,:)={id,'complete',report.includedEntries,report.excludedEntries,pathTrials,dest};
 catch ME
  rows(end+1,:)={id,['FAILED: ' ME.message],0,0,0,dest};fprintf('%s\n',getReport(ME));
 end
 R=cell2table(rows,'VariableNames',{'trial','status','includedGroupEntries','excludedGroupEntries','flatPathTrials','output'});writetable(R,fullfile(root,'barrier','control_heat_bins_manifest.csv'));
end
R=cell2table(rows,'VariableNames',{'trial','status','includedGroupEntries','excludedGroupEntries','flatPathTrials','output'});writetable(R,fullfile(root,'barrier','control_heat_bins_manifest.csv'));
end

function P=protocol(S)
d=S.d;t=d.t;f=char(S.meta.source);P=struct;P.source=f;
if endsWith(f,'fly_data.mat')
 A=load(f,'fly_data');a=A.fly_data;q=a.dq(1);p=a.protocol;assert(p.geometryAligned&&isequal(t,double(q.t(:))));
 P.control=p.control_trial(:);P.trial=p.trial_number(:);P.active=p.barrier_active(:)>0;P.heat=double(q.heat(:));
 P.geometry=repmat("barrier",size(t));if isfield(p,'heat_geometry'),P.geometry=string(p.heat_geometry(:));end
 P.x=p.pos_x_mm(:);P.y=p.pos_y_mm(:);P.bx=p.barrier_x_mm(:);P.by=p.barrier_y_mm(:);P.nx=p.barrier_normal_x(:);P.ny=p.barrier_normal_y(:);
 P.activeDefinition='Recorded protocol barrier_active >0';
else
 folder=fileparts(f);v=dir(fullfile(folder,'*virmenData.mat'));
 if ~isempty(v)
  A=load(fullfile(v.folder,v.name),'vmData');V=A.vmData;tk=V.Seconds{1};if isduration(tk),tk=seconds(tk);end;tk=double(tk(:));
  ctl=col(V,'control_trial');type=col(V,'currTt');dist=col(V,'heat_barrier_distance');
  active=ismember(round(type),[1 2])&isfinite(dist)&abs(dist)<99999;
  assert(all(abs(type-round(type))<1e-6|~isfinite(type)),'Noninteger trial labels in legacy source');
  starts=active&[true;~active(1:end-1)|diff(type)~=0|diff(ctl)~=0];trial=cumsum(starts);trial(~active)=nan;
  P.control=interp1(tk,ctl,t,'previous',nan);P.trial=interp1(tk,trial,t,'previous',nan);P.active=interp1(tk,double(active),t,'previous',0)>0;
  P.x=interp1(tk,col(V,'pos_x'),t,'linear',nan);P.y=interp1(tk,col(V,'pos_y'),t,'linear',nan);
  P.heat=interp1(tk,col(V,'heatDAQ'),t,'linear',nan);P.source=fullfile(v.folder,v.name);
  P.activeDefinition='Legacy currTt 1/2 and finite saved barrier distance with |distance|<99999; sentinel100000 excluded';
 else
  [P.control,P.trial,P.active,P.x,P.y]=native_legacy(folder,f,t);
  A=load(f,'dq');P.heat=double(A.dq(1).heat(:));P.activeDefinition='Raw recorded trial type1/2 after 5-second placement delay, aligned via recorded DAQ iteration clock';
 end
 P.geometry=repmat("barrier",size(t));for key={'bx','by','nx','ny'},P.(key{1})=nan(size(t));end
end
assert(all(cellfun(@(key)numel(P.(key))==numel(t),{'control','trial','active','heat','x','y','geometry'})));
P.active=P.active&isfinite(P.x)&isfinite(P.y)&~ismissing(P.geometry);
end
function x=col(V,key)
x=V.(key);if iscell(x),x=x{1};end;x=double(x(:));
end

function [control,trial,active,x,y]=native_legacy(folder,file,t)
S=load(file,'dq');q=S.dq(2);R=load(fullfile(folder,'virmenData','vr_test_1.mat'),'data');a=R.data;n=numel(a.time);
% The final custom-state write can outlast the final position/time sample.
assert(numel(a.currTt)>=n&&numel(a.control_trial)>=n);type=double(a.currTt(1:n))';ctl=double(a.control_trial(1:n))';rt=double(a.time(:));
assert(all(diff(rt)>0));starts=[true;diff(type)~=0];trialRaw=cumsum(starts);beg=find(starts);elapsed=zeros(n,1);
for j=1:numel(beg),last=n;if j<numel(beg),last=beg(j+1)-1;end;ix=beg(j):last;elapsed(ix)=rt(ix)-rt(beg(j));end
placed=ismember(type,[1 2])&elapsed>5;
clock=double(q.iter(:));u=clock+10*cumsum([false;diff(clock)<-8]);idx=round(10*u)-1;offset=round((n-idx(end))/100)*100;idx=idx+offset;
assert(abs(n-idx(end))<100,'Iteration-clock alignment failed');ok=idx>=1&idx<=n;ii=idx(ok);qt=double(q.t(ok));
[ii,keep]=unique(ii,'stable');qt=qt(keep);assert(all(diff(ii)>0)&&all(diff(qt)>0),'Nonmonotonic iteration clock');
control=interp1(qt,ctl(ii),t,'previous',nan);trial=interp1(qt,trialRaw(ii),t,'previous',nan);active=interp1(qt,double(placed(ii)),t,'previous',0)>0;
x=interp1(qt,double(a.position(ii,1)),t,'linear',nan);y=interp1(qt,double(a.position(ii,2)),t,'linear',nan);
end

function [n,pages]=paths(d,P,id,dest)
t=d.t;valid=P.active&isfinite(P.control)&isfinite(P.trial)&isfinite(d.dff)&d.valid;
trialIds=unique(P.trial(valid),'stable');seg=zeros(size(t));for j=1:numel(trialIds),seg(valid&P.trial==trialIds(j))=j;end;ids=(1:numel(trialIds))';n=numel(ids);pages=ceil(n/12);
z=d.dff(valid);lim=[0 1];if ~isempty(z),lim=prctile(z,[5 95]);if diff(lim)==0,lim=lim+[-.1 .1];end,end
T=cell(0,7);for j=1:n,ix=find(seg==ids(j));T(end+1,:)={j,P.trial(ix(1)),P.control(ix(1)),t(ix(1)),t(ix(end)),numel(ix),median(P.heat(ix),'omitnan')};end
writetable(cell2table(T,'VariableNames',{'pathPanel','barrierTrial','controlTrial','startTime','endTime','samples','medianHeatV'}),fullfile(dest,'flat_path_trials.csv'));
for page=1:pages
 fig=figure('Visible','off','Color','w','Position',[10 10 1800 1400]);tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','loose');
 for j=(page-1)*12+1:min(page*12,n)
  ax=nexttile(tl);hold(ax,'on');ix=find(seg==ids(j));[~,order]=sort(d.dff(ix));draw=ix(order);scatter(ax,P.x(draw),P.y(draw),10,d.dff(draw),'filled','MarkerEdgeColor','none');
  plot(ax,P.x(ix(1)),P.y(ix(1)),'go','MarkerFaceColor','g','MarkerSize',4);
  q=ix(find(isfinite(P.bx(ix))&isfinite(P.nx(ix)),1));
  if ~isempty(q),length=max(10,max([range(P.x(ix)),range(P.y(ix))]));v=[-P.ny(q),P.nx(q)];plot(ax,P.bx(q)+[-1 1]*length*v(1),P.by(q)+[-1 1]*length*v(2),'k--','LineWidth',1);end
  axis(ax,'equal');grid(ax,'on');clim(ax,lim);cb=colorbar(ax);cb.Label.String='Original dF/F';xlabel(ax,'X (mm)');ylabel(ax,'Y (mm)');label='Heat';if P.control(ix(1))==1,label='Control';end
  title(ax,sprintf('%s barrier %g | %.1f-%.1f s',label,P.trial(ix(1)),t(ix(1)),t(ix(end))));
 end
 colormap(fig,parula(256));title(tl,strrep(char(id),'_',' '),'FontWeight','bold');subtitle(tl,'Placed-barrier intervals only | shared 5th-95th percentile dF/F color scale; all samples retained | green: first point; dashed: recorded barrier where available');
 stem=fullfile(dest,sprintf('flat_paths_dff_page_%02d',page));exportgraphics(fig,[stem '.png'],'Resolution',130);savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
end

function pizza_barrier_render(X,rel,titleText,subtitleText,stem,isNormalized,flies,eventCounts,sharedLimits)
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
exportgraphics(fig,[stem '.png'],'Resolution',130);savefig_visible_on_open(fig,[stem '.fig']);
end

