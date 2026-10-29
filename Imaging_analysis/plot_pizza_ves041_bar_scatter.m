function result=plot_pizza_ves041_bar_scatter(root)
% All valid calcium samples against bar position, split by bar velocity sign.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
out=fullfile(root,'bar_tuning','bar_velocity_scatter');if ~isfolder(out),mkdir(out);end
M=readtable(fullfile(root,'manifest.csv'),'TextType','string');M=M(M.status=="complete"&M.category=="bar",:);
R=struct([]);audit=cell(0,8);
for k=1:height(M)
 name=M.trial(k);fprintf('BAR SCATTER %s\n',name);
 D=load(fullfile(root,'analysis_data',name+"_events.mat"),'d');d=D.d;t=d.t(:);
 try
  recovered=false;quality=NaN;
  if endsWith(M.source(k),'.mat')
   S=load(M.source(k),'fly_data');a=S.fly_data;q=a.dq(1);
   assert(isequal(double(q.t(:)),t));vh=double(q.vh(:));source=M.source(k);
   jump=false(size(t));
   if isfield(a,'jump_detected')&&numel(a.jump_detected)==numel(t),jump=logical(a.jump_detected(:));end
   assert(~contains(lower(name),'jump')||isfield(a,'jump_detected'),'Jump labels missing');
  else
   [vh,quality,source]=recover_legacy(M.trialFolder(k),d);recovered=true;jump=false(size(t));
  end
  first=find(isfinite(vh),1);offset=vh(first);
  latest=startsWith(name,'202610');
  if latest,angle=wrapdeg(vh);reference="Recorded zero (front)";
  else,angle=wrapdeg(vh-offset);reference="Start-relative zero (provisional front)";end
  width=max(1,round(.2/median(diff(t))));velocity=movmedian(rad2deg(gradient(unwrap(vh),t)),width,'omitnan');
  away=true(size(t));for j=find(jump)',away(abs(t-t(j))<=.25)=false;end
  gaps=find(diff(t)>2*median(diff(t)));for j=gaps',away(t>=t(j)-.25&t<=t(j+1)+.25)=false;end
  valid=d.valid(:)&isfinite(vh)&isfinite(d.dff)&isfinite(velocity)&away;
  category=zeros(size(t));category(velocity>1)=1;category(velocity< -1)=-1;
  r=struct('trial',name,'source',source,'reference',reference,'recordedInitialDeg',rad2deg(offset), ...
   't',t,'dff',d.dff,'barPositionDeg',angle,'barVelocityDegS',velocity,'category',category,'valid',valid, ...
   'recoveredFromDaq',recovered,'recoveryYawCorrelation',quality);
  if isempty(R),R=r;else,R(end+1)=r;end %#ok<AGROW>
  audit(end+1,:)={name,"plotted",source,nnz(valid),nnz(valid&category==1),nnz(valid&category==-1),nnz(valid&category==0),quality}; %#ok<AGROW>
  fig=figure('Visible','off','Color','w','Position',[20 20 1550 500]);tl=tiledlayout(fig,1,3,'TileSpacing','compact');
  draw(nexttile(tl),r,0);title('All valid samples');
  draw(nexttile(tl),r,1);title('Positive bar velocity (>1 deg/s)');
  draw(nexttile(tl),r,-1);title('Negative bar velocity (<-1 deg/s)');
  axs=findall(fig,'Type','axes');ylimAll=[min(d.dff(valid)),max(d.dff(valid))];
  if diff(ylimAll)==0,ylimAll=ylimAll+[-.5 .5];end
  for ax=axs',ylim(ax,ylimAll+[-1 1]*.025*diff(ylimAll));end
  title(tl,strrep(name,'_','\_'));subtitle(tl,reference+" | Every valid sample; no bin averaging, subsampling, or dF/F clipping | Blue: +; orange: -; gray: stationary");
  saveplots(fig,fullfile(out,name+"_all_points_by_bar_velocity"));
  save(fullfile(out,name+"_all_points_by_bar_velocity.mat"),'r');
 catch ME
  audit(end+1,:)={name,"unavailable: "+string(ME.message),M.source(k),0,0,0,0,NaN}; %#ok<AGROW>
  warning('barScatter:Missing','%s: %s',name,ME.message);
 end
end
A=cell2table(audit,'VariableNames',{'trial','status','source','plottedSamples','positiveSamples','negativeSamples','stationarySamples','recoveryYawCorrelation'});
writetable(A,fullfile(out,'scatter_audit.csv'));assert(~isempty(R));
fig=figure('Visible','off','Color','w','Position',[20 20 1650 1100]);tl=tiledlayout(fig,ceil(numel(R)/3),3,'TileSpacing','compact');
for k=1:numel(R),ax=nexttile(tl);draw(ax,R(k),0);title(ax,{strrep(R(k).trial,'_','\_'),char(R(k).reference)},'FontSize',10);end
title(tl,'VES041 versus bar position: all valid samples, colored by bar velocity');
subtitle(tl,'Blue: >+1 deg/s | Orange: <-1 deg/s | Gray: within +/-1 deg/s | No bin averaging or dF/F clipping; peri-jump/gap samples excluded');
saveplots(fig,fullfile(out,'ves041_all_bar_trials_all_points'));
methods=struct('selection','All accepted bar trials in current pizza manifest; user-excluded recordings stay excluded', ...
 'scatter','One point per valid imaging sample; chronological draw order, alpha 0.25; no subsampling, bin averaging, or dF/F clipping', ...
 'reference','October recorded angle unchanged; older recordings relative to their first saved angle, provisionally front; no left/right claim', ...
 'velocity','Derivative of unwrapped synchronized bar angle, then 0.2 s median; positive >1 deg/s, negative <-1 deg/s, stationary +/-1 deg/s', ...
 'exclusions','Invalid behavior/calcium, nonfinite bar data, +/-0.25 s around saved jump labels and large time gaps', ...
 'activity','Original ROI dF/F, same-time bar signal; no z scoring, baseline subtraction, or calcium lag', ...
 'legacyRecovery','March DAQ bar voltage converted with berg2 0-10 V calibration; circular mean within saved imaging time bins; accept only if independently reconstructed fly yaw correlates >=0.9 with saved yaw');
save(fullfile(out,'bar_scatter_summary.mat'),'R','A','methods');
fid=fopen(fullfile(out,'methods.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('folder',out,'trials',numel(R),'unavailable',height(M)-numel(R));disp(result);
end
function [vh,quality,source]=recover_legacy(folder,d)
% Never infer bar angle from fly heading. Read the independently logged channel.
P=dir(fullfile(folder,'*daqData*.mat'));assert(isscalar(P),'Expected one raw DAQ file');source=string(fullfile(P.folder,P.name));
S=load(source,'trialData');Q=S.trialData;clear S
first=find(Q.frameClock~=0,1);last=find(Q.frameClock~=0,1,'last');assert(~isempty(first));
ix=(first:10:last)'; % 1 kHz input to imaging-bin means (imaging ~12 Hz).
tr=seconds(Q.Time(ix)-Q.Time(first));bar=double(Q.g4panels(ix))/10*2*pi+pi;
fly=double(Q.ficTracYaw(ix))/10*2*pi+pi;clear Q
t=d.t(:);bounds=[t(1);(t(1:end-1)+t(2:end))/2;t(end)+(t(end)-t(end-1))/2];
bin=discretize(tr,bounds);keep=isfinite(bin);bin=bin(keep);bar=bar(keep);fly=fly(keep);
coverage=numel(unique(bin))/numel(t);
fprintf('  DAQ covers %.3f%% of imaging samples (%.3f to %.3f s)\n',100*coverage,tr(1),tr(end));
assert(coverage>=.95,'DAQ covers less than 95 percent of saved calcium trace');
bs=accumarray(bin,sin(bar),[numel(t),1],@mean,NaN);bc=accumarray(bin,cos(bar),[numel(t),1],@mean,NaN);
vh=atan2(bs,bc);
fs=accumarray(bin,sin(fly),[numel(t),1],@mean,NaN);fc=accumarray(bin,cos(fly),[numel(t),1],@mean,NaN);
flyangle=atan2(fs,fc);recoveredYaw=rad2deg(gradient(unwrap(flyangle),t));
width=max(3,round(.25/median(diff(t))));recoveredYaw=movmean(recoveredYaw,width,'omitnan');savedYaw=movmean(d.yaw,width,'omitnan');
v=d.valid(:)&isfinite(recoveredYaw)&isfinite(savedYaw);quality=corr(recoveredYaw(v),savedYaw(v));
fprintf('  DAQ recovery independent yaw correlation %.4f\n',quality);
assert(quality>=.9,'DAQ timing verification failed: yaw correlation %.3f',quality);
end
function x=wrapdeg(v),x=rad2deg(atan2(sin(v),cos(v)));end
function draw(ax,r,which)
valid=r.valid;if which~=0,valid=valid&r.category==which;end
color=repmat([.65 .65 .65],nnz(valid),1);cat=r.category(valid);
color(cat==1,:)=repmat([0 .447 .741],nnz(cat==1),1);color(cat==-1,:)=repmat([.85 .325 .098],nnz(cat==-1),1);
h=scatter(ax,r.barPositionDeg(valid),r.dff(valid),8,color,'filled','MarkerFaceAlpha',.25,'MarkerEdgeAlpha',.25);
assert(numel(h.XData)==nnz(valid),'Not all requested points were plotted');
xlim(ax,[-180 180]);xticks(ax,-180:90:180);xlabel(ax,'Bar position (deg)');ylabel(ax,'VES041 \DeltaF/F');grid(ax,'on');
end
function saveplots(fig,stem)
exportgraphics(fig,stem+".png",'Resolution',150);savefig_visible_on_open(fig,stem+".fig");close(fig);
end
