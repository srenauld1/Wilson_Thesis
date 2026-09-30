function summary=compare_ves041_peak_amplitude_pivot_all_trials()
% Compare merged VES041 event amplitude for pivot versus other outcomes.
% Includes all good-behavior trials and excludes heat-on, moving-grating,
% and +/-2 s bar-jump epochs. Pivot = mean forward speed <=2 mm/s and net
% heading change >=15 deg during +0.25 to +0.75 s after the event. Low forward
% speed is defined here as <=5 mm/s.

base='/Users/sophiarenauld/stacks';out=fullfile(base,'ves041_behavior_quantification','peak_amplitude_pivot_all_trials_strict_15deg_5mms_post025_075');
if ~isfolder(out),mkdir(out);end
standard=["20250812-20_ves041_fly3_clbarjumps";"20250812-21_ves041_fly3_olgrating";"20250812-22_ves041_fly3_bbarjumps2"; ...
 "20250813-3_ves041_fly1_olgrating";"20250819-10_ves041_fly1_olgrating5";"20251001-21_ves041_fly3_clbbar"; ...
 "20260213-2_ves041_clbar";"20260213-4_ves041_fly2_clbar";"20260217-2_ves041_vt_clbbar1"; ...
 "20260217-3_ves041_vt_clbbar2_actuallydark";"20260217-4_ves041_vt_clbbar3"];
heat=["20260427-2_ves041_fly2_heat";"20260427-3_ves041_fly2_heat_20mmbarrier"; ...
 "20260427-4_ves041_fly2_heat_20mmbarrier_2";"20260427-6_ves041_fly3_heat_80mm";"20260427-7_ves041_fly3_heat_20mm"];
rows=cell(0,11);eventsRows=cell(0,13);normOther=[];normPivot=[];
for trial=[standard;heat]'
    if any(heat==trial),[t,dff,fwd,rot,heading,stim]=load_heat(fullfile(base,trial));
    else,[t,dff,fwd,rot,heading,stim]=load_standard(fullfile(base,trial),trial);end
    n=min([numel(t),numel(dff),numel(fwd),numel(rot),numel(heading),numel(stim)]);t=t(1:n);dff=dff(1:n);fwd=fwd(1:n);rot=rot(1:n);heading=heading(1:n);stim=stim(1:n);
    fs=1/median(diff(t));[peaks,amp,absolutePeak,threshold]=events(dff,fs);
    preN=round(2*fs);postN=round(fs);forbidden=movmax(double(stim),[postN preN])>0;
    keep=peaks>preN & peaks+postN<=n & ~forbidden(peaks);peaks=peaks(keep);amp=amp(keep);absolutePeak=absolutePeak(keep);
    postStart=round(.25*fs);postEnd=round(.75*fs);
    postFwd=arrayfun(@(i)mean(fwd(i+postStart:i+postEnd),'omitnan'),peaks);
    preFwd=arrayfun(@(i)mean(fwd(i-preN:i-round(fs)),'omitnan'),peaks);
    preRot=arrayfun(@(i)mean(abs(rot(i-preN:i-round(fs))),'omitnan'),peaks);
    postRot=arrayfun(@(i)mean(abs(rot(i+postStart:i+postEnd)),'omitnan'),peaks);
    turn=arrayfun(@(i)abs(rad2deg(atan2(sin(heading(i+postEnd)-heading(i+postStart)),cos(heading(i+postEnd)-heading(i+postStart))))),peaks);
    pivot=postFwd<=5 & turn>=15;other=~pivot;
    medP=median(amp(pivot),'omitnan');medO=median(amp(other),'omitnan');meanP=mean(amp(pivot),'omitnan');meanO=mean(amp(other),'omitnan');
    ratio=medP/medO;rows(end+1,:)={trial,numel(peaks),sum(pivot),meanP,meanO,medP,medO,ratio,mean(absolutePeak(pivot),'omitnan'),mean(absolutePeak(other),'omitnan'),threshold}; %#ok<AGROW>
    scale=medO;if ~isfinite(scale)||scale<=0,scale=median(amp,'omitnan');end
    normPivot=[normPivot;amp(pivot)/scale];normOther=[normOther;amp(other)/scale]; %#ok<AGROW>
    for i=1:numel(peaks),eventsRows(end+1,:)={trial,t(peaks(i)),amp(i),absolutePeak(i),preFwd(i),postFwd(i),postFwd(i)-preFwd(i),preRot(i),postRot(i),turn(i),pivot(i),amp(i)/scale,threshold};end %#ok<AGROW>
end
summary=cell2table(rows,'VariableNames',{'trial','nEvents','nPivots','meanAmplitudePivot','meanAmplitudeOther', ...
 'medianAmplitudePivot','medianAmplitudeOther','pivotToOtherMedianRatio','meanAbsoluteDffPivot','meanAbsoluteDffOther','eventThreshold'});
summary.trial=string(summary.trial);writetable(summary,fullfile(out,'trial_summary.csv'));
E=cell2table(eventsRows,'VariableNames',{'trial','eventTime','amplitudeAboveBaseline','absoluteSmoothedDff','preFwd','postFwd', ...
 'changeFwd','preAbsRot','postAbsRot','postHeadingChangeDeg','isLowSpeedPivot','amplitudeNormalizedToTrialOtherMedian','eventThreshold'});E.trial=string(E.trial);writetable(E,fullfile(out,'events.csv'));

fig=figure('Visible','off','Color','w','Position',[80 80 1450 650]);tl=tiledlayout(1,2,'TileSpacing','compact');
ax=nexttile(tl);hold(ax,'on');valid=summary.nPivots>0;
for i=find(valid)',plot(ax,[1 2],[summary.medianAmplitudeOther(i) summary.medianAmplitudePivot(i)],'-o','LineWidth',1.3,'DisplayName',short(summary.trial(i)));end
xticks(ax,[1 2]);xticklabels(ax,{'Other outcome','Low-speed pivot'});ylabel(ax,'Median amplitude above baseline (dF/F)');title(ax,'Within-trial medians');grid(ax,'on');
ax=nexttile(tl);hold(ax,'on');scatter(ax,ones(size(normOther))+.12*(rand(size(normOther))-.5),normOther,12,[.6 .6 .6],'filled','MarkerFaceAlpha',.25);
scatter(ax,2*ones(size(normPivot))+.12*(rand(size(normPivot))-.5),normPivot,18,[.82 .12 .1],'filled','MarkerFaceAlpha',.4);
plot(ax,[.8 1.2],repmat(median(normOther),1,2),'k','LineWidth',3);plot(ax,[1.8 2.2],repmat(median(normPivot),1,2),'Color',[.8 .1 .1],'LineWidth',3);
xticks(ax,[1 2]);xticklabels(ax,{sprintf('Other (n=%d)',numel(normOther)),sprintf('Pivot (n=%d)',numel(normPivot))});
ylabel(ax,'Amplitude / trial non-pivot median');title(ax,'Pooled after within-trial normalization');grid(ax,'on');xlim(ax,[.5 2.5]);
sgtitle(tl,'VES041 event amplitude and low-speed pivot outcome across good trials');
savefig(fig,fullfile(out,'all_trials_peak_amplitude_pivot.fig'));exportgraphics(fig,fullfile(out,'all_trials_peak_amplitude_pivot.png'),'Resolution',250);close(fig);disp(summary);
end

function [t,dff,fwd,rot,heading,stim]=load_standard(folder,trial)
q=load(fullfile(folder,'smr_analysis','fly_data.mat'));x=q.fly_data;t=x.dq(1).t(:);dff=x.roi(1).ts(:);fwd=x.dq(1).bvf(:);rot=x.dq(1).bvy_deg(:);heading=x.dq(1).bh(:);stim=false(size(t));
if contains(lower(trial),'grating'),vv=x.dq(1).vvy(:);m=min(numel(stim),numel(vv));stim(1:m)=abs(vv(1:m))>.02;end
if isfield(x,'jump_detected'),j=find(logical(x.jump_detected(:)));j=j(j<=numel(t));if isfield(x.dq(1),'vh'),vh=x.dq(1).vh(:);dv=atan2(sin([vh(1);diff(vh)]),cos([vh(1);diff(vh)]));j=j(abs(rad2deg(dv(j)))>=45);end;fs=1/median(diff(t));for k=j(:)',stim(max(1,k-round(2*fs)):min(numel(stim),k+round(2*fs)))=true;end,end
end

function [t,dff,fwd,rot,heading,stim]=load_heat(folder)
f=dir(fullfile(folder,'*_virmenData.mat'));v=load(fullfile(f(1).folder,f(1).name));f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));q=load(fullfile(f(1).folder,f(1).name));
t=seconds(v.vmData.Seconds{1,1});t=t(:);fwd=smooth(v.vmData.vel_for{1,1});rot=rad2deg(smooth(v.vmData.vel_yaw{1,1}));heading=v.vmData.heading{1,1};heading=heading(:);h=v.vmData.heat;if iscell(h),h=cell2mat(h);end;stim=h(:)>1e-4;
td=v.daqVolumeTime.daqVolumeTime{1,1};dff=interp1(td(:),q.roi(1).ts(:),t,'linear',NaN);
end

function [peaks,amps,absolutePeak,threshold]=events(dff,fs)
% Strict, recording-adaptive event detector. An event must rise at least
% 1 dF/F and four robust noise SD above its 30-s local baseline. It must
% also contain a peak with that much prominence, remain elevated for at
% least 0.25 s, and is merged across interruptions of at most 0.5 s.
sm=smoothdata(dff,'gaussian',max(3,round(.25*fs)));
base=movmedian(sm,max(3,round(30*fs)),'omitmissing');
signal=sm-base;
noise=mad(diff(signal),1)/(sqrt(2)*.6745);
threshold=max(1,4*noise);
active=signal>=threshold;d=diff([false;active(:);false]);s=find(d==1);e=find(d==-1)-1;if isempty(s),peaks=[];amps=[];absolutePeak=[];return,end
gap=round(.5*fs);ms=s(1);me=[];ce=e(1);for i=2:numel(s),if s(i)-ce-1<=gap,ce=e(i);else,me(end+1,1)=ce;ms(end+1,1)=s(i);ce=e(i);end,end %#ok<AGROW>
me(end+1,1)=ce;keep=(me-ms+1)>=max(1,round(.25*fs));ms=ms(keep);me=me(keep);
sigForPeaks=signal;sigForPeaks(~isfinite(sigForPeaks))=0;
[~,candidateLoc]=findpeaks(sigForPeaks,'MinPeakHeight',threshold,'MinPeakProminence',threshold);
peaks=[];amps=[];
for i=1:numel(ms)
    inEvent=candidateLoc(candidateLoc>=ms(i) & candidateLoc<=me(i));
    if isempty(inEvent),continue,end
    [a,j]=max(signal(inEvent));peaks(end+1,1)=inEvent(j);amps(end+1,1)=a; %#ok<AGROW>
end
absolutePeak=sm(peaks);
end

function s=short(name),s=extractBefore(name+"_","_ves041");end
