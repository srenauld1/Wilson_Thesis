function summary=plot_ves041_detected_peaks_all_trials(minHeight,minDistanceSec,outputName,detectorMode)
% Plot every detected VES041 peak in each good-behavior imaging trial.
% Detector matches the current analyses: >=3 noise-sigma prominence,
% >=2 s peak separation, and >=0.25 s width. No upper dF/F cap is used.
% All peaks are displayed, including peaks during sensory stimulation.

if nargin<1 || isempty(minHeight),minHeight=-Inf;end
if nargin<2 || isempty(minDistanceSec),minDistanceSec=2;end
if nargin<3 || isempty(outputName),outputName='detected_peaks_all_good_trials';end
if nargin<4 || isempty(detectorMode),detectorMode='raw';end
base='/Users/sophiarenauld/stacks';
out=fullfile(base,'ves041_behavior_quantification',char(outputName));
if ~isfolder(out),mkdir(out);end
standard=["20250812-20_ves041_fly3_clbarjumps";"20250812-21_ves041_fly3_olgrating"; ...
 "20250812-22_ves041_fly3_bbarjumps2";"20250813-3_ves041_fly1_olgrating"; ...
 "20250819-10_ves041_fly1_olgrating5";"20251001-21_ves041_fly3_clbbar"; ...
 "20260213-2_ves041_clbar";"20260213-4_ves041_fly2_clbar"; ...
 "20260217-2_ves041_vt_clbbar1";"20260217-3_ves041_vt_clbbar2_actuallydark"; ...
 "20260217-4_ves041_vt_clbbar3"];
heat=["20260427-2_ves041_fly2_heat";"20260427-3_ves041_fly2_heat_20mmbarrier"; ...
 "20260427-4_ves041_fly2_heat_20mmbarrier_2";"20260427-6_ves041_fly3_heat_80mm"; ...
 "20260427-7_ves041_fly3_heat_20mm"];
rows=cell(0,6);
for trial=[standard;heat]'
    isHeat=any(heat==trial);jumpTimes=[];
    if isHeat
        [t,dff,stim]=load_heat(fullfile(base,trial));stimLabel='Heat on';
    else
        f=fullfile(base,trial,'smr_analysis','fly_data.mat');
        if ~isfile(f),warning('Missing %s',f);continue,end
        q=load(f);x=q.fly_data;t=x.dq(1).t(:);dff=x.roi(1).ts(:);stim=false(size(t));stimLabel='';
        if contains(lower(trial),'grating')
            vv=x.dq(1).vvy(:);m=min(numel(stim),numel(vv));stim(1:m)=abs(vv(1:m))>.02;stimLabel='Optomotor moving';
        end
        if isfield(x,'jump_detected')
            j=find(logical(x.jump_detected(:)));j=j(j<=numel(t));
            if isfield(x.dq(1),'vh')
                vh=x.dq(1).vh(:);dv=atan2(sin([vh(1);diff(vh)]),cos([vh(1);diff(vh)]));
                j=j(abs(rad2deg(dv(j)))>=45);
            end
            jumpTimes=t(j);
        end
    end
    n=min(numel(t),numel(dff));t=t(1:n);dff=dff(1:n);stim=stim(1:min(n,numel(stim)));
    if numel(stim)<n,stim(end+1:n)=false;end
    fs=1/median(diff(t));smoothDff=dff;baseline=nan(size(dff));
    if strcmpi(detectorMode,'local') || strcmpi(detectorMode,'merged_local')
        smoothN=max(3,round(.25*fs));smoothDff=smoothdata(dff,'gaussian',smoothN);
        baseline=movmedian(smoothDff,max(3,round(30*fs)),'omitmissing');
        detectionSignal=smoothDff-baseline;
        noise=mad(diff(detectionSignal),1)/(sqrt(2)*.6745);prom=3*noise;
        if strcmpi(detectorMode,'merged_local')
            eventThreshold=max(.5,prom);
            [locs,proms]=merged_local_events(detectionSignal,eventThreshold, ...
                round(.5*fs),max(1,round(.25*fs)));
        else
            peakArgs={'MinPeakProminence',prom,'MinPeakHeight',.5, ...
                'MinPeakDistance',max(1,round(minDistanceSec*fs)), ...
                'MinPeakWidth',max(1,round(.25*fs))};
            [~,locs,~,proms]=findpeaks(detectionSignal,peakArgs{:});
            if isfinite(minHeight),keep=smoothDff(locs)>=minHeight;locs=locs(keep);proms=proms(keep);end
        end
        pks=smoothDff(locs);
    else
        noise=mad(diff(dff),1)/(sqrt(2)*.6745);prom=3*noise;
        peakArgs={'MinPeakProminence',prom,'MinPeakDistance',max(1,round(minDistanceSec*fs)), ...
            'MinPeakWidth',max(1,round(.25*fs))};
        if isfinite(minHeight),peakArgs=[peakArgs {'MinPeakHeight',minHeight}];end
        [pks,locs,~,proms]=findpeaks(dff,peakArgs{:});
    end
    fig=figure('Visible','off','Color','w','Position',[60 100 1500 520]);ax=axes(fig);hold(ax,'on');
    add_epochs(ax,t,stim,stimLabel,dff);
    plot(ax,t,dff,'Color',[.1 .25 .75],'LineWidth',1,'DisplayName','dF/F');
    if strcmpi(detectorMode,'local') || strcmpi(detectorMode,'merged_local')
        plot(ax,t,smoothDff,'Color',[.05 .05 .05],'LineWidth',1.15,'DisplayName','Smoothed dF/F');
        plot(ax,t,baseline,'Color',[.35 .65 .35],'LineWidth',1,'DisplayName','30 s local baseline');
    end
    scatter(ax,t(locs),pks,42,[.85 .12 .1],'v','filled','DisplayName','Detected peak');
    for jt=jumpTimes(:)',h=xline(ax,jt,'g-','LineWidth',1.1);h.HandleVisibility='off';end
    xlabel(ax,'Time (s)');ylabel(ax,'dF/F');grid(ax,'on');xlim(ax,[t(1) t(end)]);
    if isfinite(minHeight),heightText=sprintf(' | height >= %.2g',minHeight);else,heightText='';end
    if strcmpi(detectorMode,'merged_local')
        modeText=sprintf(' | merged local events >= %.3g',max(.5,prom));heightText='';
    elseif strcmpi(detectorMode,'local'),modeText=' | local amplitude >= 0.5';else,modeText='';end
    title(ax,sprintf('%s | %d peaks | prominence >= %.3g (3 sigma)%s%s | distance >= %.2g s', ...
        strrep(trial,'_','\_'),numel(locs),prom,heightText,modeText,minDistanceSec));legend(ax,'Location','best');
    savefig(fig,fullfile(out,trial+'_detected_peaks.fig'));
    exportgraphics(fig,fullfile(out,trial+'_detected_peaks.png'),'Resolution',220);close(fig);
    rows(end+1,:)={trial,numel(locs),noise,prom,median(pks,'omitnan'),median(proms,'omitnan')}; %#ok<AGROW>
end
summary=cell2table(rows,'VariableNames',{'trial','nPeaks','noiseSigma','prominenceThreshold', ...
    'medianPeakDff','medianPeakProminence'});summary.trial=string(summary.trial);
writetable(summary,fullfile(out,'detected_peak_summary.csv'));disp(summary);
end

function [locs,amps]=merged_local_events(signal,threshold,maxGap,minDuration)
% Bridge brief dips below threshold, then retain one maximum per event.
active=signal>=threshold;d=diff([false;active(:);false]);s=find(d==1);e=find(d==-1)-1;
if isempty(s),locs=[];amps=[];return,end
mergedS=s(1);mergedE=[];currentE=e(1);
for i=2:numel(s)
    if s(i)-currentE-1<=maxGap,currentE=e(i);else,mergedE(end+1,1)=currentE;mergedS(end+1,1)=s(i);currentE=e(i);end %#ok<AGROW>
end
mergedE(end+1,1)=currentE;keep=(mergedE-mergedS+1)>=minDuration;
mergedS=mergedS(keep);mergedE=mergedE(keep);locs=zeros(numel(mergedS),1);amps=locs;
for i=1:numel(mergedS)
    [amps(i),j]=max(signal(mergedS(i):mergedE(i)));locs(i)=mergedS(i)+j-1;
end
end

function add_epochs(ax,t,mask,label,dff)
if isempty(label)||~any(mask),return,end
d=diff([false;mask(:);false]);s=find(d==1);e=find(d==-1)-1;
yl=[min(dff,[],'omitnan') max(dff,[],'omitnan')];
for i=1:numel(s)
    h=patch(ax,[t(s(i)) t(e(i)) t(e(i)) t(s(i))],[yl(1) yl(1) yl(2) yl(2)], ...
        [.8 .8 .8],'FaceAlpha',.22,'EdgeColor','none');
    if i==1,h.DisplayName=label;else,h.HandleVisibility='off';end
end
end

function [t,dff,stim]=load_heat(folder)
f=dir(fullfile(folder,'*_virmenData.mat'));v=load(fullfile(f(1).folder,f(1).name));
f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));q=load(fullfile(f(1).folder,f(1).name));
t=seconds(v.vmData.Seconds{1,1});t=t(:);h=v.vmData.heat;if iscell(h),h=cell2mat(h);end;stim=h(:)>1e-4;
td=v.daqVolumeTime.daqVolumeTime{1,1};dff=interp1(td(:),q.roi(1).ts(:),t,'linear',NaN);
end
