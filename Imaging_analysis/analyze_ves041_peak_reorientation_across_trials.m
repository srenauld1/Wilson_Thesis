%% Test VES041 low-speed reorientation enrichment across good trials
% Uses the same definitions as analyze_ves041_peak_reorientation_0812_22.m
% for standard fly_data trials. Each trial receives its own circular-shift
% null distribution, so conclusions can be evaluated across flies.

clearvars
rng(544);
baseFolder='/Users/sophiarenauld/stacks';
outputFolder=fullfile(baseFolder,'ves041_behavior_quantification','peak_reorientation');
if ~isfolder(outputFolder),mkdir(outputFolder);end

entries=dir(fullfile(baseFolder,'*ves041*')); entries=entries([entries.isdir]);
rows=cell(0,17); nPerm=2000;
for q=1:numel(entries)
    folder=fullfile(entries(q).folder,entries(q).name);
    [status,tags]=system(sprintf('/usr/bin/mdls -raw -name kMDItemUserTags "%s"',folder));
    if status~=0 || ~contains(tags,'good behavior'),continue,end
    flyFile=fullfile(folder,'smr_analysis','fly_data.mat');
    if ~isfile(flyFile),continue,end
    S=load(flyFile,'fly_data'); x=S.fly_data;
    t=x.dq(1).t(:); z=x.roi.ts(:); f=x.dq(1).bvf(:); y=x.dq(1).bvy_deg(:);
    n=min([numel(t),numel(z),numel(f),numel(y)]);
    t=t(1:n);z=z(1:n);f=f(1:n);y=y(1:n);
    fs=1/(max(t)/length(t)); edgeN=round(2*fs); eventN=round(.5*fs);
    % SNR-adjusted peak detector. Estimate trial-specific noise from the
    % residual after removing a 30 s moving-median baseline. Peaks are
    % retained based on prominence relative to this noise, not absolute
    % dF/F amplitude. The raw dF/F trace is not clipped.
    baselineN=max(3,round(30*fs));
    baseline=movmedian(z,baselineN,'omitmissing');
    residual=z-baseline;
    noiseSigma=1.4826*mad(residual,1);
    prominenceThreshold=3*noiseSigma;
    [~,peaks]=findpeaks(z,'MinPeakProminence',prominenceThreshold, ...
        'MinPeakDistance',round(2*fs),'MinPeakWidth',round(0.25*fs));
    peaks=peaks(peaks>edgeN & peaks+edgeN<=n);
    if isempty(peaks),continue,end
    obs=metrics(peaks,t,f,y,eventN);
    allowed=n-2*edgeN; null=nan(nPerm,6);
    for p=1:nPerm
        shift=randi([round(5*fs),allowed-round(5*fs)]);
        shifted=edgeN+mod((peaks-edgeN-1)+shift,allowed)+1;
        null(p,:)=metrics(shifted,t,f,y,eventN);
    end
    pHigh=(1+sum(null>=obs,1))/(nPerm+1);
    pLow=(1+sum(null<=obs,1))/(nPerm+1);
    [~,name]=fileparts(folder);
    rows(end+1,:)={name,infer_context(lower(name)),numel(peaks), ...
        noiseSigma,prominenceThreshold, ...
        obs(1),mean(null(:,1)),pHigh(1),obs(2),mean(null(:,2)),pLow(2), ...
        obs(3),mean(null(:,3)),pHigh(3),obs(4),obs(5),obs(6)}; %#ok<SAGROW>
end

results=cell2table(rows,'VariableNames',{'trial','context','nPeaks', ...
    'noiseSigma','prominenceThreshold', ...
    'fractionLowSpeedPivot','nullLowSpeedPivot','pLowSpeedPivotEnriched', ...
    'medianSpeed','nullMedianSpeed','pSpeedLower','meanPeakYaw','nullPeakYaw', ...
    'pYawHigher','fractionSlowdown','fractionSpeedup','fractionStable'});
results.trial=string(results.trial);results.context=string(results.context);
writetable(results,fullfile(outputFolder,'peak_reorientation_across_trials_snr_adjusted.csv'));
disp(results);

function out=metrics(idx,t,f,y,w)
pre=arrayfun(@(i)median(f(i-w:i-1),'omitnan'),idx);
post=arrayfun(@(i)median(f(i:i+w),'omitnan'),idx);
change=post-pre;
speed=arrayfun(@(i)median(f(i-w:i+w),'omitnan'),idx);
yaw=arrayfun(@(i)max(abs(y(i-w:i+w)),[],'omitnan'),idx);
joint=speed<=2 & yaw>=50;
out=[mean(joint),median(speed),mean(yaw),mean(change<=-1), ...
    mean(change>=1),mean(abs(change)<1)];
end

function c=infer_context(name)
if contains(name,'odor'),c="odor";
elseif contains(name,'grating'),c="grating";
elseif contains(name,'jump'),c="bar jump";
elseif contains(name,'dark'),c="dark";
else,c="closed loop";
end
end
