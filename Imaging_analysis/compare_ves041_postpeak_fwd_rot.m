%% Compare locomotion before versus after SNR-adjusted VES041 peaks
% The baseline is -2 to -1 s before each peak. The interval from -1 to 0 s
% is intentionally left unused. The test window is 0 to +1 s after peak.
% This script reports descriptive changes only; it performs no hypothesis tests.

clearvars
baseFolder='/Users/sophiarenauld/stacks';
outputFolder=fullfile(baseFolder,'ves041_behavior_quantification','peak_navigation_metrics');
if ~isfolder(outputFolder),mkdir(outputFolder);end

entries=dir(fullfile(baseFolder,'*ves041*')); entries=entries([entries.isdir]);
rows=cell(0,10);
for q=1:numel(entries)
    folder=fullfile(entries(q).folder,entries(q).name);
    [status,tags]=system(sprintf('/usr/bin/mdls -raw -name kMDItemUserTags "%s"',folder));
    if status~=0 || ~contains(tags,'good behavior'),continue,end
    flyFile=fullfile(folder,'smr_analysis','fly_data.mat');
    if ~isfile(flyFile),continue,end

    S=load(flyFile,'fly_data'); x=S.fly_data;
    t=x.dq(1).t(:); dff=x.roi.ts(:);
    fwd=x.dq(1).bvf(:); rot=x.dq(1).bvy_deg(:);
    n=min([numel(t),numel(dff),numel(fwd),numel(rot)]);
    t=t(1:n); dff=dff(1:n); fwd=fwd(1:n); rot=rot(1:n);
    fs=1/(max(t)/length(t)); edgeN=round(2*fs);

    baseline=movmedian(dff,max(3,round(30*fs)),'omitmissing');
    noiseSigma=1.4826*mad(dff-baseline,1);
    [~,peaks]=findpeaks(dff,'MinPeakProminence',3*noiseSigma, ...
        'MinPeakDistance',round(2*fs),'MinPeakWidth',round(0.25*fs));
    peaks=peaks(peaks>edgeN & peaks+round(fs)<=n);
    if isempty(peaks),continue,end

    preFwd=arrayfun(@(i)mean(fwd(i-round(2*fs):i-round(fs)),'omitnan'),peaks);
    postFwd=arrayfun(@(i)mean(fwd(i:i+round(fs)),'omitnan'),peaks);
    preRot=arrayfun(@(i)mean(abs(rot(i-round(2*fs):i-round(fs))),'omitnan'),peaks);
    postRot=arrayfun(@(i)mean(abs(rot(i:i+round(fs))),'omitnan'),peaks);

    [~,name]=fileparts(folder);
    rows(end+1,:)={name,infer_context(lower(name)),numel(peaks), ... %#ok<SAGROW>
        mean(preFwd,'omitnan'),mean(postFwd,'omitnan'), ...
        mean(postFwd-preFwd,'omitnan'), ...
        mean(preRot,'omitnan'),mean(postRot,'omitnan'), ...
        mean(postRot-preRot,'omitnan'),noiseSigma};
end

results=cell2table(rows,'VariableNames',{'trial','context','nPeaks', ...
    'preFwd','postFwd','changeFwd','preAbsRot','postAbsRot','changeAbsRot','noiseSigma'});
results.trial=string(results.trial); results.context=string(results.context);
writetable(results,fullfile(outputFolder,'postpeak_fwd_rot_descriptive.csv'));
disp(results)

function c=infer_context(name)
if contains(name,'odor'),c="odor";
elseif contains(name,'grating'),c="grating";
elseif contains(name,'jump'),c="bar jump";
elseif contains(name,'dark'),c="dark";
else,c="closed loop";
end
end
