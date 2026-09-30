%% Describe VES041 peaks across visual-stimulus contexts
% Uses the agreed SNR-adjusted peak detector. Optomotor is treated as a
% trial-wide context. Bar-jump trials are also split into stored +/-2 s
% jump windows versus the remainder of that same trial. No hypothesis tests
% are performed; output is descriptive.

clearvars
baseFolder='/Users/sophiarenauld/stacks';
outputFolder=fullfile(baseFolder,'ves041_behavior_quantification','visual_contexts');
if ~isfolder(outputFolder),mkdir(outputFolder);end

pseudoFile=fullfile(baseFolder,'ves041_behavior_quantification', ...
    'peak_navigation_metrics','postpeak_vs_pseudo_1s_descriptive.csv');
pseudo=readtable(pseudoFile,'TextType','string');
entries=dir(fullfile(baseFolder,'*ves041*')); entries=entries([entries.isdir]);
trialRows=cell(0,11); jumpRows=cell(0,9);

for q=1:numel(entries)
    folder=fullfile(entries(q).folder,entries(q).name);
    [status,tags]=system(sprintf('/usr/bin/mdls -raw -name kMDItemUserTags "%s"',folder));
    if status~=0 || ~contains(tags,'good behavior'),continue,end
    name=string(entries(q).name); lname=lower(name);
    if contains(lname,'heat') || contains(lname,'odor'),continue,end
    flyFile=fullfile(folder,'smr_analysis','fly_data.mat');
    if ~isfile(flyFile),continue,end

    S=load(flyFile,'fly_data'); x=S.fly_data;
    t=x.dq(1).t(:); dff=x.roi.ts(:);
    n=min(numel(t),numel(dff)); t=t(1:n); dff=dff(1:n);
    fs=1/(max(t)/length(t));
    baseline=movmedian(dff,max(3,round(30*fs)),'omitmissing');
    noiseSigma=1.4826*mad(dff-baseline,1);
    [pks,peaks,~,prominences]=findpeaks(dff, ...
        'MinPeakProminence',3*noiseSigma, ...
        'MinPeakDistance',round(2*fs),'MinPeakWidth',round(0.25*fs));

    if contains(lname,'grating'),context="optomotor";
    elseif contains(lname,'jump'),context="bar jump";
    elseif contains(lname,'dark'),context="dark";
    else,context="closed loop";
    end

    k=find(pseudo.trial==name,1);
    if isempty(k),ef=NaN; es=NaN; eh=NaN;
    else
        ef=pseudo.effectFwd(k); es=pseudo.effectStraightness(k);
        eh=pseudo.effectNetHeadingDeg(k);
    end
    trialRows(end+1,:)={name,context,numel(peaks), ... %#ok<SAGROW>
        60*numel(peaks)/(t(end)-t(1)),median(pks,'omitnan'), ...
        median(prominences,'omitnan'),ef,es,eh,noiseSigma,3*noiseSigma};

    if context=="bar jump" && isfield(x,'omit_jump_blocks')
        mask=logical(x.omit_jump_blocks(:)); mask=mask(1:min(n,numel(mask)));
        peaks=peaks(peaks<=numel(mask)); pks=pks(1:numel(peaks));
        in=mask(peaks); secIn=sum(mask)/fs; secOut=(numel(mask)-sum(mask))/fs;
        rateIn=sum(in)/secIn*60; rateOut=sum(~in)/secOut*60;
        jumpRows(end+1,:)={name,sum(in),sum(~in),secIn,secOut,rateIn, ... %#ok<SAGROW>
            rateOut,median(pks(in),'omitnan'),median(pks(~in),'omitnan')};
    end
end

trialTable=cell2table(trialRows,'VariableNames',{'trial','context','nPeaks', ...
    'peakRatePerMin','medianPeakDff','medianProminence','postEffectFwd', ...
    'postEffectStraightness','postEffectNetHeadingDeg','noiseSigma','threshold'});
trialTable.trial=string(trialTable.trial); trialTable.context=string(trialTable.context);
writetable(trialTable,fullfile(outputFolder,'visual_context_trials.csv'));

contexts=unique(trialTable.context,'stable'); groupRows=cell(0,9);
for i=1:numel(contexts)
    z=trialTable(trialTable.context==contexts(i),:);
    groupRows(end+1,:)={contexts(i),height(z),sum(z.nPeaks), ... %#ok<SAGROW>
        mean(z.peakRatePerMin,'omitnan'),mean(z.medianPeakDff,'omitnan'), ...
        mean(z.postEffectFwd,'omitnan'),mean(z.postEffectStraightness,'omitnan'), ...
        mean(z.postEffectNetHeadingDeg,'omitnan'),strjoin(short_trial(z.trial),', ')};
end
groupTable=cell2table(groupRows,'VariableNames',{'context','nTrials','nPeaks', ...
    'meanTrialPeakRatePerMin','meanTrialMedianPeakDff','meanPostEffectFwd', ...
    'meanPostEffectStraightness','meanPostEffectNetHeadingDeg','trials'});
groupTable.context=string(groupTable.context);groupTable.trials=string(groupTable.trials);
writetable(groupTable,fullfile(outputFolder,'visual_context_summary.csv'));

jumpTable=cell2table(jumpRows,'VariableNames',{'trial','nPeaksNearJump', ...
    'nPeaksAwayFromJump','secondsNearJump','secondsAwayFromJump', ...
    'peakRateNearJumpPerMin','peakRateAwayFromJumpPerMin', ...
    'medianPeakDffNearJump','medianPeakDffAwayFromJump'});
jumpTable.trial=string(jumpTable.trial);
writetable(jumpTable,fullfile(outputFolder,'bar_jump_peak_timing.csv'));
disp(groupTable); disp(jumpTable)

function s=short_trial(names)
s=extractBefore(names+"_","_ves041");
end
