%% VES041 timing relative to bar jumps and optomotor motion
% Uses the agreed peak detector: 3 noise-sigma prominence, >=2 s distance,
% >=0.25 s width, and no upper dF/F cutoff. Outputs are descriptive.

clearvars
baseFolder='/Users/sophiarenauld/stacks';
outFolder=fullfile(baseFolder,'ves041_behavior_quantification','stimulus_timing');
if ~isfolder(outFolder),mkdir(outFolder);end

barFolders=["20250812-20_ves041_fly3_clbarjumps"; ...
    "20250812-22_ves041_fly3_bbarjumps2"];
optoFolders=["20250812-21_ves041_fly3_olgrating"; ...
    "20250813-3_ves041_fly1_olgrating"; ...
    "20250819-10_ves041_fly1_olgrating5"];

%% Bar-jump event details
barRows=cell(0,15);
for fi=1:numel(barFolders)
    folder=fullfile(baseFolder,barFolders(fi));
    S=load(fullfile(folder,'smr_analysis','fly_data.mat'),'fly_data'); x=S.fly_data;
    [t,dff,fwd,rot,heading,peaks,pks]=signals_and_peaks(x);
    fs=1/median(diff(t)); w=round(fs);
    jumps=find(logical(x.jump_detected(:)));
    jumps=jumps(jumps<=numel(t));
    for pi=1:numel(peaks)
        [distance,j]=min(abs(jumps-peaks(pi)));
        lag=(peaks(pi)-jumps(j))/fs;
        if abs(lag)>2,continue,end
        i=peaks(pi);
        pre=i-w:i-1; post=i:i+w; centered=i-round(.5*fs):i+round(.5*fs);
        preFwd=mean(fwd(pre),'omitnan'); postFwd=mean(fwd(post),'omitnan');
        changeFwd=postFwd-preFwd;
        if changeFwd<=-1,kinClass="slowdown";
        elseif changeFwd>=1,kinClass="speedup";
        else,kinClass="stable forward speed";
        end
        postRot=mean(abs(rot(post)),'omitnan');
        peakRot=max(abs(rot(centered)),[],'omitnan');
        headingChange=abs(rad2deg(angle(exp(1i*(heading(post(end))-heading(post(1)))))));
        barRows(end+1,:)={barFolders(fi),j,t(jumps(j)),pi,t(i),lag, ... %#ok<SAGROW>
            lag<0,pks(pi),preFwd,postFwd,changeFwd,kinClass,postRot,peakRot,headingChange};
    end
end
barEvents=cell2table(barRows,'VariableNames',{'trial','jumpNumber','jumpTime', ...
    'peakNumber','peakTime','peakLagFromJump','peakBeforeJump','peakDff', ...
    'preFwd','postFwd','changeFwd','forwardClass','postAbsRotSpeed', ...
    'peakAbsRotSpeed','postNetHeadingDeg'});
barEvents.trial=string(barEvents.trial);barEvents.forwardClass=string(barEvents.forwardClass);
writetable(barEvents,fullfile(outFolder,'bar_jump_nearby_peak_events.csv'));

%% Optomotor moving versus stationary stimulus
optoRows=cell(0,15);
for fi=1:numel(optoFolders)
    folder=fullfile(baseFolder,optoFolders(fi));
    S=load(fullfile(folder,'smr_analysis','fly_data.mat'),'fly_data'); x=S.fly_data;
    [t,dff,fwd,rot,heading,peaks,pks,prom]=signals_and_peaks(x);
    fs=1/median(diff(t)); w=round(fs);

    daqFile=dir(fullfile(folder,'*daqData*_trial_001.mat'));
    D=load(fullfile(daqFile(1).folder,daqFile(1).name),'trialData');
    daqTime=seconds(D.trialData.Time);
    panelVoltage=interp1(daqTime,double(D.trialData.g4panels),t,'linear');
    panelVoltage=movmedian(panelVoltage,max(3,round(.1*fs)),'omitmissing');
    wrappedStep=mod(diff(panelVoltage)+5,10)-5;
    panelSpeed=[0;abs(wrappedStep./diff(t))];
    stimulusMoving=panelSpeed>1; % clear valley between ~0 and ~4-6 V/s modes

    for state=[false true]
        label="stationary"; if state,label="moving";end
        statePeaks=peaks(stimulusMoving(peaks)==state);
        statePks=pks(stimulusMoving(peaks)==state);
        stateProm=prom(stimulusMoving(peaks)==state);
        postFwd=nan(numel(statePeaks),1); postRot=postFwd; postHeading=postFwd;
        for pi=1:numel(statePeaks)
            i=statePeaks(pi); post=i:min(i+w,numel(t));
            postFwd(pi)=mean(fwd(post),'omitnan');
            postRot(pi)=mean(abs(rot(post)),'omitnan');
            postHeading(pi)=abs(rad2deg(angle(exp(1i*(heading(post(end))-heading(post(1)))))));
        end
        secondsInState=sum(stimulusMoving==state)/fs;
        optoRows(end+1,:)={optoFolders(fi),label,secondsInState, ... %#ok<SAGROW>
            secondsInState/(t(end)-t(1)),numel(statePeaks), ...
            60*numel(statePeaks)/secondsInState,median(dff(stimulusMoving==state),'omitnan'), ...
            median(statePks,'omitnan'),median(stateProm,'omitnan'), ...
            mean(postFwd,'omitnan'),mean(postRot,'omitnan'), ...
            mean(postHeading,'omitnan'),mean(panelSpeed(stimulusMoving==state),'omitnan'), ...
            min(panelSpeed(stimulusMoving==state),[],'omitnan'), ...
            max(panelSpeed(stimulusMoving==state),[],'omitnan')};
    end
end
optoSummary=cell2table(optoRows,'VariableNames',{'trial','stimulusState', ...
    'secondsInState','fractionOfTrial','nPeaks','peakRatePerMin','medianDff', ...
    'medianPeakDff','medianProminence','meanPostFwd','meanPostAbsRotSpeed', ...
    'meanPostNetHeadingDeg','meanPanelSpeed','minPanelSpeed','maxPanelSpeed'});
optoSummary.trial=string(optoSummary.trial);optoSummary.stimulusState=string(optoSummary.stimulusState);
writetable(optoSummary,fullfile(outFolder,'optomotor_moving_vs_stationary.csv'));

disp(barEvents); disp(optoSummary)

function [t,dff,fwd,rot,heading,peaks,pks,prom]=signals_and_peaks(x)
t=x.dq(1).t(:);dff=x.roi.ts(:);fwd=x.dq(1).bvf(:);
rot=x.dq(1).bvy_deg(:);heading=x.dq(1).bh(:);
n=min([numel(t),numel(dff),numel(fwd),numel(rot),numel(heading)]);
t=t(1:n);dff=dff(1:n);fwd=fwd(1:n);rot=rot(1:n);heading=heading(1:n);
fs=1/median(diff(t)); baseline=movmedian(dff,max(3,round(30*fs)),'omitmissing');
noiseSigma=1.4826*mad(dff-baseline,1);
[pks,peaks,~,prom]=findpeaks(dff,'MinPeakProminence',3*noiseSigma, ...
    'MinPeakDistance',round(2*fs),'MinPeakWidth',round(0.25*fs));
edge=round(2*fs); keep=peaks>edge & peaks+edge<=n;
peaks=peaks(keep);pks=pks(keep);prom=prom(keep);
end
