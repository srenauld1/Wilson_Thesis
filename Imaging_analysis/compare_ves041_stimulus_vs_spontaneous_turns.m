function [trialSummary,matchedEvents] = compare_ves041_stimulus_vs_spontaneous_turns()
%COMPARE_VES041_STIMULUS_VS_SPONTANEOUS_TURNS Match turns across context.
% This analysis asks whether VES041 is preferentially associated with turns
% in a sensory-stimulus context after controlling for turn magnitude and
% forward speed. Only VES041 trials previously classified as good behavior
% are included.
%
% Optomotor stimulus turns:
%   Direction-matched turns during a moving grating (+ pattern velocity is
%   matched to clockwise/positive yaw; - pattern velocity to negative yaw).
% Optomotor internally generated controls:
%   Direction-incongruent turns during the same moving-grating periods.
%   These experience the same visual context but oppose, rather than follow,
%   the optomotor drive and are therefore the cleanest within-trial control.
% Heat stimulus turns:
%   Turns while the heat signal is on.
% Heat spontaneous controls:
%   Turns while heat is off.
%
% Reorientations are peaks in absolute rotational velocity >=50 deg/s,
% separated by >=0.5 s, with mean forward speed <=2 mm/s and >=30 degrees
% net heading displacement from the turn through +1 s. Each stimulus event
% is greedily matched without replacement to a control event from the SAME
% TRIAL using peak rotational speed, forward speed, and heading displacement.
% A match must be within one robust within-trial scale on every variable.
% VES041 peaks use the agreed
% trial-adaptive detector: prominence >=3 noise sigma, distance >=2 s, and
% width >=0.25 s. Peak association is reported for both -2 to +0.25 s
% (activity before/concurrent with behavior) and +/-1 s.

baseFolder='/Users/sophiarenauld/stacks';
outFolder=fullfile(baseFolder,'ves041_behavior_quantification', ...
    'stimulus_vs_spontaneous_matched_turns');
if ~isfolder(outFolder),mkdir(outFolder);end

optoTrials=["20250812-21_ves041_fly3_olgrating"; ...
    "20250813-3_ves041_fly1_olgrating"; ...
    "20250819-10_ves041_fly1_olgrating5"];
heatTrials=["20260427-2_ves041_fly2_heat"; ...
    "20260427-3_ves041_fly2_heat_20mmbarrier"; ...
    "20260427-4_ves041_fly2_heat_20mmbarrier_2"; ...
    "20260427-6_ves041_fly3_heat_80mm"; ...
    "20260427-7_ves041_fly3_heat_20mm"];

allRows=cell(0,16);summaryRows=cell(0,10);pairId=0;
for context=["optomotor","heat"]
    if context=="optomotor",trials=optoTrials;else,trials=heatTrials;end
    for ti=1:numel(trials)
        trial=trials(ti);
        if context=="optomotor"
            x=load(fullfile(baseFolder,trial,'smr_analysis','fly_data.mat'));
            x=x.fly_data; t=x.dq(2).t(:); dff=x.roi(1).ts(:);
            fwd=x.dq(1).bvf(:); rot=x.dq(1).bvy_deg(:);heading=x.dq(1).bh(:);
            pattern=get_pattern_velocity(x,numel(t));
            stimulus=abs(pattern)>.02;
        else
            [t,dff,fwd,rot,heading,stimulus]=load_heat_trial(fullfile(baseFolder,trial));
            pattern=nan(size(t));
        end
        n=min([numel(t),numel(dff),numel(fwd),numel(rot),numel(heading),numel(stimulus)]);
        t=t(1:n);dff=dff(1:n);fwd=fwd(1:n);rot=rot(1:n);
        heading=heading(1:n);stimulus=logical(stimulus(1:n));pattern=pattern(1:n);
        fs=1/median(diff(t),'omitnan'); edge=round(2*fs);
        [~,turns]=findpeaks(abs(rot),'MinPeakHeight',50, ...
            'MinPeakDistance',max(1,round(0.5*fs)));
        turns=turns(turns>edge & turns+edge<=n);
        turnRot=abs(rot(turns));
        turnFwd=arrayfun(@(i)mean(fwd(i:min(n,i+round(fs))),'omitnan'),turns);
        turnHeading=arrayfun(@(i)abs(rad2deg(atan2(sin(heading(min(n,i+round(fs)))-heading(i)), ...
            cos(heading(min(n,i+round(fs)))-heading(i))))),turns);
        isReorientation=turnFwd<=2 & turnHeading>=30;
        if context=="optomotor"
            % Only direction-matched turns count as stimulus-driven.
            stimEvent=isReorientation & stimulus(turns) & sign(rot(turns))==sign(pattern(turns));
            ctrlEvent=isReorientation & stimulus(turns) & sign(rot(turns))~=sign(pattern(turns));
        else
            stimEvent=isReorientation & stimulus(turns); ctrlEvent=isReorientation & ~stimulus(turns);
        end
        [~,vesPeaks]=ves_peaks(dff,fs);
        [pairs,costs]=match_events(turnRot,turnFwd,turnHeading,stimEvent,ctrlEvent);
        for pi=1:size(pairs,1)
            pairId=pairId+1;
            si=pairs(pi,1);ci=pairs(pi,2);
            for k=1:2
                if k==1,ei=si;label="stimulus";else,ei=ci;label="spontaneous";end
                idx=turns(ei); near=any(abs(vesPeaks-idx)<=round(fs));
                before=any(vesPeaks>=idx-round(2*fs) & vesPeaks<=idx+round(.25*fs));
                amp=mean(dff(idx:min(n,idx+round(fs))),'omitnan');
                allRows(end+1,:)={pairId,trial,context,label,t(idx), ... %#ok<AGROW>
                    rot(idx),turnRot(ei),turnFwd(ei),turnHeading(ei),near,before,amp, ...
                    costs(pi),idx,stimulus(idx),pattern(idx)};
            end
        end
        sRows=allRows(strcmp(string(allRows(:,2)),trial),:);
        stim=strcmp(string(sRows(:,4)),"stimulus");ctrl=~stim;
        summaryRows(end+1,:)={trial,context,size(pairs,1),sum(stimEvent), ... %#ok<AGROW>
            sum(ctrlEvent),mean(cell2mat(sRows(stim,10)),'omitnan'), ...
            mean(cell2mat(sRows(ctrl,10)),'omitnan'), ...
            mean(cell2mat(sRows(stim,11)),'omitnan'), ...
            mean(cell2mat(sRows(ctrl,11)),'omitnan'), ...
            median(cell2mat(sRows(:,13)),'omitnan')};
    end
end

matchedEvents=cell2table(allRows,'VariableNames',{'pairId','trial','context', ...
    'turnClass','turnTime','signedRotPeak','absRotPeak','meanFwd0to1s','netHeading0to1Deg', ...
    'vesPeakWithin1s','vesPeakBeforeOrConcurrent','meanDff0to1s','matchCost', ...
    'turnIndex','stimulusOn','patternVelocity'});
matchedEvents.trial=string(matchedEvents.trial);
matchedEvents.context=string(matchedEvents.context);
matchedEvents.turnClass=string(matchedEvents.turnClass);
trialSummary=cell2table(summaryRows,'VariableNames',{'trial','context','nMatchedPairs', ...
    'nCandidateStimulusTurns','nCandidateSpontaneousTurns','stimPeakWithin1Fraction', ...
    'spontPeakWithin1Fraction','stimPeakBeforeConcurrentFraction', ...
    'spontPeakBeforeConcurrentFraction','medianMatchCost'});
trialSummary.trial=string(trialSummary.trial);
trialSummary.context=string(trialSummary.context);
writetable(matchedEvents,fullfile(outFolder,'matched_turn_events.csv'));
writetable(trialSummary,fullfile(outFolder,'matched_turn_trial_summary.csv'));
plot_results(matchedEvents,outFolder);
disp(trialSummary);
end

function pattern=get_pattern_velocity(x,n)
if isfield(x.dq(1),'vvyo'),pattern=x.dq(1).vvyo(:);else,pattern=x.dq(1).vvy(:);end
if numel(pattern)~=n
    pattern=interp1(x.dq(1).t(:),pattern,x.dq(2).t(:),'nearest','extrap');
end
end

function [t,dff,fwd,rot,heading,heatOn]=load_heat_trial(folder)
f=dir(fullfile(folder,'*_virmenData.mat'));v=load(fullfile(f(1).folder,f(1).name));
f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));q=load(fullfile(f(1).folder,f(1).name));
t=seconds(v.vmData.Seconds{1,1});t=t(:);fwd=smooth(v.vmData.vel_for{1,1});
rot=rad2deg(smooth(v.vmData.vel_yaw{1,1}));heading=v.vmData.heading{1,1};heading=heading(:);heat=v.vmData.heat;
if iscell(heat),heat=cell2mat(heat);end
heatOn=heat(:)>1e-4;
td=v.daqVolumeTime.daqVolumeTime{1,1};
dff=interp1(td(:),q.roi(1).ts(:),t,'linear',NaN);
end

function [pks,locs]=ves_peaks(dff,fs)
noise=mad(diff(dff),1)/(sqrt(2)*0.6745);
[pks,locs]=findpeaks(dff,'MinPeakProminence',3*noise, ...
    'MinPeakDistance',max(1,round(2*fs)),'MinPeakWidth',max(1,round(.25*fs)));
end

function [pairs,costs]=match_events(r,f,h,isStim,isCtrl)
s=find(isStim);c=find(isCtrl);pairs=zeros(0,2);costs=zeros(0,1);
if isempty(s)||isempty(c),return,end
sr=max(mad(r([s;c]),1),eps);sf=max(mad(f([s;c]),1),eps);sh=max(mad(h([s;c]),1),eps);
available=true(size(c));
% Match the hardest stimulus observations first to reduce greedy bias.
C=sqrt(((r(s)-r(c)')./sr).^2+((f(s)-f(c)')./sf).^2+((h(s)-h(c)')./sh).^2);
[~,order]=sort(min(C,[],2),'descend');
for ii=order(:)'
    valid=find(available & abs(r(c)-r(s(ii)))<=sr & abs(f(c)-f(s(ii)))<=sf & abs(h(c)-h(s(ii)))<=sh);
    if isempty(valid),continue,end
    [cost,j]=min(C(ii,valid));pick=valid(j);available(pick)=false;
    pairs(end+1,:)=[s(ii),c(pick)];costs(end+1,1)=cost; %#ok<AGROW>
end
end

function plot_results(E,outFolder)
fig=figure('Color','w','Position',[80 80 1450 760]);tl=tiledlayout(2,3,'TileSpacing','compact');
for ci=1:2
    context=["optomotor","heat"];z=E(E.context==context(ci),:);
    stim=z(z.turnClass=="stimulus",:);ctrl=z(z.turnClass=="spontaneous",:);
    nexttile(tl,(ci-1)*3+1);scatter(ctrl.absRotPeak,stim.absRotPeak,35,'filled');
    hold on;lims=[min([ctrl.absRotPeak;stim.absRotPeak]) max([ctrl.absRotPeak;stim.absRotPeak])];plot(lims,lims,'k--');
    axis square;grid on;xlabel('Spontaneous |rot| peak');ylabel('Stimulus |rot| peak');title(context(ci)+' matched rotation');
    nexttile(tl,(ci-1)*3+2);scatter(ctrl.meanFwd0to1s,stim.meanFwd0to1s,35,'filled');
    hold on;lims=[min([ctrl.meanFwd0to1s;stim.meanFwd0to1s]) max([ctrl.meanFwd0to1s;stim.meanFwd0to1s])];plot(lims,lims,'k--');
    axis square;grid on;xlabel('Spontaneous forward speed');ylabel('Stimulus forward speed');title(context(ci)+' matched forward speed');
    nexttile(tl,(ci-1)*3+3);vals=100*[mean(stim.vesPeakBeforeOrConcurrent),mean(ctrl.vesPeakBeforeOrConcurrent); ...
        mean(stim.vesPeakWithin1s),mean(ctrl.vesPeakWithin1s)];
    bar(vals);xticklabels({'Peak -2 to +0.25 s','Peak within +/-1 s'});ylabel('Turns with VES041 peak (%)');
    legend({'Stimulus','Spontaneous'},'Location','best');grid on;title(context(ci)+' peak probability (matched)');
end
sgtitle(tl,'VES041: stimulus-context versus kinematically matched spontaneous turns');
savefig(fig,fullfile(outFolder,'matched_turn_summary.fig'));
exportgraphics(fig,fullfile(outFolder,'matched_turn_summary.png'),'Resolution',250);close(fig);
end
