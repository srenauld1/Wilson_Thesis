function report=run_smp544_oct26_antipeaks(outputFolder)
%RUN_SMP544_OCT26_ANTIPEAKS Apply the strict local-event detector downward.
% Uses only the saved full hand-drawn ROI traces from the five Oct26 trials.
% Inputs and previous analysis figures are never overwritten.
% Gaussian smooth 0.25 s; local median baseline 30 s. Detection signal is
% baseline-smoothedDff, i.e. the sign-reversed local fluorescence residual.
% Noise sigma=MAD(diff(signal))/(sqrt(2)*0.6745); threshold=max(1,4*sigma).
% Require threshold height AND prominence; merge interruptions <=0.5 s;
% merged event span >=0.25 s; keep one deepest qualifying minimum per event.
% Duration includes merged interruptions, matching the existing strict VES
% detector. NaNs and timestamp gaps split blocks and are never bridged.
% Detection uses fluorescence only, including all stimulus periods. Behavior
% is plotted for inspection, not used to select events. No upper dF/F cap.
% All inputs retain registration-QC holds, so detections are exploratory.
% Must run MATLAB -nodisplay; saved FIGs nevertheless open visibly.
if nargin<1
    outputFolder=fullfile(imaging_analysis_artifact_root(),'diagnostics', ...
        'smp544_oct26_antipeaks_20261028');
end
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
if ~isfolder(outputFolder),mkdir(outputFolder);end
root='/Users/sophiarenauld/stacks/2026-10-26/a2p_exploratory_20261027';
files=dir(fullfile(root,'*','smp544_manual_analysis_20261028','fly_data.mat'));
assert(numel(files)==5,'Expected five processed trials.');
rows=cell(0,8);
for k=1:numel(files)
    [parent,~]=fileparts(files(k).folder);[~,trial]=fileparts(parent);
    try
        q=load(fullfile(files(k).folder,files(k).name),'fly_data');a=q.fly_data;
        assert(isequaln(a.roi(1).ts,a.roi(1).pxfilter.manualts),'Expected manual ROI');
        t=double(a.dq(1).t(:));y=double(a.roi(1).ts(:));
        assert(numel(t)==numel(y)&&all(diff(t)>0),'Invalid time/trace lengths');
        fs=1/median(diff(t));[sm,base,signal,noise,threshold,events]=detect(t,y,fs);
        E=table(repmat(string(trial),size(events,1),1),t(events(:,1)), ...
            t(events(:,2)),t(events(:,3)),y(events(:,2)),sm(events(:,2)), ...
            signal(events(:,2)),events(:,4), ...
            'VariableNames',{'trial','startSec','troughSec','endSec', ...
            'rawDff','smoothedDff','depthBelowLocalBaseline','prominence'});
        writetable(E,fullfile(outputFolder,[trial '_antipeaks.csv']));
        R=struct('t',t,'rawDff',y,'smoothedDff',sm,'localBaseline',base, ...
            'detectionSignal',signal,'noiseSigma',noise,'threshold',threshold, ...
            'events',E,'input',fullfile(files(k).folder,files(k).name), ...
            'registrationQc',a.analysisRegistrationQc);
        save(fullfile(outputFolder,[trial '_antipeaks.mat']),'R');
        draw_trial(a,R,events,trial,outputFolder);
        rows(end+1,:)={string(trial),"success",height(E),noise,threshold, ...
            median(E.depthBelowLocalBaseline,'omitnan'),"registration QC held",""}; %#ok<AGROW>
    catch ME
        rows(end+1,:)={string(trial),"error",NaN,NaN,NaN,NaN,"",string(getReport(ME,'extended','hyperlinks','off'))}; %#ok<AGROW>
    end
    report=cell2table(rows,'VariableNames',{'trial','status','nAntipeaks', ...
        'noiseSigma','thresholdDff','medianDepthDff','qcNote','error'});
    writetable(report,fullfile(outputFolder,'antipeak_report.csv'));
end
disp(report);assert(~any(report.status=="error"),'Some trials failed; see report.');
end

function [sm,base,sig,noise,threshold,events]=detect(t,y,fs)
valid=isfinite(y);breaks=[true;diff(t)>1.5/fs];
starts=find(valid&([true;~valid(1:end-1)]|breaks));
ends=find(valid&([~valid(2:end);true]|[breaks(2:end);true]));
assert(numel(starts)==numel(ends),'Invalid contiguous-block construction');
sm=nan(size(y));base=sm;sig=sm;differences=[];
for b=1:numel(starts)
    ix=starts(b):ends(b);
    sm(ix)=smoothdata(y(ix),'gaussian',max(3,round(.25*fs)));
    base(ix)=movmedian(sm(ix),max(3,round(30*fs)),'omitmissing');
    sig(ix)=base(ix)-sm(ix);differences=[differences;diff(sig(ix))]; %#ok<AGROW>
end
noise=mad(differences,1)/(sqrt(2)*.6745);threshold=max(1,4*noise);
events=zeros(0,4);
for b=1:numel(starts)
    ix=starts(b):ends(b);x=sig(ix);if numel(x)<3,continue,end
    d=diff([false;x>=threshold;false]);s=find(d==1);e=find(d==-1)-1;
    if isempty(s),continue,end
    ms=s(1);me=[];ce=e(1);
    for j=2:numel(s)
        if s(j)-ce-1<=round(.5*fs),ce=e(j);
        else,me(end+1,1)=ce;ms(end+1,1)=s(j);ce=e(j);end %#ok<AGROW>
    end
    me(end+1,1)=ce;
    [~,loc,~,prom]=findpeaks(x,'MinPeakHeight',threshold,'MinPeakProminence',threshold);
    for j=1:numel(ms)
        if (me(j)-ms(j)+1)/fs<.25,continue,end
        candidates=find(loc>=ms(j)&loc<=me(j));if isempty(candidates),continue,end
        [~,best]=max(x(loc(candidates)));c=candidates(best);
        events(end+1,:)=[ix(ms(j)),ix(loc(c)),ix(me(j)),prom(c)]; %#ok<AGROW>
    end
end
end

function draw_trial(a,R,events,trial,out)
t=R.t;fig=figure('Visible','off','Color','w','Position',[20 20 1600 1000]);
tl=tiledlayout(fig,5,1,'TileSpacing','compact','Padding','compact');axs=gobjects(5,1);
for v=1:5,axs(v)=nexttile(tl);hold(axs(v),'on');grid(axs(v),'on');end
plot(axs(1),t,R.rawDff,'Color',[.7 .78 .9],'DisplayName','Raw manual ROI');
plot(axs(1),t,R.smoothedDff,'Color',[.05 .25 .7],'DisplayName','0.25 s smoothed');
plot(axs(1),t,R.localBaseline,'Color',[.2 .6 .25],'LineWidth',1.2,'DisplayName','30 s local median');
scatter(axs(1),t(events(:,2)),R.smoothedDff(events(:,2)),45,[.85 .1 .15],'v','filled','DisplayName','Detected trough');
ylabel(axs(1),'dF/F');legend(axs(1),'Location','northeast','NumColumns',2);
plot(axs(2),t,R.detectionSignal,'Color',[.05 .25 .7]);
yline(axs(2),R.threshold,'--r',sprintf('threshold %.2f',R.threshold));
scatter(axs(2),t(events(:,2)),R.detectionSignal(events(:,2)),30,[.85 .1 .15],'v','filled');
ylabel(axs(2),{'Local trough','depth (dF/F)'});
f=double(a.dq(1).bvf(:));r=double(a.dq(1).bvy_deg(:));
if isfield(a.dq(1),'bergii')&&isfield(a.dq(1).bergii,'behaviorVolumeValid')
    valid=logical(a.dq(1).bergii.behaviorVolumeValid(:));f(~valid)=NaN;r(~valid)=NaN;
end
plot(axs(3),t,f,'k');ylabel(axs(3),{'Forward','velocity (mm/s)'});
plot(axs(4),t,abs(r),'Color',[.85 .15 .1]);ylabel(axs(4),{'Rotational','speed (deg/s)'});
isBarJump=isfield(a,'protocol')&&string(a.protocol.kind)=="bar-jump";
if ~isBarJump&&isfield(a.dq(1),'heat')&&~isempty(a.dq(1).heat)
    plot(axs(5),t,double(a.dq(1).heat(:)),'Color',[.85 .4 .05]);
    ylabel(axs(5),'Heat monitor (V)');
else
    if isfield(a,'jump_detected')
        jump=find(logical(a.jump_detected(:)));
        for j=jump(:)',xline(axs(5),t(j),'-g');end
    end
    ylim(axs(5),[0 1]);ylabel(axs(5),'Bar jumps');
end
for v=1:5
    for j=1:size(events,1),xline(axs(v),t(events(j,2)),':','Color',[.8 .4 .45],'HandleVisibility','off');end
end
linkaxes(axs,'x');xlim(axs(1),[t(1) t(end)]);xlabel(axs(5),'Time (s)');
title(tl,sprintf('%s | %d antipeaks | depth/prominence >= %.2f dF/F',trial,size(events,1),R.threshold),'Interpreter','none');
subtitle(tl,'Exploratory: registration QC held | full manual ROI | all stimulus contexts included');
exportgraphics(fig,fullfile(out,[trial '_antipeaks.png']),'Resolution',150);
savefig_visible_on_open(fig,fullfile(out,[trial '_antipeaks.fig']));close(fig);
end
