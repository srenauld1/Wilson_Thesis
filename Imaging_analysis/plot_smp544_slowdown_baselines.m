function results = plot_smp544_slowdown_baselines(trialFolders, outputFolder)
%PLOT_SMP544_SLOWDOWN_BASELINES Separate local falls from early-baseline elevation.
% Events are selected ONLY from forward velocity, never from fluorescence.
% trialFolders must be an explicit list of previously behavior-screened trials.
% For legacy inputs use roi.ts; for processed trials use saved fly_data.roi.ts.
% No ROI source or fluorescence baseline is changed in the input files.
%
% Fluorescence: Gaussian 0.25-s smoothing, 30-s moving-median detrending,
% separate robust z normalization per ROI, then mean across ROIs. This matches
% compare_smp544_ves041_reorientation (not the heavily smoothed trough detector).
% Slowdowns: peaks of negative acceleration, >=1.5 s apart, prominence
% max(0.5,2*MAD_sigma(acceleration)); forward-speed drop >=1 mm/s between
% [-1,-0.25] and [0,0.75] s. Time zero is maximum deceleration.
% Exclude full +/-3-s windows containing stimuli, missing data or invalid
% behavior volumes. Heat-active flags are preferred; legacy fallback is
% heat > median(heat)+2*MAD(heat), as in the earlier July comparison.
% Turn split: >=15 degrees net wrapped heading change over [-0.25,0.75] s.
%
% TWO DISTINCT measurements, with identical events and post window:
%   earlyChangeZ = mean(+0.25:+0.75) - mean(-3:-2)
%   localChangeZ = mean(+0.25:+0.75) - mean(-1:-0.25)
% Negative localChangeZ is a local decrease. Do NOT use a selected maximum or
% minimum: those would bias the measured change downward.
% Per-trial means weight events equally. Summary means weight recordings
% equally; recordings are NOT necessarily independent flies. No tests run.
% Run with MATLAB -nodisplay so saved .fig files open visibly without popups.

assert(nargin>=2,'Supply explicit trial folders and an output folder.');
assert(~feature('ShowFigureWindows'),'Run this function with MATLAB -nodisplay.');
trialFolders=string(trialFolders(:));
if ~isfolder(outputFolder),mkdir(outputFolder);end
rel=(-3:.05:3)'; D=cell(numel(trialFolders),1); eventTables=cell(size(D));
for k=1:numel(D)
    D{k}=read_trial(trialFolders(k),rel);
    eventTables{k}=D{k}.events;
    render({D{k}},rel,fullfile(outputFolder,D{k}.name),D{k}.name,false);
end
events=vertcat(eventTables{:});
writetable(events,fullfile(outputFolder,'slowdown_event_metrics.csv'));
render(D,rel,fullfile(outputFolder,'all_good_behavior_slowdowns'), ...
    'SMP544 behavior-selected slowdowns | equal recording weight',true);
rows=cell(0,9);
for k=1:numel(D)
    for g=1:3
        T=D{k}.events; keep=group_mask(T,g);T=T(keep,:);
        rows(end+1,:)={D{k}.name,string(group_name(g)),height(T), ...
            mean(T.earlyChangeZ,'omitnan'),mean(T.localChangeZ,'omitnan'), ...
            mean(T.preElevationZ,'omitnan'),mean(T.localChangeZ<0), ...
            mean(T.earlyChangeZ>0 & T.localChangeZ<0),D{k}.source}; %#ok<AGROW>
    end
end
summary=cell2table(rows,'VariableNames',{'trial','group','nEvents', ...
    'meanEarlyChangeZ','meanLocalChangeZ','meanPreElevationZ', ...
    'fractionLocalDecrease','fractionElevatedButFalling','roiSource'});
writetable(summary,fullfile(outputFolder,'slowdown_trial_summary.csv'));
results=struct('events',events,'summary',summary,'trials',{D},'relativeTime',rel);
save(fullfile(outputFolder,'slowdown_baseline_results.mat'),'results','-v7.3');
disp(summary);
end

function D=read_trial(folder,rel)
[~,name]=fileparts(folder);
processed=fullfile(folder,'fly_data.mat');
if isfile(processed)
    q=load(processed,'fly_data');a=q.fly_data;source="saved roi.ts";
    if isfield(a,'analysisRoiSource'),source=string(a.analysisRoiSource);end
else
    files=dir(fullfile(folder,'*_or_z0z1_s_.mat'));
    assert(numel(files)==1,'Expected one input MAT in %s',folder);
    a=load(fullfile(files(1).folder,files(1).name),'dq','roi');source="legacy roi.ts";
end
dq=a.dq(1);t=double(dq.t(:));f=double(dq.bvf(:));
if isfield(dq,'bvy_deg'),r=double(dq.bvy_deg(:));else,r=rad2deg(double(dq.bvy(:)));end
h=double(dq.bh(:));x=double(a.roi(1).ts);if isvector(x),x=x(:)';end
assert(size(x,2)==numel(t),'Fluorescence/behavior length mismatch: %s',folder);
fs=1/median(diff(t));assert(all(diff(t)>0),'Nonmonotonic timestamps');
rawValid=all(isfinite(x),1)'&isfinite(f)&isfinite(r)&isfinite(h);
z=nan(size(x));
for j=1:size(x,1)
    s=smoothdata(x(j,:)','gaussian',max(3,round(.25*fs)));
    s=s-movmedian(s,max(3,round(30*fs)),'omitmissing');
    scale=1.4826*mad(s,1);assert(isfinite(scale)&&scale>0,'Invalid fluorescence scale');
    z(j,:)=s/scale;
end
z=mean(z,1)';f=smoothdata(f,'gaussian',max(3,round(.2*fs)));
r=abs(smoothdata(r,'gaussian',max(3,round(.12*fs))));
stim=false(size(t));
if isfield(dq,'heat_active')&&~isempty(dq.heat_active)
    stim=logical(dq.heat_active(:));
elseif isfield(dq,'heat')&&~isempty(dq.heat)
    heat=double(dq.heat(:));stim=heat>median(heat,'omitnan')+2*mad(heat,1);
end
valid=rawValid&isfinite(z)&isfinite(f)&isfinite(r)&isfinite(h);
if isfield(dq,'behaviorVolumeValid'),valid=valid&logical(dq.behaviorVolumeValid(:));end
if isfield(a,'jump_detected'),stim=stim|logical(a.jump_detected(:));end
if contains(lower(name),'grating')&&isfield(dq,'vvy'),stim=stim|abs(dq.vvy(:))>.02;end
% Fill only for locating candidates; every accepted window must be finite.
acc=[0;diff(f)]*fs; finiteAcc=acc(isfinite(acc));
prom=max(.5,2*1.4826*mad(finiteAcc,1));acc(~isfinite(acc))=0;
[~,candidate]=findpeaks(-acc,'MinPeakDistance',round(1.5*fs),'MinPeakProminence',prom);
E=cell(0,12);X=cell(1,4);for v=1:4,X{v}=zeros(numel(rel),0);end
for p=candidate(:)'
    window=t>=t(p)-3&t<=t(p)+3;
    if t(p)-3<t(1)||t(p)+3>t(end)||any(stim(window))||any(~valid(window)),continue,end
    tr=t-t(p);pre=tr>=-1&tr<=-.25;post=tr>=0&tr<=.75;
    before=mean(f(pre));after=mean(f(post));if before-after<1,continue,end
    % Require continuous acquisition; do not interpolate across timestamp gaps.
    if any(diff(t(window))>1.5/fs),continue,end
    % Circular interpolation avoids false heading changes at the -pi/pi seam.
    endpoints=t(p)+[-.25 .75];
    heading=atan2(interp1(t,sin(h),endpoints),interp1(t,cos(h),endpoints));
    hd=abs(rad2deg(atan2(sin(diff(heading)),cos(diff(heading)))));
    y=interp1(tr,z,rel);early=mean(y(rel>=-3&rel<=-2));
    local=mean(y(rel>=-1&rel<=-.25));response=mean(y(rel>=.25&rel<=.75));
    E(end+1,:)={string(name),t(p),hd>=15,hd,before,after,before-after, ...
        response-early,response-local,local-early,early,local}; %#ok<AGROW>
    X{1}(:,end+1)=y-early;X{2}(:,end+1)=y-local;
    X{3}(:,end+1)=interp1(tr,f,rel);X{4}(:,end+1)=interp1(tr,r,rel);
end
T=cell2table(E,'VariableNames',{'trial','timeSec','hasTurn','headingChangeDeg', ...
    'preFwd','postFwd','speedDrop','earlyChangeZ','localChangeZ','preElevationZ', ...
    'earlyBaselineZ','localBaselineZ'});
if isempty(E)
    T=table(strings(0,1),zeros(0,1),false(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
        zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
        'VariableNames',T.Properties.VariableNames);
end
D=struct('name',string(name),'source',source,'events',T,'traces',{X});
end

function keep=group_mask(T,g)
if g==1,keep=true(height(T),1);elseif g==2,keep=~T.hasTurn;else,keep=logical(T.hasTurn);end
end
function s=group_name(g)
names={'All slowdowns','Without reorientation','With reorientation (>=15 deg)'};s=names{g};
end
function render(D,t,prefix,heading,isSummary)
prefix=string(prefix);
fig=figure('Visible','off','Position',[20 20 1500 1050],'Color','w');
tl=tiledlayout(fig,4,3,'TileSpacing','compact','Padding','compact');
labels={'dF/F change vs -3:-2 s (z)','dF/F change vs -1:-0.25 s (z)', ...
    'Forward velocity (mm/s)','Rotational speed (deg/s)'};
for v=1:4
    axs=gobjects(1,3);
    for g=1:3
        ax=nexttile(tl);axs(g)=ax;hold(ax,'on');M=[];n=0;
        for k=1:numel(D)
            keep=group_mask(D{k}.events,g);Y=D{k}.traces{v}(:,keep);
            if isempty(Y),continue,end
            plot(ax,t,Y,'Color',[.85 .87 .90],'LineWidth',.4,'HandleVisibility','off');
            M(:,end+1)=mean(Y,2);n=n+size(Y,2); %#ok<AGROW>
        end
        if ~isempty(M)
            if isSummary,plot(ax,t,M,'Color',[.55 .67 .82],'LineWidth',.9,'HandleVisibility','off');end
            plot(ax,t,mean(M,2),'Color',[.05 .25 .65],'LineWidth',2.5);
        end
        xline(ax,0,':k');yline(ax,0,':','Color',[.6 .6 .6]);grid(ax,'on');xlim(ax,[-3 3]);
        if v<=2
            if v==1,limits=[-3 -2];else,limits=[-1 -.25];end
            yl=ylim(ax);patch(ax,[limits fliplr(limits)],[yl(1) yl(1) yl(2) yl(2)], ...
                [.3 .7 .3],'FaceAlpha',.08,'EdgeColor','none');
            xline(ax,.25,':','Color',[.8 .4 .1]);xline(ax,.75,':','Color',[.8 .4 .1]);
        end
        if g==1,ylabel(ax,labels{v});end
        if v==1,title(ax,sprintf('%s | n=%d, recordings=%d',group_name(g),n,size(M,2)));end
        if v==4,xlabel(ax,'Time from maximum deceleration (s)');end
    end
    linkaxes(axs,'y');
end
title(tl,heading,'Interpreter','none');
subtitle(tl,'Gray: individual events | blue: mean | green: reference window | orange: +0.25:+0.75 s response');
exportgraphics(fig,prefix+'.png','Resolution',180);
savefig_visible_on_open(fig,prefix+'.fig');close(fig);
end
