%% VES041 and compensatory heading responses after bar jumps
% Retains visual jumps >=45 deg (removing startup transients). A jump is
% classified as VES041-peak-associated when an agreed-definition peak falls
% from -0.25 to +2 s around the jump. Positive direction-aligned heading
% means that the fly turned opposite the visual jump (compensatory).

clearvars
baseFolder='/Users/sophiarenauld/stacks';
outFolder=fullfile(baseFolder,'ves041_behavior_quantification','bar_jump_responses');
if ~isfolder(outFolder),mkdir(outFolder);end
folders=["20250812-20_ves041_fly3_clbarjumps"; ...
    "20250812-22_ves041_fly3_bbarjumps2"];

tGrid=(-2:0.02:4)';
dffMat=[];fwdMat=[];rotMat=[];headingMat=[];hasPeak=[];
eventRows=cell(0,12);

for fi=1:numel(folders)
    S=load(fullfile(baseFolder,folders(fi),'smr_analysis','fly_data.mat'),'fly_data');
    x=S.fly_data;
    t=x.dq(1).t(:);dff=x.roi.ts(:);fwd=x.dq(1).bvf(:);
    rot=x.dq(1).bvy_deg(:);heading=x.dq(1).bh(:);visualHeading=x.dq(1).vh(:);
    n=min([numel(t),numel(dff),numel(fwd),numel(rot),numel(heading),numel(visualHeading)]);
    t=t(1:n);dff=dff(1:n);fwd=fwd(1:n);rot=rot(1:n);
    heading=heading(1:n);visualHeading=visualHeading(1:n);
    fs=1/median(diff(t));

    baseline=movmedian(dff,max(3,round(30*fs)),'omitmissing');
    noiseSigma=1.4826*mad(dff-baseline,1);
    [~,peaks]=findpeaks(dff,'MinPeakProminence',3*noiseSigma, ...
        'MinPeakDistance',round(2*fs),'MinPeakWidth',round(0.25*fs));

    jumpIdx=find(logical(x.jump_detected(:)));jumpIdx=jumpIdx(jumpIdx<=n);
    signedVisualStep=atan2(sin([visualHeading(1);diff(visualHeading)]), ...
        cos([visualHeading(1);diff(visualHeading)]));
    jumpAngle=rad2deg(signedVisualStep(jumpIdx));
    valid=abs(jumpAngle)>=45 & t(jumpIdx)>=2 & t(jumpIdx)<=t(end)-4;
    jumpIdx=jumpIdx(valid);jumpAngle=jumpAngle(valid);

    for ji=1:numel(jumpIdx)
        j=jumpIdx(ji);tj=t(j);rel=t-tj;
        associated=peaks(t(peaks)>=tj-0.25 & t(peaks)<=tj+2);
        peakFlag=~isempty(associated);

        hu=unwrap(heading);
        preMask=rel>=-0.5 & rel<0;
        heading0=mean(hu(preMask),'omitnan');
        % Opposite the visual jump is plotted positive.
        compensated=-sign(jumpAngle(ji))*rad2deg(hu-heading0);

        dffEpoch=interp1(rel,dff,tGrid,'linear',NaN);
        fwdEpoch=interp1(rel,fwd,tGrid,'linear',NaN);
        rotEpoch=interp1(rel,abs(rot),tGrid,'linear',NaN);
        headingEpoch=interp1(rel,compensated,tGrid,'linear',NaN);
        dffMat(:,end+1)=dffEpoch; %#ok<SAGROW>
        fwdMat(:,end+1)=fwdEpoch; %#ok<SAGROW>
        rotMat(:,end+1)=rotEpoch; %#ok<SAGROW>
        headingMat(:,end+1)=headingEpoch; %#ok<SAGROW>
        hasPeak(end+1)=peakFlag; %#ok<SAGROW>

        at1=interp1(tGrid,headingEpoch,1);at2=interp1(tGrid,headingEpoch,2);
        post=t>=tj & t<=tj+2;
        peakTimes=t(associated)-tj;
        eventRows(end+1,:)={folders(fi),ji,tj,jumpAngle(ji),peakFlag, ... %#ok<SAGROW>
            numel(associated),strjoin(string(round(peakTimes,3)),';'),at1,at2, ...
            mean(fwd(post),'omitnan'),mean(abs(rot(post)),'omitnan'), ...
            max(abs(rot(post)),[],'omitnan')};
    end
end

events=cell2table(eventRows,'VariableNames',{'trial','jumpNumber','jumpTime', ...
    'visualJumpDeg','hasAssociatedVes041Peak','nAssociatedPeaks', ...
    'peakTimesRelativeToJump','compensatoryHeadingAt1sDeg', ...
    'compensatoryHeadingAt2sDeg','meanFwd0to2s','meanRotSpeed0to2s', ...
    'peakRotSpeed0to2s'});
events.trial=string(events.trial);events.peakTimesRelativeToJump=string(events.peakTimesRelativeToJump);
writetable(events,fullfile(outFolder,'bar_jump_response_events.csv'));

%% Plot aligned responses
colors.noPeak=[.60 .60 .60];colors.peak=[.80 .18 .16];
fig=figure('Color','w','Position',[80 80 1400 850]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

plot_epochs(nexttile(tl,1),tGrid,headingMat,hasPeak,colors, ...
    'Direction-aligned heading','Compensatory heading change (deg)');
plot_epochs(nexttile(tl,2),tGrid,dffMat,hasPeak,colors, ...
    'VES041 activity','dF/F');
plot_epochs(nexttile(tl,3),tGrid,fwdMat,hasPeak,colors, ...
    'Forward velocity','Forward velocity (mm/s)');
plot_epochs(nexttile(tl,4),tGrid,rotMat,hasPeak,colors, ...
    'Rotational speed','Rotational speed (deg/s)');

sgtitle(tl,sprintf('Bar-jump responses: VES041 peak-associated (n=%d) versus no peak (n=%d)', ...
    sum(hasPeak),sum(~hasPeak)),'FontWeight','bold');
savefig(fig,fullfile(outFolder,'ves041_bar_jump_responses.fig'));
exportgraphics(fig,fullfile(outFolder,'ves041_bar_jump_responses.png'),'Resolution',250);
disp(events)

function plot_epochs(ax,t,m,hasPeak,colors,ttl,ylab)
hasPeak=logical(hasPeak);
hold(ax,'on');
for i=find(~hasPeak),plot(ax,t,m(:,i),'Color',[colors.noPeak .25],'HandleVisibility','off');end
for i=find(hasPeak),plot(ax,t,m(:,i),'Color',[colors.peak .25],'HandleVisibility','off');end
plot(ax,t,mean(m(:,~hasPeak),2,'omitnan'),'Color',colors.noPeak,'LineWidth',2.5, ...
    'DisplayName',sprintf('No VES041 peak (n=%d)',sum(~hasPeak)));
plot(ax,t,mean(m(:,hasPeak),2,'omitnan'),'Color',colors.peak,'LineWidth',2.5, ...
    'DisplayName',sprintf('VES041 peak (n=%d)',sum(hasPeak)));
hJump=xline(ax,0,'k--','Bar jump','LabelVerticalAlignment','bottom');
hJump.HandleVisibility='off';
xlim(ax,[t(1) t(end)]);xlabel(ax,'Time from bar jump (s)');ylabel(ax,ylab);
title(ax,ttl);grid(ax,'on');legend(ax,'Location','best');
end
