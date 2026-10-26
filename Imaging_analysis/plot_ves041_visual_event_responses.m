function summary = plot_ves041_visual_event_responses(a2p_data,savepath,jump)
% Event-aligned VES041 plots for bar-jump and optomotor trials.
% Bar trials: unwrapped fly heading (bh), bar heading (vh), and dF/F.
% Grating trials: direction-matched turns, plus a stringent subset requiring
% low forward speed (mean <=2 mm/s from the turn through +1 s).

summary=table();trial=string(get_trial_name(savepath));
t=a2p_data.dq(1).t(:);dff=a2p_data.roi.ts(1,:)';
n=min(numel(t),numel(dff));t=t(1:n);dff=dff(1:n);
fs=1/median(diff(t));

if jump && isfield(a2p_data,'jump_detected')
    bh=a2p_data.dq(1).bh(:);vh=a2p_data.dq(1).vh(:);
    n=min([n,numel(bh),numel(vh)]);t=t(1:n);dff=dff(1:n);bh=bh(1:n);vh=vh(1:n);
    idx=find(logical(a2p_data.jump_detected(:)));idx=idx(idx<=n);
    dvh=atan2(sin([vh(1);diff(vh)]),cos([vh(1);diff(vh)]));jumpDeg=rad2deg(dvh(idx));
    if isfield(a2p_data,'jump_source') && ~strcmp(a2p_data.jump_source,'visual heading detector')
        keep=t(idx)>=2 & t(idx)<=t(end)-2;
    else
        keep=abs(jumpDeg)>=45 & t(idx)>=2 & t(idx)<=t(end)-2;
    end
    idx=idx(keep);jumpDeg=jumpDeg(keep);tRel=(-2:0.02:2)';
    B=nan(numel(tRel),numel(idx));V=B;D=B;ubh=unwrap(bh);uvh=unwrap(vh);
    for k=1:numel(idx)
        rel=t-t(idx(k));pre=rel>=-.25 & rel<0;
        B(:,k)=interp1(rel,rad2deg(ubh-mean(ubh(pre),'omitnan')),tRel,'linear',NaN);
        V(:,k)=interp1(rel,rad2deg(uvh-mean(uvh(pre),'omitnan')),tRel,'linear',NaN);
        D(:,k)=interp1(rel,dff,tRel,'linear',NaN);
    end
    fig=figure('Color','w','Position',[100 100 1100 800]);
    tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');cc=lines(numel(idx));
    ax=nexttile(tl);hold(ax,'on');
    for k=1:numel(idx)
        plot(ax,tRel,B(:,k),'Color',cc(k,:),'DisplayName',sprintf('Jump %d: fly',k));
        plot(ax,tRel,V(:,k),'--','Color',cc(k,:),'LineWidth',1.3, ...
            'DisplayName',sprintf('Jump %d: bar (%+.0f deg)',k,jumpDeg(k)));
    end
    plot(ax,tRel,mean(B,2,'omitnan'),'k','LineWidth',2.5,'DisplayName','Mean fly');
    plot(ax,tRel,mean(V,2,'omitnan'),'k--','LineWidth',2.5,'DisplayName','Mean bar');
    event_line(ax,'Bar jump');ylabel(ax,'Unwrapped heading change (deg)');
    title(ax,'Fly heading (solid) and bar heading, vh (dashed)');grid(ax,'on');
    legend(ax,'Location','eastoutside');
    ax=nexttile(tl);hold(ax,'on');
    for k=1:numel(idx),plot(ax,tRel,D(:,k),'Color',cc(k,:),'HandleVisibility','off');end
    plot(ax,tRel,mean(D,2,'omitnan'),'k','LineWidth',2.5,'DisplayName','Mean');
    event_line(ax,'Bar jump');xlabel(ax,'Time from bar jump (s)');ylabel(ax,'VES041 dF/F');
    title(ax,'VES041 activity following bar jump');grid(ax,'on');legend(ax,'Location','best');
    sgtitle(tl,strrep(trial,'_','\_'));
    savefig(fig,fullfile(savepath,'bar_jump_heading_vh_dff.fig'));
    exportgraphics(fig,fullfile(savepath,'bar_jump_heading_vh_dff.png'),'Resolution',200);close(fig);
    summary=table(trial,numel(idx),'VariableNames',{'trial','nBarJumps'});
end

if contains(lower(savepath),'grating')
    yaw=a2p_data.dq(1).bvy_deg(:);fwd=a2p_data.dq(1).bvf(:);vv=a2p_data.dq(1).vvy(:);
    n=min([n,numel(yaw),numel(fwd),numel(vv)]);t=t(1:n);dff=dff(1:n);
    yaw=yaw(1:n);fwd=fwd(1:n);vv=vv(1:n);tRel=(-2:0.02:2)';
    baseline=movmedian(dff,max(3,round(30*fs)),'omitmissing');noise=1.4826*mad(dff-baseline,1);
    [~,dffPeaks]=findpeaks(dff,'MinPeakProminence',3*noise, ...
        'MinPeakDistance',round(2*fs),'MinPeakWidth',round(.25*fs));
    dirs=[1 -1];Y=[];F=[];D=[];lowSpeed=[];nearPeak=[];peakBefore=[];rows=cell(0,11);
    for direction=dirs
        mask=clean_runs(direction*vv>.02,round(1.5*fs));[starts,stops]=runs(mask);
        for b=1:numel(starts)
            [turnSpeed,ii]=max(direction*yaw(starts(b):stops(b)),[],'omitnan');
            c=starts(b)+ii-1;if turnSpeed<50 || t(c)<2 || t(c)>t(end)-2,continue,end
            rel=t-t(c);pre=rel>=-1 & rel<0;post=rel>=0 & rel<=1;
            preF=mean(fwd(pre),'omitnan');postF=mean(fwd(post),'omitnan');isLowSpeed=postF<=2;
            lags=t(dffPeaks)-t(c);isNear=any(lags>=-1 & lags<=1);isBefore=any(lags>=-2 & lags<=.25);
            Y(:,end+1)=interp1(rel,direction*yaw,tRel,'linear',NaN); %#ok<AGROW>
            F(:,end+1)=interp1(rel,fwd,tRel,'linear',NaN);D(:,end+1)=interp1(rel,dff,tRel,'linear',NaN); %#ok<AGROW>
            lowSpeed(end+1)=isLowSpeed;nearPeak(end+1)=isNear;peakBefore(end+1)=isBefore; %#ok<AGROW>
            label="CW";if direction<0,label="CCW";end
            rows(end+1,:)={trial,label,b,t(c),turnSpeed,preF,postF,postF-preF, ...
                isLowSpeed,isNear,isBefore}; %#ok<AGROW>
        end
    end
    fig=figure('Color','w','Position',[80 50 1400 900]);tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
    plot_subset(nexttile(tl),tRel,Y,true(size(lowSpeed)),nearPeak,'All corrective turns','Matched-direction rotational velocity (deg/s)');
    plot_subset(nexttile(tl),tRel,Y,logical(lowSpeed),nearPeak,'Corrective turns at low forward speed','Matched-direction rotational velocity (deg/s)');
    plot_subset(nexttile(tl),tRel,F,true(size(lowSpeed)),nearPeak,'Forward velocity','Forward velocity (mm/s)');
    plot_subset(nexttile(tl),tRel,F,logical(lowSpeed),nearPeak,'Forward velocity: low-speed subset','Forward velocity (mm/s)');
    plot_subset(nexttile(tl),tRel,D,true(size(lowSpeed)),nearPeak,'VES041 activity','dF/F');
    plot_subset(nexttile(tl),tRel,D,logical(lowSpeed),nearPeak,'VES041 activity: low-speed subset','dF/F');
    sgtitle(tl,sprintf('%s | corrective optomotor turns',strrep(trial,'_','\_')));
    savefig(fig,fullfile(savepath,'optomotor_corrective_turns.fig'));
    exportgraphics(fig,fullfile(savepath,'optomotor_corrective_turns.png'),'Resolution',200);close(fig);
    eventTable=cell2table(rows,'VariableNames',{'trial','direction','stimulusBout','turnTime', ...
        'matchedTurnSpeed','preFwd','postFwd','fwdChange','isLowForwardSpeed', ...
        'vesPeakWithin1s','vesPeakBeforeOrConcurrent'});
    eventTable.trial=string(eventTable.trial);eventTable.direction=string(eventTable.direction);
    writetable(eventTable,fullfile(savepath,'optomotor_corrective_turn_events.csv'));
    summary=table(trial,numel(lowSpeed),sum(lowSpeed),sum(nearPeak),sum(nearPeak & lowSpeed), ...
        sum(peakBefore),sum(peakBefore & lowSpeed),'VariableNames',{'trial','nCorrectiveTurns', ...
        'nCorrectiveTurnsAtLowSpeed','nWithVesPeakWithin1s','nLowSpeedWithVesPeakWithin1s', ...
        'nWithVesPeakBeforeOrConcurrent','nLowSpeedWithVesPeakBeforeOrConcurrent'});
end
end

function plot_subset(ax,t,m,use,nearPeak,ttl,ylab)
use=logical(use);hold(ax,'on');
for i=find(use),plot(ax,t,m(:,i),'Color',[.65 .65 .65 .15],'HandleVisibility','off');end
plot(ax,t,mean(m(:,use),2,'omitnan'),'k','LineWidth',2.5, ...
    'DisplayName',sprintf('Mean (n=%d; peak-associated=%d)',sum(use),sum(nearPeak & use)));
h=xline(ax,0,'k--','Turn peak','LabelVerticalAlignment','bottom');h.HandleVisibility='off';
xlim(ax,[-2 2]);xlabel(ax,'Time from turn peak (s)');ylabel(ax,ylab);title(ax,ttl);grid(ax,'on');legend(ax,'Location','best');
end

function event_line(ax,label)
h=xline(ax,0,'g--',label,'LineWidth',1.4,'LabelVerticalAlignment','bottom');h.HandleVisibility='off';xlim(ax,[-2 2]);
end
function mask=clean_runs(mask,minRun)
mask=logical(mask(:));[s,e]=runs(mask);for i=1:numel(s),if e(i)-s(i)+1<minRun,mask(s(i):e(i))=false;end,end
end
function [s,e]=runs(mask),d=diff([false;logical(mask(:));false]);s=find(d==1);e=find(d==-1)-1;end
function name=get_trial_name(savepath),[folder,~]=fileparts(savepath);[~,name]=fileparts(folder);end
