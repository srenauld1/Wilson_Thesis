%% Bar-jump heading responses and direction-matched optomotor turns
% Bar jumps: overlay every valid jump from -2 to +2 s for each trial.
% Optomotor: + visual velocity is CW and - visual velocity is CCW. Within
% each stimulus bout, retain the strongest fly turn in the matching signed
% direction when it reaches 50 deg/s, then align dF/F to that turn.

clearvars
baseFolder='/Users/sophiarenauld/stacks';
outFolder=fullfile(baseFolder,'ves041_behavior_quantification','visual_turn_alignment');
if ~isfolder(outFolder),mkdir(outFolder);end
barFolders=["20250812-20_ves041_fly3_clbarjumps"; ...
    "20250812-22_ves041_fly3_bbarjumps2"];
optoFolders=["20250812-21_ves041_fly3_olgrating"; ...
    "20250813-3_ves041_fly1_olgrating"; ...
    "20250819-10_ves041_fly1_olgrating5"];
tGrid=(-2:0.02:2)';

%% One bar-jump overlay figure per trial
barRows=cell(0,7);
for fi=1:numel(barFolders)
    S=load(fullfile(baseFolder,barFolders(fi),'smr_analysis','fly_data.mat'),'fly_data');x=S.fly_data;
    t=x.dq(1).t(:);dff=x.roi.ts(:);heading=x.dq(1).bh(:);vh=x.dq(1).vh(:);
    n=min([numel(t),numel(dff),numel(heading),numel(vh)]);
    t=t(1:n);dff=dff(1:n);heading=heading(1:n);vh=vh(1:n);
    jumpIdx=find(logical(x.jump_detected(:)));jumpIdx=jumpIdx(jumpIdx<=n);
    visualStep=atan2(sin([vh(1);diff(vh)]),cos([vh(1);diff(vh)]));
    jumpDeg=rad2deg(visualStep(jumpIdx));
    valid=abs(jumpDeg)>=45 & t(jumpIdx)>=2 & t(jumpIdx)<=t(end)-2;
    jumpIdx=jumpIdx(valid);jumpDeg=jumpDeg(valid);
    hMat=nan(numel(tGrid),numel(jumpIdx));vMat=hMat;dMat=hMat;
    hu=unwrap(heading);vhu=unwrap(vh);
    for j=1:numel(jumpIdx)
        rel=t-t(jumpIdx(j));base=mean(hu(rel>=-.25 & rel<0),'omitnan');
        visualBase=mean(vhu(rel>=-.25 & rel<0),'omitnan');
        hMat(:,j)=interp1(rel,rad2deg(hu-base),tGrid,'linear',NaN);
        vMat(:,j)=interp1(rel,rad2deg(vhu-visualBase),tGrid,'linear',NaN);
        dMat(:,j)=interp1(rel,dff,tGrid,'linear',NaN);
        barRows(end+1,:)={barFolders(fi),j,t(jumpIdx(j)),jumpDeg(j), ... %#ok<SAGROW>
            interp1(tGrid,hMat(:,j),.5),interp1(tGrid,hMat(:,j),1), ...
            interp1(tGrid,hMat(:,j),2)};
    end

    fig=figure('Color','w','Position',[100 100 1100 800]);
    tl=tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
    ax=nexttile(tl);hold(ax,'on');
    cc=lines(numel(jumpIdx));
    for j=1:numel(jumpIdx)
        plot(ax,tGrid,hMat(:,j),'Color',cc(j,:),'LineWidth',1.2, ...
            'DisplayName',sprintf('Jump %d: fly',j));
        plot(ax,tGrid,vMat(:,j),'--','Color',cc(j,:),'LineWidth',1.4, ...
            'DisplayName',sprintf('Jump %d: bar (%+.0f deg)',j,jumpDeg(j)));
    end
    plot(ax,tGrid,mean(hMat,2,'omitnan'),'k','LineWidth',2.7,'DisplayName','Mean fly');
    plot(ax,tGrid,mean(vMat,2,'omitnan'),'k--','LineWidth',2.7,'DisplayName','Mean bar');
    add_jump_line(ax);ylabel(ax,'Unwrapped heading change (deg)');
    title(ax,'Fly heading (solid) and bar heading, vh (dashed)');
    grid(ax,'on');legend(ax,'Location','eastoutside');

    ax=nexttile(tl);hold(ax,'on');
    for j=1:numel(jumpIdx)
        plot(ax,tGrid,dMat(:,j),'Color',cc(j,:),'LineWidth',1.2, ...
            'HandleVisibility','off');
    end
    plot(ax,tGrid,mean(dMat,2,'omitnan'),'k','LineWidth',2.7,'DisplayName','Mean');
    add_jump_line(ax);xlabel(ax,'Time from bar jump (s)');ylabel(ax,'VES041 dF/F');
    title(ax,'VES041 activity following bar jump');grid(ax,'on');legend(ax,'Location','best');
    sgtitle(tl,strrep(barFolders(fi),'_','\_'));
    stem=fullfile(outFolder,barFolders(fi)+"_bar_jump_heading_dff");
    savefig(fig,stem+".fig");exportgraphics(fig,stem+".png",'Resolution',250);close(fig);
end
barTable=cell2table(barRows,'VariableNames',{'trial','jumpNumber','jumpTime', ...
    'visualJumpDeg','headingChangeAt0p5s','headingChangeAt1s','headingChangeAt2s'});
barTable.trial=string(barTable.trial);
writetable(barTable,fullfile(outFolder,'bar_jump_unwrapped_heading.csv'));

%% Pool matched-direction optomotor turns across trials and directions
yawMat=[];dffMat=[];fwdMat=[];eventDirection=strings(0,1);eventTrial=strings(0,1);
optoRows=cell(0,8);
for fi=1:numel(optoFolders)
    S=load(fullfile(baseFolder,optoFolders(fi),'smr_analysis','fly_data.mat'),'fly_data');x=S.fly_data;
    t=x.dq(1).t(:);dff=x.roi.ts(:);yaw=x.dq(1).bvy_deg(:);
    fwd=x.dq(1).bvf(:);visualVelocity=x.dq(1).vvy(:);
    n=min([numel(t),numel(dff),numel(yaw),numel(fwd),numel(visualVelocity)]);
    t=t(1:n);dff=dff(1:n);yaw=yaw(1:n);fwd=fwd(1:n);visualVelocity=visualVelocity(1:n);
    fs=1/median(diff(t));minRun=round(1.5*fs);
    cw=clean_runs(visualVelocity>0.02,minRun);
    ccw=clean_runs(visualVelocity<-0.02,minRun);

    for direction=[1 -1]
        mask=cw;if direction<0,mask=ccw;end
        [starts,stops]=runs(mask);
        for b=1:numel(starts)
            alignedYaw=direction*yaw(starts(b):stops(b));
            [turnSpeed,localIdx]=max(alignedYaw,[],'omitnan');
            idx=starts(b)+localIdx-1;
            if turnSpeed<50 || t(idx)<2 || t(idx)>t(end)-2,continue,end
            rel=t-t(idx);
            yawMat(:,end+1)=interp1(rel,direction*yaw,tGrid,'linear',NaN); %#ok<SAGROW>
            dffMat(:,end+1)=interp1(rel,dff,tGrid,'linear',NaN); %#ok<SAGROW>
            fwdMat(:,end+1)=interp1(rel,fwd,tGrid,'linear',NaN); %#ok<SAGROW>
            label="CW";if direction<0,label="CCW";end
            eventDirection(end+1)=label;eventTrial(end+1)=optoFolders(fi); %#ok<SAGROW>
            optoRows(end+1,:)={optoFolders(fi),label,b,t(idx),turnSpeed, ... %#ok<SAGROW>
                dff(idx),mean(dff(rel>=0 & rel<=1),'omitnan'), ...
                mean(fwd(rel>=0 & rel<=1),'omitnan')};
        end
    end
end

fig=figure('Color','w','Position',[100 100 1150 900]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
plot_turn_average(nexttile(tl),tGrid,yawMat,eventDirection, ...
    'Direction-aligned optomotor turn','Matched-direction rotational velocity (deg/s)');
plot_turn_average(nexttile(tl),tGrid,dffMat,eventDirection, ...
    'VES041 activity during matched optomotor turns','dF/F');
plot_turn_average(nexttile(tl),tGrid,fwdMat,eventDirection, ...
    'Forward velocity during matched optomotor turns','Forward velocity (mm/s)');
xlabel(nexttile(tl,3),'Time from peak matched-direction turn (s)');
sgtitle(tl,sprintf('Optomotor turns pooled across directions and trials (n=%d)',size(yawMat,2)), ...
    'FontWeight','bold');
savefig(fig,fullfile(outFolder,'ves041_optomotor_turn_average.fig'));
exportgraphics(fig,fullfile(outFolder,'ves041_optomotor_turn_average.png'),'Resolution',250);
close(fig);

optoTable=cell2table(optoRows,'VariableNames',{'trial','direction','stimulusBout', ...
    'turnTime','peakMatchedTurnSpeed','dffAtTurn','meanDff0to1s','meanFwd0to1s'});
optoTable.trial=string(optoTable.trial);optoTable.direction=string(optoTable.direction);
writetable(optoTable,fullfile(outFolder,'optomotor_matched_turn_events.csv'));
disp(barTable);disp(optoTable)

function add_jump_line(ax)
h=xline(ax,0,'g--','Bar jump','LineWidth',1.5,'LabelVerticalAlignment','bottom');
h.HandleVisibility='off';xlim(ax,[-2 2]);
end

function mask=clean_runs(mask,minRun)
mask=logical(mask(:));[s,e]=runs(mask);
for i=1:numel(s),if e(i)-s(i)+1<minRun,mask(s(i):e(i))=false;end,end
end

function [s,e]=runs(mask)
d=diff([false;logical(mask(:));false]);s=find(d==1);e=find(d==-1)-1;
end

function plot_turn_average(ax,t,m,direction,ttl,ylab)
hold(ax,'on');cw=direction=="CW";ccw=direction=="CCW";
for i=find(cw)',plot(ax,t,m(:,i),'Color',[.20 .50 .90 .12],'HandleVisibility','off');end
for i=find(ccw)',plot(ax,t,m(:,i),'Color',[.90 .40 .15 .12],'HandleVisibility','off');end
plot(ax,t,mean(m(:,cw),2,'omitnan'),'Color',[.20 .50 .90],'LineWidth',1.8, ...
    'DisplayName',sprintf('CW mean (n=%d)',sum(cw)));
plot(ax,t,mean(m(:,ccw),2,'omitnan'),'Color',[.90 .40 .15],'LineWidth',1.8, ...
    'DisplayName',sprintf('CCW mean (n=%d)',sum(ccw)));
plot(ax,t,mean(m,2,'omitnan'),'k','LineWidth',2.7,'DisplayName','Direction-pooled mean');
h=xline(ax,0,'k--','Turn peak','LabelVerticalAlignment','bottom');h.HandleVisibility='off';
xlim(ax,[-2 2]);ylabel(ax,ylab);title(ax,ttl);grid(ax,'on');legend(ax,'Location','best');
end
