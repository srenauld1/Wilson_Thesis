function [pairs,summary,sensitivity]=match_ves041_fly1_heat_slowdowns(out,analysisSubfolder,files)
% Within-recording, behavior-only heat/no-heat matching for 10/20 fly 1.
% Preserve the existing behavior-selected event catalog and refresh dF/F
% outcomes from current dff005015 inputs. No fluorescence enters matching.
% Primary: slowdowns with >=30 deg integrated absolute yaw in +/-0.5 s.
% Also report slowdowns without that turn, matched in a separate stratum.
% Match pre-event median fwd (-.75:-.25 s), speed drop to +.25:+.75 s,
% and integrated absolute yaw. Primary per-pair calipers: 1 mm/s, 1 mm/s,
% 15 deg. Optimal assignment maximizes count then minimizes scaled distance;
% matches are one-to-one without replacement and never cross recordings.
% Outcome: mean dF/F 0:+1 s minus mean -2:-1 s. No hypothesis tests.
oldVisible=get(groot,'defaultFigureVisible');c=onCleanup(@()set(groot,'defaultFigureVisible',oldVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
root='/Users/sophiarenauld/stacks/2026-10-20';
if nargin<1,out='/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis/15s_ave_f/fly1_matched_heat';end
if nargin<2,analysisSubfolder='smr_analysis';end
if ~isfolder(out),mkdir(out);end
if nargin<3,files=dir(fullfile(root,'*fly1*','*_or_z0z1_s_.mat'));end
catalog=table;data=struct;rel=(-2:.05:3)';
for k=1:numel(files)
    folder=files(k).folder;a=load(fullfile(folder,files(k).name),'dq','roi');a=import_ves041_protocol(a,folder);
    if ~ismember(a.protocol.kind,{'heat-barrier','heat-stripes'}),continue,end
    assert(strcmp(a.roi.opt.nrm.nrmstr,'dff005015'),'Expected new 15 s local-baseline dF/F');
    id=string(erase(files(k).name,'_or_z0z1_s_.mat'));
    source=fullfile(folder,analysisSubfolder,'heat_exploration','slowdown_events.csv');
    E=readtable(source,'Delimiter',',','ReadVariableNames',true,'TextType','string');
    t=a.dq(1).t(:);d=double(a.roi.ts(:));f=a.dq(1).bvf(:);r=rad2deg(a.dq(1).bvy(:));h=a.dq(1).heat(:);fs=1/median(diff(t));
    smF=smoothdata(f,'movmedian',max(3,round(.25*fs)));smR=smoothdata(r,'movmedian',max(3,round(.25*fs)));
    E.recording(:)=id;E.eventId=(1:height(E))';E.commandOffCore=false(height(E),1);
    X=nan(numel(rel),height(E),4);
    for j=1:height(E)
        et=E.eventTime(j);pre=t>=et-.75&t<et-.25;post=t>=et+.25&t<=et+.75;core=t>=et-.5&t<=et+.5;
        expected=[median(smF(pre)),median(smF(pre))-median(smF(post)),trapz(t(core),abs(smR(core)))];
        assert(max(abs(expected-[E.speedBefore(j),E.speedDrop(j),E.absoluteTurnDeg(j)]))<1e-3,'Saved behavior no longer matches current input');
        if E.context(j)=="heat",assert(all(h(core)>.05));elseif E.context(j)=="no heat",assert(all(h(core)<=.05));end
        E.commandOffCore(j)=all(isfinite(h(core))&h(core)<=.05);
        b=t>=et-2&t<et-1;response=t>=et&t<=et+1;
        E.baselineDff(j)=mean(d(b));E.responseDff(j)=mean(d(response));E.deltaDff(j)=E.responseDff(j)-E.baselineDff(j);
        vars={d,f,abs(r),h};
        for v=1:4,X(:,j,v)=interp1(t,vars{v},et+rel,'linear',NaN);end
        X(:,j,1)=X(:,j,1)-E.baselineDff(j);
    end
    catalog=[catalog;E]; %#ok<AGROW>
    data.(matlab.lang.makeValidName(char(id)))=struct('E',E,'X',X);
end
writetable(catalog,fullfile(out,'candidate_events.csv'));
caliper=[1 1 15];pairs=table;summary=table;sensitivity=table;controlPairs=table;epochSummary=table;epochPairs=table;
ids=unique(catalog.recording,'stable');
for id=ids'
    E=catalog(catalog.recording==id,:);
    for turn=0:1
        [M,S]=matchOne(E,turn,caliper,false);pairs=[pairs;M];summary=[summary;S]; %#ok<AGROW>
        for factor=[.5 1 1.5 2]
            [~,S]=matchOne(E,turn,caliper*factor,false);S.caliperMultiplier=factor;S.includeControls=false;sensitivity=[sensitivity;S]; %#ok<AGROW>
        end
        [MC,S]=matchOne(E,turn,caliper,true);controlPairs=[controlPairs;MC];S.caliperMultiplier=1;S.includeControls=true;sensitivity=[sensitivity;S]; %#ok<AGROW>
        if id=="20261020_11_2"
            % Sensitivity to the previously identified ~395 s baseline step:
            % exclude +/-10 s and prohibit matches across the step. This is
            % a time/recording-quality restriction, not a dF/F-based match.
            EE=E(abs(E.eventTime-395)>10,:);EE.matchEpoch=double(EE.eventTime>395);
            for includeControl=[false true]
                [ME,SE]=matchOne(EE,turn,caliper,includeControl);SE.includeControls=includeControl;epochSummary=[epochSummary;SE]; %#ok<AGROW>
                ME.includeControls=repmat(includeControl,height(ME),1);epochPairs=[epochPairs;ME]; %#ok<AGROW>
            end
        end
    end
end
writetable(pairs,fullfile(out,'matched_pairs.csv'));writetable(summary,fullfile(out,'matched_summary.csv'));writetable(sensitivity,fullfile(out,'caliper_sensitivity.csv'));
writetable(controlPairs,fullfile(out,'control_inclusive_pairs.csv'));writetable(epochSummary,fullfile(out,'trial11_same_epoch_summary.csv'));writetable(epochPairs,fullfile(out,'trial11_same_epoch_pairs.csv'));
settings=struct('sourceEvents','Existing behavior-selected slowdown catalogs; fluorescence refreshed from current inputs', ...
    'normalization','dff005015','primaryCaliper',[1 1 15], ...
    'matchVariables',{{'speedBefore','speedDrop','absoluteTurnDeg'}},'primaryTurnThresholdDeg',30, ...
    'withoutReplacement',true,'withinRecordingOnly',true,'outcomeWindow',[0 1],'baselineWindow',[-2 -1], ...
    'contextWindow',[-.5 .5],'includeControlInPrimary',false,'fluorescenceUsedForMatching',false);
save(fullfile(out,'matched_analysis.mat'),'pairs','summary','sensitivity','controlPairs','epochSummary','epochPairs','settings','rel','data','-v7.3');
for turn=0:1
    drawComparison(pairs(pairs.reorientation==turn,:),data,rel,out,sprintf('all_trials_turn%d',turn));
    for id=ids'
        drawComparison(pairs(pairs.reorientation==turn&pairs.recording==id,:),data,rel,out,sprintf('%s_turn%d',id,turn));
    end
end
drawComparison(controlPairs(controlPairs.reorientation==1,:),data,rel,out,'control_inclusive_turn1');
drawComparison(epochPairs(epochPairs.reorientation==1 & epochPairs.includeControls,:),data,rel,out,'trial11_same_epoch_control_inclusive_turn1');
disp(summary);disp(sensitivity(sensitivity.reorientation==1,:));
disp(epochSummary);
end

function [M,S]=matchOne(E,turn,caliper,includeControls)
% All decisions in this function use behavior/context only, never dF/F.
hi=find(E.context=="heat" & E.reorientation==turn);
off=E.context=="no heat";
if includeControls,off=off|(E.context=="control"&E.commandOffCore);end
oi=find(off&E.reorientation==turn);vars={'speedBefore','speedDrop','absoluteTurnDeg'};
A=E{hi,vars};B=E{oi,vars};C=zeros(numel(hi),numel(oi));allowed=true(size(C));
for v=1:3,D=abs(A(:,v)-B(:,v)');allowed=allowed&D<=caliper(v);C=C+(D/caliper(v)).^2;end
C=sqrt(C);C(~allowed)=Inf;
if ismember('matchEpoch',E.Properties.VariableNames)
    allowed=allowed & E.matchEpoch(hi)==E.matchEpoch(oi)';C(~allowed)=Inf;
end
if isempty(C)||~any(allowed,'all'),p=zeros(0,2);else,p=matchpairs(C,1e6);end
assert(numel(unique(p(:,1)))==size(p,1)&&numel(unique(p(:,2)))==size(p,1),'Event reused');
H=hi(p(:,1));O=oi(p(:,2));n=numel(H);cost=zeros(n,1);
for j=1:n,cost(j)=C(p(j,1),p(j,2));end
M=table(repmat(E.recording(1),n,1),repmat(turn,n,1),E.eventId(H),E.eventId(O),E.eventTime(H),E.eventTime(O), ...
    E.speedBefore(H),E.speedBefore(O),E.speedDrop(H),E.speedDrop(O),E.absoluteTurnDeg(H),E.absoluteTurnDeg(O), ...
    E.deltaDff(H),E.deltaDff(O),E.deltaDff(H)-E.deltaDff(O),cost,E.context(O), ...
    'VariableNames',{'recording','reorientation','heatEventId','noHeatEventId','heatTime','noHeatTime', ...
    'heatPreSpeed','noHeatPreSpeed','heatSpeedDrop','noHeatSpeedDrop','heatTurnDeg','noHeatTurnDeg', ...
    'heatDeltaDff','noHeatDeltaDff','heatMinusNoHeat','matchCost','noHeatContext'});
if n>0
    assert(all(abs(M.heatPreSpeed-M.noHeatPreSpeed)<=caliper(1)+eps));assert(all(abs(M.heatSpeedDrop-M.noHeatSpeedDrop)<=caliper(2)+eps));assert(all(abs(M.heatTurnDeg-M.noHeatTurnDeg)<=caliper(3)+eps));
end
S=table(E.recording(1),turn,numel(hi),numel(oi),n, ...
    mean(E.deltaDff(hi),'omitnan'),mean(E.deltaDff(oi),'omitnan'),mean(M.heatDeltaDff,'omitnan'),mean(M.noHeatDeltaDff,'omitnan'),mean(M.heatMinusNoHeat,'omitnan'),median(M.heatMinusNoHeat,'omitnan'), ...
    mean(M.heatPreSpeed,'omitnan'),mean(M.noHeatPreSpeed,'omitnan'),mean(M.heatSpeedDrop,'omitnan'),mean(M.noHeatSpeedDrop,'omitnan'),mean(M.heatTurnDeg,'omitnan'),mean(M.noHeatTurnDeg,'omitnan'), ...
    'VariableNames',{'recording','reorientation','nHeatCandidates','nNoHeatCandidates','nPairs','unmatchedHeatDelta','unmatchedNoHeatDelta','matchedHeatDelta','matchedNoHeatDelta','meanPairedDifference','medianPairedDifference', ...
    'matchedHeatPreSpeed','matchedNoHeatPreSpeed','matchedHeatDrop','matchedNoHeatDrop','matchedHeatTurn','matchedNoHeatTurn'});
end

function drawComparison(M,data,rel,out,stem)
fig=figure('Visible','off','Color','w','Position',[40 40 1500 1050]);tl=tiledlayout(2,4,'TileSpacing','compact');
vars={{'noHeatPreSpeed','heatPreSpeed'},{'noHeatSpeedDrop','heatSpeedDrop'},{'noHeatTurnDeg','heatTurnDeg'}};
labs={'Pre-slowdown forward (mm/s)','Speed drop (mm/s)','Absolute turn (deg)'};ids=unique(M.recording,'stable');cc=lines(numel(ids));
for v=1:3
    ax=nexttile(tl);hold(ax,'on');
    for k=1:numel(ids),sel=M.recording==ids(k);scatter(ax,M.(vars{v}{1})(sel),M.(vars{v}{2})(sel),40,cc(k,:),'filled','DisplayName',strrep(ids(k),'2026',''));end
    lim=[0 1];if height(M)>0,vals=[M.(vars{v}{1});M.(vars{v}{2})];lim=[min([0;vals]) max([1;vals])*1.05];end
    plot(ax,lim,lim,':k','HandleVisibility','off');xlim(ax,lim);ylim(ax,lim);axis(ax,'square');grid(ax,'on');xlabel(ax,['No heat: ' labs{v}]);ylabel(ax,['Heat: ' labs{v}]);
    if v==1&&height(M)>0,legend(ax,'Interpreter','none','Location','best');end
end
ax=nexttile(tl);hold(ax,'on');
for k=1:numel(ids),sel=M.recording==ids(k);plot(ax,[1 2],[M.noHeatDeltaDff(sel) M.heatDeltaDff(sel)]','-o','Color',cc(k,:));end
xlim(ax,[.7 2.3]);xticks(ax,[1 2]);xticklabels(ax,{'No heat','Heat'});yline(ax,0,':k');ylabel(ax,'dF/F change: 0:1 minus -2:-1 s');grid(ax,'on');title(ax,sprintf('%d matched pairs',height(M)));
Xoff=nan(numel(rel),height(M),4);Xheat=Xoff;
for j=1:height(M),z=data.(matlab.lang.makeValidName(char(M.recording(j))));Xoff(:,j,:)=z.X(:,M.noHeatEventId(j),:);Xheat(:,j,:)=z.X(:,M.heatEventId(j),:);end
labs={'dF/F minus -2:-1 s baseline','Forward speed (mm/s)','Rotational speed (deg/s)','Measured heat command (V)'};
for v=1:4
    ax=nexttile(tl);hold(ax,'on');
    if height(M)>0
        plot(ax,rel,Xoff(:,:,v),'Color',[.7 .7 .7],'HandleVisibility','off');plot(ax,rel,Xheat(:,:,v),'Color',[1 .8 .7],'HandleVisibility','off');
        plot(ax,rel,mean(Xoff(:,:,v),2,'omitnan'),'Color',[.2 .2 .2],'LineWidth',2,'DisplayName','No heat');plot(ax,rel,mean(Xheat(:,:,v),2,'omitnan'),'Color',[.9 .25 .05],'LineWidth',2,'DisplayName','Heat');
    else,text(ax,.5,.5,'No eligible matches','Units','normalized','HorizontalAlignment','center');end
    xline(ax,0,':k','HandleVisibility','off');xlim(ax,[-2 3]);ylabel(ax,labs{v});xlabel(ax,'Time from maximum deceleration (s)');grid(ax,'on');if v==1&&height(M)>0,legend(ax,'Location','best');end
end
title(tl,{strrep(stem,'_',' '),'Within-trial 1:1 matching | calipers: 1 mm/s pre-speed, 1 mm/s speed drop, 15 deg turn | no dF/F selection'});
savefig_visible_on_open(fig,fullfile(out,[stem '.fig']));exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',145);close(fig);
end
