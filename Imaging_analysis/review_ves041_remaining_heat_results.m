function [states,slowdowns]=review_ves041_remaining_heat_results(out)
% Descriptive heat/no-heat review of the seven newly analyzed recordings.
% Rate denominator uses the same +/-0.5 s context and -2:+3 s complete valid
% window as the peak-aligned plots. Control episodes and command transitions
% are kept separate. Heat means measured command >0.05 V, not temperature.
% No amplitude rejection, new event thresholds, hypothesis tests or matching.
oldVisible=get(groot,'defaultFigureVisible');
c=onCleanup(@()set(groot,'defaultFigureVisible',oldVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
if nargin<1,out='/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis/15s_ave_f/remaining_heat';end
opts=delimitedTextImportOptions('NumVariables',10,'Delimiter',',');opts.DataLines=[2 Inf];
opts.VariableNames={'recording','context','acquisitionStatus','status','nPeaks','nEligiblePeaks','nSlowdowns','nPlacedBarrierEpisodes','output','error'};
opts.VariableTypes={'string','string','string','string','double','double','double','double','string','string'};
R=readtable(fullfile(out,'run_report.csv'),opts);R=R(R.status=="success",:);
P=readtable(fullfile(out,'all_detected_peaks.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');
E=readtable(fullfile(out,'all_slowdown_events.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');
contexts=["no heat","heat","control","heat entry","heat exit","multiple switches","control transition","unknown"];
rows=cell(0,11);srows=cell(0,9);
for k=1:height(R)
    a=load(fullfile(R.output(k),'fly_data.mat'),'fly_data');a=a.fly_data;
    t=a.dq(1).t(:);d=double(a.roi.ts(:));f=a.dq(1).bvf(:);r=abs(rad2deg(a.dq(1).bvy(:)));h=a.dq(1).heat(:);n=numel(t);ix=(1:n)';fs=1/median(diff(t));
    control=nan(n,1);if isfield(a.protocol,'control_trial'),control=a.protocol.control_trial(:);end
    valid=isfinite(d)&isfinite(f)&isfinite(r);
    if isfield(a.dq(1),'bergii')&&isfield(a.dq(1).bergii,'behaviorVolumeValid'),valid=valid&logical(a.dq(1).bergii.behaviorVolumeValid(:));end
    L=interp1(t,ix,t-2,'next',NaN);U=interp1(t,ix,t+3,'previous',NaN);
    LC=interp1(t,ix,t-.5,'next',NaN);UC=interp1(t,ix,t+.5,'previous',NaN);
    bad=[0;cumsum(~valid)];gaps=[0;cumsum([false;diff(t)>2/fs])];C=repmat("unknown",n,1);eligible=false(n,1);
    for q=find(isfinite(L)&isfinite(U))'
        lo=L(q);hi=U(q);
        if bad(hi+1)-bad(lo)>0||gaps(hi+1)-gaps(lo+1)>0,continue,end
        core=LC(q):UC(q);if any(~isfinite(h(core))),continue,end
        eligible(q)=true;C(q)=classify_ves041_peak_heat_context(h(core)>.05,control(core));
    end
    weight=[diff(t);median(diff(t))];peaks=P(P.recording==R.recording(k),:);
    take=logical(peaks.eligibleAlignedWindow)&eligible(peaks.sampleIndex);
    assert(all(C(peaks.sampleIndex(take))==peaks.context(take)),'Peak context disagrees with exposure classification');
    for ctx=contexts
        sel=eligible&C==ctx;num=sum(take&peaks.context==ctx);seconds=sum(weight(sel));
        rows(end+1,:)={R.recording(k),R.context(k),ctx,seconds,num,60*num/seconds,mean(d(sel),'omitnan'),median(d(sel),'omitnan'),mean(f(sel),'omitnan'),mean(r(sel),'omitnan'),mean(abs(f(sel))<.5)*100}; %#ok<AGROW>
    end
    e=E(E.recording==R.recording(k),:);
    for ctx=["no heat","heat","control","transition"]
        for turn=0:1
            sel=e.context==ctx & e.reorientation==turn;
            srows(end+1,:)={R.recording(k),ctx,turn,sum(sel),mean(e.deltaDff(sel),'omitnan'),median(e.deltaDff(sel),'omitnan'),mean(e.speedBefore(sel),'omitnan'),mean(e.speedDrop(sel),'omitnan'),mean(e.absoluteTurnDeg(sel),'omitnan')}; %#ok<AGROW>
        end
    end
end
states=cell2table(rows,'VariableNames',{'recording','protocol','context','usableSeconds','nPeaks','peaksPerMinute','meanDff','medianDff','meanForward','meanRotSpeed','stationaryPercent'});
slowdowns=cell2table(srows,'VariableNames',{'recording','context','reorientation','nEvents','meanDeltaDff','medianDeltaDff','meanPreSpeed','meanSpeedDrop','meanAbsoluteTurnDeg'});
writetable(states,fullfile(out,'heat_state_exposure_review.csv'));writetable(slowdowns,fullfile(out,'heat_slowdown_review.csv'));
disp(states(ismember(states.context,["no heat","heat","control"]),:));
fig=figure('Visible','off','Color','w','Position',[40 40 1500 1050]);tl=tiledlayout(2,2,'TileSpacing','compact');cc=lines(height(R));
labels={'Peak rate (per usable minute)','Mean 15 s-baseline dF/F','Slowdowns without reorientation: mean dF/F change','Slowdowns with reorientation: mean dF/F change'};
for panel=1:4
    ax=nexttile(tl);hold(ax,'on');
    for k=1:height(R)
        if panel<=2
            x=states(states.recording==R.recording(k),:);y=[NaN NaN];
            for j=1:2,ctx=["no heat","heat"];v=x(x.context==ctx(j),:);if panel==1,y(j)=v.peaksPerMinute;else,y(j)=v.meanDff;end;end
        else
            x=slowdowns(slowdowns.recording==R.recording(k)&slowdowns.reorientation==panel-3,:);y=[x.meanDeltaDff(x.context=="no heat"),x.meanDeltaDff(x.context=="heat")];
        end
        plot(ax,[1 2],y,'-o','Color',cc(k,:),'LineWidth',1.5,'DisplayName',strrep(R.recording(k),'2026',''));
    end
    xlim(ax,[.7 2.3]);xticks(ax,[1 2]);xticklabels(ax,{'No heat','Heat'});ylabel(ax,labels{panel});grid(ax,'on');
    if panel>=3,yline(ax,0,':k','HandleVisibility','off');end
    if panel==1,legend(ax,'Location','eastoutside','Interpreter','none');end
end
title(tl,{'VES041 heat comparison: each line is one recording, not an independent fly','Stable context within +/-0.5 s; complete -2:+3 s windows | slowdown change: 0:1 s minus -2:-1 s | no behavior matching'});
savefig_visible_on_open(fig,fullfile(out,'heat_no_heat_review.fig'));exportgraphics(fig,fullfile(out,'heat_no_heat_review.png'),'Resolution',150);close(fig);
end
