function report=analyze_smp544_oct26_5_5_peak_behavior()
%ANALYZE_SMP544_OCT26_5_5_PEAK_BEHAVIOR Compare raw peaks and troughs in 5_5.
% All events are fluorescence-selected; behavior never enters detection.
% Input fluorescence normalization and full manual ROI are unchanged.
% Raw modes use NO smoothing. Broad-trough mode uses a centered 1.5-s
% moving mean. Each mode subtracts its own 30-s moving median for detection.
% Height/prominence >= max(1 dF/F,4*noise sigma), where sigma is estimated
% from first differences of the detection residual. Merge <=0.5-s breaks;
% one extremum per event. Raw event span >=0.25 s; broad event span >=1 s.
% Broad events additionally require >=1-s half-prominence width and a
% >=1-s descent from a preceding maximum in the last 5 s, >=60% descending
% differences. These explicit shape restrictions reject short sharp dips.
% Missing imaging data/timestamp gaps are never bridged. All heat contexts
% are included. Behavior-aligned plots require every +/-3-s volume to have
% valid behavior; detected events failing this check remain in the CSV.
% Descriptive Pearson correlations use valid paired samples, not tests.
% Event changes compare [0,1] with [-2,-1] s; signed forward velocity,
% rotational SPEED, and net heading displacement are kept distinct.
% Smoothing is centered, so it must not be used to infer neural lead times.
% Registration QC is held for this input: all results are exploratory.
assert(~feature('ShowFigureWindows'),'Use MATLAB -nodisplay.');
input='/Users/sophiarenauld/stacks/2026-10-26/a2p_exploratory_20261027/20261026_5_5_smp544_fly2_heatpyramid/smp544_manual_analysis_20261028/fly_data.mat';
out=fullfile(imaging_analysis_artifact_root(),'diagnostics','smp544_20261026_5_5_peak_behavior_20261028');
if ~isfolder(out),mkdir(out);end
q=load(input,'fly_data');a=q.fly_data;
assert(isequaln(a.roi(1).ts,a.roi(1).pxfilter.manualts),'Expected full manual ROI');
t=double(a.dq(1).t(:));raw=double(a.roi(1).ts(:));fs=1/median(diff(t));
assert(numel(t)==numel(raw)&&all(diff(t)>0),'Invalid fluorescence/time arrays');
f=double(a.dq(1).bvf(:));signedRot=double(a.dq(1).bvy_deg(:));rot=abs(signedRot);
h=double(a.dq(1).bh(:));heat=double(a.dq(1).heat(:));
valid=isfinite(f)&isfinite(rot)&isfinite(h)&isfinite(raw);
if isfield(a.dq(1),'bergii')&&isfield(a.dq(1).bergii,'behaviorVolumeValid')
    valid=valid&logical(a.dq(1).bergii.behaviorVolumeValid(:));
end
heavy=block_smooth(t,raw,fs,1.5);
names=["raw peaks","raw troughs","1.5 s broad troughs"];
series={raw,raw,heavy};signs=[1 -1 -1];minimumDuration=[.25 .25 1];
rel=(-3:.05:3)';modes=cell(3,1);allTables=cell(3,1);rows=cell(0,11);
for m=1:3
    [events,base,noise,threshold]=detect(t,series{m},signs(m),fs,minimumDuration(m),m==3);
    T=table;N=size(events,1);T.mode=repmat(names(m),N,1);
    T.timeSec=t(events(:,2));T.startSec=t(events(:,1));T.endSec=t(events(:,3));
    T.amplitudeDff=events(:,4);T.prominenceDff=events(:,5);
    T.durationSec=(events(:,3)-events(:,1)+1)/fs;
    T.eligibleBehaviorWindow=false(N,1);T.heatAtEventV=heat(events(:,2));
    fields={'preFwd','postFwd','changeFwd','preRotSpeed','postRotSpeed', ...
        'changeRotSpeed','netHeadingPostDeg'};
    for field=fields,T.(field{1})=nan(N,1);end
    traces=cell(1,5);for v=1:5,traces{v}=zeros(numel(rel),0);end
    for j=1:N
        p=events(j,2);tr=t-t(p);w=tr>=-3&tr<=3;
        if t(p)-3<t(1)||t(p)+3>t(end)||any(~valid(w))||any(diff(t(w))>1.5/fs),continue,end
        T.eligibleBehaviorWindow(j)=true;pre=tr>=-2&tr<=-1;post=tr>=0&tr<=1;
        T.preFwd(j)=mean(f(pre));T.postFwd(j)=mean(f(post));
        T.changeFwd(j)=T.postFwd(j)-T.preFwd(j);
        T.preRotSpeed(j)=mean(rot(pre));T.postRotSpeed(j)=mean(rot(post));
        T.changeRotSpeed(j)=T.postRotSpeed(j)-T.preRotSpeed(j);
        hh=atan2(interp1(t,sin(h),t(p)+[0 1]),interp1(t,cos(h),t(p)+[0 1]));
        dh=diff(hh);T.netHeadingPostDeg(j)=abs(rad2deg(atan2(sin(dh),cos(dh))));
        traces{1}(:,end+1)=interp1(tr,raw-mean(raw(pre)),rel);
        traces{2}(:,end+1)=interp1(tr,heavy-mean(heavy(pre)),rel);
        traces{3}(:,end+1)=interp1(tr,f,rel);
        traces{4}(:,end+1)=interp1(tr,rot,rel);
        traces{5}(:,end+1)=interp1(tr,heat,rel);
    end
    U=T(T.eligibleBehaviorWindow,:);
    rows(end+1,:)={names(m),N,height(U),noise,threshold, ...
        mean(U.preFwd,'omitnan'),mean(U.postFwd,'omitnan'),mean(U.changeFwd,'omitnan'), ...
        mean(U.preRotSpeed,'omitnan'),mean(U.postRotSpeed,'omitnan'),mean(U.netHeadingPostDeg,'omitnan')}; %#ok<AGROW>
    modes{m}=struct('name',names(m),'baseline',base,'noiseSigma',noise, ...
        'threshold',threshold,'indices',events,'events',T,'traces',{traces});
    allTables{m}=T;
end
events=vertcat(allTables{:});writetable(events,fullfile(out,'events.csv'));
report=cell2table(rows,'VariableNames',{'mode','nDetected','nBehaviorEligible','noiseSigma', ...
    'thresholdDff','preFwd','postFwd','changeFwd','preRotSpeed','postRotSpeed','netHeadingPostDeg'});
writetable(report,fullfile(out,'event_behavior_summary.csv'));
% Report baseline-relative correlations too: raw dF/F contains slow drift
% that can covary with recording time and declining locomotion.
correlationRows=cell(0,6);
signals={raw,heavy,raw-modes{1}.baseline,heavy-modes{3}.baseline};
signalNames=["raw dff","1.5 s dff","raw local residual","1.5 s local residual"];
variables={f,rot,signedRot,heat};variableNames=["forward velocity","rotational speed","signed yaw velocity","heat monitor"];
for s=1:numel(signals)
    for v=1:4
        for context=["all","no heat","heat"]
            mask=valid&isfinite(signals{s})&isfinite(variables{v});
            if context=="no heat",mask=mask&heat<=.05;elseif context=="heat",mask=mask&heat>.05;end
            coefficient=NaN;
            if nnz(mask)>2&&std(signals{s}(mask))>0&&std(variables{v}(mask))>0
                coefficient=corr(signals{s}(mask),variables{v}(mask));
            end
            correlationRows(end+1,:)={signalNames(s),variableNames(v),context,nnz(mask),coefficient,"descriptive Pearson r"}; %#ok<AGROW>
        end
    end
end
correlations=cell2table(correlationRows,'VariableNames',{'signal','variable','context','nSamples','r','method'});
writetable(correlations,fullfile(out,'sample_correlations.csv'));
overview(t,raw,heavy,f,rot,heat,valid,modes,out);
aligned(rel,modes,out);
results=struct('input',input,'registrationQc',a.analysisRegistrationQc,'modes',{modes}, ...
    'summary',report,'correlations',correlations,'relativeTime',rel);
save(fullfile(out,'results.mat'),'results','-v7.3');disp(report);disp(correlations);
end

function [s,e]=blocks(t,y,fs)
v=isfinite(y);b=[true;diff(t)>1.5/fs];
s=find(v&([true;~v(1:end-1)]|b));e=find(v&([~v(2:end);true]|[b(2:end);true]));
end
function sm=block_smooth(t,y,fs,seconds)
[s,e]=blocks(t,y,fs);sm=nan(size(y));
for k=1:numel(s),ix=s(k):e(k);sm(ix)=movmean(y(ix),max(1,round(seconds*fs)));end
end
function [events,base,noise,threshold]=detect(t,y,polarity,fs,minDuration,broad)
[s,e]=blocks(t,y,fs);base=nan(size(y));signal=base;d=[];
for k=1:numel(s)
    ix=s(k):e(k);base(ix)=movmedian(y(ix),max(3,round(30*fs)));
    signal(ix)=polarity*(y(ix)-base(ix));d=[d;diff(signal(ix))]; %#ok<AGROW>
end
noise=mad(d,1)/(sqrt(2)*.6745);threshold=max(1,4*noise);events=zeros(0,5);
for k=1:numel(s)
    ix=s(k):e(k);x=signal(ix);if numel(x)<3,continue,end
    dd=diff([false;x>=threshold;false]);on=find(dd==1);off=find(dd==-1)-1;
    if isempty(on),continue,end
    starts=on(1);ends=[];last=off(1);
    for j=2:numel(on)
        if on(j)-last-1<=round(.5*fs),last=off(j);
        else,ends(end+1,1)=last;starts(end+1,1)=on(j);last=off(j);end %#ok<AGROW>
    end
    ends(end+1,1)=last;
    args={'MinPeakHeight',threshold,'MinPeakProminence',threshold};
    if broad,args=[args,{'MinPeakWidth',max(1,round(fs))}];end
    [~,loc,~,prom]=findpeaks(x,args{:});
    for j=1:numel(starts)
        if (ends(j)-starts(j)+1)/fs<minDuration,continue,end
        c=find(loc>=starts(j)&loc<=ends(j));if isempty(c),continue,end
        [~,best]=max(x(loc(c)));c=c(best);p=ix(loc(c));
        if broad
            pre=(p-round(5*fs)):(p-round(fs));
            if isempty(pre)||pre(1)<s(k),continue,end
            [~,highest]=max(y(pre));ramp=pre(highest):p;
            if (numel(ramp)-1)/fs<1||mean(diff(y(ramp))<=0)<.6,continue,end
        end
        events(end+1,:)=[ix(starts(j)),p,ix(ends(j)),x(loc(c)),prom(c)]; %#ok<AGROW>
    end
end
end

function overview(t,raw,heavy,f,rot,heat,valid,M,out)
fig=figure('Visible','off','Color','w','Position',[20 20 1550 1000]);
tl=tiledlayout(fig,4,1,'TileSpacing','compact');ax=gobjects(4,1);
for v=1:4,ax(v)=nexttile(tl);hold(ax(v),'on');grid(ax(v),'on');end
plot(ax(1),t,raw,'Color',[.72 .8 .9],'DisplayName','Unsmoothed dF/F');
plot(ax(1),t,heavy,'Color',[.05 .2 .6],'LineWidth',1.4,'DisplayName','1.5 s moving mean');
colors={[.8 .2 .1],[.55 .2 .7],[.05 .6 .35]};markers={'^','v','o'};
for m=1:3
    p=M{m}.indices(:,2);y=raw(p);if m==3,y=heavy(p);end
    scatter(ax(1),t(p),y,48,colors{m},markers{m},'filled', ...
        'DisplayName',sprintf('%s (n=%d)',M{m}.name,numel(p)));
end
ylabel(ax(1),'dF/F');legend(ax(1),'Location','northeast','NumColumns',2);
f(~valid)=NaN;rot(~valid)=NaN;
plot(ax(2),t,f,'k');ylabel(ax(2),'Forward (mm/s)');
plot(ax(3),t,rot,'Color',[.85 .15 .1]);ylabel(ax(3),'Rot. speed (deg/s)');
plot(ax(4),t,heat,'Color',[.85 .45 .1]);ylabel(ax(4),'Heat monitor (V)');xlabel(ax(4),'Time (s)');
for v=2:4,for m=1:3,for p=M{m}.indices(:,2)',xline(ax(v),t(p),':','Color',colors{m});end,end,end
linkaxes(ax,'x');xlim(ax(1),[t(1) t(end)]);
title(tl,'20261026 5_5 | raw peaks, raw troughs, and broad smoothed troughs');
subtitle(tl,'Full manual ROI | fluorescence-only detection | registration QC held');
exportgraphics(fig,fullfile(out,'detection_overview.png'),'Resolution',150);
savefig_visible_on_open(fig,fullfile(out,'detection_overview.fig'));close(fig);
end
function aligned(t,M,out)
fig=figure('Visible','off','Color','w','Position',[20 20 1550 1050]);
tl=tiledlayout(fig,4,3,'TileSpacing','compact');labels={'dF/F vs -2:-1 s','Forward (mm/s)','Rot. speed (deg/s)','Heat monitor (V)'};
for v=1:4
    axesRow=gobjects(1,3);
    for m=1:3
        ax=nexttile(tl);axesRow(m)=ax;hold(ax,'on');grid(ax,'on');
        index=v;if v>1,index=v+1;end;X=M{m}.traces{index};
        if isempty(X),text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');
        else
            plot(ax,t,X,'Color',[.85 .85 .88],'LineWidth',.5);
            plot(ax,t,mean(X,2),'Color',[.05 .25 .65],'LineWidth',2.4);
            if v==1,plot(ax,t,mean(M{m}.traces{2},2),'Color',[.05 .6 .3],'LineWidth',2);end
        end
        xline(ax,0,':k');xlim(ax,[-3 3]);
        if v==1,title(ax,sprintf('%s | n=%d/%d eligible',M{m}.name,size(X,2),height(M{m}.events)));end
        if m==1,ylabel(ax,labels{v});end
        if v==4,xlabel(ax,'Time from fluorescence extremum (s)');end
    end
    linkaxes(axesRow,'y');
end
title(tl,'20261026 5_5 | behavior around fluorescence-selected events');
subtitle(tl,'Gray: individual traces; blue: raw mean; green: 1.5 s fluorescence mean | centered filtering blurs timing');
exportgraphics(fig,fullfile(out,'event_aligned_behavior.png'),'Resolution',170);
savefig_visible_on_open(fig,fullfile(out,'event_aligned_behavior.fig'));close(fig);
end
