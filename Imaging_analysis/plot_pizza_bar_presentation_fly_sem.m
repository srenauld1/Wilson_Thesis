function result=plot_pizza_bar_presentation_fly_sem(baselineMode,eventWindow)
%PLOT_PIZZA_BAR_PRESENTATION_FLY_SEM Four-column equal-fly bar-only summary.
% Same groups as the preceding individual-trial figure: original columns
% 1,4,5,6 (VES peaks; slowing+rotation; rotation peaks split by forward drop).
% Does NOT apply the separate <2 mm/s forward-speed restriction.
% Load saved event arrays, pool events within fly across included bar trials,
% average events within fly, then calculate equal-weight mean and sample
% SD/sqrt(N) across contributing FLY MEANS. Never treat events as flies.
% Missing event groups are NaN, not zero. Original dF/F units, source traces,
% event windows and eligibility are unchanged. No new smoothing or detection.
% The existing map includes three inferred fly identities; these remain
% flagged in the saved methods and mapping. No excluded or dark trials added.
% All outputs are new companions. Use MATLAB -nodisplay to prevent popups.
% Optional 'trial_mean' mode subtracts each metric's full-valid-recording
% mean separately in every trial BEFORE pooling events within a fly. This is
% not pre-event baseline subtraction and is not a z-score (no SD division).
% 'fly_zscore' uses the existing full-recording bar-only within-fly metric
% mean and sample SD, without an additional trial-mean subtraction.
if nargin<1,baselineMode='none';end
if nargin<2,eventWindow=[];end
baselineMode=string(baselineMode);assert(ismember(baselineMode,["none","trial_mean","fly_zscore"]));
centerTrials=baselineMode=="trial_mean";
zScore=baselineMode=="fly_zscore";
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';
folder=fullfile(root,'bar');
P=load(fullfile(folder,'ves041_bar_pooled_equal_fly_peaks_slowdowns.mat'));
F=P.F;flies=P.flies;nf=numel(flies);rel=P.rel(:);columns=[1 4 5 6];
expandedWindow=~isempty(eventWindow);
if expandedWindow
    assert(isequal(eventWindow,[-3 3]),'Supported extension is [-3 3].');
    % linspace includes both endpoints despite floating-point source spacing.
    rel=linspace(eventWindow(1),eventWindow(2), ...
        round(diff(eventWindow)/median(diff(P.rel)))+1)';
else,eventWindow=[-2 3];end
normalizationCenter=[];normalizationScale=[];
if zScore
    N=load(fullfile(folder,'ves041_bar_pooled_equal_fly_peaks_slowdowns_normalized.mat'));
    assert(isequal(N.flies,flies));normalizationCenter=N.center;normalizationScale=N.scale;
    assert(all(isfinite(normalizationScale(:)))&&all(normalizationScale(:)>0));
end
pooled=cell(nf,4);counts=zeros(nf,4);trialCounts=zeros(nf,4);trialRows=cell(0,6);
originalCounts=counts;windowRows=cell(0,6);
trialMeans=zeros(height(F),3);baselineRows=cell(0,5);
for k=1:height(F)
    source=fullfile(root,'analysis_data',F.trial(k)+"_events_with_rotation_peaks_and_speedups.mat");
    S=load(source,'X','rel','groups');assert(isequal(S.rel(:),P.rel(:))&&numel(S.X)==8);
    if centerTrials||expandedWindow
        A=load(fullfile(root,'analysis_data',F.trial(k)+"_events.mat"),'d');d=A.d;
    end
    if centerTrials
        valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
        variables={d.dff,d.fwd,abs(d.yaw)};metricNames=["dff","forwardVelocity","rotationalSpeed"];
        assert(nnz(valid)>1,'Insufficient full-recording baseline samples');
        for v=1:3
            reference=variables{v}(valid);trialMeans(k,v)=mean(reference);
            assert(abs(mean(reference-trialMeans(k,v)))<1e-9);
            baselineRows(end+1,:)={F.trial(k),F.flyId(k),metricNames(v),trialMeans(k,v),nnz(valid)}; %#ok<AGROW>
        end
    end
    f=find(flies==F.flyId(k));assert(isscalar(f));n=zeros(1,4);
    for c=1:4
        X=S.X{columns(c)};assert(size(X,1)==numel(P.rel)&&size(X,3)==3);
        assert(size(X,2)==numel(S.groups{columns(c)})&&all(isfinite(X(:))));
        originalCounts(f,c)=originalCounts(f,c)+size(X,2);
        if expandedWindow
            % Extend the EXISTING event selections, never redetect events.
            % Require complete valid bracketing samples for the whole -3:+3
            % window; never pad, extrapolate or interpolate across gaps.
            ids=S.groups{columns(c)};Y=nan(numel(rel),numel(ids),3);accepted=false(numel(ids),1);
            valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
            vars={d.dff,d.fwd,abs(d.yaw)};dt=median(diff(d.t));
            for j=1:numel(ids)
                t0=d.t(ids(j));lo=find(d.t<=t0+eventWindow(1),1,'last');
                hi=find(d.t>=t0+eventWindow(2),1,'first');reason="retained";
                if isempty(lo)||isempty(hi),reason="recording edge";
                elseif any(~valid(lo:hi)),reason="invalid data in expanded window";
                elseif any(diff(d.t(lo:hi))>2*dt),reason="timestamp gap";
                else
                    for v=1:3,Y(:,j,v)=interp1(d.t,vars{v},t0+rel,'linear',NaN);end
                    accepted(j)=all(isfinite(Y(:,j,:)),'all');
                    if ~accepted(j),reason="interpolation invalid";end
                end
                windowRows(end+1,:)={F.trial(k),F.flyId(k),c,t0,accepted(j),reason}; %#ok<AGROW>
            end
            X=Y(:,accepted,:);assert(all(isfinite(X(:))));
        end
        if centerTrials
            original=X;X=X-reshape(trialMeans(k,:),1,1,3);
            % Subtraction is constant within each metric, never changes
            % event membership, timing or the within-event waveform.
            recovered=X+reshape(trialMeans(k,:),1,1,3);
            assert(all(abs(recovered(:)-original(:))<1e-10));
        end
        if zScore
            original=X;
            X=(X-reshape(normalizationCenter(f,:),1,1,3))./reshape(normalizationScale(f,:),1,1,3);
            recovered=X.*reshape(normalizationScale(f,:),1,1,3)+reshape(normalizationCenter(f,:),1,1,3);
            assert(all(abs(recovered(:)-original(:))<1e-10));
        end
        pooled{f,c}=cat(2,pooled{f,c},X);n(c)=size(X,2);
        counts(f,c)=counts(f,c)+n(c);trialCounts(f,c)=trialCounts(f,c)+(n(c)>0);
    end
    trialRows(end+1,:)={F.trial(k),F.flyId(k),n(1),n(2),n(3),n(4)}; %#ok<AGROW>
end
flyMeans=nan(numel(rel),nf,3,4);
for f=1:nf,for c=1:4
    if counts(f,c)>0,flyMeans(:,f,:,c)=mean(pooled{f,c},2);end
end,end
% The first two displayed columns must exactly match the old bar-only means.
assert(all(counts<=originalCounts,'all'));
if ~expandedWindow
assert(isequal(counts(:,1:2),P.counts(:,[1 4])));
for c=1:2
    old=P.flyMeans(:,:,:,columns(c));current=flyMeans(:,:,:,c);
    assert(isequal(isnan(old),isnan(current)));
    if zScore
        expected=N.flyMeans(:,:,:,columns(c));difference=abs(expected-current);
        assert(isequal(isnan(expected),isnan(current))&&all(difference(isfinite(difference))<1e-10));
    elseif ~centerTrials
        difference=abs(old-current);assert(all(difference(isfinite(difference))<1e-10));
    end
end
end
groupMean=nan(numel(rel),3,4);sem=groupMean;nPerTime=groupMean;
for c=1:4,for v=1:3
    X=flyMeans(:,:,v,c);n=sum(isfinite(X),2);
    groupMean(:,v,c)=mean(X,2,'omitnan');
    se=std(X,0,2,'omitnan')./sqrt(n);se(n<2)=NaN;
    sem(:,v,c)=se;nPerTime(:,v,c)=n;
end,end
% Verify equal fly weighting independently of event counts.
for c=1:4
    contributing=find(counts(:,c)>0);expected=zeros(numel(rel),3);
    for f=contributing',expected=expected+squeeze(mean(pooled{f,c},2))/numel(contributing);end
    difference=abs(expected-groupMean(:,:,c));assert(all(difference(:)<1e-10));
end
titles={'VES041 peaks','Slowdowns + rotation', ...
    {'Rotation peaks','Speed drop <1 mm/s'}, ...
    {'Rotation peaks','Speed drop \geq1 mm/s'}};
labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};
if centerTrials
    labels={'dF/F - trial mean','Forward - trial mean (mm/s)','Rot. speed - trial mean (deg/s)'};
end
if zScore,labels={'dF/F (z-score)','Forward velocity (z-score)','Rotational speed (z-score)'};end
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};
thinColors={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};
fig=figure('Visible','off','Color','w','Position',[20 20 1800 1050]);
tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','loose');
limits=zeros(3,2);contributingFlies=sum(counts>0,1);
sharedZLimits=[];
if zScore
    q=[flyMeans(:);groupMean(:)-sem(:);groupMean(:)+sem(:)];q=q(isfinite(q));
    bound=max(1,ceil(2*max(abs(q)))/2);sharedZLimits=[-bound bound];
end
for row=1:3
    values=flyMeans(:,:,row,:);values=values(isfinite(values));
    lowBand=groupMean(:,row,:)-sem(:,row,:);upBand=groupMean(:,row,:)+sem(:,row,:);
    values=[values;lowBand(isfinite(lowBand));upBand(isfinite(upBand))];
    limits(row,:)=[min(values) max(values)];span=diff(limits(row,:));
    if span==0,span=1;end;limits(row,:)=limits(row,:)+[-.06 .06]*span;
    if zScore,limits(row,:)=sharedZLimits;end
    for c=1:4
        ax=nexttile(tl);hold(ax,'on');X=flyMeans(:,:,row,c);
        for f=find(counts(:,c)>0)'
            plot(ax,rel,X(:,f),'Color',thinColors{row},'LineWidth',1.1, ...
                'Tag','flyMean','UserData',flies(f));
        end
        mu=groupMean(:,row,c);se=sem(:,row,c);
        fill(ax,[rel;flipud(rel)],[mu-se;flipud(mu+se)],colors{row}, ...
            'FaceAlpha',.20,'EdgeColor','none','Tag','flySem');
        plot(ax,rel,mu,'Color',colors{row},'LineWidth',3,'Tag','equalFlyMean');
        xline(ax,0,':','Color',[.45 .45 .45],'LineWidth',1);
        if centerTrials||zScore,yline(ax,0,':','Color',[.6 .6 .6]);end
        xlim(ax,eventWindow);ylim(ax,limits(row,:));ax.FontSize=18;ax.LineWidth=1.1;
        ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=eventWindow(1):1:eventWindow(2);
        if c==1,ylabel(ax,labels{row},'FontSize',21);end
        if row==3,xlabel(ax,'Time (s)','FontSize',21);else,ax.XTickLabel=[];end
        if row==1,title(ax,titles{c},'FontSize',23,'FontWeight','bold');end
    end
end
caption=sprintf('Thin: fly means | Thick: equal-fly mean | Shading: SEM across flies | %d recordings, %d estimated flies',height(F),nf);
if centerTrials,caption={'Each metric centered by its full-recording trial mean',caption};end
if zScore,caption={'Within-fly z-scores; full-recording reference',caption};end
subtitle(tl,caption,'FontSize',18);
assert(numel(findall(fig,'Type','axes'))==12&&numel(findall(fig,'Tag','flySem'))==12);
assert(numel(findall(fig,'Tag','flyMean'))==3*sum(contributingFlies));
stem=fullfile(folder,'ves041_bar_pooled_presentation_four_selections_fly_sem');
if centerTrials,stem=[stem '_trial_mean_subtracted'];end
if zScore,stem=[stem '_zscored'];end
if expandedWindow,stem=[stem '_minus3_plus3'];end
exportgraphics(fig,[stem '.png'],'Resolution',200);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
weightRows=cell(0,6);
for c=1:4,for f=1:nf
    weight=double(counts(f,c)>0)/contributingFlies(c);
    weightRows(end+1,:)={flies(f),c,columns(c),counts(f,c),trialCounts(f,c),weight}; %#ok<AGROW>
end,end
weights=cell2table(weightRows,'VariableNames',{'flyId','displayColumn','sourceColumn','events','contributingTrials','weight'});
trials=cell2table(trialRows,'VariableNames',{'trial','flyId','peaks','slowdownsWithRotation','rotationDropBelow1','rotationDropAtLeast1'});
writetable(weights,[stem '_fly_weights.csv']);writetable(trials,[stem '_trial_counts.csv']);
methods=struct('withinFly','Pool all eligible events across included bar recordings, then mean events', ...
    'betweenFly','Equal-weight mean of contributing fly means', ...
    'sem','Sample SD of contributing fly means divided by sqrt(number of contributing flies), pointwise', ...
    'notEventSem',true,'sourceColumns',columns,'additionalLowForwardSpeedGate',false, ...
    'normalization','Original dF/F and physical kinematic units; no new normalization or smoothing', ...
    'inferredFlyIdentities',nnz(contains(lower(flies),'inferred')));
methods.baselineMode=baselineMode;
methods.eventWindowSec=eventWindow;
methods.expandedWindowReextracted=expandedWindow;
if centerTrials
    methods.normalization='Subtract each metric full-valid-recording trial mean before pooling events; no SD scaling; not pre-event subtraction';
end
if zScore
    methods.normalization='Each metric uses its existing full-recording within-fly bar-only mean and sample SD: (value-mean)/SD; no additional trial centering';
    methods.sharedZLimits=sharedZLimits;
    writetable(N.normalization,[stem '_normalization_parameters.csv']);
end
baselines=cell2table(baselineRows,'VariableNames',{'trial','flyId','metric','subtractedTrialMean','validSamples'});
windowAudit=cell2table(windowRows,'VariableNames',{'trial','flyId','displayColumn','eventTime','retained','reason'});
if expandedWindow,writetable(windowAudit,[stem '_window_eligibility.csv']);end
if centerTrials,writetable(baselines,[stem '_trial_baselines.csv']);end
save([stem '.mat'],'F','flies','rel','flyMeans','groupMean','sem','nPerTime','counts', ...
    'trialCounts','contributingFlies','columns','methods','weights','trials','baselineMode','trialMeans','baselines', ...
    'normalizationCenter','normalizationScale','sharedZLimits','originalCounts','eventWindow','windowAudit');
fid=fopen([stem '_methods.txt'],'w');assert(fid>=0);
fprintf(fid,'Bar-only cohort: %d recordings grouped as %d estimated flies; %d identities inferred in the existing mapping.\n',height(F),nf,methods.inferredFlyIdentities);
fprintf(fid,'One thin line per contributing fly. Events pooled within each fly across its bar recordings. Thick line is the equal-weight mean across contributing flies. Shaded SEM is sample SD/sqrt(N flies), NOT event-level SEM.\n');
fprintf(fid,'Same four groups as the individual presentation figure: VES041 peaks; slowdowns with rotation; rotation peaks with forward drop <1 or >=1 mm/s. No additional forward-speed <2 filter. Original saved fluorescence, kinematics, event windows and eligibility unchanged. Full cohort mapping and counts are saved with the data.\n');
if expandedWindow
    fprintf(fid,'Window extension overrides the preceding window statement: existing event centers are re-extracted over -3:+3 s, requiring all source/bracketing samples valid with no gap >2 native sampling intervals. Incomplete events are excluded from the entire average, never partially padded. Counts before and after extension are saved. No event redetection.\n');
end
if centerTrials
    fprintf(fid,'All THREE metrics are centered by their own full-valid-recording mean separately within each source trial, before events are pooled by fly. No division by SD and no event-specific pre-window baseline. Original uncentered outputs preserved.\n');
end
if zScore
    fprintf(fid,'All THREE metrics are z-scored using existing full-recording within-fly bar-only references. Absolute yaw is computed before normalization. No additional per-trial centering. All panels share the same z-score limits, retaining all fly-mean values. Means and SEM are calculated after normalization, across contributing fly means.\n');
end
fclose(fid);
result=struct('stem',stem,'recordings',height(F),'estimatedFlies',nf, ...
    'contributingFlies',contributingFlies,'eventCounts',sum(counts,1));disp(result);
if expandedWindow,fprintf('Prior counts: ');disp(sum(originalCounts,1));end
end
