function report=run_pizza_bar_no_accel_comparison()
% Compare an example trial and equal-fly crunch without an acceleration gate.
% Columns 1/3/4 preserve the saved VES peak and rotation-peak selections.
% Column 2 is regenerated using pizza_ves041_slowdowns_no_accel_threshold:
% velocity before>=1, drop>=1 mm/s, accumulated yaw>=15 deg over [0,.5] s;
% NO acceleration threshold. Max deceleration is used only for alignment.
% All groups require complete valid -3:+3 s windows. Original units/traces;
% no baseline subtraction, z-scoring or new smoothing in displayed traces.
% Example figure: one recording, lines=events, SEM across events.
% Crunch: pool events within fly, lines=fly means, SEM across fly means.
% Same screened bar-only 9-recording/7-estimated-fly cohort; no dark/heat.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';
out=fullfile(root,'bar','no_acceleration_threshold_20261102');
if ~isfolder(out),mkdir(out);end
P=load(fullfile(root,'bar','ves041_bar_pooled_equal_fly_peaks_slowdowns.mat'));
F=P.F;flies=P.flies;nf=numel(flies);rel=linspace(-3,3,121)';columns=[1 4 5 6];
pooled=cell(nf,4);counts=zeros(nf,4);oldSlowCounts=zeros(nf,1);
trialRows=cell(0,7);eventRows=cell(0,9);exampleX=[];
exampleTrial="20250812-22_ves041_fly3_bbarjumps2";
for k=1:height(F)
    trial=F.trial(k);f=find(flies==F.flyId(k));
    S=load(fullfile(root,'analysis_data',trial+"_events_with_rotation_peaks_and_speedups.mat"),'groups');
    A=load(fullfile(root,'analysis_data',trial+"_events.mat"),'d');d=A.d;
    [newIdx,angle,before,after,detector]=pizza_ves041_slowdowns_no_accel_threshold(d);
    withRotation=angle>=15;newIdx=newIdx(withRotation);angle=angle(withRotation);
    before=before(withRotation);after=after(withRotation);
    valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
    fs=1/median(diff(d.t));vars={d.dff,d.fwd,abs(d.yaw)};X=cell(1,4);n=zeros(1,4);
    oldIdx=S.groups{4};oldValid=0;
    for j=oldIdx(:)'
        if complete(d,valid,j,fs),oldValid=oldValid+1;end
    end
    oldSlowCounts(f)=oldSlowCounts(f)+oldValid;
    for c=1:4
        idx=S.groups{columns(c)};if c==2,idx=newIdx;end
        keep=false(numel(idx),1);Y=nan(numel(rel),numel(idx),3);
        for j=1:numel(idx)
            i=idx(j);if ~complete(d,valid,i,fs),continue,end
            for v=1:3,Y(:,j,v)=interp1(d.t,vars{v},d.t(i)+rel,'linear',NaN);end
            keep(j)=all(isfinite(Y(:,j,:)),'all');
        end
        X{c}=Y(:,keep,:);n(c)=nnz(keep);counts(f,c)=counts(f,c)+n(c);
        pooled{f,c}=cat(2,pooled{f,c},X{c});
        if c==2
            assert(all(keep),'New detector should ensure extended-window validity');
            for j=1:numel(idx)
                eventRows(end+1,:)={trial,flies(f),idx(j),d.t(idx(j)),before(j),after(j), ...
                    before(j)-after(j),angle(j),ismember(idx(j),oldIdx)}; %#ok<AGROW>
            end
        end
    end
    trialRows(end+1,:)={trial,flies(f),n(1),oldValid,n(2),n(3),n(4)}; %#ok<AGROW>
    if trial==exampleTrial,exampleX=X;end
    fprintf('%s | prior %d, no-acceleration-gate %d slowing+rotation\n',trial,oldValid,n(2));
end
assert(~isempty(exampleX));
flyMeans=nan(numel(rel),nf,3,4);
for f=1:nf,for c=1:4
    if counts(f,c)>0,flyMeans(:,f,:,c)=mean(pooled{f,c},2);end
end,end
% Verify that the three other groups and their fly means are untouched.
Q=load(fullfile(root,'bar','ves041_bar_pooled_presentation_four_selections_fly_sem_minus3_plus3.mat'));
assert(isequal(counts(:,[1 3 4]),Q.counts(:,[1 3 4])));
delta=abs(flyMeans(:,:,:,[1 3 4])-Q.flyMeans(:,:,:,[1 3 4]));assert(all(delta(isfinite(delta))<1e-10));
cohortX=cell(1,4);for c=1:4,cohortX{c}=flyMeans(:,:,:,c);end
exampleStats=render(exampleX,rel,fullfile(out,'20250812_22_example_no_accel'),false);
cohortStats=render(cohortX,rel,fullfile(out,'bar_equal_fly_no_accel'),true);
report=cell2table(trialRows,'VariableNames',{'trial','flyId','vesPeaks','priorSlowingRotation', ...
    'noAccelSlowingRotation','rotationDropBelow1','rotationDropAtLeast1'});
writetable(report,fullfile(out,'trial_counts.csv'));
events=cell2table(eventRows,'VariableNames',{'trial','flyId','sampleIndex','eventTimeSec', ...
    'velocityBefore','velocityAfter','velocityDrop','accumulatedRotationDeg','sameCenterAsPrior'});
writetable(events,fullfile(out,'no_accel_slowing_rotation_events.csv'));
methods=struct('detector',detector,'onlyRedetectedColumn',2,'sourceColumns',columns, ...
    'example','Individual events in 20250812-22; SEM across events', ...
    'crunch','Pool events within fly, then mean and SEM across equally weighted contributing flies', ...
    'normalization','none; original saved units','inferredFlyIdentities',3);
save(fullfile(out,'comparison.mat'),'F','flies','rel','counts','oldSlowCounts','flyMeans', ...
    'exampleX','exampleStats','cohortStats','report','events','methods');
fid=fopen(fullfile(out,'methods.txt'),'w');assert(fid>=0);
fprintf(fid,'No acceleration cutoff. Candidates are contiguous runs with pre median forward velocity >=1 mm/s and forward drop >=1 mm/s. Time zero is max deceleration within each velocity-drop run, regardless of magnitude. Centers >=2 s apart. 0.25-s median smoothing only for behavioral selection.\n');
fprintf(fid,'Before: -0.75:-0.25 s; after: +0.25:+0.75 s. Keep accumulated absolute yaw >=15 deg over 0:+0.5 s. Full valid -3:+3 s windows required. No fluorescence-based amplitude selection. Other three event selections unchanged.\n');
fprintf(fid,'Velocity-drop runs replace acceleration-threshold runs; episode boundaries can change, so this is not necessarily a strict superset of previous event centers. No baseline subtraction, z-scoring, or new display smoothing.\n');
fprintf(fid,'Example: 20250812-22, lines and SEM across events. Crunch: 9 screened bar recordings grouped into 7 estimated flies, three identities inferred. Events pooled within fly, then equal fly weight; SEM across fly means. Older outputs preserved.\n');fclose(fid);
disp(report);fprintf('Old slowing+rotation count: %d; new: %d\n',sum(oldSlowCounts),sum(counts(:,2)));
end

function yes=complete(d,valid,i,fs)
lo=find(d.t<=d.t(i)-3,1,'last');hi=find(d.t>=d.t(i)+3,1,'first');
yes=~isempty(lo)&&~isempty(hi);
if yes,yes=all(valid(lo:hi))&&all(diff(d.t(lo:hi))<=2/fs);end
end

function stats=render(X,t,stem,pooled)
titles={'VES041 peaks','Slowdowns + rotation',{'Rotation peaks','Speed drop <1 mm/s'}, ...
    {'Rotation peaks','Speed drop \geq1 mm/s'}};
labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};
thin={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};
fig=figure('Visible','off','Color','w','Position',[20 20 1800 1050]);
tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','loose');stats=cell(3,4);
for v=1:3
    values=[];
    for c=1:4
        A=X{c}(:,:,v);n=sum(isfinite(A),2);mu=mean(A,2,'omitnan');se=std(A,0,2,'omitnan')./sqrt(n);se(n<2)=NaN;
        stats{v,c}=struct('mean',mu,'sem',se,'n',n);
        q=[A(:);mu-se;mu+se];values=[values;q(isfinite(q))]; %#ok<AGROW>
    end
    limits=[min(values) max(values)];span=max(diff(limits),eps);limits=limits+[-.06 .06]*span;
    for c=1:4
        ax=nexttile(tl);hold(ax,'on');A=X{c}(:,:,v);plot(ax,t,A,'Color',thin{v},'LineWidth',.8);
        mu=stats{v,c}.mean;se=stats{v,c}.sem;
        fill(ax,[t;flipud(t)],[mu-se;flipud(mu+se)],colors{v},'FaceAlpha',.2,'EdgeColor','none');
        plot(ax,t,mu,'Color',colors{v},'LineWidth',3);xline(ax,0,':','Color',[.45 .45 .45]);
        ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;
        xlim(ax,[-3 3]);ylim(ax,limits);ax.XTick=-3:3;
        if c==1,ylabel(ax,labels{v},'FontSize',21);end
        if v==3,xlabel(ax,'Time (s)','FontSize',21);else,ax.XTickLabel=[];end
        if v==1,title(ax,titles{c},'FontSize',23);end
    end
end
if pooled
    caption='Bar trials | Thin: fly means | Thick: equal-fly mean +/- SEM | No acceleration threshold';
else,caption='20250812-22 | Thin: individual events | Thick: mean +/- SEM across events | No acceleration threshold';end
subtitle(tl,caption,'FontSize',18);
assert(numel(findall(fig,'Type','axes'))==12);
exportgraphics(fig,[stem '.png'],'Resolution',200);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
