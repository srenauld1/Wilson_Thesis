function result=plot_pizza_bar_dff_event_baseline()
% Baseline-subtract ONLY dF/F in the existing -3:+3 bar-cohort figure.
% Baseline is [-3,-2] s for each event. For the same event set, averaging
% (event minus event baseline) is algebraically identical to subtracting
% the baseline of its fly-mean trace, so use saved fly means without any
% redetection, changed eligibility, altered weighting or interpolation.
% Recompute SEM across centered FLY means. Kinematics remain unchanged.
% Save full 4-column figure AND the attached slowdown+rotation single panel.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
folder='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging/bar';
source=fullfile(folder,'ves041_bar_pooled_presentation_four_selections_fly_sem_minus3_plus3');
S=load([source '.mat']);rel=S.rel(:);mask=rel>=-3&rel<=-2;
assert(any(mask)&&rel(1)==-3&&rel(end)==3);
flyMeans=S.flyMeans;groupMean=S.groupMean;sem=S.sem;baselineOffsets=nan(numel(S.flies),4);
for c=1:4
    for f=1:numel(S.flies)
        baselineOffsets(f,c)=mean(flyMeans(mask,f,1,c),'omitnan');
        flyMeans(:,f,1,c)=flyMeans(:,f,1,c)-baselineOffsets(f,c);
        if isfinite(baselineOffsets(f,c)),assert(abs(mean(flyMeans(mask,f,1,c)))<1e-10);end
    end
    X=flyMeans(:,:,1,c);n=sum(isfinite(X),2);
    groupMean(:,1,c)=mean(X,2,'omitnan');sem(:,1,c)=std(X,0,2,'omitnan')./sqrt(n);
end
assert(isequaln(flyMeans(:,:,2:3,:),S.flyMeans(:,:,2:3,:)));
assert(isequaln(groupMean(:,2:3,:),S.groupMean(:,2:3,:))&&isequaln(sem(:,2:3,:),S.sem(:,2:3,:)));
fig=openfig([source '.fig'],'invisible');axs=findall(fig,'Type','axes');ordered=gobjects(12,1);
for ax=axs',ordered(ax.Layout.Tile)=ax;end
values=[reshape(flyMeans(:,:,1,:),[],1);reshape(groupMean(:,1,:)-sem(:,1,:),[],1);reshape(groupMean(:,1,:)+sem(:,1,:),[],1)];
values=values(isfinite(values));limits=[min(values) max(values)];limits=limits+[-.06 .06]*diff(limits);
for c=1:4
    ax=ordered(c);lines=findall(ax,'Tag','flyMean');
    for line=lines'
        f=find(S.flies==string(line.UserData));assert(isscalar(f));line.YData=flyMeans(:,f,1,c);
    end
    m=findall(ax,'Tag','equalFlyMean');assert(isscalar(m));m.YData=groupMean(:,1,c);
    band=findall(ax,'Tag','flySem');assert(isscalar(band));
    band.YData=[groupMean(:,1,c)-sem(:,1,c);flipud(groupMean(:,1,c)+sem(:,1,c))];
    yline(ax,0,':','Color',[.45 .45 .45],'LineWidth',1);ylim(ax,limits);
    if c==1,ylabel(ax,'Baseline-subtracted dF/F','FontSize',21);end
end
tl=findall(fig,'Type','tiledlayout');
subtitle(tl,{'dF/F baseline: -3 to -2 s | Forward velocity and rotation unchanged', ...
    'Thin: fly means | Thick: equal-fly mean | Shading: SEM across flies'},'FontSize',18);
stem=[source '_dff_event_baseline'];export_plot(fig,stem);
panel=figure('Visible','off','Color','w','Position',[20 20 600 1000]);
layout=tiledlayout(panel,3,1,'TileSpacing','compact','Padding','loose');
for row=1:3
    old=ordered((row-1)*4+2);ax=copyobj(old,panel);ax.Parent=layout;ax.Layout.Tile=row;
    ax.FontSize=19;
    % The standalone panel need not retain the VES-peak column's larger
    % shared limits; include all its fly traces and SEM without clipping.
    vals=[reshape(flyMeans(:,:,row,2),[],1);groupMean(:,row,2)-sem(:,row,2);groupMean(:,row,2)+sem(:,row,2)];
    vals=vals(isfinite(vals));bounds=[min(vals) max(vals)];
    ylim(ax,bounds+[-.08 .08]*diff(bounds));
    if row==1,ylabel(ax,{'Baseline-subtracted','dF/F'},'FontSize',22);
    elseif row==2,ylabel(ax,{'Forward velocity','(mm/s)'},'FontSize',22);
    else,ylabel(ax,{'Rotational speed','(deg/s)'},'FontSize',22);end
end
panelStem=[stem '_slowdowns_rotation'];export_plot(panel,panelStem);
close(panel);close(fig);
methods=S.methods;
methods.normalization='dF/F only: subtract mean -3:-2 s per event (equivalent fly-mean baseline subtraction); kinematics unchanged';
methods.dffBaselineWindowSec=[-3 -2];
flies=S.flies;counts=S.counts;F=S.F;
save([stem '.mat'],'flyMeans','groupMean','sem','baselineOffsets','rel','flies','counts','F','methods');
result=struct('fullFigure',[stem '.fig'],'panelFigure',[panelStem '.fig'], ...
    'dffBaselineWindow',[-3 -2],'kinematicsUnchanged',true);disp(result);
end
function export_plot(fig,stem)
exportgraphics(fig,[stem '.png'],'Resolution',200);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig_visible_on_open(fig,[stem '.fig']);
end
