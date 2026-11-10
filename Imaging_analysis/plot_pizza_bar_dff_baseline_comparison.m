function plot_pizza_bar_dff_baseline_comparison()
% Same all-bar slowdown+rotation events: original vs pre-event-centered dF/F.
% One thin line per fly, equal-fly mean +/- SEM, matching physical y-scale.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
folder='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging/bar';
base=fullfile(folder,'ves041_bar_pooled_presentation_four_selections_fly_sem_minus3_plus3');
A=load([base '.mat']);B=load([base '_dff_event_baseline.mat']);
assert(isequal(A.counts,B.counts)&&isequal(A.rel,B.rel));
t=A.rel(:);S={A,B};titles={'Original dF/F','Baseline-subtracted dF/F'};
values=[];
for k=1:2
    q=S{k};v=[reshape(q.flyMeans(:,:,1,2),[],1);q.groupMean(:,1,2)-q.sem(:,1,2);q.groupMean(:,1,2)+q.sem(:,1,2)];
    values=[values;v(isfinite(v))]; %#ok<AGROW>
end
limits=[min(values) max(values)];limits=limits+[-.07 .07]*diff(limits);
fig=figure('Visible','off','Color','w','Position',[20 20 1250 520]);
tl=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','loose');
for k=1:2
    q=S{k};ax=nexttile(tl);hold(ax,'on');X=q.flyMeans(:,:,1,2);
    plot(ax,t,X,'Color',[.9 .6 .6],'LineWidth',1.1);
    mu=q.groupMean(:,1,2);se=q.sem(:,1,2);
    fill(ax,[t;flipud(t)],[mu-se;flipud(mu+se)],[.85 .05 .05],'FaceAlpha',.2,'EdgeColor','none');
    plot(ax,t,mu,'Color',[.85 .05 .05],'LineWidth',3);
    xline(ax,0,':','Color',[.5 .5 .5]);yline(ax,0,':','Color',[.5 .5 .5]);
    xlim(ax,[-3 3]);ylim(ax,limits);grid(ax,'on');ax.GridAlpha=.12;ax.FontSize=18;ax.Box='off';
    xlabel(ax,'Time (s)','FontSize',21);ylabel(ax,'dF/F','FontSize',21);title(ax,titles{k},'FontSize',23);
end
title(tl,'All bar trials | slowdowns + rotation','FontSize',24);
subtitle(tl,'Same events and flies | Thin: fly means; thick: mean +/- SEM | Baseline: -3 to -2 s','FontSize',16);
stem=fullfile(folder,'ves041_bar_slowdowns_rotation_original_vs_baseline_dff');
exportgraphics(fig,[stem '.png'],'Resolution',180);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
