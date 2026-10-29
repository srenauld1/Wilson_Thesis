function plot_ves041_rectangular_heat_summary(out)
% Plot saved descriptive results without re-detecting events or changing data.
% Bars are trial-specific; recordings are not assumed independent flies.
% Peak fractions are not adjusted for chance or matched for movement.
if nargin<1,out=fullfile(imaging_analysis_artifact_root(),'diagnostics','rectangular_heat_metrics','hand_drawn');end
S=readtable(fullfile(out,'event_metrics.csv'),'TextType','string');
S=S(S.reorientation==1,:);
old=get(groot,'defaultFigureVisible');
c=onCleanup(@()set(groot,'defaultFigureVisible',old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
fig=figure('Visible','off','Color','w','Position',[30 30 1550 1050]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','loose');
contexts=["heat","control","no heat"];
colors=[.88 .27 .16;.20 .49 .76;.50 .50 .50];
ids={ ["20261019_5_1","20261019_7_1","20261019_8_2"], ...
      ["20261020_11_2","20261020_12_1","20261020_13_1"] };
labels={{'5-1: barrier','7-1: stripes','8-2: stripes'}, ...
        {'11-2: barrier','12-1: stripes','13-1: stripes'}};
for day=1:2
    counts=zeros(3,3);hits=counts;fraction=nan(3,3);delta=fraction;
    for j=1:3
        for k=1:3
            q=S(S.recording==ids{day}(j)&S.context==contexts(k),:);
            assert(height(q)==1);
            counts(j,k)=q.nEvents;hits(j,k)=q.nPeakWithin1s;
            fraction(j,k)=q.percentPeakWithin1s;delta(j,k)=q.meanDeltaDff;
        end
    end
    ax=nexttile(tl,day);hold(ax,'on');b=bar(ax,fraction,'grouped');
    for k=1:3
        b(k).FaceColor=colors(k,:);b(k).EdgeColor='none';
    end
    for k=1:3
        for j=1:3
            if counts(j,k)>0
                text(ax,b(k).XEndPoints(j),fraction(j,k)+3, ...
                    sprintf('%d/%d',hits(j,k),counts(j,k)), ...
                    'HorizontalAlignment','center','FontSize',11);
            end
        end
    end
    set(ax,'XTick',1:3,'XTickLabel',labels{day},'FontSize',12,'Box','off');
    ylim(ax,[0 72]);yticks(ax,0:20:60);ylabel(ax,'Events with a VES041 peak (%)');
    title(ax,sprintf('10/%d | peak within +/-1 s of slowdown',18+day));
    if day==1,legend(ax,b,{'Heat','Control barrier','No heat'},'Location','northoutside','Orientation','horizontal');end
    grid(ax,'on');ax.XGrid='off';
    ax=nexttile(tl,day+2);hold(ax,'on');b=bar(ax,delta,'grouped');
    for k=1:3
        b(k).FaceColor=colors(k,:);b(k).EdgeColor='none';
        for j=1:3
            if counts(j,k)>0
                y=delta(j,k);offset=.035;if y<0,offset=-.045;end
                text(ax,b(k).XEndPoints(j),y+offset,sprintf('%+.2f',y), ...
                    'HorizontalAlignment','center','FontSize',11);
            end
        end
    end
    yline(ax,0,'k-','HandleVisibility','off');ylim(ax,[-.58 .88]);
    set(ax,'XTick',1:3,'XTickLabel',labels{day},'FontSize',12,'Box','off');
    ylabel(ax,'Mean event-related dF/F change');
    title(ax,'Mean 0 to 1 s minus mean -2 to -1 s');
    grid(ax,'on');ax.XGrid='off';
end
title(tl,{'VES041 during slowdowns with turning: heat versus control', ...
    'Full hand-drawn ROIs | rectangular recordings | 15 s local-baseline dF/F'},'FontSize',19);
xlabel(tl,{'Turn: >=30 deg accumulated absolute yaw within +/-0.5 s of slowdown; counts above bars = peak-associated / all events.', ...
    'Peak and slowdown must share the same context. Missing bars = no qualifying events. Descriptive, not behavior-matched or chance-adjusted.'},'FontSize',11);
savefig_visible_on_open(fig,fullfile(out,'ves041_rectangular_heat_summary.fig'));
exportgraphics(fig,fullfile(out,'ves041_rectangular_heat_summary.png'),'Resolution',180);
close(fig);
end
