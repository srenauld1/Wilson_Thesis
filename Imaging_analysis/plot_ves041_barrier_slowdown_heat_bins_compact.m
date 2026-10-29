function plot_ves041_barrier_slowdown_heat_bins_compact(out)
% Read saved event segments; make a compact physiology/kinematics display.
% No event selection or normalization is recomputed here.
if nargin<1,out=fullfile(imaging_analysis_artifact_root(),'diagnostics','barrier_slowdown_heat_bins_20261027');end
old=get(groot,'defaultFigureVisible');c=onCleanup(@()set(groot,'defaultFigureVisible',old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
R=readtable(fullfile(out,'run_report.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');
colors=[.25 .45 .65;.75 .60 .12;.95 .45 .10;.78 .12 .12];
names={'Control barrier','Heat trial: 0 to <1 V','Heat trial: 1 to <2.5 V','Heat trial: >=2.5 V'};
labs={'Baseline-subtracted dF/F','Forward speed (mm/s)','Rotational speed (deg/s)','Measured heat command (V)'};
for k=find(R.status=="success")'
    dest=char(R.output(k));a=load(fullfile(dest,'segmented_events.mat'));E=a.E;X=a.X;rel=a.rel;
    for geometry=unique(E.geometry(E.included),'stable')'
        take=E.included&E.geometry==geometry;
        fig=figure('Visible','off','Color','w','Position',[20 20 1900 1200]);tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
        for row=1:4
            v=[2 3 4 6];v=v(row);values=X(:,take,v);values=values(isfinite(values));lim=[0 1];
            if ~isempty(values),lim=[min(values) max(values)];if diff(lim)==0,lim=lim+[-.1 .1];end;lim=lim+[-1 1]*.05*diff(lim);end
            for g=1:4
                ax=nexttile(tl,(row-1)*4+g);hold(ax,'on');sel=take&E.group==g;
                if any(sel)
                    plot(ax,rel,X(:,sel,v),'Color',.8+.2*colors(g,:),'LineWidth',.6);
                    plot(ax,rel,mean(X(:,sel,v),2,'omitnan'),'Color',colors(g,:),'LineWidth',2.5);
                else,text(ax,.5,.5,'No eligible slowdowns','Units','normalized','HorizontalAlignment','center');end
                xline(ax,0,':k');xlim(ax,[-2 3]);ylim(ax,lim);grid(ax,'on');ylabel(ax,labs{row});
                if row==1,title(ax,{names{g},sprintf('n=%d | %d with turn, %d without',sum(sel),sum(E.reorientation(sel)),sum(~E.reorientation(sel)))});end
                if row==4,xlabel(ax,'Time from strongest deceleration (s)');end
            end
        end
        title(tl,{R.recording(k)+" | "+R.protocol(k)+" | "+geometry, ...
            "Hand-drawn ROI | "+string(a.settings.normalization)+" | thin: individual slowdowns; thick: mean"},'Interpreter','none','FontSize',17);
        xlabel(tl,{'Slowdowns with and without turns pooled. Invalid behavior and windows crossing barrier trials/geometries excluded.', ...
            'Bins use median measured command within +/-0.5 s; not measured temperature. No behavior matching across columns.'},'FontSize',11);
        tl.Units='normalized';tl.OuterPosition=[.01 .035 .98 .93];drawnow;
        stem=[matlab.lang.makeValidName(char(geometry)) '_slowdown_compact'];
        savefig_visible_on_open(fig,fullfile(dest,[stem '.fig']));
        drawnow;
        exportgraphics(fig,fullfile(dest,[stem '.png']),'Resolution',155);close(fig);
    end
end
end
