function plot_pizza_0812_22_presentation()
% Presentation copy of four selections from the existing trace-only figure.
% No analysis or trace data are changed. Copy original axes and limits.
% Keep original columns 1,4,5,6 and every individual-event trace. Add pointwise
% mean +/- SEM across events (sample SD/sqrt(finite N)), NOT across flies.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
folder='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging/bar';
base='20250812-22_ves041_fly3_bbarjumps2_peaks_slowdowns';
src=openfig(fullfile(folder,[base '_traces_only_red_black_blue.fig']),'invisible');
cleanup=onCleanup(@()close(src)); %#ok<NASGU>
allAxes=findall(src,'Type','axes');old=gobjects(24,1);
for ax=allAxes',old(ax.Layout.Tile)=ax;end
columns=[1 4 5 6];
titles={'VES041 peaks','Slowdowns + rotation', ...
    {'Rotation peaks','Speed drop <1 mm/s'}, ...
    {'Rotation peaks','Speed drop >=1 mm/s'}};
labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};
fig=figure('Visible','off','Color','w','Position',[20 20 1800 1050]);
tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','loose');
semData=cell(3,4);eventCounts=zeros(3,4);expectedLines=0;
for row=1:3
    for c=1:4
        previous=old((row-1)*8+columns(c));
        ax=copyobj(previous,fig);ax.Parent=tl;ax.Layout.Tile=(row-1)*4+c;
        lines=findall(ax,'Type','line');
        % The saved mean is the thick line; retain its exact X/Y data. Every
        % other line is one original event. Verify before adding the SEM band.
        keep=[lines.LineWidth]>=1.5;assert(nnz(keep)==1);
        meanLine=lines(keep);eventLines=lines(~keep);meanLine.LineWidth=3;
        n=numel(eventLines);assert(n>1);eventCounts(row,c)=n;
        xx=meanLine.XData(:);X=nan(numel(xx),n);
        for j=1:n
            assert(isequal(eventLines(j).XData(:),xx));
            X(:,j)=eventLines(j).YData(:);eventLines(j).LineWidth=.5;
        end
        mu=mean(X,2,'omitnan');count=sum(isfinite(X),2);
        se=std(X,0,2,'omitnan')./sqrt(count);
        delta=abs(mu-meanLine.YData(:));assert(all(delta(isfinite(delta))<1e-10));
        hold(ax,'on');
        band=fill(ax,[xx;flipud(xx)],[mu-se;flipud(mu+se)],meanLine.Color, ...
            'FaceAlpha',.23,'EdgeColor','none','HandleVisibility','off');
        uistack(band,'top');uistack(meanLine,'top');
        semData{row,c}=struct('time',xx,'mean',mu,'sem',se,'n',count);
        expectedLines=expectedLines+n+1;
        ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';
        ax.GridAlpha=.12;ax.XMinorGrid='off';ax.YMinorGrid='off';
        ax.XTick=-2:1:3;
        if c==1,ylabel(ax,labels{row},'FontSize',21);else,ylabel(ax,'');end
        if row==3,xlabel(ax,'Time (s)','FontSize',21);else,xlabel(ax,'');ax.XTickLabel=[];end
        if row==1
            title(ax,titles{c},'FontSize',23,'FontWeight','bold');
        else,title(ax,'');end
        guides=findall(ax,'Type','constantline');
        for j=1:numel(guides),guides(j).Color=[.45 .45 .45];guides(j).LineWidth=1;end
        assert(isequal(ax.XLim,previous.XLim)&&isequal(ax.YLim,previous.YLim));
    end
end
% Keep metadata off the slide figure; it remains in the filename and source.
subtitle(tl,'Thin: individual events | Thick: mean | Shading: +/- SEM across events','FontSize',18);
assert(numel(findall(fig,'Type','axes'))==12);
assert(numel(findall(fig,'Type','line'))==expectedLines);
assert(numel(findall(fig,'Type','patch'))==12);
assert(isequal(eventCounts(1,:),eventCounts(2,:))&&isequal(eventCounts(1,:),eventCounts(3,:)));
stem=fullfile(folder,[base '_presentation_four_selections_with_events_sem']);
exportgraphics(fig,[stem '.png'],'Resolution',200);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
save([stem '_sem.mat'],'semData','eventCounts','columns');
fprintf('Saved %s\n',stem);
fprintf('Events per selection: ');disp(eventCounts(1,:));
end
