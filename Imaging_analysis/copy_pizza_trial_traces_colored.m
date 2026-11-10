function result=copy_pizza_trial_traces_colored(sourceStem)
%COPY_PIZZA_TRIAL_TRACES_COLORED Copy the existing plot without heatmap row.
% Copies all trace axes and their original X/Y data, limits and event groups.
% Dark trace colors: dF/F red, forward black, rotational speed blue.
% Individual events use lighter shades. No analysis, filtering or averaging
% is rerun. Original files unchanged. Run MATLAB -nodisplay.
if nargin<1
    sourceStem='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging/bar/20250812-22_ves041_fly3_bbarjumps2_peaks_slowdowns';
end
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
source=openfig([sourceStem '.fig'],'invisible');cleanup=onCleanup(@()close(source)); %#ok<NASGU>
axesOld=findall(source,'Type','axes');assert(numel(axesOld)==32,'Expected original 4x8 plot');
old=gobjects(32,1);
for a=axesOld',old(a.Layout.Tile)=a;end
fig=figure('Visible','off','Color','w','Position',[10 10 3300 1100]);
tl=tiledlayout(fig,3,8,'TileSpacing','compact','Padding','loose');
meanColors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};
thinColors={[.92 .73 .73],[.75 .75 .75],[.72 .77 .94]};
for row=1:3
    for col=1:8
        previous=old(row*8+col);
        ax=copyobj(previous,fig);ax.Parent=tl;ax.Layout.Tile=(row-1)*8+col;
        lines=findall(ax,'Type','line');originalLines=findall(previous,'Type','line');
        assert(numel(lines)==numel(originalLines));
        for k=1:numel(lines)
            assert(isequaln(lines(k).XData,originalLines(k).XData)&&isequaln(lines(k).YData,originalLines(k).YData));
            if lines(k).LineWidth>=1.5,lines(k).Color=meanColors{row};
            else,lines(k).Color=thinColors{row};end
        end
        % Keep the event-zero guide neutral so red specifically denotes dF/F.
        guides=findall(ax,'Type','constantline');
        for k=1:numel(guides),guides(k).Color=[.5 .5 .5];end
        if row==1,title(ax,old(col).Title.String,'Interpreter',old(col).Title.Interpreter);end
        assert(isequal(ax.XLim,previous.XLim)&&isequal(ax.YLim,previous.YLim));
    end
end
layout=findall(source,'Type','tiledlayout');
if ~isempty(layout)
    title(tl,layout(1).Title.String,'Interpreter',layout(1).Title.Interpreter);
end
subtitle(tl,'Same saved events and axes | dF/F: red; forward velocity: black; rotational speed: blue | heatmaps omitted');
assert(numel(findall(fig,'Type','axes'))==24);
assert(isempty(findall(fig,'Type','image'))&&isempty(findall(fig,'Type','colorbar')));
target=[sourceStem '_traces_only_red_black_blue'];
exportgraphics(fig,[target '.png'],'Resolution',170);
savefig_visible_on_open(fig,[target '.fig']);close(fig);
result=struct('png',[target '.png'],'fig',[target '.fig'],'dataUnchanged',true,'panels',24);disp(result);
end
