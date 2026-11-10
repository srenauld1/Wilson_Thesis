function export_chrimson_flatpaths_matched_scale()
% Identical square plotting areas and limits; original trajectory coordinates.
folders={ ...
 '/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers/fly019/cell01_processed/trial01Astim2', ...
 '/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers/fly05/cell01_processed/trial01Astim2'};
names={'Chrimson_fly019_VT_csChrimson_20260313_trial01A_stim2_flatpath', ...
 'Chrimson_fly005_VT_GFP_20260306_trial01A_stim2_flatpath'};
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples';
preview=fullfile(pwd,'outputs','chrimson_example_svg_exports_20261103');
figs=gobjects(1,2);xs=[];ys=[];lines=cell(1,2);oldX=cell(1,2);oldY=cell(1,2);
for i=1:2
    figs(i)=openfig(fullfile(folders{i},'flatpath_vs_opto.fig'),'invisible');
    lines{i}=findall(figs(i),'Type','line');
    oldX{i}=get(lines{i},'XData');oldY{i}=get(lines{i},'YData');
    for j=1:numel(lines{i})
        xs=[xs,reshape(lines{i}(j).XData,1,[])]; %#ok<AGROW>
        ys=[ys,reshape(lines{i}(j).YData,1,[])]; %#ok<AGROW>
    end
end
xs=xs(isfinite(xs));ys=ys(isfinite(ys));
span=ceil((max(range(xs),range(ys))+60)/100)*100;
cx=round((min(xs)+max(xs))/100)*50;
cy=round((min(ys)+max(ys))/100)*50;
xl=cx+[-.5 .5]*span;yl=cy+[-.5 .5]*span;
assert(min(xs)>=xl(1)&&max(xs)<=xl(2)&&min(ys)>=yl(1)&&max(ys)<=yl(2));
for i=1:2
    fig=figs(i);
    set(fig,'Units','pixels','Position',[100 100 1000 1000], ...
        'PaperUnits','inches','PaperPosition',[0 0 10 10], ...
        'PaperSize',[10 10],'PaperPositionMode','manual','Color','white');
    ax=findobj(fig,'Type','axes');assert(numel(ax)==1);
    set(ax,'Units','normalized','Position',[.14 .14 .78 .78], ...
        'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1], ...
        'XLim',xl,'YLim',yl,'XTick',xl(1):100:xl(2),'YTick',yl(1):100:yl(2));
    lg=findobj(fig,'Type','legend');set(lg,'Location','northeast');
    drawnow;
    assert(isequaln(oldX{i},get(lines{i},'XData'))&&isequaln(oldY{i},get(lines{i},'YData')));
    dest=fullfile(out,[names{i} '_matched_scale.svg']);
    assert(~isfile(dest),'Destination already exists: %s',dest);
    print(fig,dest,'-dsvg','-painters');
    print(fig,fullfile(preview,[names{i} '_matched_scale.png']),'-dpng','-r120');
    assert(~contains(fileread(dest),'<image'));
    fprintf('SAVED %s\n',dest);
    close(fig);
end
fprintf('Both plots: X=[%g %g], Y=[%g %g], span=%g; identical 10-inch square canvas.\n',xl,yl,span);
end
