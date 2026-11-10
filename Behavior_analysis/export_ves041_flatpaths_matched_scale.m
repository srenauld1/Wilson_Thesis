function export_ves041_flatpaths_common_axes()
% Identical square plotting areas and limits; original trajectory coordinates.
root=fullfile(pwd,'outputs','ves041_fly8_fly45_master_rerun_20261102');
sources={fullfile(root,'fly8_full_recording_clean.fig'), ...
    fullfile(root,'fly45_full_recording_clean.fig'), ...
 '/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers/fly019/cell01_processed/trial01Astim2/flatpath_vs_opto.fig', ...
 '/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers/fly05/cell01_processed/trial01Astim2/flatpath_vs_opto.fig'};
names={'VES041_fly008_ACR_20260108_trial01B_stim5_flatpath', ...
    'VES041_fly045_GFP_20260224_trial01B_stim5_flatpath', ...
 'Chrimson_fly019_VT_csChrimson_20260313_trial01A_stim2_flatpath', ...
 'Chrimson_fly005_VT_GFP_20260306_trial01A_stim2_flatpath'};
out=fullfile(pwd,'outputs','four_flatpaths_common_axes_20261103');
preview=out;
if ~isfolder(out),mkdir(out);end
figs=gobjects(1,4);xs=[];ys=[];lines=cell(1,4);oldX=cell(1,4);oldY=cell(1,4);
for i=1:4
    figs(i)=openfig(sources{i},'invisible');
    lines{i}=findall(figs(i),'Type','line');
    oldX{i}=get(lines{i},'XData');oldY{i}=get(lines{i},'YData');
    for j=1:numel(lines{i})
        xs=[xs,reshape(lines{i}(j).XData,1,[])]; %#ok<AGROW>
        ys=[ys,reshape(lines{i}(j).YData,1,[])]; %#ok<AGROW>
    end
end
xs=xs(isfinite(xs));ys=ys(isfinite(ys));
% Round outward and use the same square bounds for all four recordings.
xl=[floor((min(xs)-30)/100)*100 ceil((max(xs)+30)/100)*100];
yl=[floor((min(ys)-30)/100)*100 ceil((max(ys)+30)/100)*100];
span=max(diff(xl),diff(yl));
xl(2)=xl(1)+span;yl(2)=yl(1)+span;
assert(min(xs)>=xl(1)&&max(xs)<=xl(2)&&min(ys)>=yl(1)&&max(ys)<=yl(2));
for i=1:4
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
    assert(isequal(ax.XLim,xl)&&isequal(ax.YLim,yl));
    assert(isequal(ax.Position,[.14 .14 .78 .78]));
    assert(isequaln(oldX{i},get(lines{i},'XData'))&&isequaln(oldY{i},get(lines{i},'YData')));
    dest=fullfile(out,[names{i} '_common_axes.svg']);
    assert(~isfile(dest),'Destination already exists: %s',dest);
    print(fig,dest,'-dsvg','-painters');
    print(fig,fullfile(preview,[names{i} '_common_axes.png']),'-dpng','-r120');
    assert(~contains(fileread(dest),'<image'));
    fprintf('SAVED %s\n',dest);
    close(fig);
end
fprintf('All four plots: X=[%g %g], Y=[%g %g], span=%g; identical 10-inch square canvas.\n',xl,yl,span);
end
