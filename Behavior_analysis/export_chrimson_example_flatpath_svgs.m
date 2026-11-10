function export_chrimson_example_flatpath_svgs()
% Export the original two example FIGs without changing data or plot styling.
sources={ ...
 '/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers/fly019/cell01_processed/trial01Astim2/flatpath_vs_opto.fig', ...
 '/Volumes/Neurobio/wilsonlab/Sophia/behavior/Chrimson Data/keepers/fly05/cell01_processed/trial01Astim2/flatpath_vs_opto.fig'};
names={ ...
 'Chrimson_fly019_VT_csChrimson_20260313_trial01A_stim2_flatpath', ...
 'Chrimson_fly005_VT_GFP_20260306_trial01A_stim2_flatpath'};
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples';
preview=fullfile(pwd,'outputs','chrimson_example_svg_exports_20261103');
assert(isfolder(out),'Requested destination does not exist.');
if ~isfolder(preview),mkdir(preview);end
for i=1:numel(sources)
    destination=fullfile(out,[names{i} '.svg']);
    assert(~isfile(destination),'Refusing to overwrite an existing SVG: %s',destination);
    fig=openfig(sources{i},'invisible');
    cleanup=onCleanup(@()close(fig));
    lines=findall(fig,'Type','line');
    oldX=get(lines,'XData');oldY=get(lines,'YData');
    set(fig,'PaperPositionMode','auto');
    print(fig,destination,'-dsvg','-painters');
    print(fig,fullfile(preview,[names{i} '.png']),'-dpng','-r120');
    assert(isequaln(oldX,get(lines,'XData')) && isequaln(oldY,get(lines,'YData')),'Trajectory changed during export.');
    xml=fileread(destination);
    assert(~contains(xml,'<image'),'SVG contains an embedded raster image.');
    fprintf('VECTOR SVG SAVED: %s\n',destination);
    clear cleanup;
end
end
