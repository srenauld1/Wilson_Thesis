function export_ves041_example_flatpaths(root)
% Presentation exports of master-generated figures; no trajectory edits.
trials = {'2026_01_08_fly02_cell01_trial01B','2026_02_24_fly03_cell01_trial01B'};
labels = {'Fly 8 — VES041_1 x ACR','Fly 45 — VES041_1 x GFP'};
ids = [8 45];
for k = 1:2
    folder = fullfile(root,trials{k},'trial01Bstim5');
    for aligned = [false true]
        if aligned
            name = 'flatpath_aligned_2s_rawY_vs_opto';
            suffix = 'aligned';
        else
            name = 'flatpath_vs_opto';
            suffix = 'full_recording';
        end
        fig = openfig(fullfile(folder,[name '.fig']),'invisible');
        delete(findall(fig,'Type','textboxshape'));
        ax = findobj(fig,'Type','axes');
        delete(findobj(ax,'Type','text')); % removes heading-stat annotation only
        set(fig,'Color','white','Position',[100 100 1000 850]);
        set(ax,'FontName','Arial','FontSize',16,'Box','off');
        title(ax,labels{k},'Interpreter','none','FontSize',20);
        if aligned
            xlabel(ax,'Lateral position (mm)'); ylabel(ax,'Forward position (mm)');
            subtitle(ax,'All 40 presentations; 2 s pre / 5 s opto / 2 s post; no mirroring');
        else
            xlabel(ax,'X position (mm)'); ylabel(ax,'Y position (mm)');
            subtitle(ax,'Full recording; black: opto-off; red: opto-on');
        end
        axis(ax,'equal'); grid(ax,'off');
        set(ax,'Units','normalized','Position',[0.15 0.14 0.78 0.68]);
        lg = findobj(fig,'Type','legend');
        if ~isempty(lg), set(lg,'Location','northeast','FontSize',12); end
        file = fullfile(root,sprintf('fly%d_%s_clean',ids(k),suffix));
        exportgraphics(fig,[file '.png'],'Resolution',220);
        print(fig,[file '.svg'],'-dsvg');
        savefig(fig,[file '.fig']);
        close(fig);
    end
end
end
