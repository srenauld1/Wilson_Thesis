function plot_ves041_sync_audit(out)
% Render read-only raw CSV/DAQ verification results without showing windows.
% audit_ves041_sync_raw.py checks every marker; the plots show raw start/end
% excerpts, full-recording timing, gaps, and independent motion validation.
if nargin<1,out=fullfile(imaging_analysis_artifact_root(),'diagnostics','sync_audit_20261026');end
old=get(groot,'defaultFigureVisible');c=onCleanup(@()set(groot,'defaultFigureVisible',old)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
files=dir(fullfile(out,'*_sync_audit.mat'));
overview=figure('Visible','off','Color','w','Position',[30 30 1700 1200]);
ov=tiledlayout(overview,4,4,'TileSpacing','compact','Padding','loose');
gapRows=cell(0,4);
for k=1:numel(files)
    a=load(fullfile(files(k).folder,files(k).name));
    for j=find(a.updateIntervals(:)>.15)'
        gapRows(end+1,:)={string(a.recording),a.intervalTimes(j)-a.updateIntervals(j),a.intervalTimes(j),a.updateIntervals(j)}; %#ok<AGROW>
    end
    fig=figure('Visible','off','Color','w','Position',[30 30 1700 1350]);
    tl=tiledlayout(fig,4,2,'TileSpacing','compact','Padding','loose');
    for col=1:2
        if col==1,prefix='start';else,prefix='end';end
        for bit=1:2
            if bit==1,signal='Fast';else,signal='Slow';end
            ax=nexttile(tl,(bit-1)*2+col);hold(ax,'on');
            stairs(ax,a.([prefix 'RawTime']),a.([prefix 'Raw' signal]),'k','LineWidth',1.2,'DisplayName','Raw DAQ');
            plot(ax,a.([prefix 'CsvTime']),a.([prefix 'Csv' signal]),'o','Color',[.85 .2 .1],'MarkerSize',4,'DisplayName','Matched CSV row');
            ylim(ax,[-.2 1.2]);yticks(ax,[0 1]);grid(ax,'on');xlabel(ax,'Time from first imaging frame (s)');
            ylabel(ax,[signal ' marker bit']);title(ax,[prefix ' of recording']);
            if col==1&&bit==1,legend(ax,'Location','best');end
        end
    end
    ax=nexttile(tl,5);plot(ax,a.intervalTimes,a.updateIntervals,'k.','MarkerSize',3);hold(ax,'on');
    yline(ax,.15,':r','Importer gap threshold');xlabel(ax,'Time (s)');ylabel(ax,'Update interval (s)');grid(ax,'on');title(ax,'All tracking updates; long intervals are missing behavior');
    ax=nexttile(tl,6);plot(ax,a.updateTimes,a.timingError*1e6,'b');grid(ax,'on');
    xlabel(ax,'Time (s)');ylabel(ax,'Raw minus saved time (microseconds)');title(ax,'Raw reconstruction versus A2P mapping');
    ax=nexttile(tl,7);plot(ax,a.forwardErrorTimes,a.forwardErrorRad,'k.','MarkerSize',3);hold(ax,'on');
    yline(ax,.05,':r');yline(ax,-.05,':r');grid(ax,'on');xlabel(ax,'Time (s)');ylabel(ax,'Wrapped forward-position error (rad)');title(ax,'Independent DAQ motion voltage versus CSV');
    ax=nexttile(tl,8);plot(ax,a.offsetShifts,100*a.offsetMismatchFraction,'-o','MarkerSize',3);grid(ax,'on');
    xlabel(ax,'Shift from recovered CSV offset (rows)');ylabel(ax,'Updates with a marker mismatch (%)');title(ax,'Local offset check; exact matcher searches all offsets');
    title(tl,{strrep(a.recording,'_',' '),sprintf('%d matched updates | fast mismatches %d | slow mismatches %d',a.matchedUpdates,a.fastMismatches,a.slowMismatches)},'FontSize',17);
    note=a.notes;if isempty(note),note='No coverage flags';end
    xlabel(tl,{note,'Markers timestamp processed tracking updates, not camera exposure. Heat voltage is not temperature.'},'Interpreter','none','FontSize',10);
    savefig_visible_on_open(fig,fullfile(out,[a.recording '_sync_alignment.fig']));
    exportgraphics(fig,fullfile(out,[a.recording '_sync_alignment.png']),'Resolution',150);close(fig);
    ax=nexttile(ov,k);hold(ax,'on');
    stairs(ax,a.startRawTime,a.startRawSlow,'k','LineWidth',1);
    plot(ax,a.startCsvTime,a.startCsvSlow,'o','Color',[.85 .2 .1],'MarkerSize',3);
    ylim(ax,[-.2 1.2]);yticks(ax,[0 1]);grid(ax,'on');
    color=[0 0 0];if contains(string(a.notes),'gap')||contains(string(a.notes),'failed'),color=[.7 .35 0];end
    title(ax,{strrep(a.recording,'_',' '),sprintf('Bit errors: %d/%d | gaps: %d',a.fastMismatches,a.slowMismatches,sum(a.updateIntervals>.15))},'FontSize',11,'Color',color);
    xlabel(ax,'DAQ time (s)');ylabel(ax,'Slow bit');
end
title(ov,{'VES041 Oct 19/20: raw DAQ versus aligned CSV markers', ...
    'Black: DAQ slow code | red circles: matched CSV code | excerpts shown; every update checked'},'FontSize',18);
xlabel(ov,'Coverage gaps and incomplete acquisitions are reported separately; a marker match does not certify camera-exposure latency.','FontSize',11);
savefig_visible_on_open(overview,fullfile(out,'all_trials_sync_alignment.fig'));
exportgraphics(overview,fullfile(out,'all_trials_sync_alignment.png'),'Resolution',160);close(overview);
G=cell2table(gapRows,'VariableNames',{'recording','gapStartSec','gapEndSec','gapDurationSec'});
writetable(G,fullfile(out,'tracking_gap_intervals.csv'));
end
