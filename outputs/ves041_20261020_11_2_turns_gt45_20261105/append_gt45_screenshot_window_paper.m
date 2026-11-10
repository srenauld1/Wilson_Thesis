function append_gt45_screenshot_window_paper(out,lims)
% Second paper page: original screenshot context and +/-0.5 s turn window.
source='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/20261020_11_2/presentation_control_vs_heat_20261102/original_criteria_turns/comparison_data.mat';
S=load(source);rel=S.rel(:);selected=S.selected;X=cell(1,2);counts=zeros(1,2);mu=nan(numel(rel),3,2);sem=mu;names=["control","heat"];
for phase=1:2
 E=selected(selected.context==names(phase),:);keep=E.absoluteTurnDeg>45;
 X{phase}=S.X{phase}(:,keep,:);counts(phase)=nnz(keep);mu(:,:,phase)=squeeze(mean(X{phase},2));sem(:,:,phase)=squeeze(std(X{phase},0,2))/sqrt(counts(phase));
end
assert(isequal(counts,[7 12]));selected=selected(selected.absoluteTurnDeg>45,:);
writetable(selected,fullfile(out,'screenshot_window_selected_turns_gt45.csv'));
fig=figure('Visible','off','Color','white','Units','inches','Position',[1 1 7.2 8.2]);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 7.2 8.2],'PaperSize',[7.2 8.2],'PaperPositionMode','manual');
positions=[.14 .70 .35 .20;.60 .70 .35 .20;.14 .43 .35 .20;.60 .43 .35 .20;.14 .16 .35 .20;.60 .16 .35 .20];
colors={[.86 .06 .06],[.08 .08 .08],[.05 .2 .8]};labels={'\DeltaF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};phaseNames={'Placed control','Heat'};
for row=1:3
 for phase=1:2
  n=(row-1)*2+phase;ax=axes(fig,'Units','normalized','Position',positions(n,:),'FontName','Arial','FontSize',10,'Box','off','TickDir','out','LineWidth',.6,'XLim',[-2 3],'YLim',lims{row},'XTick',-2:1:3);hold(ax,'on');
  m=mu(:,row,phase);se=sem(:,row,phase);fill(ax,[rel;flipud(rel)],[m-se;flipud(m+se)],colors{row},'FaceAlpha',.18,'EdgeColor','none');
  plot(ax,rel,m,'Color',colors{row},'LineWidth',1.5);xline(ax,0,':','Color',[.5 .5 .5],'LineWidth',.7);
  if phase==1,ylabel(ax,labels{row},'FontSize',11);end
  if row==1,title(ax,sprintf('%s (n = %d)',phaseNames{phase},counts(phase)),'FontSize',12,'FontWeight','bold');end
  if row==3,xlabel(ax,'Time from slowdown (s)','FontSize',10);else,ax.XTickLabel=[];end
 end
 annotation(fig,'textbox',[.02 positions((row-1)*2+1,2)+.19 .05 .035],'String',char('A'+row-1),'FontName','Arial','FontSize',13,'FontWeight','bold','EdgeColor','none','Margin',0);
end
annotation(fig,'textbox',[.14 .955 .82 .025],'String','VES041 | 2026-10-20, fly 1 | turns >45 deg','FontName','Arial','FontSize',12,'FontWeight','bold','EdgeColor','none','Margin',0);
annotation(fig,'textbox',[.14 .02 .82 .055],'String',{'Screenshot selection: accumulated yaw over -0.5 to +0.5 s.','Original control/heat contexts; mean +/- event SEM; one fly.'},'FontName','Arial','FontSize',8,'EdgeColor','none','Margin',0);
file=fullfile(tempdir,'ves041_gt45_screenshot_page.pdf');print(fig,file,'-dpdf','-painters');
python='/Users/sophiarenauld/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3';
merger=fullfile(fileparts(mfilename('fullpath')),'merge_gt45_paper.py');
[status,message]=system(sprintf('"%s" "%s" "%s" "%s"',python,merger,fullfile(out,'20261020_11_2_turns_gt45_paper.pdf'),file));
assert(status==0,'Could not combine paper pages: %s',message);
stem=fullfile(out,'20261020_11_2_turns_gt45_screenshot_window_paper');print(fig,[stem '.svg'],'-dsvg','-painters');print(fig,[stem '.png'],'-dpng','-r600');
set(fig,'Visible','on');savefig(fig,[stem '.fig']);close(fig);
save(fullfile(out,'analysis_gt45_screenshot_window.mat'),'X','counts','mu','sem','selected','rel','lims','source','-v7.3');
fprintf('Screenshot-window >45 deg: %d control, %d heat events.\n',counts);
end
