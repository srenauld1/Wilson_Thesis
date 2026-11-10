function plot_20261020_11_2_turns_gt45_paper(out)
% Strict >45-degree sensitivity of the linked PDF's original pizza rule.
% Retain all original slowdown candidates, timing, heat labels and traces.
if nargin<1,out=fullfile(pwd,'outputs','ves041_20261020_11_2_turns_gt45_20261105');end
if ~isfolder(out),mkdir(out);end
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier';
source=fullfile(root,'20261020_11_2','analysis.mat');A=load(source);
idx=A.groups{2}(:);angle=pizza_ves041_turn_angle(A.d,idx,[0 .5]);
old=angle>=15;selected=angle>45;
assert(isequal(idx(old),A.groups{4}(:)),'Original 15-degree group failed reconstruction.');
heat=A.d.heat(idx)>A.meta.threshold;countsOld=[sum(old&~heat),sum(old&heat)];
assert(isequal(countsOld,[20 19]),'Original PDF cohort changed.');
counts=[sum(selected&~heat),sum(selected&heat)];rel=A.rel(:);
X=cell(1,2);oldX=X;mu=nan(numel(rel),3,2);sem=mu;oldMu=mu;
for phase=1:2
 take=selected&(heat==(phase==2));X{phase}=A.X{2}(:,take,:);
 oldX{phase}=A.X{2}(:,old&(heat==(phase==2)),:);
 if counts(phase)>0,mu(:,:,phase)=squeeze(mean(X{phase},2));end
 if counts(phase)>1,sem(:,:,phase)=squeeze(std(X{phase},0,2))/sqrt(counts(phase));end
 oldMu(:,:,phase)=squeeze(mean(oldX{phase},2));
 assert(all(isfinite(X{phase}(:))));
end
% Audit all saved candidate slowdowns, not only the retained turns.
E=table(idx,A.d.t(idx),angle,logical(heat),logical(old),logical(selected), ...
 'VariableNames',{'sampleIndex','eventTime','absoluteTurnDeg_0_to_0p5s','heatOn','passesOld15deg','passesStrict45deg'});
writetable(E,fullfile(out,'all_slowdown_candidates.csv'));
writetable(E(selected,:),fullfile(out,'selected_turn_events_gt45.csv'));
phaseNames=["Heat off";"Heat on"];
T=table(phaseNames,countsOld(:),counts(:), ...
 'VariableNames',{'condition','originalEvents_ge15deg','newEvents_gt45deg'});
for phase=1:2
 if counts(phase)>0
  rot=X{phase}(:,:,3);core=rel>=0&rel<=.5;
  T.peakMeanRotation_0_to_0p5s(phase)=max(mean(rot(core,:),2));
  T.medianIndividualPeakRotation_0_to_0p5s(phase)=median(max(rot(core,:),[],1));
  T.medianTurnAngleDeg(phase)=median(angle(selected&(heat==(phase==2))));
 else
  T.peakMeanRotation_0_to_0p5s(phase)=NaN;T.medianIndividualPeakRotation_0_to_0p5s(phase)=NaN;T.medianTurnAngleDeg(phase)=NaN;
 end
end
writetable(T,fullfile(out,'selection_summary.csv'));
% Include the screenshot-window sensitivity when choosing shared page limits.
O=load(fullfile(root,'20261020_11_2','presentation_control_vs_heat_20261102','original_criteria_turns','comparison_data.mat'));
otherMu=mu;otherSem=sem;
for phase=1:2
 label="control";if phase==2,label="heat";end
 q=O.selected(O.selected.context==label,:);take=q.absoluteTurnDeg>45;v=O.X{phase}(:,take,:);
 otherMu(:,:,phase)=squeeze(mean(v,2));otherSem(:,:,phase)=squeeze(std(v,0,2))/sqrt(size(v,2));
end
% Shared paper axes include both definitions and all mean/SEM values.
lims=cell(1,3);steps=[.5 2 50];
for row=1:3
 v=[reshape(mu(:,row,:),[],1);reshape(mu(:,row,:)-sem(:,row,:),[],1);reshape(mu(:,row,:)+sem(:,row,:),[],1)];v=[v;reshape(otherMu(:,row,:)-otherSem(:,row,:),[],1);reshape(otherMu(:,row,:)+otherSem(:,row,:),[],1)];v=v(isfinite(v));
 if isempty(v),lims{row}=[0 1];else
  lims{row}=[floor(min([0;v])/steps(row))*steps(row),ceil(max(v)/steps(row))*steps(row)];
  if diff(lims{row})<=0,lims{row}=[0 1];end
 end
end
colors={[.86 .06 .06],[.08 .08 .08],[.05 .2 .8]};
labels={'\DeltaF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};
fig=figure('Visible','off','Color','white','Units','inches','Position',[1 1 7.2 8.2]);
set(fig,'PaperUnits','inches','PaperPosition',[0 0 7.2 8.2],'PaperSize',[7.2 8.2],'PaperPositionMode','manual');
positions=[.14 .70 .35 .20;.60 .70 .35 .20;.14 .43 .35 .20;.60 .43 .35 .20;.14 .16 .35 .20;.60 .16 .35 .20];
for row=1:3
 for phase=1:2
  n=(row-1)*2+phase;ax=axes(fig,'Units','normalized','Position',positions(n,:), ...
   'FontName','Arial','FontSize',10,'Box','off','TickDir','out','LineWidth',.6, ...
   'XLim',[-2 3],'YLim',lims{row},'XTick',-2:1:3);hold(ax,'on');
  if counts(phase)>1
   m=mu(:,row,phase);se=sem(:,row,phase);
   fill(ax,[rel;flipud(rel)],[m-se;flipud(m+se)],colors{row},'FaceAlpha',.18,'EdgeColor','none','Tag','event_sem');
  end
  if counts(phase)>0
   plot(ax,rel,mu(:,row,phase),'Color',colors{row},'LineWidth',1.5,'Tag','event_mean');
  else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
  xline(ax,0,':','Color',[.5 .5 .5],'LineWidth',.7,'HandleVisibility','off');
  if phase==1,ylabel(ax,labels{row},'FontSize',11);end
  if row==1,title(ax,sprintf('%s (n = %d)',phaseNames(phase),counts(phase)),'FontSize',12,'FontWeight','bold');end
  if row==3,xlabel(ax,'Time from slowdown (s)','FontSize',10);else,ax.XTickLabel=[];end
 end
 annotation(fig,'textbox',[.02 positions((row-1)*2+1,2)+.19 .05 .035], ...
  'String',char('A'+row-1),'FontName','Arial','FontSize',13,'FontWeight','bold','EdgeColor','none','Margin',0);
end
annotation(fig,'textbox',[.14 .955 .82 .025],'String','VES041 | 2026-10-20, fly 1 | turns >45 deg', ...
 'FontName','Arial','FontSize',12,'FontWeight','bold','EdgeColor','none','Margin',0);
annotation(fig,'textbox',[.14 .02 .82 .055],'String', ...
 {'Mean +/- SEM across events; one fly.','Turn: accumulated absolute yaw in 0-0.5 s after slowdown.'}, ...
 'FontName','Arial','FontSize',8,'EdgeColor','none','Margin',0);
stem=fullfile(out,'20261020_11_2_turns_gt45_paper');
print(fig,[stem '.pdf'],'-dpdf','-painters');print(fig,[stem '.svg'],'-dsvg','-painters');print(fig,[stem '.png'],'-dpng','-r600');
set(fig,'Visible','on');savefig(fig,[stem '.fig']);close(fig);
% Separate matched-axis rotation comparison, with no selection changes.
fig=figure('Visible','off','Color','white','Position',[50 50 1200 520]);tl=tiledlayout(fig,1,2,'TileSpacing','compact');
for phase=1:2
 ax=nexttile(tl);hold(ax,'on');plot(ax,rel,oldMu(:,3,phase),'Color',[.55 .55 .55],'LineWidth',1.5,'DisplayName',sprintf('Original >=15 deg (n=%d)',countsOld(phase)));
 plot(ax,rel,mu(:,3,phase),'Color',colors{3},'LineWidth',2,'DisplayName',sprintf('Strict >45 deg (n=%d)',counts(phase)));
 xline(ax,0,':','Color',[.5 .5 .5]);set(ax,'FontName','Arial','FontSize',13,'Box','off','XLim',[-2 3]);
 title(ax,phaseNames(phase));xlabel(ax,'Time from strongest deceleration (s)');ylabel(ax,'Mean rotational speed (deg/s)');legend(ax,'Location','best','Box','off');
end
axs=findall(fig,'Type','axes');values=[oldMu(:,3,:);mu(:,3,:)];ymax=ceil(max(values,[],'all')/50)*50;set(axs,'YLim',[0 ymax]);
print(fig,fullfile(out,'rotation_threshold_comparison.svg'),'-dsvg','-painters');exportgraphics(fig,fullfile(out,'rotation_threshold_comparison.png'),'Resolution',180);
set(fig,'Visible','on');savefig(fig,fullfile(out,'rotation_threshold_comparison.fig'));close(fig);
methods=struct('recording','20261020_11_2','source',source,'strictThresholdDeg',45, ...
 'comparison','Original linked PDF >=15 deg vs new strict >45 deg, all other rules held fixed', ...
 'turnMeasurement','Accumulated absolute 0.25-second median-smoothed yaw, exact integration over 0:+0.5 s after strongest deceleration', ...
 'slowdownCandidates','Exact saved original group2; acceleration <=-2 mm/s^2, speed before >=1 mm/s, drop >=1 mm/s, >=2 s accepted spacing and full-valid -2:+3 s window', ...
 'split','Original measured heat command at event sample >0.05 V; transitions elsewhere in window retained', ...
 'alignment','Original strongest forward deceleration', ...
 'display','Original saved dF/F and physical kinematics, no extra smoothing or normalization; rotational speed is absolute before averaging', ...
 'uncertainty','Sample SD/sqrt(N events); single-fly descriptive SEM', ...
 'validation','Original >=15 group reconstructed exactly and reproduces linked PDF counts 20/19; candidate selection is unchanged');
save(fullfile(out,'analysis_gt45.mat'),'idx','angle','selected','heat','countsOld','counts','rel','X','oldX','mu','sem','oldMu','E','T','lims','methods','-v7.3');
fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);
append_gt45_screenshot_window_paper(out,lims);
disp(T);fprintf('Saved paper figure to %s\n',out);
end
