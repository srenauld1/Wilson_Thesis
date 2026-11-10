function plot_other_heat_flies_individually(out)
% Exact saved per-event arrays; no re-detection or extra behavior smoothing.
if nargin<1,out=fullfile(pwd,'outputs','ves041_other_heat_flies_individual_20261105');end
if ~isfolder(out),mkdir(out);end
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier';
source=fullfile(root,'pooled_four_original_pizza_rules_20261102','baseline_subtracted_dff','excluding_april27_fly3','crunch_data.mat');
S=load(source);ids=["20260427-2_ves041_fly2_heat","20261019_4_1"];
C=struct([]);allBounds=cell(1,3);rel=S.rel(:);
for k=1:2
 ix=find(string(S.ids)==ids(k));assert(isscalar(ix));
 c=struct('id',ids(k),'fly',string(S.flies(ix)),'X',{S.X(ix,:)},'counts',S.counts(ix,:),'mu',[],'sem',[]);
 for phase=1:2
  % Match the pooled screenshot's event baseline subtraction plus one.
  c.X{phase}(:,:,1)=c.X{phase}(:,:,1)+1;
  c.mu(:,:,phase)=squeeze(mean(c.X{phase},2));
  c.sem(:,:,phase)=squeeze(std(c.X{phase},0,2))/sqrt(size(c.X{phase},2));
  assert(size(c.X{phase},2)==c.counts(phase));
  for row=1:3,allBounds{row}=[allBounds{row};c.mu(:,row,phase)-c.sem(:,row,phase);c.mu(:,row,phase)+c.sem(:,row,phase)];end %#ok<AGROW>
 end
 if k==1,C=c;else,C(k)=c;end
end
assert(isequal(C(1).counts,[6 3])&&isequal(C(2).counts,[17 7]));
lims=cell(1,3);
for row=1:3
 v=allBounds{row};step=[.5 2 50];step=step(row);
 lims{row}=[floor(min([v;0])/step)*step ceil(max(v)/step)*step];
 if diff(lims{row})==0,lims{row}=[0 1];end
end
colors={[.90 .04 .03],[.1 .1 .1],[.05 .20 .8]};labels={'dF/F (baseline set to 1)','Forward velocity (mm/s)','Rotational speed (deg/s)'};
for k=1:2
 c=C(k);fig=figure('Visible','off','Color','w','Position',[50 50 1200 1150]);
 tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');
 for row=1:3
  for phase=1:2
   ax=nexttile(tl);hold(ax,'on');mu=c.mu(:,row,phase);se=c.sem(:,row,phase);
   fill(ax,[rel;flipud(rel)],[mu-se;flipud(mu+se)],colors{row},'FaceAlpha',.18,'EdgeColor','none','Tag','event_sem');
   plot(ax,rel,mu,'Color',colors{row},'LineWidth',2.4,'Tag','event_mean');
   xline(ax,0,':','Color',[.5 .5 .5],'LineWidth',1);
   set(ax,'XLim',[-2 2],'YLim',lims{row},'FontName','Arial','FontSize',17,'Box','off','TickDir','out');
   if phase==1,ylabel(ax,labels{row},'FontSize',19);end
   if row==3,xlabel(ax,'Time from strongest deceleration (s)','FontSize',16);else,ax.XTickLabel=[];end
   if row==1
    names={'Heat off','Heat on'};title(ax,sprintf('%s | n=%d events',names{phase},c.counts(phase)),'FontSize',20);
   end
  end
 end
 title(tl,strrep(char(c.id),'_',' '),'FontName','Arial','FontSize',23,'Interpreter','none');
 subtitle(tl,{'Saved slowdown + reorientation events; mean +/- event SEM', ...
  'Heat off is the source category, not a strict spontaneous-event classification'},'FontSize',13);
 dest=fullfile(out,char(c.id+"_heat_off_on_individual"));
 print(fig,[dest '.svg'],'-dsvg','-painters');exportgraphics(fig,[dest '.png'],'Resolution',240);
 set(fig,'Visible','on');savefig(fig,[dest '.fig']);close(fig);
end
% Audit the precise single-fly arrays used in the screenshot's left pair.
exampleSource=fullfile(root,'20261020_11_2','presentation_control_vs_heat_20261102','original_criteria_turns','comparison_data.mat');
O=load(exampleSource);assert(isequal(O.counts,[15 18]));
R=cell(0,9);core=rel>=-.5&rel<=.5;diagTime=(-1:.05:1)';aligned=cell(1,2);
for phase=1:2
 rotation=O.X{phase}(:,:,3);[individualPeaks,ii]=max(rotation(core,:),[],1);
 coreTime=rel(core);peakTimes=coreTime(ii);meanTrace=mean(rotation,2);meanCore=meanTrace(core);
 aligned{phase}=nan(numel(diagTime),size(rotation,2));
 for j=1:size(rotation,2),aligned{phase}(:,j)=interp1(rel,rotation(:,j),diagTime+peakTimes(j),'linear');end
 R(end+1,:)={phase,size(rotation,2),median(individualPeaks),mean(individualPeaks),max(meanCore),min(peakTimes),max(peakTimes),max(meanTrace),max(mean(aligned{phase},2))}; %#ok<AGROW>
end
rotationAudit=cell2table(R,'VariableNames',{'condition','nEvents','medianIndividualPeakWithinHalfSecond','meanIndividualPeakWithinHalfSecond','maxMeanWithinHalfSecond','earliestPeakTime','latestPeakTime','maxMeanFullWindow','maxMeanYawPeakAligned'});
writetable(rotationAudit,fullfile(out,'example_rotation_averaging_audit.csv'));
fig=figure('Visible','off','Color','w','Position',[50 50 1250 650]);tl=tiledlayout(fig,1,2,'TileSpacing','compact');
for phase=1:2
 ax=nexttile(tl);hold(ax,'on');rot=O.X{phase}(:,:,3);take=rel>=-1&rel<=1;
 plot(ax,rel(take),rot(take,:),'Color',[.78 .82 .93],'LineWidth',.6);
 plot(ax,rel(take),mean(rot(take,:),2),'Color',[.05 .2 .8],'LineWidth',2.5,'DisplayName','Aligned to deceleration');
 plot(ax,diagTime,mean(aligned{phase},2),'Color',[.7 .2 .6],'LineWidth',2.5,'DisplayName','Same events aligned to yaw peak');
 xline(ax,0,':','Color',[.5 .5 .5]);set(ax,'FontName','Arial','FontSize',16,'Box','off','XLim',[-1 1]);
 names={'Placed control barriers','Heat'};title(ax,sprintf('%s | n=%d',names{phase},O.counts(phase)),'FontSize',20);
 xlabel(ax,'Time from alignment point (s)');ylabel(ax,'Rotational speed (deg/s)');
end
sgtitle(tl,'Example fly: event peaks are larger than the deceleration-aligned mean','FontSize',21);
subtitle(tl,'Thin: individual events | Blue: original alignment | Purple: diagnostic yaw-peak alignment','FontSize',13);
print(fig,fullfile(out,'example_rotation_alignment_diagnostic.svg'),'-dsvg','-painters');
exportgraphics(fig,fullfile(out,'example_rotation_alignment_diagnostic.png'),'Resolution',200);
set(fig,'Visible','on');savefig(fig,fullfile(out,'example_rotation_alignment_diagnostic.fig'));close(fig);
counts=table(string({C.id})',vertcat(C.counts),'VariableNames',{'recording','heatOffOnCounts'});
writetable(counts,fullfile(out,'event_counts.csv'));
methods=struct('source',source,'exampleSource',exampleSource,'split','Exact original heat-state at event, not a new spontaneous classification', ...
 'displayDff','Existing per-event -2:-1 s baseline subtraction plus 1, matching the screenshot pooled panel', ...
 'kinematics','Existing physical forward and absolute yaw arrays; absolute value before event averaging; no new smoothing', ...
 'selection','Exact original saved slowdown+rotation group4; no re-detection or event filtering', ...
 'sem','Sample SD across events divided by sqrt(N events), separately for each fly and condition', ...
 'alignment','Strongest forward deceleration; yaw-peak re-alignment is a separate diagnostic only', ...
 'originalExampleSelection','The single-fly screenshot uses original-criteria >=30 deg in +/-0.5 s control/heat contexts; pooled other-fly source uses original pizza >=15 deg in 0:+0.5 s and heat at event time. These sources are not identical selections.');
save(fullfile(out,'individual_fly_data.mat'),'C','rel','lims','methods','O','rotationAudit','diagTime','aligned','-v7.3');
fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);
disp(rotationAudit);fprintf('Saved individual fly plots to %s\n',out);
end
