function result=plot_pizza_1020_11_control_heat_presentation()
% Match Imaging Analysis four-selection presentation for one heat recording.
% Preserve saved event selection; controls are placed protocol control trials.
% Heat pools prior measured-voltage bins 2:4, including heat-trial zero voltage.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
folder='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/20261020_11_2';
out=fullfile(folder,'presentation_control_vs_heat_20261102');if ~isfolder(out),mkdir(out);end
source=fullfile(folder,'control_and_heat_bins','binned_analysis.mat');A=load(source);B=load(fullfile(folder,'analysis.mat'));
rel=A.rel(:);columns=[1 4 5 6];D=cell(2,4);counts=zeros(2,4);barrierCounts=counts;condition=["Control","Heat"];
E=A.E(ismember(A.E.eventGroup,columns),:);E.condition=repmat("Excluded",height(E),1);E.condition(E.heatBin==1)="Control";E.condition(E.heatBin>1)="Heat";
for c=1:4
 D{1,c}=A.C{1,columns(c)};D{2,c}=cat(2,A.C{2:4,columns(c)});
 % Rebuild every array from saved original groups in chronological order.
 for s=1:2
  q=E.eventGroup==columns(c)&E.condition==condition(s);idx=E.sampleIndex(q);[tf,ix]=ismember(idx,B.groups{columns(c)});assert(all(tf));
  rebuilt=B.X{columns(c)}(:,ix,:);original=D{s,c};
  assert(size(original,2)==numel(idx));
  % Pooling voltage bins reorders heat events; use the original chronology.
  if s==1,sourceIndices=E.sampleIndex(q);else,sourceIndices=[];for b=2:4,sourceIndices=[sourceIndices;E.sampleIndex(q&E.heatBin==b)];end,end
  [found,permutation]=ismember(idx,sourceIndices);assert(all(found)&&isequaln(original(:,permutation,:),rebuilt));
  D{s,c}=rebuilt;counts(s,c)=size(rebuilt,2);barrierCounts(s,c)=numel(unique(E.barrierTrial(q)));
  assert(all(isfinite(rebuilt(:))));
  for j=idx'
   t=B.d.t;lo=find(t<=t(j)-2,1,'last');hi=find(t>=t(j)+3,1,'first');win=lo:hi;
   assert(all(B.d.valid(win))&&all(A.P.active(win))&&all(A.P.control(win)==(s==1)));
   assert(numel(unique(A.P.trial(win)))==1&&numel(unique(A.P.geometry(win)))==1);
  end
 end
end
mu=nan(numel(rel),3,2,4);sem=mu;
for s=1:2,for c=1:4
 X=D{s,c};if isempty(X),continue,end
 mu(:,:,s,c)=squeeze(mean(X,2));if size(X,2)>1,sem(:,:,s,c)=squeeze(std(X,0,2))/sqrt(size(X,2));end
 % Independent sum-of-squared-deviation SEM check.
 if size(X,2)>1
  centered=X-mean(X,2);expected=squeeze(sqrt(sum(centered.^2,2)/(size(X,2)-1))/sqrt(size(X,2)));
  delta=abs(expected-sem(:,:,s,c));assert(all(delta(:)<1e-10));
 end
end,end
rowLimits=zeros(3,2);
for v=1:3
 z=[];for s=1:2,for c=1:4,q=D{s,c}(:,:,v);z=[z;q(:)];end,end
 band=[reshape(mu(:,v,:,:)-sem(:,v,:,:),[],1);reshape(mu(:,v,:,:)+sem(:,v,:,:),[],1)];z=[z;band];z=z(isfinite(z));pad=.06*(max(z)-min(z));rowLimits(v,:)=[min(z)-pad max(z)+pad];
end
titles={'VES041 peaks','Slowdowns + rotation',{'Rotation peaks','Speed drop <1 mm/s'},{'Rotation peaks','Speed drop \geq1 mm/s'}};
labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};thinColors={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};
files=strings(3,1);
for mode=1:3
 fig=figure('Visible','off','Color','w','Position',[20 20 1800 1050]);tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','loose');
 for row=1:3,for c=1:4
  ax=nexttile(tl);hold(ax,'on');states=1:2;if mode<3,states=mode;end
  for s=states
   style='-';if mode==3&&s==1,style='--';end
   X=D{s,c}(:,:,row);light=thinColors{row};if mode==3&&s==1,light=1-.6*(1-light);end
   if ~isempty(X),plot(ax,rel,X,'Color',light,'LineWidth',.5,'LineStyle',style,'Tag','eventTrace');end
   m=mu(:,row,s,c);se=sem(:,row,s,c);fill(ax,[rel;flipud(rel)],[m-se;flipud(m+se)],colors{row},'FaceAlpha',.17,'EdgeColor','none','Tag','eventSem');
  end
  for s=states
   style='-';if mode==3&&s==1,style='--';end
   plot(ax,rel,mu(:,row,s,c),'Color',colors{row},'LineWidth',3,'LineStyle',style,'Tag','conditionMean','UserData',condition(s));
  end
  xline(ax,0,':','Color',[.45 .45 .45],'LineWidth',1);xlim(ax,[-2 3]);ylim(ax,rowLimits(row,:));
  ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=-2:1:3;
  if c==1,ylabel(ax,labels{row},'FontSize',21);end
  if row==3,xlabel(ax,'Time (s)','FontSize',21);else,ax.XTickLabel=[];end
  if row==1
   title(ax,titles{c},'FontSize',23,'FontWeight','bold');
   if mode<3,subtitle(ax,sprintf('n=%d events | %d barrier trials',counts(mode,c),barrierCounts(mode,c)),'FontSize',12);
   else,subtitle(ax,sprintf('Control n=%d | Heat n=%d',counts(1,c),counts(2,c)),'FontSize',12);end
  end
 end,end
 if mode<3
  caption={sprintf('2026-10-20 trial 11-2 | %s barriers',condition(mode)),'Thin: individual events | Thick: mean | Shading: SEM across events | One recording, one fly'};
  stem=fullfile(out,['20261020_11_2_' lower(char(condition(mode))) '_four_selections_events_sem']);
 else
  caption={'2026-10-20 trial 11-2 | Solid: heat barriers | Dashed: control barriers','Thin: individual events | Thick: condition means | Shading: SEM across events | One recording, one fly'};
  stem=fullfile(out,'20261020_11_2_control_vs_heat_four_selections_events_sem');
 end
 subtitle(tl,caption,'FontSize',16);assert(numel(findall(fig,'Type','axes'))==12);
 expected=3*sum(counts(states,:),'all');assert(numel(findall(fig,'Tag','eventTrace'))==expected);
 exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);files(mode)=stem;
end
rows=cell(0,5);for s=1:2,for c=1:4,rows(end+1,:)={condition(s),c,columns(c),counts(s,c),barrierCounts(s,c)};end,end
T=cell2table(rows,'VariableNames',{'condition','displayColumn','sourceColumn','events','barrierTrials'});writetable(T,fullfile(out,'selection_counts.csv'));writetable(E,fullfile(out,'selected_and_excluded_events.csv'));
methods=struct('recording','20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo','source',source, ...
 'referenceChat','Imaging Analysis','referenceChatId','019f6872-db96-7780-ac0f-9e899008b83e', ...
 'referenceFormat','plot_pizza_bar_presentation_fly_sem.m original-units mode', ...
 'sourceColumns',columns,'lowForwardSpeedLt2Filter',false,'conditionDefinition','Placed control_trial=1 versus placed control_trial=0; heat pools prior voltage bins2:4, including zero-voltage heat-trial periods', ...
 'fullWindow','-2:+3 seconds valid, same placed barrier trial, geometry and control flag; original bin exclusions retained', ...
 'normalization','Original saved dF/F and physical kinematics; no new baseline subtraction, z-score, smoothing, or detector', ...
 'sem','Pointwise sample SD/sqrt(number of events) within condition; descriptive single-fly SEM, not between-fly SEM', ...
 'turnDefinition','Accumulated absolute 0.25-second median-smoothed yaw >=15 deg within 0:+0.5 s after strongest deceleration', ...
 'rotationDropDefinition','0.25-second median-smoothed forward; median(-0.75:-0.25 s) minus median(+0.25:+0.75 s); split <1 and >=1 mm/s');
methods.referenceUserTurns={'01a19d16-7220-7f72-a298-bbf271f1366c','01a19d0c-47cc-7f42-9bdd-012566936a38','01a19d10-19bd-7291-a100-aa5d852441ed'};
methods.peakDetector=B.peakDetector;methods.rotationDetector=B.rotationDetector;
methods.peakRule='Original dF/F findpeaks: prominence>=3*1.4826*MAD(dF/F minus 30-second moving median), spacing>=2 s, width>=0.25 s';
methods.slowdownDetector=struct('smoothSec',.25,'maximumAccelerationMmS2',-2,'minimumForwardBeforeMmS',1,'minimumForwardDropMmS',1,'minimumSpacingSec',2,'beforeWindowSec',[-.75 -.25],'afterWindowSec',[.25 .75]);
save(fullfile(out,'comparison_data.mat'),'D','mu','sem','counts','barrierCounts','rel','columns','rowLimits','methods','T','-v7.3');
fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('folder',out,'controlCounts',counts(1,:),'heatCounts',counts(2,:),'figures',files);disp(result);
end
