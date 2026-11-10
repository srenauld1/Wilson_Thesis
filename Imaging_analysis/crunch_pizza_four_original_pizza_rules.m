function result=crunch_pizza_four_original_pizza_rules()
% Exact original pizza-talk event groups and original heat-at-event split.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier';
ids=["20260427-2_ves041_fly2_heat","20260427-7_ves041_fly3_heat_20mm","20261019_4_1","20261020_11_2"];
flies=["20260427_fly2","20260427_fly3","20261019_fly1","20261020_fly1"];
out=fullfile(root,'pooled_four_original_pizza_rules_20261102');if ~isfolder(out),mkdir(out);end
individual=fullfile(out,'individual_trials');if ~isfolder(individual),mkdir(individual);end
X=cell(4,2);counts=zeros(4,2);rows=cell(0,7);allE=table;flyMeans=[];
for k=1:4
 file=fullfile(root,ids(k),'analysis.mat');A=load(file);E=readtable(fullfile(root,ids(k),'events.csv'),'Delimiter',',','TextType','string');E=E(E.group==4,:);
 idx=A.groups{4};[found,order]=ismember(idx,E.sampleIndex);assert(all(found)&&height(E)==numel(idx));E=E(order,:);assert(isequal(E.sampleIndex,idx));on=A.d.heat(idx)>A.meta.threshold;assert(isequal(logical(E.heatOn),on));
 if k==1,rel=A.rel(:);flyMeans=nan(numel(rel),4,3,2);else,assert(isequal(rel,A.rel(:)));end
 for s=1:2
  take=on==(s==2);X{k,s}=A.X{4}(:,take,:);counts(k,s)=nnz(take);assert(size(X{k,s},2)==counts(k,s)&&all(isfinite(X{k,s}(:))));
  flyMeans(:,k,:,s)=mean(X{k,s},2);
 end
 E.recording=repmat(ids(k),height(E),1);E.flyId=repmat(flies(k),height(E),1);allE=[allE;E];
 rows(end+1,:)={ids(k),flies(k),A.meta.source,A.meta.heatSource,A.meta.threshold,counts(k,1),counts(k,2)};
 fprintf('%s: %d heat-off / %d heat-on turns\n',ids(k),counts(k,:));
end
assert(isequal(counts,[6 3;43 3;17 7;20 19]));assert(numel(unique(flies))==4);
meanTrace=nan(numel(rel),3,2);sem=meanTrace;
for s=1:2
 meanTrace(:,:,s)=squeeze(mean(flyMeans(:,:,:,s),2));sem(:,:,s)=squeeze(std(flyMeans(:,:,:,s),0,2))/sqrt(4);
 expected=zeros(numel(rel),3);for k=1:4,expected=expected+squeeze(mean(X{k,s},2))/4;end
 assert(max(abs(expected-meanTrace(:,:,s)),[],'all')<1e-10);
end
stem=fullfile(out,'ves041_four_trials_original_pizza_rules_heat_off_on_turns');draw(flyMeans,meanTrace,sem,rel,sum(counts,1),stem,'4 requested recordings | Original pizza-talk rules',true);
for k=1:4
 values=cell(1,2);mu=nan(numel(rel),3,2);se=mu;
 for s=1:2,values{s}=X{k,s};mu(:,:,s)=squeeze(mean(values{s},2));se(:,:,s)=squeeze(std(values{s},0,2))/sqrt(counts(k,s));end
 draw(values,mu,se,rel,counts(k,:),fullfile(individual,char(ids(k)+"_original_pizza_rules_heat_off_on_turns")),strrep(char(ids(k)),'_',' '),false);
end
F=cell2table(rows,'VariableNames',{'recording','flyId','traceSource','heatSource','heatThreshold','heatOffTurns','heatOnTurns'});writetable(F,fullfile(out,'cohort_and_counts.csv'));writetable(allE,fullfile(out,'selected_events.csv'));
methods=struct('selection','Exact saved original analysis.mat group4 (slowdowns with rotation); no re-detection, exclusions or amplitude filtering', ...
 'turn','Accumulated absolute .25-second median-smoothed yaw >=15 degrees in 0:+0.5 s after strongest deceleration', ...
 'slowdown','Acceleration<=-2 mm/s^2 after .25-second median smoothing, speed-before>=1 mm/s, forward drop>=1 mm/s, >=2 s accepted-event spacing', ...
 'speedWindows','Median forward -0.75:-0.25 versus +0.25:+0.75 s','eventWindow','Original full-valid -2:+3 s windows', ...
 'heat','Original source channel and threshold at event sample only: legacy Virmen command>0.0001; October DAQ command>0.05 V', ...
 'heatOffMeaning','Heat off is not synonymous with control barrier. Accommodation, lead-in, and protocol heat-trial zero-command periods remain included when present in the original selection. Heat transitions elsewhere in the event window are retained.', ...
 'pooledMean','Equal weight per fly, after averaging that fly events within each condition; all four flies contribute to both conditions', ...
 'pooledSEM','Pointwise sample SD of four fly means divided by sqrt(4)', ...
 'individualSEM','Pointwise sample SD of events divided by sqrt(N events); descriptive single-recording uncertainty', ...
 'normalization','Original saved dF/F and physical kinematics; no added baseline subtraction, z-score or display smoothing', ...
 'validation','Every selected sample index and heat-state label matches original events.csv; saved arrays reused without alteration');
save(fullfile(out,'crunch_data.mat'),'X','flyMeans','meanTrace','sem','counts','ids','flies','rel','F','allE','methods','-v7.3');fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('output',out,'heatOffTurns',sum(counts(:,1)),'heatOnTurns',sum(counts(:,2)),'fliesPerCondition',[4 4]);disp(result);
end
function draw(values,mu,se,rel,counts,stem,caption,pooled)
fig=figure('Visible','off','Color','w','Position',[20 20 1120 1050]);tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};light={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};names={'Heat-off turns','Heat-on turns'};
for row=1:3
 z=[];for s=1:2,if pooled,q=values(:,:,row,s);else,q=values{s}(:,:,row);end;z=[z;q(:);mu(:,row,s)-se(:,row,s);mu(:,row,s)+se(:,row,s)];end;z=z(isfinite(z));pad=max(.01,.06*range(z));lim=[min(z)-pad max(z)+pad];
 for s=1:2
  ax=nexttile(tl);hold(ax,'on');if pooled,q=values(:,:,row,s);else,q=values{s}(:,:,row);end
  plot(ax,rel,q,'Color',light{row},'LineWidth',.9,'Tag','thinTrace');m=mu(:,row,s);e=se(:,row,s);fill(ax,[rel;flipud(rel)],[m-e;flipud(m+e)],colors{row},'FaceAlpha',.2,'EdgeColor','none');plot(ax,rel,m,'Color',colors{row},'LineWidth',3);
  xline(ax,0,':','Color',[.45 .45 .45]);xlim(ax,[-2 3]);ylim(ax,lim);ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=-2:1:3;
  if s==1,ylabel(ax,labels{row},'FontSize',21);end
  if row==3,xlabel(ax,'Time from strongest deceleration (s)','FontSize',17);else,ax.XTickLabel=[];end
  if row==1,title(ax,names{s},'FontSize',23,'FontWeight','bold');label=sprintf('n=%d events',counts(s));if pooled,label=sprintf('4 flies | %d events',counts(s));end;subtitle(ax,label,'FontSize',15);end
 end
end
text='Thin: events | Thick: mean | Shading: SEM across events';if pooled,text='Thin: fly means | Thick: equal-fly mean | Shading: SEM across flies';end
subtitle(tl,{caption,text},'FontSize',15);assert(numel(findall(fig,'Type','axes'))==6);expected=3*sum(counts);if pooled,expected=24;end;assert(numel(findall(fig,'Tag','thinTrace'))==expected);
exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
