function result=exclude_april27_fly3_pizza_baseline()
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/pooled_four_original_pizza_rules_20261102/baseline_subtracted_dff';
S=load(fullfile(root,'crunch_data.mat'));excluded="20260427-7_ves041_fly3_heat_20mm";keep=S.ids~=excluded;assert(nnz(keep)==3);
out=fullfile(root,'excluding_april27_fly3');if ~isfolder(out),mkdir(out);end
S.X=S.X(keep,:);S.counts=S.counts(keep,:);S.flyMeans=S.flyMeans(:,keep,:,:);S.ids=S.ids(keep);S.flies=S.flies(keep);S.F=S.F(keep,:);S.baselines=S.baselines(keep,:);S.allE=S.allE(S.allE.recording~=excluded,:);
for s=1:2
 S.meanTrace(:,:,s)=squeeze(mean(S.flyMeans(:,:,:,s),2));S.sem(:,:,s)=squeeze(std(S.flyMeans(:,:,:,s),0,2))/sqrt(3);
 for row=1:3
  Q=squeeze(S.flyMeans(:,:,row,s));mu=sum(Q,2)/3;se=sqrt(sum((Q-mu).^2,2)/2)/sqrt(3);
  assert(max(abs(mu-S.meanTrace(:,row,s)))<1e-10&&max(abs(se-S.sem(:,row,s)))<1e-10);
 end
end
assert(isequal(sum(S.counts,1),[43 29]));assert(max(abs(mean(S.meanTrace(S.baselineMask,1,:),1)),[],'all')<1e-10);
stem=fullfile(out,'ves041_three_flies_heat_off_on_baseline_subtracted_dff');
draw(S.flyMeans,S.meanTrace,S.sem,S.rel,sum(S.counts,1),stem,'3 recordings | April 27 fly 3 excluded | dF/F baseline: -2 to -1 s',true);
S.methods.excludedRecording=excluded;S.methods.pooledMean='Equal mean of remaining three fly means';S.methods.pooledSEM='Sample SD of remaining three fly means divided by sqrt(3)';
S.methods.validation='Retained event traces and baselines unchanged; three-fly means and SEM independently verified';
save(fullfile(out,'crunch_data.mat'),'-struct','S','-v7.3');writetable(S.F,fullfile(out,'cohort_and_counts.csv'));fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(S.methods,PrettyPrint=true));fclose(fid);
result=struct('output',out,'remainingFlies',S.flies,'heatOffEvents',43,'heatOnEvents',29);disp(result);
end

function draw(values,mu,se,rel,counts,stem,caption,pooled)
fig=figure('Visible','off','Color','w','Position',[20 20 1120 1050]);tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};light={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};labels={'Baseline-subtracted dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};names={'Heat-off turns','Heat-on turns'};
for row=1:3
 z=[];for s=1:2,if pooled,q=values(:,:,row,s);else,q=values{s}(:,:,row);end;z=[z;q(:);mu(:,row,s)-se(:,row,s);mu(:,row,s)+se(:,row,s)];end;z=z(isfinite(z));pad=max(.01,.06*range(z));lim=[min(z)-pad max(z)+pad];
 for s=1:2
  ax=nexttile(tl);hold(ax,'on');if pooled,q=values(:,:,row,s);else,q=values{s}(:,:,row);end
  plot(ax,rel,q,'Color',light{row},'LineWidth',.9,'Tag','thinTrace');m=mu(:,row,s);e=se(:,row,s);fill(ax,[rel;flipud(rel)],[m-e;flipud(m+e)],colors{row},'FaceAlpha',.2,'EdgeColor','none');plot(ax,rel,m,'Color',colors{row},'LineWidth',3);
  xline(ax,0,':','Color',[.45 .45 .45]);if row==1,yline(ax,0,':','Color',[.6 .6 .6]);end;xlim(ax,[-2 3]);ylim(ax,lim);ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=-2:1:3;
  if s==1,ylabel(ax,labels{row},'FontSize',21);end
  if row==3,xlabel(ax,'Time from strongest deceleration (s)','FontSize',17);else,ax.XTickLabel=[];end
  if row==1,title(ax,names{s},'FontSize',23,'FontWeight','bold');label=sprintf('n=%d events',counts(s));if pooled,label=sprintf('3 flies | %d events',counts(s));end;subtitle(ax,label,'FontSize',15);end
 end
end
text='Thin: events | Thick: mean | Shading: SEM across events';if pooled,text='Thin: fly means | Thick: equal-fly mean | Shading: SEM across flies';end
subtitle(tl,{caption,text},'FontSize',15);assert(numel(findall(fig,'Type','axes'))==6);expected=3*sum(counts);if pooled,expected=18;end;assert(numel(findall(fig,'Tag','thinTrace'))==expected);
exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
