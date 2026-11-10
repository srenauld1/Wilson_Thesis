function result=baseline_pizza_four_original_rules_dff()
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/pooled_four_original_pizza_rules_20261102';
S=load(fullfile(root,'crunch_data.mat'));original=S;out=fullfile(root,'baseline_subtracted_dff');if ~isfolder(out),mkdir(out);end
individual=fullfile(out,'individual_trials');if ~isfolder(individual),mkdir(individual);end
baselineMask=S.rel>=-2&S.rel< -1;assert(nnz(baselineMask)==20);baselines=cell(4,2);rows=cell(0,6);
for k=1:4,for s=1:2
 X=S.X{k,s};base=mean(X(baselineMask,:,1),1);baselines{k,s}=base;X(:,:,1)=X(:,:,1)-base;
 assert(max(abs(mean(X(baselineMask,:,1),1)))<1e-10);assert(isequaln(X(:,:,2:3),original.X{k,s}(:,:,2:3)));
 reconstructed=X(:,:,1)+base;assert(max(abs(reconstructed-original.X{k,s}(:,:,1)),[],'all')<1e-10);
 S.X{k,s}=X;S.flyMeans(:,k,:,s)=mean(X,2);
 events=S.allE(S.allE.recording==S.ids(k)&S.allE.heatOn==(s==2),:);assert(height(events)==numel(base));
 for j=1:numel(base),rows(end+1,:)={S.ids(k),S.flies(k),s==2,events.sampleIndex(j),events.eventTime(j),base(j)};end
end,end
for s=1:2
 S.meanTrace(:,:,s)=squeeze(mean(S.flyMeans(:,:,:,s),2));S.sem(:,:,s)=squeeze(std(S.flyMeans(:,:,:,s),0,2))/sqrt(4);
end
assert(isequaln(S.meanTrace(:,2:3,:),original.meanTrace(:,2:3,:))&&isequaln(S.sem(:,2:3,:),original.sem(:,2:3,:)));
assert(max(abs(mean(S.meanTrace(baselineMask,1,:),1)),[],'all')<1e-10);
% Independently verify the mean and SEM after event-wise subtraction.
for s=1:2,for v=1:3
 Q=squeeze(S.flyMeans(:,:,v,s));mu=sum(Q,2)/4;se=sqrt(sum((Q-mu).^2,2)/3)/2;
 assert(max(abs(mu-S.meanTrace(:,v,s)))<1e-10&&max(abs(se-S.sem(:,v,s)))<1e-10);
end,end
stem=fullfile(out,'ves041_four_trials_original_pizza_rules_heat_off_on_turns_baseline_subtracted_dff');
draw(S.flyMeans,S.meanTrace,S.sem,S.rel,sum(S.counts,1),stem,'4 requested recordings | dF/F baseline: -2 to -1 s',true);
for k=1:4
 values=cell(1,2);mu=nan(numel(S.rel),3,2);se=mu;
 for s=1:2,values{s}=S.X{k,s};mu(:,:,s)=squeeze(mean(values{s},2));se(:,:,s)=squeeze(std(values{s},0,2))/sqrt(S.counts(k,s));end
 draw(values,mu,se,S.rel,S.counts(k,:),fullfile(individual,char(S.ids(k)+"_original_pizza_rules_baseline_subtracted_dff")),[strrep(char(S.ids(k)),'_',' ') ' | dF/F baseline: -2 to -1 s'],false);
end
B=cell2table(rows,'VariableNames',{'recording','flyId','heatOn','sampleIndex','eventTime','subtractedDffBaseline'});writetable(B,fullfile(out,'event_dff_baselines.csv'));
S.baselines=baselines;S.baselineMask=baselineMask;S.methods.normalization='Subtract each event mean dF/F over -2 <= time < -1 s before averaging events within fly, then compute equal-fly mean and across-fly SEM. Forward velocity and yaw are unchanged.';
S.methods.baselineWindowSec=[-2 -1];S.methods.baselineEndExclusive=true;S.methods.validation='Each event pre-baseline mean is zero; original dF/F recoverable by adding baseline; all event counts and kinematic arrays/means/SEM unchanged; between-fly SEM independently verified';
save(fullfile(out,'crunch_data.mat'),'-struct','S','-v7.3');fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(S.methods,PrettyPrint=true));fclose(fid);
result=struct('output',out,'baselineWindowSec',[-2 -1],'baselineSubtractedEvents',height(B),'countsUnchanged',isequal(S.counts,original.counts),'kinematicsUnchanged',true);disp(result);
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
  if row==1,title(ax,names{s},'FontSize',23,'FontWeight','bold');label=sprintf('n=%d events',counts(s));if pooled,label=sprintf('4 flies | %d events',counts(s));end;subtitle(ax,label,'FontSize',15);end
 end
end
text='Thin: events | Thick: mean | Shading: SEM across events';if pooled,text='Thin: fly means | Thick: equal-fly mean | Shading: SEM across flies';end
subtitle(tl,{caption,text},'FontSize',15);assert(numel(findall(fig,'Type','axes'))==6);expected=3*sum(counts);if pooled,expected=24;end;assert(numel(findall(fig,'Tag','thinTrace'))==expected);
exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
