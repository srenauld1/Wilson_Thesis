function plot_other_heat_flies_calcium_and_turns(out)
% Separate true neural-peak averages from slowdown/turn averages per fly.
if nargin<1,out=fullfile(pwd,'outputs','ves041_other_heat_flies_individual_20261105');end
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier';
S=load(fullfile(out,'individual_fly_data.mat'),'C');C=S.C;P=struct([]);bounds=cell(1,3);rows=cell(0,6);
for k=1:2
 A=load(fullfile(root,char(C(k).id),'analysis.mat'));rel=A.rel(:);base=rel>=-2&rel< -1;
 p=struct('id',C(k).id,'X',{{}},'mu',[],'sem',[],'counts',zeros(1,4));
 for c=1:4
  if c<=2
   idx=A.groups{1};on=A.d.heat(idx)>A.meta.threshold;X=A.X{1}(:,on==(c==2),:);
   X(:,:,1)=X(:,:,1)-mean(X(base,:,1),1)+1;
  else,X=C(k).X{c-2};end
  p.X{c}=X;p.counts(c)=size(X,2);assert(p.counts(c)>1);
  p.mu(:,:,c)=squeeze(mean(X,2));p.sem(:,:,c)=squeeze(std(X,0,2))/sqrt(size(X,2));
  for v=1:3,bounds{v}=[bounds{v};p.mu(:,v,c)-p.sem(:,v,c);p.mu(:,v,c)+p.sem(:,v,c)];end %#ok<AGROW>
  on=mod(c,2)==0;kind="Calcium peak";if c>2,kind="Slowdown + turn";end
  rot=X(:,:,3);core=rel>=-.5&rel<=.5;
  rows(end+1,:)={p.id,kind,on,p.counts(c),max(mean(rot,2)),median(max(rot(core,:),[],1))}; %#ok<AGROW>
 end
 if k==1,P=p;else,P(k)=p;end
end
lims=cell(1,3);steps=[.5 2 50];
for v=1:3,lims{v}=[floor(min([bounds{v};0])/steps(v))*steps(v),ceil(max(bounds{v})/steps(v))*steps(v)];end
colors={[.9 .04 .03],[.1 .1 .1],[.05 .2 .8]};labels={'dF/F (baseline set to 1)','Forward velocity (mm/s)','Rotational speed (deg/s)'};
for k=1:2
 p=P(k);fig=figure('Visible','off','Color','white','Position',[50 50 1800 1050]);tl=tiledlayout(fig,3,4,'TileSpacing','compact','Padding','loose');
 for v=1:3
  for c=1:4
   ax=nexttile(tl);hold(ax,'on');mu=p.mu(:,v,c);se=p.sem(:,v,c);
   fill(ax,[rel;flipud(rel)],[mu-se;flipud(mu+se)],colors{v},'FaceAlpha',.18,'EdgeColor','none');
   plot(ax,rel,mu,'Color',colors{v},'LineWidth',2.4);xline(ax,0,':','Color',[.5 .5 .5]);
   set(ax,'XLim',[-2 2],'YLim',lims{v},'FontName','Arial','FontSize',15,'Box','off','TickDir','out');
   if c==1,ylabel(ax,labels{v},'FontSize',18);end
   if v==3
    if c<=2,xlabel(ax,'Time from calcium peak (s)','FontSize',14);else,xlabel(ax,'Time from deceleration (s)','FontSize',14);end
   else,ax.XTickLabel=[];end
   if v==1
    kind='Calcium peaks';if c>2,kind='Slowdowns + turns';end
    condition='Heat off';if mod(c,2)==0,condition='Heat on';end
    title(ax,{kind,sprintf('%s | n=%d',condition,p.counts(c))},'FontSize',18);
   end
  end
 end
 names={'April 27, fly 2','October 19, fly 1'};
 title(tl,names{k},'FontName','Arial','FontSize',24);
 subtitle(tl,{'Each fly separately | original saved event selections | mean +/- event SEM', ...
  'Neural baseline -2:-1 s set to 1; heat-state split at event time'},'FontSize',13);
 dest=fullfile(out,char(p.id+"_calcium_peaks_and_turns_by_heat_state"));
 print(fig,[dest '.svg'],'-dsvg','-painters');exportgraphics(fig,[dest '.png'],'Resolution',200);
 set(fig,'Visible','on');savefig(fig,[dest '.fig']);close(fig);
end
T=cell2table(rows,'VariableNames',{'recording','eventType','heatOn','nEvents','maxMeanRotationalSpeed','medianIndividualCoreRotationPeak'});
writetable(T,fullfile(out,'calcium_and_turn_event_counts.csv'));save(fullfile(out,'calcium_and_turn_data.mat'),'P','rel','lims','T','-v7.3');disp(T);
end
