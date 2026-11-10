function plot_smp544_troughs_per_fly()
% Replot saved event means without changing the original event selection.
src='/Users/sophiarenauld/stacks/smp544_behavior_quantification/ramp_turn_go_bidirectional';
out=fullfile(fileparts(fileparts(mfilename('fullpath'))),'outputs','smp544_troughs_per_fly_20261102');
if ~isfolder(out),mkdir(out);end
T=readtable(fullfile(src,'smp544_ramp_turn_go_events.csv'),'TextType','string');
trials=unique(T.trial,'stable');flyIds=[1 1 2 2 2 3];
assert(numel(trials)==numel(flyIds),'Unexpected recording inventory');
t=(-5:.05:5)';sums=zeros(numel(t),4,3);counts=zeros(3,1);
for i=1:numel(trials)
 n=sum(T.trial==trials(i)&T.selectionArm=="right-physiology");
 f=openfig(fullfile(src,trials(i)+'_bidirectional.fig'),'invisible');
 cleanup=onCleanup(@()close(f));
 ax=findall(f,'Type','axes');pos=vertcat(ax.Position);left=ax(pos(:,1)<min(pos(:,1))+.03);
 p=vertcat(left.Position);[~,order]=sort(p(:,2),'descend');left=left(order);
 assert(numel(left)==4,'Expected four left-column panels');
 for v=1:4
  lines=findall(left(v),'Type','line');chosen=[];
  for k=1:numel(lines)
   if startsWith(string(lines(k).DisplayName),'Events'),chosen=lines(k);break;end
  end
  assert(~isempty(chosen),'Missing event mean');
  y=interp1(chosen.XData,chosen.YData,t,'linear',NaN);
  sums(:,v,flyIds(i))=sums(:,v,flyIds(i))+n*y;
 end
 counts(flyIds(i))=counts(flyIds(i))+n;
 clear cleanup
end
means=sums./reshape(counts,1,1,3);colors=[.05 .32 .75;.05 .32 .75;0 0 0;.85 .12 .10];
labels={'Jul 20, fly 2','Jul 21, fly 1','Jul 22, fly 1'};
ylabels={'Right SMP544 (robust z)','Left SMP544 (robust z)','Forward speed (mm/s)','Rotational speed (deg/s)'};
fig=figure('Visible','off','Color','w','Position',[100 100 730 1150]);
tl=tiledlayout(fig,4,1,'TileSpacing','compact','Padding','compact');
for v=1:4
 a=nexttile(tl);hold(a,'on');
 c=colors(v,:);m=mean(means(:,v,:),3,'omitnan');se=std(means(:,v,:),0,3,'omitnan')/sqrt(3);
 fill(a,[t;flipud(t)],[m-se;flipud(m+se)],c,'FaceAlpha',.18,'EdgeColor','none','DisplayName','SEM across flies');
 for j=1:3,plot(a,t,means(:,v,j),'-','Color',.65+.35*c,'LineWidth',.8,'HandleVisibility','off');end
 plot(a,t,m,'-','Color',c,'LineWidth',2.8,'DisplayName','Mean (3 flies)');
 xlim(a,[-5 5]);ylabel(a,ylabels{v});set(a,'FontSize',14,'LineWidth',1,'Box','off');grid(a,'off');
 if v==1,title(a,'Right dF/F-selected troughs','FontSize',18);legend(a,'Location','northoutside','FontSize',12,'Box','off');end
 if v==4,xlabel(a,'Time from right SMP544 trough (s)');end
end
savefig(fig,fullfile(out,'smp544_right_troughs_mean_sem.fig'));
exportgraphics(fig,fullfile(out,'smp544_right_troughs_mean_sem.png'),'Resolution',200);
exportgraphics(fig,fullfile(out,'smp544_right_troughs_mean_sem.pdf'),'ContentType','vector');
save(fullfile(out,'smp544_right_troughs_per_fly.mat'),'t','means','counts','labels','trials','flyIds');
assert(all(arrayfun(@(h)strcmp(h.LineStyle,'-'),findall(fig,'Type','line'))),'Unexpected non-solid line');
close(fig);disp(table(string(labels'),counts,'VariableNames',{'Fly','RightTroughEvents'}));
end
