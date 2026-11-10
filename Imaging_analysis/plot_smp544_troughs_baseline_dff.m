function plot_smp544_troughs_baseline_dff()
% Preserve saved trough selection; display dF/F with event baseline subtraction.
repo=fileparts(fileparts(mfilename('fullpath')));out=fullfile(repo,'outputs','smp544_troughs_per_fly_20261102');
load(fullfile(out,'smp544_right_troughs_per_fly.mat'),'t','means','counts','labels','trials','flyIds');
T=readtable('/Users/sophiarenauld/stacks/smp544_behavior_quantification/ramp_turn_go_bidirectional/smp544_ramp_turn_go_events.csv','TextType','string');
X=cell(3,2);baseline=t>=-3&t<=-2;
for i=1:numel(trials)
 files=dir(fullfile('/Users/sophiarenauld/stacks',trials(i),'*_or_z0z1_s_.mat'));q=load(fullfile(files(1).folder,files(1).name),'dq','roi');
 time=q.dq(1).t(:);d=q.roi.ts;fs=1/median(diff(time));
 events=T.time(T.trial==trials(i)&T.selectionArm=="right-physiology");
 for v=1:2
  y=smoothdata(d(v,:)','gaussian',max(3,round(.25*fs)));y=smoothdata(y,'gaussian',max(3,round(1.5*fs)));
  for k=1:numel(events)
   z=interp1(time,y,events(k)+t,'linear',NaN);z=z-mean(z(baseline),'omitnan');X{flyIds(i),v}(:,end+1)=z;
  end
 end
end
for j=1:3
 assert(size(X{j,1},2)==counts(j));
 for v=1:2,means(:,v,j)=mean(X{j,v},2,'omitnan');end
end
colors=[.05 .32 .75;.05 .32 .75;0 0 0;.85 .12 .10];
labels={'Jul 20, fly 2','Jul 21, fly 1','Jul 22, fly 1'};
ylabels={'Right SMP544 (dF/F change)','Left SMP544 (dF/F change)','Forward speed (mm/s)','Rotational speed (deg/s)'};
fig=figure('Visible','off','Color','w','Position',[100 100 730 900]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
for v=[1 3 4]
 a=nexttile(tl);hold(a,'on');
 c=colors(v,:);m=mean(means(:,v,:),3,'omitnan');se=std(means(:,v,:),0,3,'omitnan')/sqrt(3);
 fill(a,[t;flipud(t)],[m-se;flipud(m+se)],c,'FaceAlpha',.18,'EdgeColor','none','DisplayName','SEM across flies');
 for j=1:3,plot(a,t,means(:,v,j),'-','Color',.65+.35*c,'LineWidth',.8,'HandleVisibility','off');end
 plot(a,t,m,'-','Color',c,'LineWidth',2.8,'DisplayName','Mean (3 flies)');
 xlim(a,[-5 5]);ylabel(a,ylabels{v});set(a,'FontSize',14,'LineWidth',1,'Box','off');grid(a,'off');
 if v==1,title(a,'Right dF/F-selected troughs','FontSize',18);legend(a,'Location','northoutside','FontSize',12,'Box','off');end
 if v==4,xlabel(a,'Time from right SMP544 trough (s)');end
end
savefig(fig,fullfile(out,'smp544_right_only_troughs_baseline_dff_mean_sem.fig'));
exportgraphics(fig,fullfile(out,'smp544_right_only_troughs_baseline_dff_mean_sem.png'),'Resolution',200);
exportgraphics(fig,fullfile(out,'smp544_right_only_troughs_baseline_dff_mean_sem.pdf'),'ContentType','vector');
save(fullfile(out,'smp544_right_only_troughs_baseline_dff_mean_sem.mat'),'t','means','counts','labels','trials','flyIds');
assert(all(arrayfun(@(h)strcmp(h.LineStyle,'-'),findall(fig,'Type','line'))),'Unexpected non-solid line');
close(fig);disp(table(string(labels'),counts,'VariableNames',{'Fly','RightTroughEvents'}));
end
