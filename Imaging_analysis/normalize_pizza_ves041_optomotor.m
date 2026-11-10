function result=normalize_pizza_ves041_optomotor(root)
% Presentation summaries: whole-recording neural z scores, equal fly mean/SEM.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
out=fullfile(root,'gratings','optomotor_exploration');
S=load(fullfile(out,'optomotor_summary.mat'),'R','rel','settings');R=S.R;rel=S.rel;
assert(numel(R)==3,'This cohort expects three distinct-date flies');
dates=extractBefore(string({R.trial}),9);assert(numel(unique(dates))==numel(R),'Repeated flies require pooling before normalization');
normalization=table();
for k=1:numel(R)
 d=R(k).d;valid=R(k).valid&isfinite(d.dff);
 mu=mean(d.dff(valid));sigma=std(d.dff(valid),0);assert(isfinite(sigma)&&sigma>0);
 % Apply the identical full-recording transform to every event, not per event.
 R(k).onset(:,:,4)=(R(k).onset(:,:,4)-mu)/sigma;
 R(k).turn(:,:,4)=(R(k).turn(:,:,4)-mu)/sigma;
 normalization=[normalization;table(R(k).trial,mu,sigma,nnz(valid),'VariableNames',{'trial','fullRecordingMeanDff','fullRecordingSdDff','validSamples'})]; %#ok<AGROW>
 z=(d.dff(valid)-mu)/sigma;assert(abs(mean(z))<1e-10&&abs(std(z,0)-1)<1e-10);
end
methods=struct('normalization','Neural dF/F only: (dF/F - full valid recording mean) / full valid recording sample SD, within each fly; motion and stationary periods included', ...
 'eventBaseline','No per-event baseline subtraction; z=0 is that fly full-recording mean', ...
 'kinematics','Original physical units; stimulus and yaw signs aligned to stimulus direction', ...
 'mean','First average eligible bouts within each fly, then equal-weight average contributing fly means', ...
 'error','Shaded +/-1 SEM across fly means: sample SD(fly means)/sqrt(number of finite contributing flies); not across bouts; SEM undefined for fewer than 2 flies', ...
 'selection','Reuse saved motion bouts and turn labels exactly; no redetection', ...
 'display','Four rows with shared row limits across columns; colored fly means, black equal-fly mean, shaded between-fly SEM; no temporal smoothing', ...
 'cohort','Three distinct recording dates treated as three flies, as in original optomotor summary');
summaries=cell(1,2);
for mode=1:2
 byDirection=mode==2;
 titles={'All stimulus motion','Motion + matched turn','Motion + weak/no matched turn','Matched turn peaks'};
 alignment={'Stimulus onset','Stimulus onset','Stimulus onset','Turn peak'};
 stem=fullfile(out,'ves041_optomotor_equal_fly_summary_zscored');
 if byDirection
  titles={'Positive stimulus motion','Negative stimulus motion','Positive-stimulus turns','Negative-stimulus turns'};
  alignment={'Stimulus onset','Stimulus onset','Turn peak','Turn peak'};
  stem=fullfile(out,'ves041_optomotor_by_stimulus_direction_zscored');
 end
 flyMeans=nan(numel(rel),numel(R),4,4);counts=zeros(numel(R),4);
 for c=1:4
  for k=1:numel(R)
   r=R(k);select=r.E.onsetWindowValid;turnAligned=c==4;
   if byDirection
    turnAligned=c>=3;if turnAligned,select=r.E.turnWindowValid;end
    direction=1;if mod(c,2)==0,direction=-1;end
    select=select&r.E.stimulusSign==direction;
   else
    if c==2,select=select&r.E.compensatory;elseif c==3,select=select&~r.E.compensatory;elseif c==4,select=r.E.turnWindowValid;end
   end
   counts(k,c)=nnz(select);
   if any(select)
    if turnAligned,X=r.turn(:,select,1:4);else,X=r.onset(:,select,1:4);end
    flyMeans(:,k,:,c)=mean(X,2,'omitnan');
   end
  end
 end
 n=sum(isfinite(flyMeans),2);groupMean=mean(flyMeans,2,'omitnan');sem=std(flyMeans,0,2,'omitnan')./sqrt(n);sem(n<2)=NaN;
 assert(all(sum(counts,1)>0));
 colors=[.15 .43 .70;.86 .34 .17;.48 .31 .65];
 fig=figure('Visible','off','Color','w','Position',[20 20 1700 1050]);
 tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
 labels={'Stimulus speed\newline(sign aligned, deg/s)','Fly yaw\newline(sign aligned, deg/s)','Forward speed\newline(mm/s)','VES041 activity\newline(z-scored \DeltaF/F)'};
 legendHandles=[];
 for v=1:4
  lo=flyMeans(:,:,v,:);a=groupMean(:,:,v,:)-sem(:,:,v,:);b=groupMean(:,:,v,:)+sem(:,:,v,:);
  limits=[min([lo(:);a(:)],[],'omitnan'),max([lo(:);b(:)],[],'omitnan')];
  if v==4,limits=[min(limits(1),0),max(limits(2),0)];end
  if diff(limits)==0,limits=limits+[-.5 .5];end
  limits=limits+[-1 1]*.08*diff(limits);
  for c=1:4
   ax=nexttile(tl);hold(ax,'on');
   if v==4,ax.Color=[.98 .975 .95];end
   m=groupMean(:,1,v,c);e=sem(:,1,v,c);ok=isfinite(m)&isfinite(e);
   % SEM is based only on independent fly means, never individual event samples.
   if any(ok),fill(ax,[rel(ok);flipud(rel(ok))],[m(ok)-e(ok);flipud(m(ok)+e(ok))],[.25 .25 .25],'FaceAlpha',.13,'EdgeColor','none','HandleVisibility','off');end
   h=gobjects(1,numel(R)+1);
   for k=1:numel(R),h(k)=plot(ax,rel,flyMeans(:,k,v,c),'Color',colors(k,:),'LineWidth',1.1);end
   h(end)=plot(ax,rel,m,'Color',[.08 .08 .08],'LineWidth',2.4);
   if v==1&&c==1,legendHandles=h;end
   xline(ax,0,':','Color',[.15 .45 .4],'LineWidth',1.2,'HandleVisibility','off');
   if v==4,yline(ax,0,':','Color',[.6 .6 .6],'HandleVisibility','off');end
   xlim(ax,[-2 4]);ylim(ax,limits);xticks(ax,-2:2:4);
   set(ax,'Box','off','TickDir','out','FontSize',10,'YGrid','on','GridAlpha',.10,'LineWidth',.7);
   if c==1,ylabel(ax,labels{v},'FontSize',11);end
   if v==1,title(ax,{titles{c},sprintf('%d bouts | n per fly: %d / %d / %d',sum(counts(:,c)),counts(:,c))},'FontSize',11);end
   if v==4,xlabel(ax,['Time from ' lower(alignment{c}) ' (s)']);else,ax.XTickLabel=[];end
  end
 end
 flyLabels={'Aug 12, fly 3','Aug 13, fly 1','Aug 19, fly 1','Equal-fly mean +/- SEM'};
 lg=legend(legendHandles,flyLabels,'Orientation','horizontal','Box','off','FontSize',11);lg.Layout.Tile='south';
 title(tl,'VES041 during optomotor stimulation | 3 flies','FontSize',17,'FontWeight','bold');
 subtitle(tl,{'Each fly normalized using its full recording; kinematics retain physical units. Shading = SEM across fly means.', ...
  'Matched turn: >=50 deg/s in stimulus direction. Weak/no matched turn can include smaller or opposite turns. No event-baseline subtraction.'},'FontSize',11);
 exportgraphics(fig,[stem '.png'],'Resolution',180);savefig_visible_on_open(fig,[stem '.fig']);close(fig);
 summary=struct('flyMeans',flyMeans,'mean',groupMean,'sem',sem,'contributingFlies',n,'boutCounts',counts,'columnTitles',{titles},'alignment',{alignment});
 save([stem '.mat'],'summary','normalization','rel','methods');summaries{mode}=summary;
end
writetable(normalization,fullfile(out,'optomotor_zscore_normalization.csv'));
fid=fopen(fullfile(out,'optomotor_zscore_methods.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('folder',out,'flies',numel(R),'mainBoutCounts',sum(summaries{1}.boutCounts,1));disp(result);
end
