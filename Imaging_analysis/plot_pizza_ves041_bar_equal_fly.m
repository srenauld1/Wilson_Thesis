function result=plot_pizza_ves041_bar_equal_fly(out,includeDark)
% Pool saved events within fly, then average fly means with equal weights.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
if nargin<2,includeDark=false;end
F=readtable(fullfile(out,'bar_fly_mapping.csv'),'Delimiter',',','TextType','string');
M=readtable(fullfile(out,'manifest.csv'),'Delimiter',',','TextType','string');
cohortLabel='bar';mappingName='bar_fly_mapping.csv';stem=fullfile(out,'bar','ves041_bar_pooled_equal_fly_peaks_slowdowns');
if includeDark
 dark=M(M.status=="complete"&M.category=="dark",:);
 for k=1:height(dark)
  name=dark.trial(k);date=extractBefore(name,9);
  if date=="20260217",id="20260217_fly1_inferred";basis="Prior date/fly mapping; same fly as February 17 bar recordings";
  elseif date=="20261020"
   q=regexp(char(dark.trialFolder(k)),'fly(\d+)','tokens','once');assert(~isempty(q));id=date+"_fly"+string(q{1});basis="User-confirmed October 20 grouping and explicit fly number";
  else,error('No dark-trial fly mapping established');end
  row=table(name,id,string(basis),dark.trialFolder(k),'VariableNames',F.Properties.VariableNames);F=[F;row];
 end
 cohortLabel='dark, closed-loop bar, and bar-jump';mappingName='dark_bar_barjump_fly_mapping.csv';
 stem=fullfile(out,'ves041_dark_bar_barjump_pooled_equal_fly_peaks_slowdowns');
 writetable(F,fullfile(out,mappingName));
end
flies=sort(unique(F.flyId));nf=numel(flies);assert(height(F)>0&&nf>0,'Empty source cohort');
inferred=nnz(contains(lower(flies),'inferred'));
identityText=sprintf('%d estimated flies; %d inferred identities; see %s',nf,inferred,mappingName);
labels={'VES041 peaks','All slowdowns','Slowdowns without turns','Slowdowns with turns'};
X=cell(nf,4);counts=zeros(nf,4);trialCounts=zeros(nf,4);rel=[];sourceSettings=[];
trialRows=cell(0,6);
for k=1:height(F)
 name=char(F.trial(k));assert(any(M.trial==string(name)&M.status=="complete"&(M.category=="bar"|(includeDark&M.category=="dark"))));
 S=load(fullfile(out,'analysis_data',[name '_events.mat']),'X','rel','groups','settings');
 if isempty(rel),rel=S.rel;sourceSettings=S.settings;else,assert(isequal(rel,S.rel));assert(isequal(sourceSettings,S.settings));end
 assert(isequal(sort([S.groups{3};S.groups{4}]),S.groups{2}));
 f=find(flies==F.flyId(k));
 for c=1:4
  if isempty(X{f,c}),X{f,c}=S.X{c};else,X{f,c}=cat(2,X{f,c},S.X{c});end
  counts(f,c)=counts(f,c)+size(S.X{c},2);trialCounts(f,c)=trialCounts(f,c)+(size(S.X{c},2)>0);
 end
 trialRows(end+1,:)={string(name),F.flyId(k),size(S.X{1},2),size(S.X{2},2),size(S.X{3},2),size(S.X{4},2)}; %#ok<AGROW>
end
assert(all(counts(:,2)==counts(:,3)+counts(:,4)));
flyMeans=nan(numel(rel),nf,3,4);groupMean=nan(numel(rel),3,4);contributingFlies=zeros(1,4);
for c=1:4
 for f=1:nf
  if counts(f,c)>0,flyMeans(:,f,:,c)=mean(X{f,c},2);end
 end
 contributingFlies(c)=nnz(counts(:,c)>0);groupMean(:,:,c)=squeeze(mean(flyMeans(:,:,:,c),2,'omitnan'));
end
% Explicitly check the equal-fly mean; an absent group never becomes a zero.
for c=1:4
 ids=find(counts(:,c)>0);expected=zeros(numel(rel),3);
 for f=ids',expected=expected+squeeze(mean(X{f,c},2))/numel(ids);end
 q=abs(expected-groupMean(:,:,c));assert(max(q(:))<1e-10);
end
fig=figure('Visible','off','Color','w','Position',[10 10 1800 1300]);cleanup=onCleanup(@()close(fig));tl=tiledlayout(fig,4,4,'TileSpacing','compact','Padding','loose');
lims=cell(1,3);
for v=1:3,z=flyMeans(:,:,v,:);z=z(isfinite(z));lims{v}=[min(z) max(z)];if diff(lims{v})==0,lims{v}=lims{v}+[-.5 .5];else,lims{v}=lims{v}+[-1 1]*.025*diff(lims{v});end,end
flyLabel=replace(flies,'_fly',' F');flyLabel=replace(flyLabel,'_unlabeled_inferred',' *');flyLabel=replace(flyLabel,'_inferred','*');
for row=1:4
 for c=1:4
  ax=nexttile(tl);hold(ax,'on');
  if row==1
   im=imagesc(ax,rel,1:nf,flyMeans(:,:,1,c)');im.AlphaData=isfinite(im.CData);set(ax,'YDir','normal','Color',[.9 .9 .9]);ylim(ax,[.5 nf+.5]);
   yticks(ax,1:nf);yticklabels(ax,flyLabel);ax.YAxis.FontSize=8;caxis(ax,lims{1});colormap(ax,parula(256));cb=colorbar(ax);cb.Label.String='Fly mean \DeltaF/F';
   title(ax,{labels{c},sprintf('%d flies | %d events',contributingFlies(c),sum(counts(:,c)))});ylabel(ax,'Fly (* inferred identity)');
  else
   v=row-1;plot(ax,rel,flyMeans(:,:,v,c),'Color',[.7 .7 .7],'LineWidth',.9);plot(ax,rel,groupMean(:,v,c),'k','LineWidth',2.3);
   names={'Original \DeltaF/F','Forward velocity (mm/s)','Rotational speed |yaw| (deg/s)'};ylabel(ax,names{v});ylim(ax,lims{v});grid(ax,'on');
  end
  xline(ax,0,':','Color',[.8 .2 .2]);xlim(ax,[-2 3]);if row==4,xlabel(ax,'Time from event (s)');end
 end
end
title(tl,['VES041 ' cohortLabel ' trials | equal-fly pooled peaks and slowdowns'],'FontWeight','bold');
subtitle(tl,{sprintf('%d recordings | %d estimated flies (%d inferred identities) | Original dF/F | Turn: accumulated |yaw| >=15 deg over 0:+0.5 s after slowdown',height(F),nf,inferred), ...
 'Events pooled within each fly; gray: fly means; black: equal-weight mean across contributing flies | Columns 3 + 4 partition the slowdown events'},'FontSize',10);
exportgraphics(fig,[stem '.png'],'Resolution',170);savefig_visible_on_open(fig,[stem '.fig']);
rows=cell(0,6);
for c=1:4,for f=1:nf,weight=double(counts(f,c)>0)/contributingFlies(c);rows(end+1,:)={flies(f),string(labels{c}),counts(f,c),trialCounts(f,c),counts(f,c)>0,weight};end,end %#ok<AGROW>
W=cell2table(rows,'VariableNames',{'flyId','eventGroup','events','contributingTrials','included','flyWeight'});
T=cell2table(trialRows,'VariableNames',{'trial','flyId','peaks','slowdowns','slowdownsWithoutTurn','slowdownsWithTurn'});
writetable(W,[stem '_fly_weights.csv']);writetable(T,[stem '_trial_counts.csv']);
methods=struct('source','Saved event arrays; no event redetection or normalization change', ...
 'withinFly',['Pool all eligible events from that fly across its included ' cohortLabel ' recordings, then mean events'], ...
 'betweenFly','Arithmetic mean of contributing fly means, weight 1/N per fly per column', ...
 'absentEvents','Excluded from that column, never replaced by zeros', ...
 'heatmap','One chronological-date-ordered row per fly; each row is that fly event mean', ...
 'partition','Columns 3 and 4 partition events in column 2; their equally weighted group means need not sum to the column 2 mean', ...
 'identities',identityText);
save([stem '.mat'],'flies','flyMeans','groupMean','counts','trialCounts','contributingFlies','rel','F','W','T','methods','sourceSettings','cohortLabel','-v7.3');
fid=fopen([stem '_methods.txt'],'w');assert(fid>=0);
fprintf(fid,'Equal-fly pooled bar figure\n\n%s\n%s\n%s\n%s\n%s\n%s\n\nOriginal source dF/F is preserved, without baseline subtraction or z-scoring. Source event definitions and -2 to +3 s windows are unchanged. Event counts differ among flies; each contributing fly nevertheless has exactly equal weight in the black trace. Within a fly, recordings contribute in proportion to their eligible events. The all-slowdown mean is computed from that fly''s full slowdown set, not from an unweighted average of its turn/no-turn category means.\n',methods.withinFly,methods.betweenFly,methods.absentEvents,methods.heatmap,methods.partition,methods.identities);fclose(fid);
% Verify the saved editable figure contains exactly the requested 16 panels.
assert(numel(findall(fig,'Type','axes'))==16);
result=struct('stem',stem,'flies',nf,'recordings',height(F),'contributingFlies',contributingFlies,'events',sum(counts,1));disp(result);
end
