function result=plot_pizza_ves041_rotation_columns(out)
% Add behavior-selected rotation peaks to the combined equal-fly figure.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
sourceStem=fullfile(out,'ves041_dark_bar_barjump_pooled_equal_fly_peaks_slowdowns');
S=load([sourceStem '.mat']);N=load([sourceStem '_normalized.mat']);assert(isequal(S.flies,N.flies));
flies=S.flies;nf=numel(flies);rel=S.rel;F=S.F;
counts=[S.counts zeros(nf,2)];flyMeans=nan(numel(rel),nf,3,6);flyMeans(:,:,:,1:4)=S.flyMeans;
rotationArrays=cell(nf,2);events=cell(0,9);detectors=cell(0,4);
for k=1:height(F)
 trial=F.trial(k);f=find(flies==F.flyId(k));A=load(fullfile(out,'analysis_data',trial+"_events.mat"),'d');d=A.d;
 [idx,drop,preSpeed,postSpeed,info]=pizza_ves041_rotation_peaks(d);
 detectors(end+1,:)={trial,info.prominence,info.minimumHeight,info.minimumDistanceSec}; %#ok<AGROW>
 for g=1:2
  selected=idx((drop>=1)==(g==2));X=nan(numel(rel),numel(selected),3);vars={d.dff,d.fwd,abs(d.yaw)};
  for j=1:numel(selected),for v=1:3,X(:,j,v)=interp1(d.t,vars{v},d.t(selected(j))+rel,'linear',NaN);end,end
  assert(all(isfinite(X(:))));
  rotationArrays{f,g}=cat(2,rotationArrays{f,g},X);counts(f,g+4)=counts(f,g+4)+numel(selected);
 end
 for j=1:numel(idx)
  events(end+1,:)={trial,flies(f),idx(j),d.t(idx(j)),preSpeed(j),postSpeed(j),drop(j),drop(j)>=1,abs(d.yaw(idx(j)))}; %#ok<AGROW>
 end
end
for f=1:nf,for g=1:2,if counts(f,g+4)>0,flyMeans(:,f,:,g+4)=mean(rotationArrays{f,g},2);end,end,end
contributingFlies=sum(counts>0,1);groupMean=nan(numel(rel),3,6);
for c=1:6,groupMean(:,:,c)=squeeze(mean(flyMeans(:,:,:,c),2,'omitnan'));end
assert(isequaln(flyMeans(:,:,:,1:4),S.flyMeans));assert(isequaln(groupMean(:,:,1:4),S.groupMean));
assert(sum(sum(counts(:,5:6)))==size(events,1));assert(all(counts(:,2)==counts(:,3)+counts(:,4)));
normalizedFlyMeans=flyMeans;
for f=1:nf,for v=1:3,normalizedFlyMeans(:,f,v,:)=(flyMeans(:,f,v,:)-N.center(f,v))/N.scale(f,v);end,end
q=abs(normalizedFlyMeans(:,:,:,1:4)-N.flyMeans);assert(all(q(isfinite(q))<1e-10));
normalizedGroupMean=nan(size(groupMean));for c=1:6,normalizedGroupMean(:,:,c)=squeeze(mean(normalizedFlyMeans(:,:,:,c),2,'omitnan'));end
E=cell2table(events,'VariableNames',{'trial','flyId','sampleIndex','eventTime','forwardBefore','forwardAfter','forwardDrop','forwardDropAtLeast1','rawAbsoluteYawAtEvent'});
D=cell2table(detectors,'VariableNames',{'trial','rotPeakProminenceDegS','rotPeakMinimumHeightDegS','rotPeakSpacingSec'});
stem=[sourceStem '_with_rotation_peaks'];
render(flyMeans,groupMean,false,[stem '.png'],[stem '.fig']);
render(normalizedFlyMeans,normalizedGroupMean,true,[stem '_normalized.png'],[stem '_normalized.fig']);
writetable(E,[stem '_rotation_events.csv']);writetable(D,[stem '_rotation_detector_settings.csv']);
weightRows=cell(0,4);
for c=1:6,for f=1:nf,weightRows(end+1,:)={flies(f),c,counts(f,c),double(counts(f,c)>0)/contributingFlies(c)};end,end %#ok<AGROW>
W=cell2table(weightRows,'VariableNames',{'flyId','column','events','flyWeight'});writetable(W,[stem '_fly_weights.csv']);
methods=struct('rotationDetector','Absolute yaw after 0.12-second Gaussian smoothing; findpeaks height >=20 deg/s, prominence >=max(20,2*1.4826*MAD(abs(yaw)-5-second moving median)), minimum spacing 1.5 s', ...
 'forwardSplit','0.25-second median-smoothed forward speed; median(-0.75:-0.25 s) minus median(+0.25:+0.75 s); column 5 <1 mm/s, column 6 >=1 mm/s', ...
 'selection','Rotational peaks selected from behavior only, without calcium-amplitude, pre-speed, net-turn, or post-speed requirements; full valid -2:+3 s windows required', ...
 'averaging','Events pooled across included dark/bar/bar-jump recordings within fly, then contributing fly means averaged equally; flies without events excluded from that column', ...
 'normalization','Same full-recording within-fly means and sample SDs as the existing combined normalized 4x4 figure; first four columns preserve exact event selections', ...
 'labels','Column 5 may contain small decreases (<1 mm/s), unchanged speed, or acceleration. Rotational peaks are aligned to peak absolute yaw, not maximum deceleration.');
save([stem '.mat'],'flies','F','rel','counts','contributingFlies','flyMeans','groupMean','normalizedFlyMeans','normalizedGroupMean','E','D','W','methods','sourceStem','-v7.3');
fid=fopen([stem '_methods.txt'],'w');assert(fid>=0);keys=fieldnames(methods);for j=1:numel(keys),fprintf(fid,'%s: %s\n',keys{j},methods.(keys{j}));end;fprintf(fid,'\nVerified: first four raw columns identical to source; first four normalized columns agree to 1e-10; columns 5/6 partition every retained rotation peak; all 24 data axes present in each figure.\n');fclose(fid);
result=struct('stem',stem,'recordings',height(F),'flies',nf,'contributingFlies',contributingFlies,'events',sum(counts,1));disp(result);

 function render(Y,G,isNormalized,pngPath,figPath)
  fig=figure('Visible','off','Color','w','Position',[10 10 2500 1350]);cleanup=onCleanup(@()close(fig));tl=tiledlayout(fig,4,6,'TileSpacing','compact','Padding','loose');
  labels={{'VES041 peaks'},{'All slowdowns'},{'Slowdowns without turns'},{'Slowdowns with turns'},{'Rotation peaks','Forward drop <1 mm/s'},{'Rotation peaks','Forward drop >=1 mm/s'}};
  limits=cell(1,3);
  if isNormalized,z=Y(isfinite(Y));bound=max(1,ceil(2*max(abs(z)))/2);for v=1:3,limits{v}=[-bound bound];end
  else,for v=1:3,z=Y(:,:,v,:);z=z(isfinite(z));limits{v}=[min(z) max(z)];if diff(limits{v})==0,limits{v}=limits{v}+[-.5 .5];else,limits{v}=limits{v}+[-1 1]*.025*diff(limits{v});end,end,end
  flyLabel=replace(flies,'_fly',' F');flyLabel=replace(flyLabel,'_unlabeled_inferred',' *');flyLabel=replace(flyLabel,'_inferred','*');
  for row=1:4,for c=1:6
   ax=nexttile(tl);hold(ax,'on');
   if row==1
    im=imagesc(ax,rel,1:nf,Y(:,:,1,c)');im.AlphaData=isfinite(im.CData);set(ax,'YDir','normal','Color',[.9 .9 .9]);ylim(ax,[.5 nf+.5]);yticks(ax,1:nf);yticklabels(ax,flyLabel);ax.YAxis.FontSize=8;
    caxis(ax,limits{1});cb=colorbar(ax);cb.Label.String='Fly mean \DeltaF/F';if isNormalized,cb.Label.String='Fly mean \DeltaF/F (z-score)';end
    title(ax,[labels{c} {sprintf('%d flies | %d events',contributingFlies(c),sum(counts(:,c)))}]);ylabel(ax,'Fly (* inferred identity)');
   else
    v=row-1;plot(ax,rel,Y(:,:,v,c),'Color',[.7 .7 .7],'LineWidth',.9);plot(ax,rel,G(:,v,c),'k','LineWidth',2.3);ylim(ax,limits{v});grid(ax,'on');
    if isNormalized,names={'\DeltaF/F (z-score)','Forward velocity (z-score)','|Yaw| speed (z-score)'};yline(ax,0,':','Color',[.45 .45 .45]);else,names={'Original \DeltaF/F','Forward velocity (mm/s)','|Yaw| speed (deg/s)'};end
    ylabel(ax,names{v});
   end
   xline(ax,0,':','Color',[.8 .2 .2]);xlim(ax,[-2 3]);if row==4,xlabel(ax,'Time from event (s)');end
  end,end
  if isNormalized
   cmap=[linspace(.15,1,128)' linspace(.3,1,128)' ones(128,1);ones(128,1) linspace(1,.2,128)' linspace(1,.15,128)'];colormap(fig,cmap);mode='within-fly z-scores';
  else,colormap(fig,parula(256));mode='original units';end
  title(tl,['VES041 dark, bar, and bar jumps | equal fly weight | ' mode],'FontWeight','bold');
  subtitle(tl,{'12 recordings | 8 estimated flies (3 inferred) | Gray: fly means; black: equal-fly mean | Slowdown turn: >=15 deg over 0:+0.5 s', ...
   'Columns 5/6: behavior-selected absolute-yaw peaks, split by forward drop before versus after peak (1 mm/s threshold); no calcium-based selection'},'FontSize',10);
  assert(numel(findall(fig,'Type','axes'))==24);exportgraphics(fig,pngPath,'Resolution',170);savefig_visible_on_open(fig,figPath);
 end
end
