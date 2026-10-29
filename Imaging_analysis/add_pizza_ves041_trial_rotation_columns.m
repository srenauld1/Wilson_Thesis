function report=add_pizza_ves041_trial_rotation_columns(out)
% Extend every included no-heat trial to 4x6 using the pooled rotation detector.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
M=readtable(fullfile(out,'manifest.csv'),'Delimiter',',','TextType','string');
included=find(M.status=="complete");assert(numel(included)==17,'Unexpected included cohort');
backup=fullfile(out,'analysis_data','original_4x4_figures');if ~isfolder(backup),mkdir(backup);end
M.rotationPeaksDropBelow1=nan(height(M),1);M.rotationPeaksDropAtLeast1=nan(height(M),1);
pooled=readtable(fullfile(out,'ves041_dark_bar_barjump_pooled_equal_fly_peaks_slowdowns_with_rotation_peaks_rotation_events.csv'),'Delimiter',',','TextType','string');
rows=cell(0,5);
for k=included'
 trial=M.trial(k);S=load(fullfile(out,'analysis_data',trial+"_events.mat"),'d','X','groups','rel');
 [idx,drop,before,after,detector]=pizza_ves041_rotation_peaks(S.d);
 X=S.X(1:4);groups=S.groups(1:4);d=S.d;rel=S.rel;
 for g=1:2
  groups{g+4}=idx((drop>=1)==(g==2));X{g+4}=nan(numel(rel),numel(groups{g+4}),3);vars={d.dff,d.fwd,abs(d.yaw)};
  for j=1:numel(groups{g+4}),for v=1:3,X{g+4}(:,j,v)=interp1(d.t,vars{v},d.t(groups{g+4}(j))+rel,'linear',NaN);end,end
  assert(all(isfinite(X{g+4}(:))));
 end
 for c=1:4,assert(isequal(X{c},S.X{c}));assert(isequal(groups{c},S.groups{c}));end
 assert(isequal(sort([groups{5};groups{6}]),idx));assert(isempty(intersect(groups{5},groups{6})));
 % Match the previously plotted pooled events exactly for its 12 source trials.
 p=pooled(pooled.trial==trial,:);
 if ~isempty(p)
  assert(isequal(p.sampleIndex,idx));assert(max(abs(p.forwardDrop-drop))<1e-10);
 end
 stem=char(M.eventFigureStem(k));[~,name]=fileparts(stem);
 for ext={'.png','.fig'}
  dest=fullfile(backup,[name ext{1}]);if ~isfile(dest),copyfile([stem ext{1}],dest);end
 end
 render(X,rel,char(trial),stem);
 events=table(idx,d.t(idx),before,after,drop,drop>=1,'VariableNames',{'sampleIndex','eventTime','forwardBefore','forwardAfter','forwardDrop','forwardDropAtLeast1'});
 writetable(events,fullfile(out,'analysis_data',trial+"_rotation_events.csv"));
 save(fullfile(out,'analysis_data',trial+"_events_with_rotation_peaks.mat"),'X','groups','rel','events','detector','-v7.3');
 M.rotationPeaksDropBelow1(k)=numel(groups{5});M.rotationPeaksDropAtLeast1(k)=numel(groups{6});
 rows(end+1,:)={trial,numel(groups{5}),numel(groups{6}),height(events),string(stem)}; %#ok<AGROW>
 fprintf('%s: rotation peaks %d below 1 + %d at least 1 mm/s\n',trial,numel(groups{5}),numel(groups{6}));
end
writetable(M,fullfile(out,'manifest.csv'));
report=cell2table(rows,'VariableNames',{'trial','rotationPeaksDropBelow1','rotationPeaksDropAtLeast1','allRotationPeaks','figureStem'});
writetable(report,fullfile(out,'individual_rotation_columns_report.csv'));
settings=jsondecode(fileread(fullfile(out,'analysis_settings.json')));
settings.individualFigureLayout='4 rows x 6 columns; original four groups followed by rotation peaks with forward drop <1 and >=1 mm/s';
settings.rotationDetector='Absolute yaw after 0.12 s Gaussian smoothing; height >=20 deg/s, prominence >=max(20,2*1.4826*MAD(abs(yaw)-5 s moving median)), spacing >=1.5 s';
settings.rotationForwardSplit='0.25 s median-smoothed forward; median(-0.75:-0.25 s) minus median(+0.25:+0.75 s), threshold 1 mm/s';
fid=fopen(fullfile(out,'analysis_settings.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(settings,PrettyPrint=true));fclose(fid);
fid=fopen(fullfile(out,'README.txt'),'a');assert(fid>=0);fprintf(fid,'\nIndividual figure extension: all 17 included no-heat trial peak/slowdown figures are now 4x6 at their original filenames and folder locations. Columns 5 and 6 align to rotational-speed peaks and split forward-speed drop at 1 mm/s. The same detector and complete-window rules as the pooled 4x6 figure are used. Rows remain original dF/F heatmap, original dF/F traces, forward velocity, and absolute yaw speed. The first four groups are unchanged. Original 4x4 figure files are preserved under analysis_data/original_4x4_figures. Rotational events and six-column arrays are saved separately; the original four-column event data and flat paths are unchanged.\n');fclose(fid);
disp(report);
end

function render(X,rel,trial,stem)
fig=figure('Visible','off','Color','w','Position',[10 10 2500 1350]);cleanup=onCleanup(@()close(fig));tl=tiledlayout(fig,4,6,'TileSpacing','compact','Padding','loose');
labels={{'VES041 peaks'},{'All slowdowns'},{'Slowdowns without turns'},{'Slowdowns with turns'},{'Rotation peaks','Forward drop <1 mm/s'},{'Rotation peaks','Forward drop >=1 mm/s'}};
lims=cell(1,3);
for v=1:3
 z=[];for c=1:6,q=X{c}(:,:,v);z=[z;q(:)];end;z=z(isfinite(z)); %#ok<AGROW>
 if isempty(z),lims{v}=[0 1];else,lims{v}=[min(z) max(z)];if diff(lims{v})==0,lims{v}=lims{v}+[-.5 .5];else,lims{v}=lims{v}+[-1 1]*.025*diff(lims{v});end,end
end
for row=1:4,for c=1:6
 ax=nexttile(tl);hold(ax,'on');n=size(X{c},2);
 if row==1
  if n>0,imagesc(ax,rel,1:n,X{c}(:,:,1)');set(ax,'YDir','normal');ylim(ax,[.5 n+.5]);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');ylim(ax,[.5 1.5]);end
  caxis(ax,lims{1});cb=colorbar(ax);cb.Label.String='Original \DeltaF/F';ylabel(ax,'Event (chronological)');title(ax,[labels{c} {sprintf('n=%d',n)}]);
 else
  v=row-1;if n>0,plot(ax,rel,X{c}(:,:,v),'Color',[.8 .8 .8],'LineWidth',.4);plot(ax,rel,mean(X{c}(:,:,v),2),'k','LineWidth',2);else,text(ax,.5,.5,'No eligible events','Units','normalized','HorizontalAlignment','center');end
  names={'Original \DeltaF/F','Forward velocity (mm/s)','|Yaw| speed (deg/s)'};ylabel(ax,names{v});ylim(ax,lims{v});grid(ax,'on');
 end
 xline(ax,0,':','Color',[.8 .2 .2]);xlim(ax,[-2 3]);if row==4,xlabel(ax,'Time from event (s)');end
end,end
colormap(fig,parula(256));title(tl,strrep(trial,'_',' '),'FontWeight','bold');
subtitle(tl,{'Original dF/F | Thin: events; black: mean | Columns 1-4 unchanged | Turn groups: accumulated |yaw| >=15 deg within +/-0.5 s', ...
 'Columns 5/6: behavior-selected absolute-yaw peaks; forward-speed drop split at 1 mm/s | Visual and odor epochs retained'},'FontSize',10);
assert(numel(findall(fig,'Type','axes'))==24);exportgraphics(fig,[stem '.png'],'Resolution',170);savefig_visible_on_open(fig,[stem '.fig']);
end
