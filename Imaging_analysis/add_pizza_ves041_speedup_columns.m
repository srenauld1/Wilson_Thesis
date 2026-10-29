function report=add_pizza_ves041_speedup_columns(out)
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
M=readtable(fullfile(out,'manifest.csv'),'Delimiter',',','TextType','string');included=find(M.status=="complete");assert(numel(included)==17);
backup=fullfile(out,'analysis_data','original_4x6_figures');if ~isfolder(backup),mkdir(backup);end
M.speedupsWithoutRotation=nan(height(M),1);M.speedupsWithRotation=nan(height(M),1);rows=cell(0,5);
for k=included'
 trial=M.trial(k);A=load(fullfile(out,'analysis_data',trial+"_events.mat"),'d');S=load(fullfile(out,'analysis_data',trial+"_events_with_rotation_peaks.mat"));d=A.d;rel=S.rel;X=S.X;groups=S.groups;
 [idx,angle,before,after,detector]=pizza_ves041_speedups(d);hasTurn=angle>=15;
 for g=1:2
  groups{g+6}=idx(hasTurn==(g==2));X{g+6}=nan(numel(rel),numel(groups{g+6}),3);vars={d.dff,d.fwd,abs(d.yaw)};
  for j=1:numel(groups{g+6}),for v=1:3,X{g+6}(:,j,v)=interp1(d.t,vars{v},d.t(groups{g+6}(j))+rel,'linear',NaN);end,end
  assert(all(isfinite(X{g+6}(:))));
 end
 assert(isequal(sort([groups{7};groups{8}]),idx)&&isempty(intersect(groups{7},groups{8})));
 for c=1:6,assert(isequal(X{c},S.X{c})&&isequal(groups{c},S.groups{c}));end
 events=table(idx,d.t(idx),before,after,after-before,angle,hasTurn,'VariableNames',{'sampleIndex','eventTime','forwardBefore','forwardAfter','forwardGain','absoluteTurnDeg','withRotation'});
 assert(all(events.forwardGain>=1)&all(events.forwardAfter>=1));
 stem=char(M.eventFigureStem(k));[~,name]=fileparts(stem);for ext={'.png','.fig'},dest=fullfile(backup,[name ext{1}]);if ~isfile(dest),copyfile([stem ext{1}],dest);end,end
 pizza_ves041_render_eight_columns(X,rel,strrep(char(trial),'_',' '), ...
  {'Original dF/F | Thin: events; black: mean | First six groups unchanged', ...
  'Speed-ups: maximum positive acceleration; forward gain >=1 mm/s; after-speed >=1 mm/s | Turn split: slowdown 0:+0.5 s; speed-up -0.5:0 s; accumulated |yaw| >=15 deg'},stem,false);
 writetable(events,fullfile(out,'analysis_data',trial+"_speedup_events.csv"));
 save(fullfile(out,'analysis_data',trial+"_events_with_rotation_peaks_and_speedups.mat"),'X','groups','rel','events','detector','-v7.3');
 M.speedupsWithoutRotation(k)=nnz(~hasTurn);M.speedupsWithRotation(k)=nnz(hasTurn);
 rows(end+1,:)={trial,nnz(~hasTurn),nnz(hasTurn),numel(idx),string(stem)}; %#ok<AGROW>
 fprintf('%s: %d speed-ups without rotation + %d with rotation\n',trial,nnz(~hasTurn),nnz(hasTurn));
end
writetable(M,fullfile(out,'manifest.csv'));report=cell2table(rows,'VariableNames',{'trial','withoutRotation','withRotation','allSpeedups','figureStem'});writetable(report,fullfile(out,'individual_speedup_columns_report.csv'));
% Extend the active combined pooled figure, preserving its exact first six columns.
base=fullfile(out,'ves041_dark_bar_barjump_pooled_equal_fly_peaks_slowdowns');P=load([base '_with_rotation_peaks.mat']);N=load([base '_normalized.mat']);assert(isequal(P.flies,N.flies));
nf=numel(P.flies);counts=[P.counts zeros(nf,2)];means=nan(numel(P.rel),nf,3,8);means(:,:,:,1:6)=P.flyMeans;pooled=cell(nf,2);
for k=1:height(P.F)
 f=find(P.flies==P.F.flyId(k));S=load(fullfile(out,'analysis_data',P.F.trial(k)+"_events_with_rotation_peaks_and_speedups.mat"),'X');
 for g=1:2,pooled{f,g}=cat(2,pooled{f,g},S.X{g+6});counts(f,g+6)=counts(f,g+6)+size(S.X{g+6},2);end
end
for f=1:nf,for g=1:2,if counts(f,g+6)>0,means(:,f,:,g+6)=mean(pooled{f,g},2);end,end,end
normalized=means;for f=1:nf,for v=1:3,normalized(:,f,v,:)=(means(:,f,v,:)-N.center(f,v))/N.scale(f,v);end,end
assert(isequaln(means(:,:,:,1:6),P.flyMeans));q=abs(normalized(:,:,:,1:6)-P.normalizedFlyMeans);assert(all(q(isfinite(q))<1e-10));
rel=P.rel;flies=P.flies;F=P.F;stem=[base '_with_rotation_peaks_and_speedups'];
for mode=1:2
 X=cell(1,8);for c=1:8,if mode==1,X{c}=means(:,:,:,c);else,X{c}=normalized(:,:,:,c);end,end
 dest=stem;if mode==2,dest=[dest '_normalized'];end
 pizza_ves041_render_eight_columns(X,rel,'VES041 dark, bar, and bar jumps | equal fly weight', ...
  {'12 recordings | 8 estimated flies (3 inferred) | Thin: fly means; black: equal-fly mean | First six groups unchanged', ...
   'Columns 7/8: behavior-selected speed-ups, split by accumulated |yaw| <15 or >=15 deg over -0.5:0 s before acceleration peak'},dest,mode==2,flies,sum(counts,1));
end
contributingFlies=sum(counts>0,1);weights=cell(0,4);for c=1:8,for f=1:nf,weights(end+1,:)={flies(f),c,counts(f,c),double(counts(f,c)>0)/contributingFlies(c)};end,end %#ok<AGROW>
W=cell2table(weights,'VariableNames',{'flyId','column','events','flyWeight'});writetable(W,[stem '_fly_weights.csv']);
methods=struct('speedupDetector',detector,'source','Saved first six event groups unchanged; new speed-ups are behavior-selected only', ...
 'averaging','Pool events within fly across included recordings, then average contributing fly means equally', ...
 'normalization','Existing full-recording within-fly mean/sample-SD references retained for every metric', ...
 'split','All retained speed-ups partition exactly into accumulated absolute yaw <15 and >=15 deg over -0.5:0 s before acceleration peak; slowdowns use 0:+0.5 s after deceleration peak; no fluorescence amplitude selection');
save([stem '.mat'],'flies','F','rel','counts','contributingFlies','means','normalized','W','methods','-v7.3');
settings=jsondecode(fileread(fullfile(out,'analysis_settings.json')));settings.turnWindowSec=[0 .5];settings.speedupTurnWindowSec=[-.5 0];settings.individualFigureLayout='4 rows x 8 columns: peaks; slowdowns; slowdowns without/with turns; rotation peaks below/above 1 mm/s forward drop; speed-ups without/with turns';settings.speedupDetector=detector;
fid=fopen(fullfile(out,'analysis_settings.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(settings,PrettyPrint=true));fclose(fid);
fid=fopen(fullfile(out,'README.txt'),'a');assert(fid>=0);fprintf(fid,'\nSpeed-up extension: all 17 included individual no-heat figures are now 4x8 at their original filenames, with columns 7/8 showing speed-ups without/with rotation. Events align to maximum positive acceleration >=2 mm/s^2 after 0.25 s median smoothing; median forward gain >=1 mm/s and after-speed >=1 mm/s are required. Pre/post windows are -0.75:-0.25 and +0.25:+0.75 s, spacing >=2 s, full valid windows -2:+3 s, and speed-up rotation is accumulated absolute yaw >=15 deg over -0.5:0 s before acceleration peak; slowdown rotation uses 0:+0.5 s after deceleration peak. This includes starts from rest and excludes merely becoming less negative while still reversing. The original six event groups and full-recording normalization references are unchanged. Pooled raw/normalized 4x8 companions have suffix _with_rotation_peaks_and_speedups. Previous individual 4x6 figures are under analysis_data/original_4x6_figures.\n');fclose(fid);
disp(report);fprintf('Pooled eight-column counts: ');disp(sum(counts,1));fprintf('Contributing flies: ');disp(contributingFlies);
end
