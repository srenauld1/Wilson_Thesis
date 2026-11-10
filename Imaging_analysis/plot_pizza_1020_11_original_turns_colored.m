function result=plot_pizza_1020_11_original_turns_colored()
% Exact original-folder selections and traces, redrawn in presentation colors.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
trial='/Users/sophiarenauld/stacks/2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo';
source=fullfile(trial,'hand_drawn_15s_ave_f_20261026','heat_exploration');
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/20261020_11_2/presentation_control_vs_heat_20261102/original_criteria_turns';if ~isfolder(out),mkdir(out);end
E=readtable(fullfile(source,'slowdown_events.csv'),'Delimiter',',','TextType','string');
dataFile=fullfile(fileparts(source),'fly_data.mat');A=load(dataFile,'fly_data');a=A.fly_data;t=a.dq(1).t(:);d=a.roi.ts(1,:)';fraw=a.dq(1).bvf(:);h=a.dq(1).heat(:);p=a.protocol;
if isfield(a.dq(1),'bvy_deg'),rraw=a.dq(1).bvy_deg(:);else,rraw=rad2deg(a.dq(1).bvy(:));end
n=max(3,round(.25/median(diff(t))));f=smoothdata(fraw,'movmedian',n);r=smoothdata(rraw,'movmedian',n);acc=gradient(f,t);fs=1/median(diff(t));
valid=isfinite(f)&isfinite(r)&isfinite(h)&isfinite(d)&isfinite(p.pos_x_mm(:));if isfield(a.dq(1),'bergii')&&isfield(a.dq(1).bergii,'behaviorVolumeValid'),valid=valid&logical(a.dq(1).bergii.behaviorVolumeValid(:));end
candidate=valid&acc<=-2;beg=find(diff([false;candidate])==1);last=find(diff([candidate;false])==-1);times=[];contexts=strings(0,1);turns=[];lastAccepted=-Inf;
for k=1:numel(beg)
 [~,j]=min(acc(beg(k):last(k)));q=beg(k)+j-1;et=t(q);if et-lastAccepted<2||et<t(1)+2||et>t(end)-3,continue,end
 win=t>=et-2&t<=et+3;pre=t>=et-.75&t<et-.25;post=t>=et+.25&t<=et+.75;core=t>=et-.5&t<=et+.5;
 if ~all(valid(win))||any(diff(t(win))>2/fs),continue,end
 before=median(f(pre));drop=before-median(f(post));if before<1||drop<1,continue,end
 angle=trapz(t(core),abs(r(core)));on=h(core)>.05;
 if all(on),ctx="heat";elseif all(~on),ctx="no heat";else,ctx="transition";end
 if any(~isfinite(p.control_trial(core)))||numel(unique(p.trial_number(core)))~=1,ctx="transition";
 elseif all(p.control_trial(core)==1),ctx="control";end
 times(end+1,1)=et;contexts(end+1,1)=ctx;turns(end+1,1)=angle>=30;lastAccepted=et;
end
assert(numel(times)==height(E)&&max(abs(times-E.eventTime))<1e-6&&isequal(contexts,E.context)&&isequal(logical(turns),logical(E.reorientation)),'Recomputed events differ from original table');
for j=1:height(E),base=t>=E.eventTime(j)-2&t<E.eventTime(j)-1;assert(abs(mean(d(base))-E.baselineDff(j))<1e-5,'Fluorescence source differs from original event table');end
selected=E(E.reorientation==1 & ismember(E.context,["control","heat"]),:);conditions=["control","heat"];
X=cell(2,1);mu=cell(2,1);sem=cell(2,1);counts=zeros(1,2);rel=(-2:.05:3)';vars={d,fraw,abs(rraw)};
for s=1:2
 events=selected.eventTime(selected.context==conditions(s));counts(s)=numel(events);X{s}=nan(numel(rel),counts(s),3);
 for j=1:counts(s),for v=1:3,X{s}(:,j,v)=interp1(t,vars{v},events(j)+rel,'linear',NaN);end,end
 assert(all(isfinite(X{s}(:))));mu{s}=squeeze(mean(X{s},2));sem{s}=squeeze(std(X{s},0,2))/sqrt(counts(s));
end
assert(isequal(counts,[15 18]));
fig=figure('Visible','off','Color','w','Position',[20 20 1120 1050]);tl=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','loose');
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};thinColors={[.88 .6 .6],[.65 .65 .65],[.57 .65 .9]};labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};names={'Control turns','Heat turns'};
for row=1:3
 z=[reshape(X{1}(:,:,row),[],1);reshape(X{2}(:,:,row),[],1)];for s=1:2,z=[z;mu{s}(:,row)-sem{s}(:,row);mu{s}(:,row)+sem{s}(:,row)];end
 pad=max(.01,.06*range(z));lim=[min(z)-pad max(z)+pad];
 for s=1:2
  ax=nexttile(tl);hold(ax,'on');plot(ax,rel,X{s}(:,:,row),'Color',thinColors{row},'LineWidth',.65,'Tag','eventTrace');
  m=mu{s}(:,row);se=sem{s}(:,row);fill(ax,[rel;flipud(rel)],[m-se;flipud(m+se)],colors{row},'FaceAlpha',.2,'EdgeColor','none','Tag','eventSem');plot(ax,rel,m,'Color',colors{row},'LineWidth',3,'Tag','eventMean');
  xline(ax,0,':','Color',[.45 .45 .45],'LineWidth',1);xlim(ax,[-2 3]);ylim(ax,lim);ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.12;ax.XTick=-2:1:3;
  if s==1,ylabel(ax,labels{row},'FontSize',21);end
  if row==3,xlabel(ax,'Time from strongest deceleration (s)','FontSize',17);else,ax.XTickLabel=[];end
  if row==1,title(ax,names{s},'FontSize',23,'FontWeight','bold');subtitle(ax,sprintf('n=%d events',counts(s)),'FontSize',15);end
 end
end
subtitle(tl,{'2026-10-20 trial 11-2 | Original turn criteria','Thin: events | Thick: mean | Shading: SEM across events'},'FontSize',15);
assert(numel(findall(fig,'Type','axes'))==6&&numel(findall(fig,'Tag','eventTrace'))==3*sum(counts));
stem=fullfile(out,'20261020_11_2_original_criteria_control_heat_turns');exportgraphics(fig,[stem '.png'],'Resolution',200);exportgraphics(fig,[stem '.pdf'],'ContentType','vector');savefig_visible_on_open(fig,[stem '.fig']);close(fig);
writetable(selected,fullfile(out,'selected_original_turn_events.csv'));
methods=jsondecode(fileread(fullfile(trial,'smr_analysis','heat_exploration','analysis_settings.json')));methods.display='Original hand-drawn ROI traces; no added display smoothing or normalization; old selection rerun and every timestamp/context/turn flag verified against original event CSV; event SEM added';
methods.sourceFolder=source;methods.traceSource=dataFile;methods.originalTableVerified=true;methods.conditions='Original context==control and context==heat, both with reorientation==1; no-heat and transition contexts excluded';
methods.controlDefinition='All control_trial==1 within -0.5:+0.5 s, finite flags and one barrier trial in that core; no extra full-window same-control or barrier-active restriction';
methods.heatDefinition='All measured heat command samples >0.05 V within -0.5:+0.5 s; core must stay in one trial and not be classified control';
methods.uncertainty='Sample SD/sqrt(N events) within condition, descriptive single-recording single-fly SEM';
save(fullfile(out,'comparison_data.mat'),'X','mu','sem','counts','rel','selected','methods','-v7.3');fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('output',out,'controlTurns',counts(1),'heatTurns',counts(2),'originalEventsRecomputedAndVerified',true);disp(result);
end
