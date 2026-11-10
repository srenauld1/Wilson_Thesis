function export_pizza_barrier_turn_examples()
% Illustrative, response-selected examples; not an independent statistical test.
assert(~feature('ShowFigureWindows'),'Use -nodisplay');
root='/Users/sophiarenauld/stacks/pizza_talk_analysis';
out=fullfile(root,'example traces','SVG examples','heat_control_turn_examples_20261104');if ~isfolder(out),mkdir(out);end
base=fullfile(root,'ves041_heat_analysis','barrier');
S=load(fullfile(base,'pooled_four_requested_trials_oldcriteria_20261102','crunch_data.mat'));
A=load(fullfile(base,'20261020_11_2','analysis.mat'));d=A.d;
folder='/Users/sophiarenauld/stacks/2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo';
bf=dir(fullfile(folder,'*_bergii_.mat'));z=load(fullfile(bf.folder,bf.name),'bergii');b=z.bergii;
assert(isequal(double(b.volumeTimes(:)),d.t(:)));
opts=detectImportOptions(b.sourceCsv,'Delimiter',',');opts.SelectedVariableNames={'protocol_state_json','kinematic_valid','forward_velocity_rad_s','yaw_velocity_rad_s'};C=readtable(b.sourceCsv,opts);
T=S.allEvents(S.allEvents.recording=="20261020_11_2"&S.allEvents.selected,:);
seeds=[545.97 58.359];names={'Control','Heat'};P=struct;selected=table;
for k=1:2
 [err,j]=min(abs(T.eventTime-seeds(k)));assert(err<.01);E=T(j,:);et=E.eventTime;
 assert(E.context==lower(string(names{k})));
 take=b.updateTimes>=et-2&b.updateTimes<=et+3;ix=find(take);rows=b.csvRowIndicesOneBased(ix);t=double(b.updateTimes(ix));
 xy=nan(numel(t),2);distance=nan(size(t));trial=distance;active=false(size(t));control=active;
 geom=nan(numel(t),5);
 for n=1:numel(t)
  v=jsondecode(char(string(C.protocol_state_json(rows(n)))));
  xy(n,:)=[v.pos_x_mm v.pos_y_mm];trial(n)=v.trial_number;active(n)=v.barrier_active;control(n)=v.control_trial;
  if ~isempty(v.distance_to_barrier_mm),distance(n)=v.distance_to_barrier_mm;end
  fields={'barrier_x_mm','barrier_y_mm','barrier_normal_x','barrier_normal_y','barrier_radius_mm'};
  for h=1:5,if ~isempty(v.(fields{h})),geom(n,h)=v.(fields{h});end,end
 end
 assert(all(active)&&numel(unique(trial))==1&&all(control==(k==1)),'Window crosses protocol contexts');
 ok=ismember(lower(string(C.kinematic_valid(rows))),["true","1"]);assert(all(ok)&&all(isfinite(xy),'all')&&all(diff(t)<.06));
 di=interp1(d.t,d.dff,t,'linear',nan);heat=interp1(d.t,d.heat,t,'linear',nan);assert(all(isfinite(di)));
 core=abs(t-et)<=.5;if k==2,assert(all(heat(core)>.05));end
 P(k).time=t;P(k).xy=xy;P(k).dff=di;P(k).heat=heat;P(k).rows=rows;P(k).geometry=geom(ceil(end/2),:);P(k).eventTime=et;P(k).trial=trial(1);P(k).distance=distance;
 P(k).fwd=double(C.forward_velocity_rad_s(rows))*double(b.ballDiameterMm)/2;P(k).yaw=-rad2deg(double(C.yaw_velocity_rad_s(rows)));
 q=d.t-et;pre=q>=-2&q< -1;post=q>=0&q<=1.5;
 E.baselineDff=mean(d.dff(pre));E.postDff=mean(d.dff(post));E.responseDff=E.postDff-E.baselineDff;E.protocolTrial=trial(1);E.windowStart=et-2;E.windowEnd=et+3;selected=[selected;E];
end
% Shared acquired-sample percentiles over the union of the two displayed windows.
mask=false(size(d.t));for k=1:2,mask=mask|(d.t>=P(k).eventTime-2&d.t<=P(k).eventTime+3);end
limits=prctile(d.dff(mask&isfinite(d.dff)&d.valid),[5 95]);
cmap=interp1([0 .5 1],[.05 .20 .8;.70 .70 .70;.9 .04 .03],linspace(0,1,256));
span=0;for k=1:2,span=max(span,max(range(P(k).xy)));end;span=max(12,span*1.3);barLength=5;
for k=1:2
 f=figure('Visible','off','Color','w','Position',[20 20 760 650]);ax=axes(f,'Position',[.08 .12 .72 .76]);drawpath(ax,P(k),names{k},span,barLength,cmap,limits);
 cb=colorbar(ax);cb.Label.String='ΔF/F';cb.FontSize=17;caxis(ax,limits);
 stem=fullfile(out,sprintf('20261020_11_2_%s_trial%d_turn_%.3fs_flat_path',lower(names{k}),P(k).trial,P(k).eventTime));saveall(f,stem);close(f);
end
f=figure('Visible','off','Color','w','Position',[20 20 1250 650]);tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','loose');
for k=1:2,ax=nexttile(tl);drawpath(ax,P(k),names{k},span,barLength,cmap,limits);end
cb=colorbar(ax);cb.Layout.Tile='east';cb.Label.String='ΔF/F';cb.FontSize=17;
subtitle(tl,'Same recording | Shared ΔF/F and spatial scales | Open circle: start; black ring: strongest deceleration','FontSize',13);
saveall(f,fullfile(out,'control_vs_heat_turn_examples_flat_paths'));close(f);
% Companion traces document what each path illustrates, using native kinematics.
f=figure('Visible','off','Color','w','Position',[20 20 1050 850]);tl=tiledlayout(f,3,2,'TileSpacing','compact');
colors={[.85 .05 .05],[0 0 0],[.05 .15 .9]};labels={'ΔF/F','Forward velocity (mm/s)','Rotational speed (°/s)'};
for r=1:3
 for k=1:2
  ax=nexttile(tl);hold(ax,'on');if r==1,y=P(k).dff;elseif r==2,y=P(k).fwd;else,y=abs(P(k).yaw);end
  plot(ax,P(k).time-P(k).eventTime,y,'Color',colors{r},'LineWidth',1.8);xline(ax,0,':');xlim(ax,[-2 3]);ax.FontSize=14;ax.Box='off';if k==1,ylabel(ax,labels{r});end;if r==1,title(ax,sprintf('%s | %.3f s',names{k},P(k).eventTime));end;if r==3,xlabel(ax,'Time from strongest deceleration (s)');end
 end
end
saveall(f,fullfile(out,'control_vs_heat_selected_event_traces'));close(f);
writetable(selected,fullfile(out,'selected_examples.csv'));
methods=struct('source',A.meta.source,'nativeBehavior',b.sourceCsv,'selection','Illustrative response-selected pair from previously qualified >=30-degree slowdown turns; similar speed drop; both windows entirely within one active protocol barrier period; April27 fly3 excluded', ...
 'response','Mean dF/F 0:+1.5 s minus mean -2:-1 s; acquired imaging samples', ...
 'path','Native 50-Hz protocol X/Y in mm; raw fluorescence linearly interpolated to coded-sync native times; original orientation preserved; positions centered for display only', ...
 'color','Shared raw dF/F limits: pooled displayed-window acquired-sample 5th/95th percentiles; same blue-gray-red palette as previous heat flat paths; saturation at limits, no point removal', ...
 'scaleBarMm',barLength,'windowSeconds',[-2 3],'interpretation','Selected examples illustrate the requested contrast, not its frequency or significance');
save(fullfile(out,'example_data.mat'),'P','selected','limits','cmap','methods');fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);
copyfile([mfilename('fullpath') '.m'],out);disp(selected(:,{'eventTime','context','protocolTrial','speedDrop','absoluteTurnDeg','baselineDff','postDff','responseDff'}));disp(out);
end
function drawpath(ax,P,name,span,barLength,cmap,limits)
hold(ax,'on');xy=P.xy;origin=mean([min(xy);max(xy)],1);xy=xy-origin;lim=[-span/2 span/2];xlim(ax,lim);ylim(ax,lim);
% Show virtual barrier band, including the unheated control counterpart.
g=P.geometry;normal=g(3:4)/norm(g(3:4));tangent=[-normal(2) normal(1)];center=g(1:2)-origin;width=g(5);edges=linspace(-width,width,101);
for j=1:100
 co=center+[edges(j)*normal-span*3*tangent;edges(j+1)*normal-span*3*tangent;edges(j+1)*normal+span*3*tangent;edges(j)*normal+span*3*tangent];
 if strcmp(name,'Heat'),z=1-abs(mean(edges(j:j+1)))/width;col=[1 1-.13*z 1-.17*z];else,col=[.97 .97 .97];end
 patch(ax,co(:,1),co(:,2),col,'EdgeColor','none');
end
line(ax,center(1)+[-3*span 3*span]*tangent(1),center(2)+[-3*span 3*span]*tangent(2),'Color',[.65 .65 .65],'LineStyle',':');
for j=1:size(xy,1)-1
 c=mean(P.dff(j:j+1));u=min(1,max(0,(c-limits(1))/diff(limits)));rgb=interp1(linspace(0,1,256),cmap,u);
 line(ax,xy(j:j+1,1),xy(j:j+1,2),'Color',rgb,'LineWidth',3);
end
plot(ax,xy(1,1),xy(1,2),'ko','MarkerSize',7,'LineWidth',1.1);
[~,e]=min(abs(P.time-P.eventTime));plot(ax,xy(e,1),xy(e,2),'ko','MarkerSize',10,'LineWidth',1.5);
plot(ax,xy(end,1),xy(end,2),'ks','MarkerSize',5,'MarkerFaceColor','k');
x=lim(1)+.08*span;y=lim(1)+.11*span;plot(ax,[x x+barLength],[y y],'k-','LineWidth',3);text(ax,x+barLength/2,y-.035*span,sprintf('%g mm',barLength),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',15);
axis(ax,'equal');xlim(ax,lim);ylim(ax,lim);axis(ax,'off');colormap(ax,cmap);caxis(ax,limits);title(ax,sprintf('%s barrier | trial %d\n%.2f–%.2f s',name,P.trial,P.eventTime-2,P.eventTime+3),'FontSize',19,'FontWeight','normal');
end
function saveall(f,stem)
exportgraphics(f,[stem '.png'],'Resolution',180);exportgraphics(f,[stem '.pdf'],'ContentType','vector');print(f,[stem '.svg'],'-dsvg','-painters');savefig(f,[stem '.fig']);
end
