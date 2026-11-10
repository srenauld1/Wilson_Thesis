function report=export_all_heat_barrier_dff_pngs(manifestFile,out)
% Export PNGs only; do not rerun event detection or alter source data.
% Modern recordings: one PNG per logged barrier/control episode plus pages
% of six episodes. Older recordings: full fictive path without invented heat
% boundaries. DFF color saturation is 5th/95th percentile within recording;
% no fluorescence amplitude exclusion. All plots have equal spatial axes.
% Hand-drawn ROI1 is preferred. Registration-fixed duplicates are resolved
% in the supplied source manifest. Invalid samples/time gaps break lines.
if nargin<1,manifestFile='/tmp/heat_barrier_sources.json';end
if nargin<2,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/heat_barrier_flat_paths_png_20261105';end
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay');
sources=jsondecode(fileread(manifestFile));if ~isfolder(out),mkdir(out);end
rows=cell(0,8);
cm=interp1([0 .5 1],[.05 .20 .8;.70 .70 .70;.9 .04 .03],linspace(0,1,256));
for k=1:numel(sources)
 s=sources(k);id=s.identity;fprintf('FLAT PATH %s [%s]\n',id,s.kind);
 try
  if strcmp(s.kind,'legacy-barrier'),loaded=load(s.file,'dq','roi');else,loaded=load(s.file,'fly_data');end
  if isfield(loaded,'fly_data'),a=loaded.fly_data;else,a=loaded;end
  if strcmp(s.kind,'legacy-barrier')&&(isempty(a.dq)||isempty(a.roi))
   if contains(lower(s.name),'smp544'),celltype='SMP544';else,celltype='VES041';end
   dest=fullfile(out,celltype,id);if ~isfolder(dest),mkdir(dest);end
   note=recover_saved_path(fileparts(s.file),dest,id,cm);
   rows(end+1,:)={id,celltype,s.kind,s.file,1,0,'recovered',note}; %#ok<AGROW>
   fprintf('  RECOVERED: %s\n',note);continue
  end
  a=select_ves041_manual_roi(a);
  t=double(a.dq(1).t(:));d=double(a.roi.ts(1,:))';
  assert(numel(t)==numel(d),'DFF/time length mismatch: %d vs %d',numel(d),numel(t));
  assert(all(diff(t)>0),'Nonmonotonic imaging clock');
  finite=isfinite(d);assert(any(finite),'No finite fluorescence');cl=prctile(d(finite),[5 95]);
  if cl(1)==cl(2),cl=cl+[-.1 .1];end
  if contains(lower(s.name),'smp544'),celltype='SMP544';else,celltype='VES041';end
  dest=fullfile(out,celltype,id);if ~isfolder(dest),mkdir(dest);end
  legacy=strcmp(s.kind,'legacy-barrier');warningText='';
  if legacy
   assert(isfield(a.dq(1),'pxb')&&isfield(a.dq(1),'pyb'),'Missing legacy ball-heading flat path');
   x=double(a.dq(1).pxb(:));y=double(a.dq(1).pyb(:));p=[];
   sets={find(isfinite(x)&isfinite(y)&finite)};episodes=NaN;
   warningText='Legacy heat recording: no aligned barrier geometry; full path only, no boundary reconstructed.';
  else
   p=a.protocol;assert(p.geometryAligned,'No aligned barrier geometry');
   assert(isequal(double(p.t(:)),t),'Protocol/image clock mismatch');
   x=double(p.pos_x_mm(:));y=double(p.pos_y_mm(:));
   valid=isfinite(x)&isfinite(y)&finite;
   if isfield(p,'coordinate_step_valid'),v=p.coordinate_step_valid(:);valid=valid&isfinite(v)&v==1;end
   if isfield(a.dq(1),'bergii')&&isfield(a.dq(1).bergii,'behaviorVolumeValid')
    v=a.dq(1).bergii.behaviorVolumeValid(:);valid=valid&isfinite(v)&v==1;
   end
   episodes=unique(p.trial_number(p.barrier_active>0&isfinite(p.trial_number)&p.trial_number>0))';
   assert(~isempty(episodes),'No logged barrier episodes');
   sets=arrayfun(@(e)find(valid&p.barrier_active>0&p.trial_number==e),episodes,'UniformOutput',false);
  end
  assert(numel(x)==numel(t)&&numel(y)==numel(t),'Coordinate/time length mismatch');
  exported=0;empty=0;
  for j=1:numel(sets)
   ix=sets{j};if numel(ix)<2,empty=empty+1;continue,end
   fig=figure('Visible','off','Color','w','Position',[30 30 1100 850]);ax=axes(fig);
   drawpath(ax,p,t,x,y,d,ix,cl,cm,episodes(j));
   title(ax,{strrep([celltype ' | ' id],'_',' '),context(p,ix,episodes(j)),sprintf('%.2f–%.2f s | ROI1: hand drawn',t(ix(1)),t(ix(end)))},'FontSize',13);
   if legacy,stem='legacy_full_path';else,stem=sprintf('episode_%03d_%s',episodes(j),state(p,ix));end
   exportgraphics(fig,fullfile(dest,[stem '.png']),'Resolution',160);close(fig);exported=exported+1;
  end
  if ~legacy
   for page=1:ceil(numel(sets)/6)
    fig=figure('Visible','off','Color','w','Position',[30 30 1650 1050]);tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
    for j=(page-1)*6+1:min(page*6,numel(sets))
     ax=nexttile(tl);ix=sets{j};
     if numel(ix)<2,title(ax,sprintf('Episode %g: insufficient valid data',episodes(j)));axis(ax,'off');continue,end
     drawpath(ax,p,t,x,y,d,ix,cl,cm,episodes(j));title(ax,context(p,ix,episodes(j)),'FontSize',10);
    end
    title(tl,{strrep([celltype ' | ' id],'_',' '),'Hand-drawn ROI1 ΔF/F | Same color limits across episodes | Heat shading = command, not temperature'},'FontSize',15);
    exportgraphics(fig,fullfile(dest,sprintf('overview_%02d.png',page)),'Resolution',140);close(fig);
   end
  end
  if empty>0,warningText=[warningText sprintf(' %d episodes had fewer than 2 valid points.',empty)];end
  rows(end+1,:)={id,celltype,s.kind,s.file,exported,empty,'complete',warningText}; %#ok<AGROW>
  fprintf('  %d paths; %d insufficient episodes; %s\n',exported,empty,warningText);
 catch err
  close(findall(groot,'Type','figure','Visible','off'));
  rows(end+1,:)={id,'unknown',s.kind,s.file,0,0,'ERROR',err.message}; %#ok<AGROW>
  fprintf(2,'  ERROR %s: %s\n',id,err.message);
 end
end
report=cell2table(rows,'VariableNames',{'recording','cellType','kind','source','pathPNGs','emptyEpisodes','status','notes'});
disp(report(:,{'recording','pathPNGs','status','notes'}));
end

function note=recover_saved_path(folder,out,id,cm)
% Recovery is explicit: do not pretend the original source structures exist.
file=fullfile(folder,'smr_analysis','heat_and_dff_colored_vs_flatpath.fig');
assert(isfile(file),'Source ROI/behavior empty and no saved flat-path figure exists');
old=openfig(file,'invisible');closer=onCleanup(@()close(old));
points=findall(old,'Type','scatter');assert(numel(points)==1,'Ambiguous saved fluorescence path');
x=double(points.XData(:));y=double(points.YData(:));d=double(points.CData(:));
assert(numel(d)==numel(x)&&numel(y)==numel(x),'Saved figure lacks scalar per-point fluorescence');
valid=isfinite(x)&isfinite(y)&isfinite(d);assert(nnz(valid)>1);
fig=figure('Visible','off','Color','w','Position',[30 30 1100 850]);ax=axes(fig);hold(ax,'on');
scatter(ax,x(valid),y(valid),10,d(valid),'filled');axis(ax,'equal');colormap(ax,cm);clim(ax,prctile(d(valid),[5 95]));
plot(ax,x(find(valid,1)),y(find(valid,1)),'gs','MarkerFaceColor','g');
plot(ax,x(find(valid,1,'last')),y(find(valid,1,'last')),'s','Color',[0 .8 .9],'MarkerFaceColor',[0 .8 .9]);
cb=colorbar(ax);cb.Label.String='Saved right-ROI ΔF/F (5th–95th percentile colors)';xlabel(ax,'X (mm)');ylabel(ax,'Y (mm)');
title(ax,{strrep(['SMP544 | ' id],'_',' '),'Recovered saved right-ROI flat path','Source MAT ROI/behavior are empty; no new processing or heat geometry'},'FontSize',12);
exportgraphics(fig,fullfile(out,'legacy_saved_right_roi_path_RECOVERED.png'),'Resolution',160);close(fig);
note='Source MAT has empty ROI and behavior; recovered right-ROI scalar fluorescence and XY from existing heat_and_dff_colored_vs_flatpath.fig. Not reprocessed; no heat boundary.';
end

function s=state(p,ix)
if isempty(p),s='legacy';return,end
v=unique(p.control_trial(ix));
if isequal(v,1),s='control';elseif isequal(v,0),s='heat';else,s='mixed';end
end
function label=context(p,ix,e)
if isempty(p),label='Legacy heat recording | geometry unavailable';return,end
label=sprintf('Episode %g | %s',e,state(p,ix));
if isfield(p,'heat_geometry')
 v=unique(string(p.heat_geometry(ix)));label=[label ' | ' strjoin(cellstr(v),', ')];
end
end
function drawpath(ax,p,t,x,y,d,ix,cl,cm,e) %#ok<INUSD>
hold(ax,'on');xx=x(ix);yy=y(ix);dd=d(ix);
bx=[min(xx),max(xx)];by=[min(yy),max(yy)];pad=max(.5,.07*max([diff(bx),diff(by)]));bx=bx+[-pad,pad];by=by+[-pad,pad];
if ~isempty(p)
 q=ix(1);keys={'barrier_x_mm','barrier_y_mm','barrier_normal_x','barrier_normal_y','barrier_radius_mm'};
 if all(isfield(p,keys))
  nx=p.barrier_normal_x(q);ny=p.barrier_normal_y(q);cx=p.barrier_x_mm(q);cy=p.barrier_y_mm(q);
  assert(all(isfinite([nx ny cx cy]))&&hypot(nx,ny)>0,'Invalid barrier geometry');
  [X,Y]=meshgrid(linspace(bx(1),bx(2),180),linspace(by(1),by(2),180));distance=(cx-X)*nx+(cy-Y)*ny;
  radius=p.barrier_radius_mm(q);Z=max(0,1-abs(distance)/radius);
  if isfield(p,'heat_profile')&&string(p.heat_profile(q))=="barrier"&&isfield(p,'ramp_distance_mm')
   Z=max(0,1-distance/p.ramp_distance_mm(q));Z(distance<0)=0;
  end
  if strcmp(state(p,ix),'control'),rgb=repmat(1-.16*Z,1,1,3);else,rgb=cat(3,ones(size(Z)),1-.20*Z,1-.27*Z);end
  image(ax,bx,by,rgb);
  if min(distance(:))<=0&&max(distance(:))>=0,contour(ax,X,Y,distance,[0 0],'Color',[.5 .25 .2],'LineWidth',1);end
 end
end
% Do not connect across invalid coordinates, dropped samples, or clock gaps.
gap=[false;diff(ix)>1|diff(t(ix))>2.5*median(diff(t))];xx(gap)=NaN;yy(gap)=NaN;dd(gap)=NaN;
surface(ax,[xx xx],[yy yy],zeros(numel(ix),2),[dd dd],'FaceColor','none','EdgeColor','interp','LineWidth',2.3);
plot(ax,x(ix(1)),y(ix(1)),'s','MarkerFaceColor',[0 .75 .1],'MarkerEdgeColor','none','MarkerSize',7);
plot(ax,x(ix(end)),y(ix(end)),'s','MarkerFaceColor',[0 .8 .9],'MarkerEdgeColor','none','MarkerSize',7);
axis(ax,'equal');xlim(ax,bx);ylim(ax,by);set(ax,'YDir','normal','FontSize',11,'Box','off');
colormap(ax,cm);clim(ax,cl);cb=colorbar(ax);cb.Label.String='ΔF/F (5th–95th percentile colors)';
xlabel(ax,'X (mm)');ylabel(ax,'Y (mm)');
end
