function export_pizza_example_svgs(scatterPaths,heatOnly)
% Presentation copies only: native behavior + linearly interpolated dF/F.
% Percentiles use acquired fluorescence samples, not interpolated samples.
% No smoothing, extrapolation, amplitude rejection, or time-axis warping.
assert(~feature('ShowFigureWindows'),'Use MATLAB -nodisplay.');
if nargin<1,scatterPaths=true;end
if nargin<2,heatOnly=false;end
root='/Users/sophiarenauld/stacks';
out=fullfile(root,'pizza_talk_analysis','example traces','SVG examples');
if ~isfolder(out),mkdir(out);end
folders={ '20250812-29_ves041_fly3_bbarjumps2 copy', ...
 '20260722-2_smp544_fly1_clbar', ...
 '2026-10-20/20261020_12_1_ves041_syt8m_fly1_stripes', ...
 '2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo'};
ids={'20250812-29','20260722-2','20261020_12_1','20261020_11_2'};
trialIndices=1:4;if heatOnly,trialIndices=3:4;end
for k=trialIndices
 folder=fullfile(root,folders{k});f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));
 a=load(fullfile(f.folder,f.name),'dq','roi');
 if k>=3
  % October source figures used the saved smr_analysis fluorescence, not
  % the current acquisition MAT's subsequently updated ROI trace.
  savedSource=fullfile(folder,'smr_analysis','fly_data.mat');
 end
 % Preserve exactly the fluorescence used by the supplied source figures.
 ti=double(a.dq(1).t(:));di=double(a.roi.ts(1,:))';
 if k>=3
  ti=read_saved(savedSource,'dq/t');di=read_saved(savedSource,'roi/ts');
 end
 assert(numel(ti)==numel(di)&&all(diff(ti)>0));
 bg=struct; nativeSource='saved dq(2)';rotationMatrix=eye(2);
 if k<=2
  b=a.dq(2);t=double(b.t(:));xy=[b.pxb(:),b.pyb(:)];
  fwd=double(b.bvf(:));rot=rad2deg(double(b.bvy(:)));
  % Independently reconstructed streams can differ in initial heading.
  % Match coordinate orientation/origin to the supplied imaging-rate path
  % by rigid rotation/translation ONLY (never stretch distances).
  old=[a.dq(1).pxb(:),a.dq(1).pyb(:)];q=interp1(t,xy,ti);
  good=all(isfinite(q),2)&all(isfinite(old),2);
  mu=mean(q(good,:),1);mv=mean(old(good,:),1);
  [U,~,V]=svd((q(good,:)-mu)'*(old(good,:)-mv));
  rotationMatrix=U*diag([1,det(U*V')])*V';
  xy=(xy-mu)*rotationMatrix+mv;
  rmsPathError=sqrt(mean(sum((interp1(t,xy,ti(good))-old(good,:)).^2,2)));
  fprintf('%s rigid path alignment RMS %.4f mm\n',ids{k},rmsPathError);
  valid=all(isfinite(xy),2)&isfinite(fwd)&isfinite(rot);
  chosen=true(size(t));imageChosen=true(size(ti));
 else
  bf=dir(fullfile(folder,'*_bergii_.mat'));z=load(fullfile(bf.folder,bf.name),'bergii');b=z.bergii;
  assert(isequal(double(b.volumeTimes(:)),ti),'Clock mismatch');
  t=double(b.updateTimes(:));rows=double(b.csvRowIndicesOneBased(:));
  cf=dir(fullfile(folder,'*_fictrac_frames.csv'));nativeSource=fullfile(cf.folder,cf.name);
  opts=detectImportOptions(nativeSource,'Delimiter',',');
  opts.SelectedVariableNames={'kinematic_valid','forward_velocity_rad_s','yaw_velocity_rad_s','protocol_state_json'};
  C=readtable(nativeSource,opts);
  fwd=double(C.forward_velocity_rad_s(rows))*double(b.ballDiameterMm)/2;
  rot=-rad2deg(double(C.yaw_velocity_rad_s(rows)));
  valid=ismember(lower(string(C.kinematic_valid(rows))),["true","1"]);
  xy=nan(numel(t),2);episode=nan(size(t));active=false(size(t));
  fields={'barrier_x_mm','barrier_y_mm','barrier_normal_x','barrier_normal_y','barrier_radius_mm'};
  geom=nan(numel(t),5);
  for j=1:numel(t)
   s=jsondecode(char(string(C.protocol_state_json(rows(j)))));
   if isfield(s,'pos_x_mm')&&~isempty(s.pos_x_mm),xy(j,:)=[s.pos_x_mm,s.pos_y_mm];end
   if k==4
    episode(j)=s.trial_number;active(j)=s.barrier_active;
    for v=1:5,if isfield(s,fields{v})&&~isempty(s.(fields{v})),geom(j,v)=s.(fields{v});end,end
   end
  end
  valid=valid&all(isfinite(xy),2)&isfinite(fwd)&isfinite(rot);
  if k==3
   chosen=t>=314&t<=325;imageChosen=ti>=314&ti<=325;bg.kind='stripe';bg.boundary=-7.5;
  else
   chosen=episode==7&active;assert(any(chosen),'No barrier episode 7');
   imageChosen=ti>=min(t(chosen))&ti<=max(t(chosen));
   bg.kind='barrier';bg.geometry=geom(find(chosen&all(isfinite(geom),2),1),:);
  end
 end
 d=interp_valid(ti,di,t);
 limits=prctile(di(imageChosen&isfinite(di)),[5 95]);
 nativeHz=1/median(diff(t));imagingHz=1/median(diff(ti));
 selection=chosen&valid&isfinite(d);idx=find(selection);
 x=xy(idx,1);y=xy(idx,2);dn=d(idx);tn=t(idx);
 gap=[false;diff(idx)>1|diff(tn)>2.5/nativeHz];
 fig=figure('Visible','off','Color','w','Position',[40 40 1200 850]);
 if k==3,fig.Position=[40 40 1400 480];end
 ax=axes(fig);hold(ax,'on');
 pad=max(.4,.04*max([range(x),range(y)]));
 xl=[min(x)-pad,max(x)+pad];yl=[min(y)-pad,max(y)+pad];
 xlim(ax,xl);ylim(ax,yl);
 if isfield(bg,'kind')
  if strcmp(bg.kind,'stripe')
   patch(ax,[xl fliplr(xl)],[yl(1) yl(1) bg.boundary bg.boundary],[1 .928 .88],'EdgeColor','none');
  else
   g=bg.geometry;normal=g(3:4);normal=normal/norm(normal);tangent=[-normal(2),normal(1)];
   span=3*max([range(xl),range(yl)]);edges=linspace(-g(5),g(5),101);
   for j=1:100
    corners=g(1:2)+[edges(j)*normal-span*tangent;edges(j+1)*normal-span*tangent;edges(j+1)*normal+span*tangent;edges(j)*normal+span*tangent];
    z=1-abs(mean(edges(j:j+1)))/g(5);
    patch(ax,corners(:,1),corners(:,2),[1,1-.24*z,1-.32*z],'EdgeColor','none');
   end
   line(ax,g(1)+[-span span]*tangent(1),g(2)+[-span span]*tangent(2),'Color',[.55 .2 .15],'LineWidth',1.2);
  end
 end
 cmap=parula(256);
 if heatOnly
  cmap=interp1([0 .5 1],[.05 .20 .8;.70 .70 .70;.9 .04 .03],linspace(0,1,256));
 end
 if scatterPaths
  scatter(ax,x,y,10,dn,'filled','MarkerEdgeColor','none');
 else
  draw_vector_path(ax,x,y,dn,limits,cmap,gap);
 end
 plot(ax,x(1),y(1),'o','Color',[0 .6 0],'MarkerFaceColor',[0 .8 0],'MarkerSize',7);
 plot(ax,x(end),y(end),'ks','MarkerFaceColor','k','MarkerSize',6);
 axis(ax,'equal');xlim(ax,xl);ylim(ax,yl);ax.FontSize=18;ax.Box='off';
 xlabel(ax,'X position (mm)','FontSize',22);ylabel(ax,'Y position (mm)','FontSize',22);
 title(ax,strrep(ids{k},'_',' '),'FontSize',24);
 colormap(ax,cmap);caxis(ax,limits);cb=colorbar(ax);cb.Label.String='dF/F';cb.FontSize=18;
 stem=fullfile(out,[ids{k},'_flat_path_dff_native_p05_p95']);
 if scatterPaths,stem=[stem '_scatter'];end
 if heatOnly
  % Final heat companions: SVG only. Keep original acquired-sample limits
  % while drawing the measured 50-Hz coordinates with interpolated dF/F.
  stem=[stem '_blue_gray_red'];
  print(fig,[stem '.svg'],'-dsvg','-painters');close(fig);
  old=fullfile(out,[ids{k} '_flat_path_dff_imaging_rate_scatter10_p05_p95_blue_gray_red.svg']);
  if isfile(old),delete(old);end
  fprintf('%s: connected native path, %d samples, %.1f Hz; limits %.4f %.4f\n',ids{k},numel(tn),nativeHz,limits);
  continue
 end
 print(fig,[stem '.svg'],'-dsvg','-painters');exportgraphics(fig,[stem '.png'],'Resolution',150);
 savefig_visible_on_open(fig,[stem '.fig']);close(fig);
 save([stem '.mat'],'tn','x','y','dn','limits','nativeHz','imagingHz','nativeSource','rotationMatrix','bg');
 fprintf('%s %d native samples, %.3f/%.3f Hz native/imaging; limits %.4f %.4f\n',ids{k},numel(tn),nativeHz,imagingHz,limits);
 if k==1
  % dq(2).bvf/bvy are moving-slope estimates, not raw frame-to-frame
  % velocities. Recompute the shortest possible difference from the saved
  % 60-Hz positions; locate each difference at its actual interval midpoint.
  % No moving fit, Gaussian smoothing, or interpolation of behavior.
  rawTime=(t(1:end-1)+t(2:end))/2;dt=diff(t);
  rawForward=diff(double(b.bf(:)))./dt;
  dh=diff(double(b.bh(:)));rawRotation=rad2deg(atan2(sin(dh),cos(dh))./dt);
  rawValid=valid(1:end-1)&valid(2:end)&dt<2.5/nativeHz;
  rawForward(~rawValid)=NaN;rawRotation(~rawValid)=NaN;
  rawDff=interp_valid(ti,di,rawTime);
  if isfield(b.opt,'dvlensec'),fprintf('Saved legacy velocity moving-slope window: %.3f s\n',b.opt.dvlensec);end
  fig=figure('Visible','off','Color','w','Position',[40 40 1550 850]);tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
  vals={rawDff,rawForward,rawRotation};colors={'r','k','b'};labs={'dF/F','Forward velocity (mm/s)','Rotational velocity (deg/s)'};
  for v=1:3
   ax=nexttile(tl);hold(ax,'on');vv=vals{v};vv(~rawValid)=NaN;
   plot(ax,rawTime,vv,colors{v},'LineWidth',1.4);xlim(ax,[176 226.5]);ylabel(ax,labs{v});grid(ax,'on');ax.FontSize=17;
   if v==3,yline(ax,0,'Color',[.3 .3 .3]);xlabel(ax,'Time (s)');end
  end
  title(tl,'20250812-29 | dF/F, forward and rotational velocity','FontSize',24);
  stem=fullfile(out,'20250812-29_time_traces_176_226_native_frame_differences');
  print(fig,[stem '.svg'],'-dsvg','-painters');exportgraphics(fig,[stem '.png'],'Resolution',150);
  savefig_visible_on_open(fig,[stem '.fig']);close(fig);
  take=rawTime>=176&rawTime<=226.5;nativeTime=rawTime(take);dff=rawDff(take);forward=rawForward(take);rotation=rawRotation(take);
  save([stem '.mat'],'nativeTime','dff','forward','rotation');
 end
end
end

function d=interp_valid(ti,di,t)
d=nan(size(t));valid=isfinite(di);cut=[true;diff(ti)>2.5*median(diff(ti))];
starts=find(valid&([true;~valid(1:end-1)]|cut));ends=find(valid&([~valid(2:end);true]|[cut(2:end);true]));
for j=1:numel(starts)
 ix=starts(j):ends(j);if numel(ix)<2,continue,end
 take=t>=ti(ix(1))&t<=ti(ix(end));d(take)=interp1(ti(ix),di(ix),t(take),'linear');
end
end

function draw_vector_path(ax,x,y,d,limits,cmap,gap)
% Quantized segment colors produce genuine vector SVG paths (no image).
c=(d(1:end-1)+d(2:end))/2;bins=max(1,min(256,1+floor(255*(c-limits(1))/diff(limits))));
ok=~gap(2:end)&isfinite(c);
for b=1:256
 ix=find(ok&bins==b);if isempty(ix),continue,end
 xx=[x(ix),x(ix+1),nan(size(ix))]';yy=[y(ix),y(ix+1),nan(size(ix))]';
 plot(ax,xx(:),yy(:),'Color',cmap(b,:),'LineWidth',2.8);
end
end

function data=read_saved(file,field)
% Read only a numeric saved-analysis field, avoiding the image stack.
f=H5F.open(file,'H5F_ACC_RDONLY','H5P_DEFAULT');cf=onCleanup(@()H5F.close(f)); %#ok<NASGU>
d=H5D.open(f,['/fly_data/' field]);cd=onCleanup(@()H5D.close(d)); %#ok<NASGU>
type=H5D.get_type(d);ct=onCleanup(@()H5T.close(type)); %#ok<NASGU>
if H5T.get_class(type)==H5ML.get_constant_value('H5T_REFERENCE')
 ref=H5D.read(d);target=H5R.dereference(d,'H5R_OBJECT',ref(:,1));
 closeTarget=onCleanup(@()H5D.close(target)); %#ok<NASGU>
 data=H5D.read(target);
else,data=H5D.read(d);end
data=double(data(:));
end
