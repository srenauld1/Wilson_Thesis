function export_pizza_final_imaging_scatter_svgs()
% SVG-only presentation export. All coordinates and fluorescence sampled on
% the acquired imaging clock, no upsampling/smoothing. Saturate colors at
% displayed-epoch acquired-sample percentiles; keep all finite samples.
assert(~feature('ShowFigureWindows'),'Run headlessly.');
root='/Users/sophiarenauld/stacks';
out=fullfile(root,'pizza_talk_analysis','example traces','SVG examples');
folders={'20250812-29_ves041_fly3_bbarjumps2 copy','20260722-2_smp544_fly1_clbar', ...
 '2026-10-20/20261020_12_1_ves041_syt8m_fly1_stripes', ...
 '2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo'};
ids={'20250812-29','20260722-2','20261020_12_1','20261020_11_2'};
cm=interp1([0 .5 1],[.05 .20 .8;.70 .70 .70;.9 .04 .03],linspace(0,1,256));
keep=cell(1,4);
for k=1:4
 folder=fullfile(root,folders{k});g=[];
 if k<=2
  f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));a=load(fullfile(f.folder,f.name),'dq','roi');
  t=double(a.dq(1).t(:));x=double(a.dq(1).pxb(:));y=double(a.dq(1).pyb(:));d=double(a.roi.ts(1,:))';
  chosen=true(size(t));
 else
  source=fullfile(folder,'smr_analysis','fly_data.mat');
  t=read_saved(source,'dq/t');d=read_saved(source,'roi/ts');
  x=read_saved(source,'protocol/pos_x_mm');y=read_saved(source,'protocol/pos_y_mm');
  if k==3
   chosen=t>=314&t<=325;
  else
   episode=read_saved(source,'protocol/trial_number');active=read_saved(source,'protocol/barrier_active');
   chosen=episode==7&active>0;
   fields={'barrier_x_mm','barrier_y_mm','barrier_normal_x','barrier_normal_y','barrier_radius_mm'};
   q=find(chosen,1);g=zeros(1,5);
   for j=1:5,v=read_saved(source,['protocol/' fields{j}]);g(j)=v(q);end
  end
 end
 assert(isequal(size(t),size(x),size(y),size(d)));
 limits=prctile(d(chosen&isfinite(d)),[5 95]);v=chosen&isfinite(x)&isfinite(y)&isfinite(d);
 x=x(v);y=y(v);d=d(v);assert(~isempty(x));
 fig=figure('Visible','off','Color','w','Position',[40 40 1200 850]);
 if k==3,fig.Position=[40 40 1400 480];end
 ax=axes(fig);hold(ax,'on');pad=max(.4,.04*max([range(x),range(y)]));
 xl=[min(x)-pad,max(x)+pad];yl=[min(y)-pad,max(y)+pad];
 xlim(ax,xl);ylim(ax,yl);
 if k==3
  patch(ax,[xl fliplr(xl)],[yl(1) yl(1) -7.5 -7.5],[1 .928 .88],'EdgeColor','none');
 elseif k==4
  normal=g(3:4)/norm(g(3:4));tangent=[-normal(2),normal(1)];
  span=3*max([range(xl),range(yl)]);edges=linspace(-g(5),g(5),101);
  for j=1:100
   corners=g(1:2)+[edges(j)*normal-span*tangent;edges(j+1)*normal-span*tangent;edges(j+1)*normal+span*tangent;edges(j)*normal+span*tangent];
   z=1-abs(mean(edges(j:j+1)))/g(5);
   patch(ax,corners(:,1),corners(:,2),[1,1-.24*z,1-.32*z],'EdgeColor','none');
  end
  line(ax,g(1)+[-span span]*tangent(1),g(2)+[-span span]*tangent(2),'Color',[.55 .2 .15],'LineWidth',1.2);
 end
 scatter(ax,x,y,10,d,'filled','MarkerEdgeColor','none');
 axis(ax,'equal');xlim(ax,xl);ylim(ax,yl);ax.FontSize=18;ax.Box='off';
 xlabel(ax,'X position (mm)','FontSize',22);ylabel(ax,'Y position (mm)','FontSize',22);
 title(ax,strrep(ids{k},'_',' '),'FontSize',24);
 colormap(ax,cm);caxis(ax,limits);cb=colorbar(ax);cb.Label.String='dF/F';cb.FontSize=18;
 keep{k}=[ids{k} '_flat_path_dff_imaging_rate_scatter10_p05_p95_blue_gray_red.svg'];
 print(fig,fullfile(out,keep{k}),'-dsvg','-painters');close(fig);
 fprintf('%s: %d acquired samples, limits %.4f %.4f\n',ids{k},numel(d),limits);
end
% Cleanup only this explicitly requested output folder, and only after all
% four SVG exports exist. Never touch source trial folders or time series.
assert(all(cellfun(@(s)isfile(fullfile(out,s)),keep)));
files=dir(out);deleted=0;preserved=0;
for j=1:numel(files)
 if files(j).isdir,continue,end
 name=files(j).name;
 if ismember(name,keep)||contains(lower(name),'time_traces')||contains(lower(name),'timeseries')
  preserved=preserved+1;continue
 end
 delete(fullfile(out,name));deleted=deleted+1;
end
fprintf('Kept %d files (4 final SVGs + all time-series files); removed %d superseded files.\n',preserved,deleted);
end

function data=read_saved(file,field)
% Resolve saved numeric struct fields without loading image stacks.
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
