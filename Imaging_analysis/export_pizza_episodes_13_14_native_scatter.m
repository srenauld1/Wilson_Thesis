function export_pizza_episodes_13_14_native_scatter()
assert(~feature('ShowFigureWindows'),'Use -nodisplay');
base='/Users/sophiarenauld/stacks/2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo';
out='/Users/sophiarenauld/stacks/pizza_talk_analysis/example traces/SVG examples/barrier_episodes_13_14_20261104';
stem=fullfile(out,'20261020_11_2_barrier_episodes_13_heat_14_control_equal_xy_blue_red');
f=openfig([stem '.fig'],'invisible');axesList=findall(f,'Type','axes');
A=load(fullfile(base,'smr_analysis/fly_data.mat'),'fly_data');p=A.fly_data.protocol;ti=double(A.fly_data.dq(1).t(:));
z=dir(fullfile(base,'*_bergii_.mat'));B=load(fullfile(z.folder,z.name),'bergii');b=B.bergii;assert(isequal(double(b.volumeTimes(:)),ti));
opts=detectImportOptions(b.sourceCsv,'Delimiter',',');opts.SelectedVariableNames={'protocol_state_json','kinematic_valid'};C=readtable(b.sourceCsv,opts);
P=struct;pointSize=15;
for ep=[13 14]
 ax=axesList(arrayfun(@(a)startsWith(string(a.Title.String),sprintf('Barrier episode %d |',ep)),axesList));assert(isscalar(ax));
 h=findall(ax,'Type','surface');assert(isscalar(h));
 mask=p.trial_number(:)==ep&p.barrier_active(:)>0&isfinite(p.pos_x_mm(:))&isfinite(p.pos_y_mm(:))&isfinite(A.fly_data.roi.ts(1,:))';
 ix=find(mask);x=h.XData(:,1);y=h.YData(:,1);di=h.CData(:,1);assert(numel(ix)==numel(di));
 good=isfinite(x)&isfinite(y);assert(max(abs(x(good)-p.pos_x_mm(ix(good))))<1e-8&&max(abs(y(good)-p.pos_y_mm(ix(good))))<1e-8);
 % Use exactly the plotted fluorescence at its original imaging timestamps,
 % preserving the source figure even if an ROI file was subsequently edited.
 imagingTime=ti(ix);cut=[false;diff(ix)>1|diff(imagingTime)>2*median(diff(ti))];di(cut)=nan;
 ni=find(b.updateTimes>=imagingTime(1)&b.updateTimes<=imagingTime(end));tn=double(b.updateTimes(ni));rows=double(b.csvRowIndicesOneBased(ni));
 xy=nan(numel(tn),2);episode=nan(size(tn));active=false(size(tn));
 for j=1:numel(tn)
  s=jsondecode(char(string(C.protocol_state_json(rows(j)))));xy(j,:)=[s.pos_x_mm s.pos_y_mm];episode(j)=s.trial_number;active(j)=s.barrier_active;
 end
 dn=interp1(imagingTime,di,tn,'linear',nan);
 valid=ismember(lower(string(C.kinematic_valid(rows))),["true","1"])&episode==ep&active&all(isfinite(xy),2)&isfinite(dn);
 assert(nnz(valid)>numel(ix));assert(all(diff(tn)<.06));
 hold(ax,'on');delete(h);scatter(ax,xy(valid,1),xy(valid,2),pointSize,dn(valid),'filled','MarkerEdgeColor','none');
 % Start/end markers remain above the denser scatter points.
 lines=findall(ax,'Type','line');for j=1:numel(lines),if numel(lines(j).XData)==1,uistack(lines(j),'top');end;end
 k=ep-12;P(k).episode=ep;P(k).imagingTime=imagingTime;P(k).originalDff=di;P(k).nativeTime=tn(valid);P(k).xy=xy(valid,:);P(k).dff=dn(valid);P(k).csvRows=rows(valid);
 fprintf('Episode %d: %d imaging samples -> %d native scatter points, %.2f Hz\n',ep,numel(ix),nnz(valid),1/median(diff(tn)));
end
f.Position=[20 20 1700 510];
for k=1:numel(axesList)
 ax=axesList(k);j=1+startsWith(string(ax.Title.String),'Barrier episode 14 |');ax.FontSize=14;ax.Title.FontSize=18;ax.Position=[.065+(j-1)*.48 .20 .385 .65];
 cb=colorbar(ax);cb.FontSize=12;cb.Label.FontSize=16;cb.Position=[.458+(j-1)*.48 .265 .012 .48];ax.Position=[.065+(j-1)*.48 .20 .385 .65];
end
stem=[stem '_scatter_native50Hz'];
exportgraphics(f,[stem '.png'],'Resolution',180);exportgraphics(f,[stem '.pdf'],'ContentType','vector');print(f,[stem '.svg'],'-dsvg','-painters');savefig(f,[stem '.fig']);close(f);
methods=struct('fluorescence','Exact saved displayed path CData, linearly interpolated from imaging timestamps to coded-sync native kinematic update timestamps; no extrapolation or gap bridging','positions','Native CSV protocol X/Y in mm at the same native timestamps; valid tracking and matching active barrier episode only','pointAreaPt2',pointSize,'preserved','Shared spatial spans, equal mm aspect, original numeric color limits, study blue-gray-red palette, backgrounds, barrier lines, 10-mm scale bars and start/end markers','sourceCsv',b.sourceCsv);
save([stem '.mat'],'P','methods');fid=fopen([stem '_methods.json'],'w');fprintf(fid,'%s',jsonencode(methods,PrettyPrint=true));fclose(fid);copyfile([mfilename('fullpath') '.m'],out);
end
