function plot_ves041_boundary_native_flatpath()
% Presentation flat path using 50-Hz native logged coordinates and dF/F
% linearly interpolated from imaging timestamps onto the same coded-sync
% DAQ update times. No smoothing, extrapolation, or extra fluorescence data.
% Keep the original interval-specific 5th/95th percentile limits calculated
% on ACQUIRED imaging samples, so interpolation does not change the scale.
% Colors outside the limits saturate; no samples are discarded for color.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
out=fullfile(fileparts(mfilename('fullpath')),'diagnostics', ...
    'ves041_20261020_12_1_boundary_turns_20261102');
S=load(fullfile(out,'boundary_native_kinematics_imaging_upsampled_314_325_signed_rotation_presentation.mat'), ...
    'nativeTime','upsampledDff','rowIds','sourceCsv','ok');
C=load(fullfile(out,'flat_path_dff_314_325_p05_p95.mat'),'limits');limits=C.limits;
opts=detectImportOptions(S.sourceCsv,'Delimiter',',');opts.SelectedVariableNames={'protocol_state_json'};
T=readtable(S.sourceCsv,opts);n=numel(S.nativeTime);x=nan(n,1);y=x;inside=false(n,1);
for k=1:n
    state=jsondecode(char(string(T.protocol_state_json(S.rowIds(k)))));
    x(k)=state.pos_x_mm;y(k)=state.pos_y_mm;inside(k)=state.in_heat_stripe;
end
assert(all(isfinite(x))&&all(isfinite(y))&&all(S.ok));
% The displayed region crosses only this one horizontal boundary. Verify
% the background against the native protocol's logged zone membership.
boundary=-7.5;assert(mean((y<boundary)==inside)>.98,'Background does not match native stripe log');
dff=S.upsampledDff(:);t=S.nativeTime(:);
assert(all(t>=314&t<=325)&&all(isfinite(dff)));
fig=figure('Visible','off','Color','w','Position',[20 20 1250 450]);
ax=axes(fig);hold(ax,'on');
heatColor=[1 .928 .88];
hOn=patch(ax,[min(x)-.5 max(x)+.5 max(x)+.5 min(x)-.5], ...
    [min(y)-.4 min(y)-.4 boundary boundary],heatColor,'EdgeColor','none');
hOff=patch(ax,NaN,NaN,[1 1 1],'EdgeColor',[.6 .6 .6]);
% Break potential missing native update intervals rather than connecting gaps.
xx=x;yy=y;gap=[false;diff(t)>.15];xx(gap)=NaN;yy(gap)=NaN;
surface(ax,[xx xx],[yy yy],zeros(n,2),[dff dff], ...
    'FaceColor','none','EdgeColor','interp','LineWidth',4);
plot(ax,x(1),y(1),'go','MarkerFaceColor','g','MarkerSize',7);
plot(ax,x(end),y(end),'ks','MarkerFaceColor','k','MarkerSize',6);
axis(ax,'equal');xlim(ax,[min(x)-.5 max(x)+.5]);ylim(ax,[min(y)-.4 max(y)+.4]);
ax.FontSize=18;ax.LineWidth=1.1;ax.Box='off';grid(ax,'on');ax.GridAlpha=.1;
xlabel(ax,'X position (mm)','FontSize',22);ylabel(ax,'Y position (mm)','FontSize',22);
colormap(ax,parula(256));caxis(ax,limits);cb=colorbar(ax);cb.FontSize=18;cb.Label.String='dF/F';cb.Label.FontSize=22;
title(ax,'VES041 | 314-325 s','FontSize',24,'FontWeight','bold');
legend(ax,[hOn hOff],{'Heat region','Cool region'},'Location','northoutside', ...
    'Orientation','horizontal','Box','off','FontSize',18);
stem=fullfile(out,'flat_path_dff_314_325_native50Hz_p05_p95_presentation');
exportgraphics(fig,[stem '.png'],'Resolution',180);
exportgraphics(fig,[stem '.pdf'],'ContentType','vector');
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
sourceCsv=S.sourceCsv;rowIds=S.rowIds;nativeTime=t;
save([stem '.mat'],'nativeTime','x','y','dff','inside','limits','sourceCsv','rowIds');
fprintf('Native path samples: %d; acquired-sample color limits: %.4f %.4f\n',n,limits);
end
