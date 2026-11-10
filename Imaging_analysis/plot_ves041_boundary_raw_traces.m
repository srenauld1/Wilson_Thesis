function plot_ves041_boundary_raw_traces(showHeat,lineWidth)
% Plot saved, imaging-aligned measurements without additional smoothing.
% Rotational speed is abs(saved signed yaw velocity), not a filtered trace.
% This does not undo upstream acquisition/resampling or dF/F normalization.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
if nargin<1,showHeat=false;end
if nargin<2,lineWidth=1;end
assert(isscalar(lineWidth)&&isfinite(lineWidth)&&lineWidth>0);
out=fullfile(fileparts(mfilename('fullpath')),'diagnostics', ...
    'ves041_20261020_12_1_boundary_turns_20261102');
S=load(fullfile(out,'boundary_turns_data.mat'),'t','dff','fwd','yaw','heat');
keep=S.t>=314&S.t<=325;
fig=figure('Visible','off','Color','w','Position',[20 20 1450 800]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
values={S.dff,S.fwd,abs(S.yaw)};
colors={[0 0 .9],[0 0 0],[.9 0 0]};
labels={'dF/F','Forward velocity (mm/s)','Rotational speed (deg/s)'};
axs=gobjects(3,1);
for k=1:3
    axs(k)=nexttile(tl);
    plot(axs(k),S.t(keep),values{k}(keep),'Color',colors{k},'LineWidth',lineWidth);
    ylabel(axs(k),labels{k});grid(axs(k),'on');xlim(axs(k),[314 325]);
    if showHeat
        hold(axs(k),'on');tt=S.t(keep);on=S.heat(keep)>.05;
        change=diff([false;on;false]);starts=find(change==1);ends=find(change==-1)-1;
        yl=ylim(axs(k));
        for j=1:numel(starts)
            right=tt(ends(j));if ends(j)<numel(tt),right=tt(ends(j)+1);end
            p=patch(axs(k),[tt(starts(j)) right right tt(starts(j))], ...
                [yl(1) yl(1) yl(2) yl(2)],[1 .55 .25], ...
                'FaceAlpha',.16,'EdgeColor','none');
            uistack(p,'bottom');
        end
        ylim(axs(k),yl);
    end
end
linkaxes(axs,'x');xlabel(axs(3),'Time from recording start (s)');
title(tl,'20261020 12_1 | 314-325 s | saved traces, no added smoothing');
stem=fullfile(out,'boundary_dff_fwd_rot_raw_314_325');
if showHeat
    subtitle(tl,'Orange: heat on (recorded monitor >0.05 V); unshaded: heat off | no added smoothing');
    stem=[stem '_heat_on_off'];
end
if lineWidth~=1,stem=[stem '_thick_lines'];end
exportgraphics(fig,[stem '.png'],'Resolution',170);
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
end
