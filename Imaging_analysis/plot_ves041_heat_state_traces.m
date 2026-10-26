function plot_ves041_heat_state_traces(a,savepath)
% Same three raw traces as fwd_rot_vs_dff_separate, colored by heat command.
% Heat on: measured DAQ heat command >0.05 V. Nonfinite command is uncolored.
% No smoothing, event selection, fluorescence normalization, or rejection.
if ~isfield(a.dq(1),'heat')||isempty(a.dq(1).heat),return,end
 t=a.dq(1).t(:);heat=a.dq(1).heat(:);d=a.roi.ts(1,:)';f=a.dq(1).bvf(:);
if isfield(a.dq(1),'bvy_deg'),r=a.dq(1).bvy_deg(:);else,r=rad2deg(a.dq(1).bvy(:));end
assert(all([numel(heat) numel(d) numel(f) numel(r)]==numel(t)),'Heat and trace lengths differ');
assert(all(isfinite(t))&&all(diff(t)>0),'Time must be finite and increasing');
if ~any(isfinite(heat)&heat>.05),return,end
if ~isfolder(savepath),mkdir(savepath);end
% Reuse the saved figure's actual signals and limits when available, so the
% requested variant changes only color, even if the MAT file was updated.
vars={d,f,r};times={t,t,t};limits=cell(1,3);
original=fullfile(savepath,'fwd_rot_vs_dff_separate.fig');
if isfile(original)
    src=openfig(original,'invisible');
    sourceAxes=findall(src,'Type','axes');
    assert(numel(sourceAxes)==3,'Expected three axes in the original figure');
    positions=cell2mat(get(sourceAxes,'Position'));[~,order]=sort(positions(:,2),'descend');
    for v=1:3
        ax=sourceAxes(order(v));traces=findall(ax,'Type','line');
        [~,j]=max(arrayfun(@(z)numel(z.XData),traces));
        times{v}=traces(j).XData(:);vars{v}=traces(j).YData(:);limits{v}=ax.YLim;
    end
    close(src);
end
fig=figure('Visible','off','Color','w','Position',[50 50 1500 950]);
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','loose');
labels={'dF/F','Forward velocity (mm/s)','Rotational velocity (deg/s)'};
colors=[.35 .35 .35;.9 .25 .08];names={'Heat off','Heat on'};axs=gobjects(3,1);
% Split every interval by its left sample's command state; include its end
% point so transitions have no artificial missing trace segment. Outages break.
for v=1:3
    ax=nexttile(tl);axs(v)=ax;hold(ax,'on');y=vars{v};tt=times{v};
    hh=interp1(t,heat,tt,'previous',NaN);state=double(hh>.05);state(~isfinite(hh))=NaN;
    gap=diff(tt)>2*median(diff(tt));
    for s=0:1
        ix=find(state(1:end-1)==s&~gap);
        xx=[tt(ix) tt(ix+1) nan(numel(ix),1)]';
        yy=[y(ix) y(ix+1) nan(numel(ix),1)]';
        plot(ax,xx(:),yy(:),'Color',colors(s+1,:),'LineWidth',1.2,'DisplayName',names{s+1});
    end
    ylabel(ax,labels{v});title(ax,labels{v});grid(ax,'on');xlim(ax,[t(1) t(end)]);xlabel(ax,'Time (s)');
    if v==2,ylim(ax,[-5 20]);end
    if ~isempty(limits{v}),ylim(ax,limits{v});end
    if v==1,legend(ax,'Location','northeast','Orientation','horizontal');end
end
linkaxes(axs,'x');[~,trial]=fileparts(fileparts(savepath));
if strcmp(trial,'smr_analysis'),[~,trial]=fileparts(fileparts(fileparts(savepath)));end
title(tl,{strrep(trial,'_',' '),'dF/F, forward, and rotational velocity | heat on: measured command >0.05 V'});
stem=fullfile(savepath,'fwd_rot_vs_dff_separate_heat_on_off');
savefig(fig,[stem '.fig']);exportgraphics(fig,[stem '.png'],'Resolution',150);close(fig);
end
