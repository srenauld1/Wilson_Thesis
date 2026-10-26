function E=plot_ves041_peak_heat_context(a,out,locs,prominence)
% Save every detected peak, plus descriptive heat-context aligned figures.
% Optional locs/prominence reuse the main routine's exact detector output.
t=a.dq(1).t(:);d=a.roi.ts(1,:)';d=d(:);f=a.dq(1).bvf(:);
if isfield(a.dq(1),'bvy_deg'),r=abs(a.dq(1).bvy_deg(:));else,r=abs(rad2deg(a.dq(1).bvy(:)));end
fs=1/median(diff(t));assert(all(diff(t)>0),'Time must increase');
if nargin<3
    baseline=movmedian(d,max(3,round(30*fs)),'omitmissing');
    residual=d-baseline;residual=residual(isfinite(residual));
    prominence=3*1.4826*mad(residual,1);
    % Split invalid samples/outages instead of bridging them for detection.
    good=isfinite(d);beg=find(good&[true;~good(1:end-1)|diff(t)>2/fs]);
    ends=find(good&[~good(2:end)|diff(t)>2/fs;true]);locs=[];
    for j=1:numel(beg)
        ix=beg(j):ends(j);if numel(ix)<3,continue,end
        [~,q]=findpeaks(d(ix),'MinPeakProminence',max(eps,prominence), ...
            'MinPeakDistance',max(1,round(2*fs)),'MinPeakWidth',max(1,round(.25*fs)));
        locs=[locs;ix(q(:))']; %#ok<AGROW>
    end
end
locs=locs(:);assert(all(locs>=1&locs<=numel(t)),'Peak indices outside trace');
if ~isfolder(out),mkdir(out);end
h=nan(size(t));if isfield(a.dq(1),'heat')&&numel(a.dq(1).heat)==numel(t),h=a.dq(1).heat(:);end
if ~isfield(a,'protocol'),a=import_ves041_protocol(a,fileparts(out));end
control=nan(size(t));
if isfield(a,'protocol')&&isfield(a.protocol,'control_trial'),control=a.protocol.control_trial(:);end
contexts=["no heat" "heat" "heat entry" "heat exit" "control" "multiple switches" "control transition" "unknown"];
colors=[.35 .35 .35;.9 .25 .08;.15 .6 .3;.6 .25 .7;.15 .5 .8;.7 .5 .1;.25 .65 .7;.1 .1 .1];
C=repmat("unknown",numel(locs),1);eligible=false(size(locs));
for j=1:numel(locs)
    q=locs(j);core=t>=t(q)-.5&t<=t(q)+.5;
    if t(q)-.5>=t(1)&&t(q)+.5<=t(end)&&all(isfinite(h(core)))&&all(diff(t(core))<=2/fs)
        on=h(core)>.05;
        C(j)=classify_ves041_peak_heat_context(on,control(core));
    end
    win=t>=t(q)-2&t<=t(q)+3;
    eligible(j)=t(q)-2>=t(1)&&t(q)+3<=t(end)&&all(isfinite(d(win))&isfinite(f(win))&isfinite(r(win)))&&all(diff(t(win))<=2/fs);
    if isfield(a.dq(1),'bergii')&&isfield(a.dq(1).bergii,'behaviorVolumeValid')
        valid=a.dq(1).bergii.behaviorVolumeValid(:);eligible(j)=eligible(j)&&all(valid(win));
    end
end
E=table(locs,t(locs),d(locs),h(locs),C,eligible,'VariableNames', ...
    {'sampleIndex','peakTime','peakDff','heatCommandAtPeakV','context','eligibleAlignedWindow'});
writetable(E,fullfile(out,'ves041_detected_peaks.csv'));
settings=struct('prominenceDff',prominence,'minimumPeakDistanceSec',2, ...
    'minimumPeakWidthSec',.25,'baselineWindowSec',30,'heatThresholdV',.05, ...
    'contextWindowSec',[-.5 .5],'alignedWindowSec',[-2 3]);
settings.contextDefinition='Exactly one off-to-on crossing = heat entry; exactly one on-to-off crossing = heat exit within +/-0.5 s of the activity peak. Multiple switches and control-only transitions remain separate in the diagnostic and CSV. Alignment stays at the activity peak.';
settings.traceSource='Supplied ROI 1 and DAQ signals';
if isfield(a,'peakTraceSource'),settings.traceSource=a.peakTraceSource;end
save(fullfile(out,'ves041_peak_detection_settings.mat'),'settings','E');
fig=figure('Visible','off','Color','w','Position',[40 40 1550 700]);tl=tiledlayout(fig,2,1,'TileSpacing','compact');
ax=nexttile(tl);hold(ax,'on');plot(ax,t,d,'Color',[.35 .35 .4],'DisplayName','dF/F');
for c=1:numel(contexts)
    sel=C==contexts(c);if any(sel),scatter(ax,t(locs(sel)),d(locs(sel)),35,colors(c,:),'v','filled','DisplayName',sprintf('%s peaks (n=%d)',contexts(c),sum(sel)));end
end
xlabel(ax,'Time (s)');ylabel(ax,'dF/F');grid(ax,'on');xlim(ax,[t(1) t(end)]);legend(ax,'Location','best');
title(ax,sprintf('All %d selected peaks, including edges | prominence >= %.3g | separation >=2 s | width >=0.25 s',numel(locs),prominence));
ax2=nexttile(tl);plot(ax2,t,h,'k');hold(ax2,'on');yline(ax2,.05,':r','0.05 V');xlabel(ax2,'Time (s)');ylabel(ax2,'Measured heat command (V)');grid(ax2,'on');xlim(ax2,[t(1) t(end)]);linkaxes([ax ax2],'x');
[~,trial]=fileparts(fileparts(out));title(tl,strrep(trial,'_',' '));
saveplot(fig,out,'ves041_detected_peaks');
% Every analysis gets the detection diagnostic above. Heat comparisons are
% produced only when the recording contains a finite, positive heat command.
if ~any(h>.05),return,end
rel=(-2:.05:3)';X=nan(numel(rel),numel(locs),5);vars={d,d,f,r,h};
for j=find(eligible)'
    pre=t>=t(locs(j))-2&t<t(locs(j))-1;b=mean(d(pre));
    for v=1:5,X(:,j,v)=interp1(t,vars{v},t(locs(j))+rel,'linear',NaN);end
    X(:,j,2)=X(:,j,2)-b;
end
for panel=1:2
if panel==1,groups=1:5;stem='ves041_peaks_heat_context_aligned';else,groups=3:4;stem='ves041_peaks_heat_entry_exit_aligned';end
ncol=numel(groups);
fig=figure('Visible','off','Color','w','Position',[30 30 420*ncol 1250]);tl=tiledlayout(fig,5,ncol,'TileSpacing','compact');
labs={'Raw dF/F','dF/F minus pre-event baseline','Forward (mm/s)','Abs. yaw (deg/s)','Heat command (V)'};
for col=1:ncol
    c=groups(col);sel=find(eligible&C==contexts(c));
    for v=1:5
        ax=nexttile(tl,(v-1)*ncol+col);hold(ax,'on');
        if ~isempty(sel),plot(ax,rel,X(:,sel,v),'Color',[.8 .8 .8]);plot(ax,rel,mean(X(:,sel,v),2,'omitnan'),'Color',colors(c,:),'LineWidth',2);else,text(ax,.5,.5,'No eligible peaks','Units','normalized','HorizontalAlignment','center');end
        xline(ax,0,':k');ylabel(ax,labs{v});xlim(ax,[-2 3]);grid(ax,'on');
        if v==1,title(ax,sprintf('%s | n=%d',contexts(c),numel(sel)));end
        if v==5,xlabel(ax,'Time from activity peak (s)');end
    end
end
for v=1:5
    axs=gobjects(ncol,1);for col=1:ncol,axs(col)=nexttile(tl,(v-1)*ncol+col);end
    yl=cell2mat(get(axs,'YLim'));set(axs,'YLim',[min(yl(:,1)) max(yl(:,2))]);
end
title(tl,{strrep(trial,'_',' '),'Activity-selected peaks | context within +/-0.5 s | thin: individual events; color: mean | pre-event baseline: -2:-1 s'});
saveplot(fig,out,stem);
end
fig=figure('Visible','off','Color','w','Position',[30 30 2100 1000]);tl=tiledlayout(fig,3,5,'TileSpacing','compact');
for row=1:3
    v=row+1;val=X(:,:,v);finite=val(isfinite(val));if isempty(finite),lim=[0 1];else,lim=prctile(finite,[2 98]);end
    if diff(lim)<=0,lim=lim+[-.1 .1];end
    for c=1:5
        ax=nexttile(tl,(row-1)*5+c);sel=find(eligible&C==contexts(c));
        if isempty(sel),text(ax,.5,.5,'No eligible peaks','Units','normalized','HorizontalAlignment','center');else,imagesc(ax,rel,1:numel(sel),X(:,sel,v)');clim(ax,lim);colorbar(ax);end
        set(ax,'YDir','normal');hold(ax,'on');xline(ax,0,':k');xlim(ax,[-2 3]);ylabel(ax,'Peak (chronological)');
        title(ax,sprintf('%s | %s | n=%d',contexts(c),labs{v},numel(sel)));
        if row==3,xlabel(ax,'Time from activity peak (s)');end
    end
end
colormap(fig,parula);title(tl,{strrep(trial,'_',' '),'Individual activity-selected events | shared color limits by row (2nd-98th percentile display only)'});
saveplot(fig,out,'ves041_peaks_heat_context_heatmaps');
fprintf('%s: %d peaks; %d complete windows\n',trial,numel(locs),sum(eligible));
for c=1:numel(contexts),fprintf('  %s: %d selected / %d eligible\n',contexts(c),sum(C==contexts(c)),sum(eligible&C==contexts(c)));end
end
function saveplot(fig,out,stem)
savefig(fig,fullfile(out,[stem '.fig']));exportgraphics(fig,fullfile(out,[stem '.png']),'Resolution',140);close(fig);
end
