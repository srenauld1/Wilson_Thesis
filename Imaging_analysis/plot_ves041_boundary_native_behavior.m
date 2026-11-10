function plot_ves041_boundary_native_behavior(signedRotation,presentation)
% Native CSV kinematics with imaging linearly interpolated to their DAQ clock.
% Uses saved coded-sync CSV-row/update-time map, NOT host/camera timestamps.
% Plot forward_velocity_rad_s * measured ball radius (4.5 mm); rotational
% SPEED = abs(yaw_velocity_rad_s)*180/pi. No differentiation, volume averaging,
% temporal smoothing, interpolation of kinematics, or artifact removal.
% Raw signed CSV yaw has a different sign convention from dq heading-derived
% yaw; this figure uses absolute rotational speed, for which sign is irrelevant.
% Fluorescence is unchanged source roi.ts, interpolated linearly only inside
% finite, valid imaging segments (no extrapolation or gap bridging). The
% interpolated line is shown without sample markers. Heat shading uses the SAME saved recorded-monitor
% intervals as the previous figure, not the CSV commanded-heat values.
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
if nargin<1,signedRotation=false;end
if nargin<2,presentation=false;end
folder='/Users/sophiarenauld/stacks/2026-10-20/20261020_12_1_ves041_syt8m_fly1_stripes';
out=fullfile(fileparts(mfilename('fullpath')),'diagnostics', ...
    'ves041_20261020_12_1_boundary_turns_20261102');
S=load(fullfile(out,'boundary_turns_data.mat'),'t','dff','heat','valid');
q=load(fullfile(folder,'20261020_12_1_bergii_.mat'),'bergii');b=q.bergii;
assert(isequal(double(b.volumeTimes(:)),S.t(:)),'Imaging/DAQ clock mismatch');
[~,name,ext]=fileparts(b.sourceCsv);sourceCsv=fullfile(folder,[name ext]);
opts=detectImportOptions(sourceCsv,'Delimiter',',');
opts.SelectedVariableNames={'frame','kinematic_valid','forward_velocity_rad_s','yaw_velocity_rad_s'};
T=readtable(sourceCsv,opts);
updateTime=double(b.updateTimes(:));rowIds=double(b.csvRowIndicesOneBased(:));
assert(numel(updateTime)==numel(rowIds)&&all(diff(updateTime)>0));
assert(all(rowIds>=1)&all(rowIds<=height(T))&all(rowIds==round(rowIds)));
take=updateTime>=314&updateTime<=325;nativeTime=updateTime(take);rowIds=rowIds(take);
ok=strcmpi(string(T.kinematic_valid(rowIds)),"true")|string(T.kinematic_valid(rowIds))=="1";
nativeForward=double(T.forward_velocity_rad_s(rowIds))*double(b.ballDiameterMm)/2;
nativeRotSpeed=abs(rad2deg(double(T.yaw_velocity_rad_s(rowIds))));
% CSV yaw has the opposite sign to the saved dq.bvy_deg heading convention.
% Negate it for a signed trace consistent with the earlier analysis plots.
nativeRotVelocity=-rad2deg(double(T.yaw_velocity_rad_s(rowIds)));
ok=ok&isfinite(nativeForward)&isfinite(nativeRotSpeed);
nativeForward(~ok)=NaN;nativeRotSpeed(~ok)=NaN;
nativeRotVelocity(~ok)=NaN;
imagingTime=S.t(:);dff=S.dff(:);imagingDt=median(diff(imagingTime));
valid=isfinite(dff)&logical(S.valid(:));cut=[true;diff(imagingTime)>2*imagingDt];
starts=find(valid&([true;~valid(1:end-1)]|cut));
ends=find(valid&([~valid(2:end);true]|[cut(2:end);true]));
upsampledDff=nan(size(nativeTime));
for k=1:numel(starts)
    ix=starts(k):ends(k);if numel(ix)<2,continue,end
    v=nativeTime>=imagingTime(ix(1))&nativeTime<=imagingTime(ix(end));
    upsampledDff(v)=interp1(imagingTime(ix),dff(ix),nativeTime(v),'linear',NaN);
end
assert(all(isfinite(upsampledDff)),'Missing imaging coverage in displayed interval');
nativeHz=1/median(diff(nativeTime));imagingHz=1/imagingDt;
fig=figure('Visible','off','Color','w','Position',[20 20 1450 850]);
if presentation,fig.Position=[20 20 950 900];end
tl=tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');axs=gobjects(3,1);
data={upsampledDff,nativeForward,nativeRotSpeed};colors={[.9 0 0],[0 0 0],[0 0 .9]};
labels={'dF/F','Native forward velocity (mm/s)','Native rotational speed (deg/s)'};
if signedRotation
    data{3}=nativeRotVelocity;labels{3}='Native rotational velocity (deg/s)';
end
if presentation
    labels={'dF/F',{'Forward velocity','(mm/s)'},{'Rotational speed','(deg/s)'}};
    if signedRotation,labels{3}={'Rotational velocity','(deg/s)'};end
end
imageTake=imagingTime>=314&imagingTime<=325&valid;
heatTake=imagingTime>=314&imagingTime<=325;heatTime=imagingTime(heatTake);
heatOn=S.heat(heatTake)>.05;transition=diff([false;heatOn;false]);on=find(transition==1);off=find(transition==-1)-1;
for k=1:3
    ax=nexttile(tl);axs(k)=ax;hold(ax,'on');
    plot(ax,nativeTime,data{k},'Color',colors{k},'LineWidth',3);
    if k==3&&signedRotation,yline(ax,0,'-','Color',[.1 .1 .1],'LineWidth',1.8);end
    ylabel(ax,labels{k});grid(ax,'on');xlim(ax,[314 325]);yl=ylim(ax);
    for j=1:numel(on)
        right=heatTime(off(j));if off(j)<numel(heatTime),right=heatTime(off(j)+1);end
        p=patch(ax,[heatTime(on(j)) right right heatTime(on(j))], ...
            [yl(1) yl(1) yl(2) yl(2)],[1 .55 .25],'FaceAlpha',.16,'EdgeColor','none');
        uistack(p,'bottom');
    end
    ylim(ax,yl);
    if presentation
        ax.FontSize=18;ax.LineWidth=1.2;ax.Box='off';ax.GridAlpha=.1;
        ylabel(ax,labels{k},'FontSize',22);ax.XTick=314:325;
        if k<3,ax.XTickLabel=[];end
    end
end
linkaxes(axs,'x');xlabel(axs(3),'Time from recording start (s)');
title(tl,'20261020 12_1 | 314-325 s | native behavior, no smoothing');
subtitle(tl,sprintf('dF/F: linear interpolation from %.1f to %.0f Hz | Behavior: native %.0f Hz | Orange: recorded heat on',imagingHz,nativeHz,nativeHz));
if presentation
    title(tl,'VES041 | heat stripes','FontSize',24,'FontWeight','bold');
    subtitle(tl,'');xlabel(axs(3),'Time (s)','FontSize',22);
    % Dummy patches give an explicit two-color key without changing data.
    hOn=patch(axs(1),NaN,NaN,[1 .55 .25],'FaceAlpha',.16,'EdgeColor',[.8 .65 .55]);
    hOff=patch(axs(1),NaN,NaN,[1 1 1],'EdgeColor',[.6 .6 .6]);
    lg=legend(axs(1),[hOn hOff],{'Heat on','Heat off'}, ...
        'Orientation','horizontal','Box','off','FontSize',20);
    lg.Layout.Tile='north';
end
stem=fullfile(out,'boundary_native_kinematics_imaging_upsampled_314_325');
if signedRotation,stem=[stem '_signed_rotation'];end
if presentation,stem=[stem '_presentation'];end
exportgraphics(fig,[stem '.png'],'Resolution',170);
if presentation,exportgraphics(fig,[stem '.pdf'],'ContentType','vector');end
savefig_visible_on_open(fig,[stem '.fig']);close(fig);
save([stem '.mat'],'nativeTime','nativeForward','nativeRotSpeed','nativeRotVelocity','signedRotation','upsampledDff', ...
    'rowIds','sourceCsv','nativeHz','imagingHz','imagingTime','dff','ok');
fprintf('Native %.3f Hz, imaging %.3f Hz; %d native samples, %d valid, %d imaging samples\n', ...
    nativeHz,imagingHz,numel(nativeTime),nnz(ok),nnz(imageTake));
end
