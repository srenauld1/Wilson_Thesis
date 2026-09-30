function summary=plot_ves041_spontaneous_reorientation_phase_space(selectedTrials,outputName)
% Plot stimulus-free behavior after VES041 peaks in speed-heading space.
% Each highlighted point is behavior during 0 to +1 s after a VES041 peak.
% Gray points are non-overlapping 1 s control windows sampled through the
% same trial. X is mean forward speed; Y is absolute net heading change.
% Dashed lines show the working low-speed reorientation definition:
% mean forward speed <=2 mm/s and net heading change >=30 degrees.
%
% Stimulus exclusion:
%   bar jumps: +/-2 s around each detected jump
%   optomotor: all moving-grating samples
%   heat: all heat-on samples
% A candidate is retained only if its complete -2 to +1 s neighborhood is
% stimulus-free. Odor trials are omitted because their odor-state encoding
% is not standardized across the current processed files.

base='/Users/sophiarenauld/stacks';
out=fullfile(base,'ves041_behavior_quantification','spontaneous_reorientation_phase_space');
if nargin>=2 && strlength(string(outputName))>0,out=fullfile(out,char(outputName));end
if ~isfolder(out),mkdir(out);end
standard=["20250812-20_ves041_fly3_clbarjumps";"20250812-21_ves041_fly3_olgrating"; ...
 "20250812-22_ves041_fly3_bbarjumps2";"20250813-3_ves041_fly1_olgrating"; ...
 "20250819-10_ves041_fly1_olgrating5";"20251001-21_ves041_fly3_clbbar"; ...
 "20260213-2_ves041_clbar";"20260213-4_ves041_fly2_clbar"; ...
 "20260217-2_ves041_vt_clbbar1";"20260217-3_ves041_vt_clbbar2_actuallydark"; ...
 "20260217-4_ves041_vt_clbbar3"];
heat=["20260427-2_ves041_fly2_heat";"20260427-3_ves041_fly2_heat_20mmbarrier"; ...
 "20260427-4_ves041_fly2_heat_20mmbarrier_2";"20260427-6_ves041_fly3_heat_80mm"; ...
 "20260427-7_ves041_fly3_heat_20mm"];
trials=[standard;heat];
if nargin>=1 && ~isempty(selectedTrials),trials=string(selectedTrials(:));end
pooledControl=[];pooledPeak=[];pooledPre=[];rows=cell(0,7);
for i=1:numel(trials)
    trial=trials(i);isHeat=any(heat==trial);
    if isHeat
        [t,dff,fwd,heading,xpos,ypos,stim]=load_heat(fullfile(base,trial));
    else
        file=fullfile(base,trial,'smr_analysis','fly_data.mat');
        if ~isfile(file),warning('Skipping missing %s',file);continue,end
        q=load(file);x=q.fly_data;t=x.dq(1).t(:);dff=x.roi(1).ts(:);
        fwd=x.dq(1).bvf(:);heading=x.dq(1).bh(:);xpos=x.dq(1).px(:);ypos=x.dq(1).py(:);stim=false(size(t));
        if contains(lower(trial),'grating')
            vv=x.dq(1).vvy(:);m=min(numel(stim),numel(vv));stim(1:m)=abs(vv(1:m))>.02;
        end
        if isfield(x,'jump_detected')
            j=find(logical(x.jump_detected(:)));j=j(j<=numel(t));
            if isfield(x.dq(1),'vh')
                vh=x.dq(1).vh(:);dv=atan2(sin([vh(1);diff(vh)]),cos([vh(1);diff(vh)]));
                j=j(abs(rad2deg(dv(j)))>=45);
            end
            fs=1/median(diff(t));for k=j(:)',stim(max(1,k-round(2*fs)):min(numel(stim),k+round(2*fs)))=true;end
        end
    end
    n=min([numel(t),numel(dff),numel(fwd),numel(heading),numel(xpos),numel(ypos),numel(stim)]);
    t=t(1:n);dff=dff(1:n);fwd=fwd(1:n);heading=heading(1:n);xpos=xpos(1:n);ypos=ypos(1:n);stim=stim(1:n);
    fs=1/median(diff(t));preN=round(2*fs);postN=round(fs);
    forbidden=movmax(double(stim),[postN preN])>0; % center requires -2:+1 clean
    peaks=merged_local_peaks(dff,fs);
    peaks=peaks(peaks>preN & peaks+postN<=n & ~forbidden(peaks));
    controls=(1+preN):round(fs):(n-postN);
    controls=controls(~forbidden(controls));
    if ~isempty(peaks),controls=controls(all(abs(controls(:)-peaks(:)')>round(2*fs),2));end
    P=metrics(peaks,fwd,heading,xpos,ypos,postN);Pre=metrics(peaks-preN,fwd,heading,xpos,ypos,postN);
    C=metrics(controls,fwd,heading,xpos,ypos,postN);
    pooledPeak=[pooledPeak;P];pooledPre=[pooledPre;Pre];pooledControl=[pooledControl;C]; %#ok<AGROW>
    rows(end+1,:)={trial,numel(peaks),size(C,1),sum(P(:,1)<=2 & P(:,2)>=30), ...
        mean(P(:,1),'omitnan'),mean(P(:,2),'omitnan'),mean(P(:,1)<=2 & P(:,2)>=30,'omitnan')}; %#ok<AGROW>
    make_plot(C,P,trial,fullfile(out,trial+'_spontaneous_phase_space'));
end
summary=cell2table(rows,'VariableNames',{'trial','nStimulusFreeVesPeaks','nControlWindows', ...
    'nPeakLowSpeedReorientations','meanPostPeakFwd','meanPostPeakHeadingChange', ...
    'fractionPeaksFollowedByLowSpeedReorientation'});
summary.trial=string(summary.trial);writetable(summary,fullfile(out,'trial_summary.csv'));
make_plot(pooledControl,pooledPeak,'All good-behavior trials (stimulus epochs removed)', ...
    fullfile(out,'ALL_TRIALS_spontaneous_phase_space'));
writematrix(pooledPeak,fullfile(out,'pooled_peak_windows_fwd_heading.csv'));
writematrix(pooledPre,fullfile(out,'pooled_prepeak_windows_fwd_heading.csv'));
writematrix(pooledControl,fullfile(out,'pooled_control_windows_fwd_heading.csv'));
disp(summary);
end

function peaks=merged_local_peaks(dff,fs)
% One peak per locally elevated calcium event.
smoothDff=smoothdata(dff,'gaussian',max(3,round(.25*fs)));
baseline=movmedian(smoothDff,max(3,round(30*fs)),'omitmissing');
signal=smoothDff-baseline;noise=mad(diff(signal),1)/(sqrt(2)*.6745);
threshold=max(.5,3*noise);active=signal>=threshold;
d=diff([false;active(:);false]);s=find(d==1);e=find(d==-1)-1;
if isempty(s),peaks=[];return,end
maxGap=round(.5*fs);ms=s(1);me=[];ce=e(1);
for i=2:numel(s)
    if s(i)-ce-1<=maxGap,ce=e(i);else,me(end+1,1)=ce;ms(end+1,1)=s(i);ce=e(i);end %#ok<AGROW>
end
me(end+1,1)=ce;keep=(me-ms+1)>=max(1,round(.25*fs));ms=ms(keep);me=me(keep);
peaks=zeros(numel(ms),1);for i=1:numel(ms),[~,j]=max(signal(ms(i):me(i)));peaks(i)=ms(i)+j-1;end
end

function M=metrics(idx,fwd,heading,x,y,n)
M=nan(numel(idx),3);
for k=1:numel(idx)
    j=idx(k);M(k,1)=mean(fwd(j:j+n),'omitnan');
    M(k,2)=abs(rad2deg(atan2(sin(heading(j+n)-heading(j)),cos(heading(j+n)-heading(j)))));
    dx=diff(x(j:j+n));dy=diff(y(j:j+n));pathLength=sum(hypot(dx,dy),'omitnan');
    if pathLength>=0.5,M(k,3)=hypot(x(j+n)-x(j),y(j+n)-y(j))/pathLength;end
end
end

function make_plot(C,P,titleText,path)
fig=figure('Visible','off','Color','w','Position',[100 100 850 680]);hold on;
scatter(C(:,1),C(:,2),18,[.7 .7 .7],'filled','MarkerFaceAlpha',.22,'DisplayName','Control 1 s windows');
scatter(P(:,1),P(:,2),65,[.82 .15 .12],'filled','MarkerEdgeColor','w','LineWidth',.5, ...
    'DisplayName','0-1 s after VES041 peak');
hx=xline(2,'k--','Low-speed boundary');hx.HandleVisibility='off';
hy=yline(30,'k--','Reorientation boundary');hy.HandleVisibility='off';
xlabel('Mean forward speed, 0 to +1 s (mm/s)');ylabel('|Net heading change|, 0 to +1 s (deg)');
title(strrep(titleText,'_','\_'));grid on;legend('Location','northeast');
xlim([max(-5,min([C(:,1);P(:,1)],[],'omitnan')) min(20,max([C(:,1);P(:,1)],[],'omitnan'))]);ylim([0 180]);
text(.02,.98,sprintf('VES peaks: %d\nLow-speed reorientations: %d (%.1f%%)',size(P,1), ...
    sum(P(:,1)<=2 & P(:,2)>=30),100*mean(P(:,1)<=2 & P(:,2)>=30,'omitnan')), ...
    'Units','normalized','VerticalAlignment','top','BackgroundColor','w');
path=char(path);savefig(fig,[path '.fig']);exportgraphics(fig,[path '.png'],'Resolution',220);close(fig);
end

function [t,dff,fwd,heading,xpos,ypos,stim]=load_heat(folder)
f=dir(fullfile(folder,'*_virmenData.mat'));v=load(fullfile(f(1).folder,f(1).name));
f=dir(fullfile(folder,'*_or_z0z1_s_.mat'));q=load(fullfile(f(1).folder,f(1).name));
t=seconds(v.vmData.Seconds{1,1});t=t(:);fwd=smooth(v.vmData.vel_for{1,1});
heading=v.vmData.heading{1,1};heading=heading(:);xpos=v.vmData.pos_x{1,1};xpos=xpos(:);ypos=v.vmData.pos_y{1,1};ypos=ypos(:);
h=v.vmData.heat;if iscell(h),h=cell2mat(h);end;stim=h(:)>1e-4;
td=v.daqVolumeTime.daqVolumeTime{1,1};dff=interp1(td(:),q.roi(1).ts(:),t,'linear',NaN);
end
