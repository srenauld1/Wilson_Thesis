function results=compare_ves041_peak_amplitude_pivot_20260213()
% Compare VES041 event amplitude for pivot versus non-pivot events.
% A pivot is mean forward speed <=2 mm/s and >=30 deg absolute net heading
% change during 0 to +1 s after a merged local VES041 event.

base='/Users/sophiarenauld/stacks';out=fullfile(base,'ves041_behavior_quantification', ...
    'peak_amplitude_pivot_20260213');if ~isfolder(out),mkdir(out);end
trials=["20260213-2_ves041_clbar";"20260213-4_ves041_fly2_clbar"];
rows=cell(0,8);eventRows=cell(0,8);
fig=figure('Visible','off','Color','w','Position',[100 100 1350 470]);tl=tiledlayout(1,3,'TileSpacing','compact');
allAmp=[];allPivot=[];
for ti=1:numel(trials)
    trial=trials(ti);q=load(fullfile(base,trial,'smr_analysis','fly_data.mat'));x=q.fly_data;
    t=x.dq(1).t(:);dff=x.roi(1).ts(:);fwd=x.dq(1).bvf(:);heading=x.dq(1).bh(:);
    n=min([numel(t),numel(dff),numel(fwd),numel(heading)]);t=t(1:n);dff=dff(1:n);fwd=fwd(1:n);heading=heading(1:n);
    fs=1/median(diff(t));[allPeaks,allEventAmp,smoothDff,threshold]=events(dff,fs);
    edge=round(2*fs);keep=allPeaks>edge & allPeaks+round(fs)<=n;
    peaks=allPeaks(keep);amp=allEventAmp(keep);
    postFwd=arrayfun(@(i)mean(fwd(i:i+round(fs)),'omitnan'),peaks);
    turn=arrayfun(@(i)abs(rad2deg(atan2(sin(heading(i+round(fs))-heading(i)), ...
        cos(heading(i+round(fs))-heading(i))))),peaks);
    pivot=postFwd<=2 & turn>=30;absolutePeak=smoothDff(peaks);
    rows(end+1,:)={trial,numel(peaks),sum(pivot),mean(amp(pivot),'omitnan'), ... %#ok<AGROW>
        mean(amp(~pivot),'omitnan'),median(amp(pivot),'omitnan'),median(amp(~pivot),'omitnan'),threshold};
    for i=1:numel(peaks),eventRows(end+1,:)={trial,t(peaks(i)),amp(i),absolutePeak(i), ... %#ok<AGROW>
            postFwd(i),turn(i),pivot(i),threshold};end
    ax=nexttile(tl);draw_groups(ax,amp,pivot);title(ax,strrep(trial,'_','\_'));
    allAmp=[allAmp;amp(:)];allPivot=[allPivot;pivot(:)]; %#ok<AGROW>
end
ax=nexttile(tl);draw_groups(ax,allAmp,allPivot);title(ax,'Pooled 20260213 trials');
sgtitle(tl,'VES041 merged-event amplitude: low-speed pivots versus other events');
savefig(fig,fullfile(out,'peak_amplitude_pivot_comparison.fig'));
exportgraphics(fig,fullfile(out,'peak_amplitude_pivot_comparison.png'),'Resolution',250);close(fig);
results=cell2table(rows,'VariableNames',{'trial','nEvents','nPivots','meanAmplitudePivot', ...
    'meanAmplitudeOther','medianAmplitudePivot','medianAmplitudeOther','eventThreshold'});
eventsTable=cell2table(eventRows,'VariableNames',{'trial','eventTime','amplitudeAboveBaseline', ...
    'absoluteSmoothedDff','postFwd','postHeadingChangeDeg','isLowSpeedPivot','eventThreshold'});
results.trial=string(results.trial);eventsTable.trial=string(eventsTable.trial);
writetable(results,fullfile(out,'peak_amplitude_summary.csv'));writetable(eventsTable,fullfile(out,'peak_amplitude_events.csv'));disp(results);
end

function [peaks,amps,smoothDff,threshold]=events(dff,fs)
smoothDff=smoothdata(dff,'gaussian',max(3,round(.25*fs)));baseline=movmedian(smoothDff,max(3,round(30*fs)),'omitmissing');
signal=smoothDff-baseline;noise=mad(diff(signal),1)/(sqrt(2)*.6745);threshold=max(.5,3*noise);
active=signal>=threshold;d=diff([false;active(:);false]);s=find(d==1);e=find(d==-1)-1;
if isempty(s),peaks=[];amps=[];return,end
gap=round(.5*fs);ms=s(1);me=[];ce=e(1);for i=2:numel(s),if s(i)-ce-1<=gap,ce=e(i);else,me(end+1,1)=ce;ms(end+1,1)=s(i);ce=e(i);end,end %#ok<AGROW>
me(end+1,1)=ce;keep=(me-ms+1)>=max(1,round(.25*fs));ms=ms(keep);me=me(keep);peaks=zeros(numel(ms),1);amps=peaks;
for i=1:numel(ms),[amps(i),j]=max(signal(ms(i):me(i)));peaks(i)=ms(i)+j-1;end
end

function draw_groups(ax,amp,pivot)
pivot=logical(pivot);hold(ax,'on');x=[ones(sum(~pivot),1);2*ones(sum(pivot),1)];y=[amp(~pivot);amp(pivot)];
scatter(ax,x+.12*(rand(size(x))-.5),y,25,[.55 .55 .55],'filled','MarkerFaceAlpha',.45);
plot(ax,[.78 1.22],repmat(median(amp(~pivot),'omitnan'),1,2),'k','LineWidth',3);
plot(ax,[1.78 2.22],repmat(median(amp(pivot),'omitnan'),1,2),'Color',[.8 .1 .1],'LineWidth',3);
xticks(ax,[1 2]);xticklabels(ax,{sprintf('Other (n=%d)',sum(~pivot)),sprintf('Pivot (n=%d)',sum(pivot))});
ylabel(ax,'Event amplitude above local baseline (dF/F)');grid(ax,'on');xlim(ax,[.5 2.5]);
end
