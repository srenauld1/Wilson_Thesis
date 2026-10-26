function [trials,events]=review_ves041_fly1_20261020()
% Descriptive within-fly review. No dF/F amplitude exclusion. Event changes
% compare 0:+1 s with -2:-1 s, keeping the immediately preceding second free.
% Heat is measured command voltage >0.05 V, not measured temperature.
root='/Users/sophiarenauld/stacks/2026-10-20';nums=[7 8 9 11 12 13];rows=cell(0,12);erows=cell(0,9);
for num=nums
    folder=dir(fullfile(root,sprintf('20261020_%d_*',num)));folder=fullfile(folder(1).folder,folder(1).name);
    ff=dir(fullfile(folder,'*_s_.mat'));a=load(fullfile(ff(1).folder,ff(1).name),'dq','roi');a=import_ves041_protocol(a,folder);
    t=a.dq(1).t(:);d=a.roi.ts(1,:)';f=a.dq(1).bvf(:);r=abs(rad2deg(a.dq(1).bvy(:)));heat=a.dq(1).heat(:);on=heat>.05;fs=1/median(diff(t));
    valid=isfinite(d)&isfinite(f)&isfinite(r)&isfinite(heat);
    baseline=movmedian(d,max(3,round(30*fs)),'omitmissing');noise=1.4826*mad(d-baseline,1);
    [~,pk]=findpeaks(d,'MinPeakProminence',3*noise,'MinPeakDistance',round(2*fs),'MinPeakWidth',max(1,round(.25*fs)));
    z=d-baseline;
    rows(end+1,:)={num,string(a.protocol.kind),t(end)-t(1),100*mean(f(valid)>=.5),100*mean(on(valid)), ...
        mean(d(valid&on),'omitnan'),mean(d(valid&~on),'omitnan'),mean(f(valid&on),'omitnan'),mean(f(valid&~on),'omitnan'), ...
        numel(pk),corr(z(valid),f(valid),'Rows','complete'),corr(z(valid),r(valid),'Rows','complete')}; %#ok<AGROW>
    if num==11
        for bounds=[0 390;400 600]'
            ix=valid&t>=bounds(1)&t<bounds(2);fprintf('Trial11 %.0f:%.0f median dff %.3f heat%% %.1f\n',bounds,median(d(ix)),100*mean(on(ix)));
        end
    end
    if strcmp(a.protocol.kind,'bar-jump')
        c=a.dq(1).jump_count(:);idx=find([false;diff(c)>0]);labels=repmat("bar jump",numel(idx),1);
    elseif strcmp(a.protocol.kind,'dark')
        idx=pk(:);labels=repmat("dff peak",numel(idx),1);
    else
        E=readtable(fullfile(folder,'smr_analysis','heat_geometry_events.csv'),'TextType','string');
        idx=arrayfun(@(v)find(abs(t-v)==min(abs(t-v)),1),E.time);labels=E.event;
    end
    for k=1:numel(idx)
        p=idx(k);pre=t>=t(p)-2&t<t(p)-1;post=t>=t(p)&t<t(p)+1;
        if t(p)<t(1)+2||t(p)>t(end)-1,continue,end
        erows(end+1,:)={num,labels(k),t(p),mean(d(pre),'omitnan'),mean(d(post),'omitnan'), ...
            mean(d(post),'omitnan')-mean(d(pre),'omitnan'),mean(f(pre),'omitnan'),mean(f(post),'omitnan'),mean(r(post),'omitnan')}; %#ok<AGROW>
    end
end
trials=cell2table(rows,'VariableNames',{'trial','context','durationS','movingPercent','heatPercent','meanDffHeat','meanDffNoHeat','meanFwdHeat','meanFwdNoHeat','nPeaks','detrendedDffFwdR','detrendedDffRotSpeedR'});
events=cell2table(erows,'VariableNames',{'trial','event','time','preDff','postDff','dffChange','preFwd','postFwd','postRotSpeed'});
out='/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis';writetable(trials,fullfile(out,'fly1_descriptive_review.csv'));writetable(events,fullfile(out,'fly1_event_review.csv'));disp(trials);
for num=nums
    E=events(events.trial==num,:);groups=unique(E.event);
    for g=1:numel(groups),ix=E.event==groups(g);fprintf('T%d %s n=%d dffChange=%.3f median=%.3f positive=%d preFwd=%.2f postFwd=%.2f rot=%.1f\n',num,groups(g),sum(ix),mean(E.dffChange(ix),'omitnan'),median(E.dffChange(ix),'omitnan'),sum(E.dffChange(ix)>0),mean(E.preFwd(ix),'omitnan'),mean(E.postFwd(ix),'omitnan'),mean(E.postRotSpeed(ix),'omitnan'));end
end
end
