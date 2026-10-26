function [summary,events]=compare_ves041_15s_ave_f(files,out)
% Compare preserved fixed-baseline outputs with the 15 s baseline rerun.
% Read only old ROI trace directly from HDF5, not the saved movie stack.
% Same detector on both: 3 MAD sigma prominence, 2 s distance, .25 s width.
% Event response: mean 0:+1 s minus mean -2:-1 s. No hypothesis tests.
oldVisible=get(groot,'defaultFigureVisible');c=onCleanup(@()set(groot,'defaultFigureVisible',oldVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');rows=cell(0,22);erows=cell(0,7);
for k=1:numel(files)
    folder=files(k).folder;[~,trial]=fileparts(folder);id=erase(files(k).name,'_or_z0z1_s_.mat');
    previous=fullfile(folder,'smr_analysis','fly_data.mat');current=fullfile(folder,'15s_ave_f','fly_data.mat');
    if ~isfile(current),warning('Missing rerun: %s',current);continue,end
    a=load(current,'fly_data');a=a.fly_data;t=a.dq(1).t(:);f=a.dq(1).bvf(:);r=abs(rad2deg(a.dq(1).bvy(:)));heat=a.dq(1).heat(:)>.05;
    old=double(h5read(previous,'/fly_data/roi/ts'));old=old(:);new=double(a.roi.ts(:));assert(numel(old)==numel(new));
    oldMask=logical(h5read(previous,'/fly_data/roi/mask'));sameMask=isequal(oldMask,a.roi.mask);
    oldRaw=double(h5read(previous,'/fly_data/roi/pxfilter/filteredtsraw'));newRaw=double(a.roi.pxfilter.filteredtsraw);
    rawMaxDifference=max(abs(oldRaw(:)-newRaw(:)));
    fs=1/median(diff(t));D={old,new};locs=cell(2,1);noise=zeros(2,1);nPeaks=noise;preF=noise;postF=noise;postR=noise;lowF=noise;rhoF=noise;rhoR=noise;
    for v=1:2
        d=D{v};z=d-movmedian(d,max(3,round(30*fs)),'omitmissing');noise(v)=1.4826*mad(z,1);
        [~,pk]=findpeaks(d,'MinPeakProminence',3*noise(v),'MinPeakDistance',round(2*fs),'MinPeakWidth',round(.25*fs));
        pk=pk(pk>round(2*fs)&pk<=numel(t)-round(2*fs));locs{v}=pk;nPeaks(v)=numel(pk);
        b=nan(numel(pk),4);
        for j=1:numel(pk)
            pre=t>=t(pk(j))-2&t<t(pk(j))-1;post=t>=t(pk(j))&t<t(pk(j))+1;
            b(j,:)=[mean(f(pre),'omitnan'),mean(f(post),'omitnan'),mean(r(post),'omitnan'),mean(f(post),'omitnan')<=5];
        end
        m=mean(b,1,'omitnan');preF(v)=m(1);postF(v)=m(2);postR(v)=m(3);lowF(v)=m(4)*100;
        rhoF(v)=corr(d,f,'Rows','complete');rhoR(v)=corr(d,r,'Rows','complete');
    end
    matched=sum(arrayfun(@(p)any(abs(t(locs{2})-t(p))<=.5),locs{1}));
    rows(end+1,:)={string(id),string(a.protocol.kind),nPeaks(1),nPeaks(2),matched,noise(1),noise(2),preF(1),postF(1),preF(2),postF(2), ...
        lowF(1),lowF(2),rhoF(1),rhoF(2),rhoR(1),rhoR(2),corr(old,new,'Rows','complete'),sameMask,rawMaxDifference,postR(1),postR(2)}; %#ok<AGROW>
    % Use behavior/stimulus-selected events, so the event set cannot change
    % as a consequence of changing the fluorescence normalization.
    ef=fullfile(folder,'smr_analysis','heat_geometry_events.csv');
    if isfile(ef)
        E=readtable(ef,'TextType','string');et=E.time;labels=E.event;
        fresh=readtable(fullfile(folder,'15s_ave_f','heat_geometry_events.csv'),'TextType','string');
        assert(isequaln(E.time,fresh.time)&&isequaln(E.event,fresh.event),'Heat-event selection changed.');
    elseif strcmp(a.protocol.kind,'bar-jump'),et=t(logical(a.jump_detected(:)));labels=repmat("bar jump",numel(et),1);
    else,et=[];labels=strings(0,1);end
    for j=1:numel(et)
        if et(j)<t(1)+2||et(j)>t(end)-1,continue,end
        pre=t>=et(j)-2&t<et(j)-1;post=t>=et(j)&t<et(j)+1;
        erows(end+1,:)={string(id),labels(j),et(j),mean(old(post))-mean(old(pre)),mean(new(post))-mean(new(pre)),mean(f(post)),mean(r(post))}; %#ok<AGROW>
    end
    fig=figure('Visible','off','Position',[50 50 1500 900]);tl=tiledlayout(4,1,'TileSpacing','compact');
    for v=1:2,nexttile;plot(t,D{v},'b');hold on;scatter(t(locs{v}),D{v}(locs{v}),12,'r','filled');ylabel('dF/F');if v==1,title('Previous: fixed baseline');else,title('New: 15 s fifth-percentile baseline');end;grid on;end
    nexttile;plot(t,f,'k');ylabel('Forward mm/s');grid on;
    nexttile;plot(t,r,'r');hold on;yyaxis right;plot(t,heat,'Color',[.8 .5 0]);ylabel('Heat command >0.05 V');yyaxis left;ylabel('|Rotation| deg/s');xlabel('s');grid on;
    axs=findall(fig,'Type','axes');linkaxes(axs,'x');xlim(axs(1),[t(1) t(end)]);
    title(tl,strrep([id ' | ' trial],'_','\_'));savefig_visible_on_open(fig,fullfile(out,[id '_old_vs_15s.fig']));exportgraphics(fig,fullfile(out,[id '_old_vs_15s.png']),'Resolution',130);close(fig);
    if contains(id,'20261020_11_2')
        pre=t>=375&t<385;post=t>=405&t<415;
        fprintf('Trial11 local step: OLD %.3f -> %.3f; NEW %.3f -> %.3f\n',median(old(pre)),median(old(post)),median(new(pre)),median(new(post)));
    end
end
summary=cell2table(rows,'VariableNames',{'trial','context','oldPeaks','newPeaks','oldPeaksMatchedWithinHalfSecond','oldNoiseSigma','newNoiseSigma','oldPreFwd','oldPostFwd','newPreFwd','newPostFwd','oldPeakLowFwdPercent','newPeakLowFwdPercent','oldDffFwdR','newDffFwdR','oldDffRotSpeedR','newDffRotSpeedR','oldNewDffR','sameRoiMask','rawFluorescenceMaxDifference','oldPeakPostRotSpeed','newPeakPostRotSpeed'});
events=cell2table(erows,'VariableNames',{'trial','event','time','oldDffChange','newDffChange','postFwd','postRotSpeed'});
writetable(summary,fullfile(out,'normalization_comparison.csv'));writetable(events,fullfile(out,'stimulus_event_comparison.csv'));disp(summary);
end
