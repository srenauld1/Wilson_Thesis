%% Compare 0-1 s after VES041 peaks with trial-wide pseudo-event windows
% For each good-behavior trial, randomly scatter the same number of 1 s
% control windows throughout the recording. Repeating this 2000 times gives
% a stable descriptive trial-wide baseline without inferential testing.

clearvars
rng(544)
baseFolder='/Users/sophiarenauld/stacks';
outputFolder=fullfile(baseFolder,'ves041_behavior_quantification','peak_navigation_metrics');
if ~isfolder(outputFolder),mkdir(outputFolder);end
nDraws=2000; minPathLength=0.5;
entries=dir(fullfile(baseFolder,'*ves041*')); entries=entries([entries.isdir]);
rows=cell(0,13);

for q=1:numel(entries)
    folder=fullfile(entries(q).folder,entries(q).name);
    [status,tags]=system(sprintf('/usr/bin/mdls -raw -name kMDItemUserTags "%s"',folder));
    if status~=0 || ~contains(tags,'good behavior'),continue,end
    flyFile=fullfile(folder,'smr_analysis','fly_data.mat');
    if ~isfile(flyFile),continue,end

    S=load(flyFile,'fly_data'); x=S.fly_data;
    t=x.dq(1).t(:); dff=x.roi.ts(:); fwd=x.dq(1).bvf(:);
    xpos=x.dq(1).px(:); ypos=x.dq(1).py(:); heading=x.dq(1).bh(:);
    n=min([numel(t),numel(dff),numel(fwd),numel(xpos),numel(ypos),numel(heading)]);
    t=t(1:n); dff=dff(1:n); fwd=fwd(1:n); xpos=xpos(1:n);
    ypos=ypos(1:n); heading=heading(1:n);
    fs=1/(max(t)/length(t)); windowN=round(fs);

    baseline=movmedian(dff,max(3,round(30*fs)),'omitmissing');
    noiseSigma=1.4826*mad(dff-baseline,1);
    [~,peaks]=findpeaks(dff,'MinPeakProminence',3*noiseSigma, ...
        'MinPeakDistance',round(2*fs),'MinPeakWidth',round(0.25*fs));
    peaks=peaks(peaks+windowN<=n);
    if isempty(peaks),continue,end

    observed=summarize_windows(peaks,fwd,xpos,ypos,heading,windowN,minPathLength);
    pseudo=nan(nDraws,3);
    for draw=1:nDraws
        pseudoStarts=randi([1,n-windowN],numel(peaks),1);
        pseudo(draw,:)=summarize_windows(pseudoStarts,fwd,xpos,ypos,heading,windowN,minPathLength);
    end
    pseudoMean=mean(pseudo,1,'omitnan');
    [~,name]=fileparts(folder);
    rows(end+1,:)={name,infer_context(lower(name)),numel(peaks), ... %#ok<SAGROW>
        observed(1),pseudoMean(1),observed(1)-pseudoMean(1), ...
        observed(2),pseudoMean(2),observed(2)-pseudoMean(2), ...
        observed(3),pseudoMean(3),observed(3)-pseudoMean(3),noiseSigma};
end

results=cell2table(rows,'VariableNames',{'trial','context','nPeaks', ...
    'postFwd','pseudoFwd','effectFwd','postStraightness','pseudoStraightness', ...
    'effectStraightness','postNetHeadingDeg','pseudoNetHeadingDeg', ...
    'effectNetHeadingDeg','noiseSigma'});
results.trial=string(results.trial); results.context=string(results.context);
writetable(results,fullfile(outputFolder,'postpeak_vs_pseudo_1s_descriptive.csv'));
disp(results)

function out=summarize_windows(starts,fwd,x,y,h,w,minPath)
n=numel(starts); f=nan(n,1); s=nan(n,1); hd=nan(n,1);
for k=1:n
    win=starts(k):starts(k)+w;
    f(k)=mean(fwd(win),'omitnan');
    pathLength=sum(hypot(diff(x(win)),diff(y(win))),'omitnan');
    displacement=hypot(x(win(end))-x(win(1)),y(win(end))-y(win(1)));
    if pathLength>=minPath,s(k)=displacement/pathLength;end
    hd(k)=abs(rad2deg(angle(exp(1i*(h(win(end))-h(win(1)))))));
end
out=[mean(f,'omitnan'),median(s,'omitnan'),mean(hd,'omitnan')];
end

function c=infer_context(name)
if contains(name,'odor'),c="odor";
elseif contains(name,'grating'),c="grating";
elseif contains(name,'jump'),c="bar jump";
elseif contains(name,'dark'),c="dark";
else,c="closed loop";
end
end
