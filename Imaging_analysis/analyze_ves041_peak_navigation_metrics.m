%% Navigation metrics surrounding SNR-adjusted VES041 peaks
% Compares path straightness, net heading displacement, and heading
% consistency (rho) around VES041 peaks with circularly shifted controls.
% Includes centered sensitivity windows, a -2 to -1 s pre-peak control,
% and a 0 to +1 s post-peak test window. The intervening -1 to 0 s is
% deliberately excluded from the direct temporal comparison.

clearvars
rng(544);
baseFolder='/Users/sophiarenauld/stacks';
outputFolder=fullfile(baseFolder,'ves041_behavior_quantification','peak_navigation_metrics');
if ~isfolder(outputFolder),mkdir(outputFolder);end

% Keep a one-second gap between the pre-peak control and the peak itself so
% behavior immediately surrounding the fluorescence peak does not leak into
% the baseline comparison.
windowDefinitions=[-0.5 0.5; -1 1; -2 2; -2 -1; 0 1];
windowLabels=["centered_1s";"centered_2s";"centered_4s"; ...
    "pre_control_1s";"post_peak_1s"];
nPerm=2000;
minPathLength=0.5; % mm required for a meaningful straightness ratio
entries=dir(fullfile(baseFolder,'*ves041*')); entries=entries([entries.isdir]);
rows=cell(0,20);

for q=1:numel(entries)
    folder=fullfile(entries(q).folder,entries(q).name);
    [status,tags]=system(sprintf('/usr/bin/mdls -raw -name kMDItemUserTags "%s"',folder));
    if status~=0 || ~contains(tags,'good behavior'),continue,end
    flyFile=fullfile(folder,'smr_analysis','fly_data.mat');
    if ~isfile(flyFile),continue,end
    S=load(flyFile,'fly_data'); x=S.fly_data;
    t=x.dq(1).t(:); dff=x.roi.ts(:); xpos=x.dq(1).px(:); ypos=x.dq(1).py(:);
    heading=x.dq(1).bh(:); % circular fly-heading signal in radians
    n=min([numel(t),numel(dff),numel(xpos),numel(ypos),numel(heading)]);
    t=t(1:n);dff=dff(1:n);xpos=xpos(1:n);ypos=ypos(1:n);heading=heading(1:n);
    fs=1/(max(t)/length(t));

    % Same retained detector: 3 sigma prominence, >=0.25 s width, and 2 s
    % minimum distance. No upper dF/F cutoff is applied.
    baseline=movmedian(dff,max(3,round(30*fs)),'omitmissing');
    noiseSigma=1.4826*mad(dff-baseline,1);
    [~,peaks]=findpeaks(dff,'MinPeakProminence',3*noiseSigma, ...
        'MinPeakDistance',round(2*fs),'MinPeakWidth',round(0.25*fs));

    edgeN=round(max(abs(windowDefinitions),[],'all')*fs);
    peaks=peaks(peaks>edgeN & peaks+edgeN<=n);
    if isempty(peaks),continue,end
    allowed=n-2*edgeN;

    for wi=1:size(windowDefinitions,1)
        windowSec=windowDefinitions(wi,:);
        startN=round(windowSec(1)*fs); endN=round(windowSec(2)*fs);
        observed=event_metrics(peaks,xpos,ypos,heading,startN,endN,minPathLength);
        obsSummary=summarize(observed);
        null=nan(nPerm,4);
        for p=1:nPerm
            shift=randi([round(5*fs),allowed-round(5*fs)]);
            shifted=edgeN+mod((peaks-edgeN-1)+shift,allowed)+1;
            null(p,:)=summarize(event_metrics(shifted,xpos,ypos,heading,startN,endN,minPathLength));
        end

        % Expected directions for reorientation: lower straightness and rho,
        % but greater net heading displacement than shifted controls.
        pStraight=(1+sum(null(:,1)<=obsSummary(1)))/(nPerm+1);
        pHeading=(1+sum(null(:,2)>=obsSummary(2)))/(nPerm+1);
        pRho=(1+sum(null(:,3)<=obsSummary(3)))/(nPerm+1);
        [~,name]=fileparts(folder);
        rows(end+1,:)={name,infer_context(lower(name)),windowLabels(wi), ...
            windowSec(1),windowSec(2),numel(peaks), ... %#ok<SAGROW>
            obsSummary(1),mean(null(:,1),'omitnan'),pStraight, ...
            obsSummary(2),mean(null(:,2),'omitnan'),pHeading, ...
            obsSummary(3),mean(null(:,3),'omitnan'),pRho, ...
            obsSummary(4),mean(null(:,4),'omitnan'),noiseSigma,3*noiseSigma,minPathLength};
    end
end

results=cell2table(rows,'VariableNames',{'trial','context','windowLabel', ...
    'windowStartSec','windowEndSec','nPeaks', ...
    'medianStraightness','nullStraightness','pLowerStraightness', ...
    'meanNetHeadingDeg','nullNetHeadingDeg','pHigherNetHeading', ...
    'meanRho','nullRho','pLowerRho','fractionStraightnessValid', ...
    'nullFractionStraightnessValid','noiseSigma','prominenceThreshold','minPathLength'});
results.trial=string(results.trial);results.context=string(results.context);
results.windowLabel=string(results.windowLabel);
writetable(results,fullfile(outputFolder,'peak_navigation_metrics_across_trials.csv'));
disp(results(results.windowLabel=="post_peak_1s",:));

function m=event_metrics(indices,x,y,h,startN,endN,minPath)
n=numel(indices); m=nan(n,4);
for k=1:n
    win=indices(k)+startN:indices(k)+endN;
    dx=diff(x(win));dy=diff(y(win));pathLength=sum(hypot(dx,dy),'omitnan');
    displacement=hypot(x(win(end))-x(win(1)),y(win(end))-y(win(1)));
    if pathLength>=minPath,m(k,1)=displacement/pathLength;end
    m(k,2)=abs(rad2deg(angle(exp(1i*(h(win(end))-h(win(1)))))));
    m(k,3)=abs(mean(exp(1i*h(win)),'omitnan'));
    m(k,4)=pathLength>=minPath;
end
end

function s=summarize(m)
s=[median(m(:,1),'omitnan'),mean(m(:,2),'omitnan'), ...
   mean(m(:,3),'omitnan'),mean(m(:,4),'omitnan')];
end

function c=infer_context(name)
if contains(name,'odor'),c="odor";
elseif contains(name,'grating'),c="grating";
elseif contains(name,'jump'),c="bar jump";
elseif contains(name,'dark'),c="dark";
else,c="closed loop";
end
end
