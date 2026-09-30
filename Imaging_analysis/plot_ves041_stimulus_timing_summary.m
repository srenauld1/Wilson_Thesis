%% Summary figure: VES041 around bar jumps and optomotor motion
% Reads the descriptive outputs from analyze_ves041_stimulus_timing.m and
% saves an editable MATLAB figure plus a high-resolution PNG.

clearvars
baseFolder='/Users/sophiarenauld/stacks/ves041_behavior_quantification/stimulus_timing';
bar=readtable(fullfile(baseFolder,'bar_jump_nearby_peak_events.csv'),'TextType','string');
opto=readtable(fullfile(baseFolder,'optomotor_moving_vs_stationary.csv'),'TextType','string');

trialColors=lines(3);
fig=figure('Color','w','Position',[80 80 1500 820]);
tl=tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

%% A: timing of VES041 peaks relative to jumps
ax=nexttile(tl,1); hold(ax,'on');
barTrials=unique(bar.trial,'stable');
for i=1:numel(barTrials)
    z=bar(bar.trial==barTrials(i),:);
    scatter(ax,z.peakLagFromJump,repmat(i,height(z),1),80,trialColors(i,:), ...
        'filled','MarkerEdgeColor','k','DisplayName',short_name(barTrials(i)));
end
xline(ax,0,'k--','Bar jump','LabelVerticalAlignment','bottom');
xlim(ax,[-2.1 2.1]); ylim(ax,[0.5 numel(barTrials)+0.5]);
yticks(ax,1:numel(barTrials));yticklabels(ax,short_name(barTrials));
xlabel(ax,'Peak time relative to jump (s)');ylabel(ax,'Trial');
title(ax,'A  Peaks near bar jumps');grid(ax,'on');

%% B: forward velocity before and after nearby peaks
ax=nexttile(tl,2); hold(ax,'on');
for i=1:height(bar)
    c=trialColors(find(barTrials==bar.trial(i),1),:);
    plot(ax,[1 2],[bar.preFwd(i) bar.postFwd(i)],'-o','Color',c, ...
        'LineWidth',1.4,'MarkerFaceColor',c);
end
xlim(ax,[0.7 2.3]);xticks(ax,[1 2]);xticklabels(ax,{'-1 to 0 s','0 to +1 s'});
ylabel(ax,'Forward velocity (mm/s)');
title(ax,'B  Forward velocity around peak');grid(ax,'on');

%% C: reorientation following bar-associated peaks
ax=nexttile(tl,3); hold(ax,'on');
for i=1:numel(barTrials)
    z=bar(bar.trial==barTrials(i),:);
    scatter(ax,z.postAbsRotSpeed,z.postNetHeadingDeg,90,trialColors(i,:), ...
        'filled','MarkerEdgeColor','k','DisplayName',short_name(barTrials(i)));
end
xlabel(ax,'Mean rotational speed, 0 to +1 s (deg/s)');
ylabel(ax,'Net heading change, 0 to +1 s (deg)');
title(ax,'C  Post-peak reorientation');grid(ax,'on');legend(ax,'Location','best');

%% Prepare paired optomotor rows
optoTrials=unique(opto.trial,'stable');
stationary=opto(opto.stimulusState=="stationary",:);
moving=opto(opto.stimulusState=="moving",:);
[~,is]=ismember(optoTrials,stationary.trial);[~,im]=ismember(optoTrials,moving.trial);
stationary=stationary(is,:);moving=moving(im,:);

%% D: peak rate by panel-motion state
ax=nexttile(tl,4); hold(ax,'on');
paired_plot(ax,stationary.peakRatePerMin,moving.peakRatePerMin,trialColors);
ylabel(ax,'VES041 peaks/min');title(ax,'D  Peak frequency');grid(ax,'on');

%% E: baseline activity and peak amplitude
ax=nexttile(tl,5); hold(ax,'on');
x=[1 2 4 5];
for i=1:numel(optoTrials)
    c=trialColors(i,:);
    plot(ax,x(1:2),[stationary.medianDff(i),moving.medianDff(i)],'-o', ...
        'Color',c,'LineWidth',1.4,'MarkerFaceColor',c);
    plot(ax,x(3:4),[stationary.medianPeakDff(i),moving.medianPeakDff(i)],'-o', ...
        'Color',c,'LineWidth',1.4,'MarkerFaceColor',c);
end
xticks(ax,x);xticklabels(ax,{'All: still','All: moving','Peaks: still','Peaks: moving'});
xline(ax,3,'Color',[.65 .65 .65]);
ylabel(ax,'Median dF/F');title(ax,'E  Activity by stimulus state');grid(ax,'on');

%% F: heading outcome after peaks
ax=nexttile(tl,6); hold(ax,'on');
paired_plot(ax,stationary.meanPostNetHeadingDeg,moving.meanPostNetHeadingDeg,trialColors);
ylabel(ax,'Net heading change, 0 to +1 s (deg)');
title(ax,'F  Post-peak heading change');grid(ax,'on');
legend(ax,arrayfun(@short_name,optoTrials,'UniformOutput',false), ...
    'Location','best','Interpreter','none');

sgtitle(tl,'VES041 activity during bar jumps and optomotor motion', ...
    'FontWeight','bold');

figPath=fullfile(baseFolder,'ves041_stimulus_timing_summary.fig');
pngPath=fullfile(baseFolder,'ves041_stimulus_timing_summary.png');
savefig(fig,figPath);
exportgraphics(fig,pngPath,'Resolution',250);
fprintf('Saved %s\nSaved %s\n',figPath,pngPath);

function paired_plot(ax,a,b,colors)
for i=1:numel(a)
    plot(ax,[1 2],[a(i) b(i)],'-o','Color',colors(i,:), ...
        'LineWidth',1.5,'MarkerSize',6,'MarkerFaceColor',colors(i,:));
end
xlim(ax,[0.7 2.3]);xticks(ax,[1 2]);xticklabels(ax,{'Stationary','Moving'});
end

function s=short_name(name)
s=extractBefore(name+"_","_ves041");
end
