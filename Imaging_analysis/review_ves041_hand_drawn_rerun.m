function review_ves041_hand_drawn_rerun()
% Final comparisons to preserved segmented-ROI outputs; no new thresholds.
root='/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis';
out=fullfile(root,'hand_drawn_15s_ave_f_20261026');
R=readtable(fullfile(out,'run_report.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');
assert(height(R)==16&&all(R.status=="success"),'Incomplete manual ROI rerun');
E=readtable(fullfile(out,'heat_review','all_slowdown_events.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');
P=readtable(fullfile(out,'all_detected_peaks.csv'),'Delimiter',',','ReadVariableNames',true,'TextType','string');
ids=["20261019_5_1","20261019_7_1","20261019_8_2", ...
    "20261020_11_2","20261020_12_1","20261020_13_1"];
rows=cell(0,9);
for id=ids
    for turn=0:1
        for context=["heat","control","no heat"]
            e=E(E.recording==id&E.reorientation==turn&E.context==context,:);
            p=P.peakTime(P.recording==id&logical(P.eligibleAlignedWindow)&P.context==context);
            hit=false(height(e),1);post=hit;
            for j=1:height(e),dt=p-e.eventTime(j);hit(j)=any(abs(dt)<=1);post(j)=any(dt>=0&dt<=1);end
            rows(end+1,:)={id,turn,context,height(e),sum(hit),100*mean(hit),sum(post),100*mean(post),mean(e.deltaDff,'omitnan')}; %#ok<AGROW>
        end
    end
end
S=cell2table(rows,'VariableNames',{'recording','reorientation','context','nEvents','nPeakWithin1s','percentPeakWithin1s','nPeakFollowing1s','percentPeakFollowing1s','meanDeltaDff'});
writetable(S,fullfile(out,'event_metrics.csv'));plot_ves041_rectangular_heat_summary(out);
rows=cell(0,10);
for file=["matched_pairs.csv","control_inclusive_pairs.csv","trial11_same_epoch_pairs.csv"]
    prior=readtable(fullfile(root,'15s_ave_f','fly1_matched_heat',file),'Delimiter',',','ReadVariableNames',true,'TextType','string');
    current=readtable(fullfile(out,'fly1_matched_heat',file),'Delimiter',',','ReadVariableNames',true,'TextType','string');
    keys={'recording','reorientation','heatTime','noHeatTime'};
    if ismember('includeControls',prior.Properties.VariableNames),keys{end+1}='includeControls';end
    prior=sortrows(prior,keys);current=sortrows(current,keys);
    assert(isequaln(prior(:,keys),current(:,keys)),'Behavior-matched pairs changed');
    for id=unique(current.recording)'
        for turn=0:1
            q=current.recording==id&current.reorientation==turn;
            if ~any(q),continue,end
            rows(end+1,:)={file,id,turn,sum(q),mean(prior.heatDeltaDff(q)),mean(current.heatDeltaDff(q)), ...
                mean(prior.noHeatDeltaDff(q)),mean(current.noHeatDeltaDff(q)), ...
                mean(prior.heatMinusNoHeat(q)),mean(current.heatMinusNoHeat(q))}; %#ok<AGROW>
        end
    end
end
C=cell2table(rows,'VariableNames',{'comparison','recording','reorientation','nPairs','priorHeatDelta','manualHeatDelta','priorOffDelta','manualOffDelta','priorPairedDifference','manualPairedDifference'});
writetable(C,fullfile(out,'matched_roi_comparison.csv'));
disp(C);disp(S(S.reorientation==1,:));
end
