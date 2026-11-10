function result=plot_pizza_bar_low_speed_rotation(out,destination)
%PLOT_PIZZA_BAR_LOW_SPEED_ROTATION Expand bar-only pooled crunches to 8 groups.
% Uses the exact established bar-only trial-to-fly map (no dark recordings).
% Columns 1:4 and 7:8 are copied from saved eight-group event arrays. Only
% columns 5:6 are filtered: retain rotation peaks with median abs(forward)
% <2 mm/s over [-0.25,+0.25] s about the absolute-yaw peak. Forward velocity
% is first median-smoothed over 0.25 s within contiguous valid data blocks,
% matching the existing rotation detector's preprocessing. This is a local
% median criterion, not a requirement that every frame be below 2 mm/s.
% Absolute forward velocity excludes fast reverse walking. Sideways speed is
% not part of this criterion. Fluorescence does not enter the low-speed gate.
% Retains the existing rotation split: pre-minus-post forward velocity drop
% <1 versus >=1 mm/s. Original detector/window eligibility stays unchanged.
% Pool events within each fly, then weight contributing flies equally.
% Z-score each metric using its EXISTING full-recording bar-only fly mean/SD,
% not the restricted event set. No source files or older figures overwritten.
% Run MATLAB -nodisplay; save visible-on-open FIG and PNG companions.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
if nargin<2,destination=fullfile(out,'bar','rotation_low_fwd_lt2_20261102');end
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
if ~isfolder(destination),mkdir(destination);end
stem=fullfile(out,'bar','ves041_bar_pooled_equal_fly_peaks_slowdowns');
P=load([stem '.mat']);N=load([stem '_normalized.mat']);
assert(isequal(P.flies,N.flies),'Normalization/cohort mismatch');
F=P.F;flies=P.flies;nf=numel(flies);rel=P.rel;pooled=cell(nf,8);
counts=zeros(nf,8);originalCounts=counts;eventRows=cell(0,8);
for k=1:height(F)
    trial=F.trial(k);f=find(flies==F.flyId(k));
    S=load(fullfile(out,'analysis_data',trial+"_events_with_rotation_peaks_and_speedups.mat"),'X','groups','rel');
    A=load(fullfile(out,'analysis_data',trial+"_events.mat"),'d');d=A.d;
    assert(numel(S.X)==8&&isequal(rel,S.rel));
    fs=1/median(diff(d.t));valid=d.valid&isfinite(d.fwd)&isfinite(d.yaw)&isfinite(d.dff);
    cut=[true;diff(d.t)>2/fs];s=find(valid&(cut|[true;~valid(1:end-1)]));
    e=find(valid&([~valid(2:end);true]|[cut(2:end);true]));filtered=nan(size(d.fwd));
    for b=1:numel(s),ix=s(b):e(b);filtered(ix)=smoothdata(d.fwd(ix),'movmedian',max(3,round(.25*fs)));end
    for c=1:8
        x=S.X{c};idx=S.groups{c};assert(size(x,2)==numel(idx));
        originalCounts(f,c)=originalCounts(f,c)+numel(idx);keep=true(numel(idx),1);
        if ismember(c,[5 6])
            for j=1:numel(idx)
                p=idx(j);window=d.t>=d.t(p)-.25&d.t<=d.t(p)+.25;
                assert(any(window)&&all(isfinite(filtered(window))),'Invalid speed gate window');
                localSpeed=median(abs(filtered(window)));keep(j)=localSpeed<2;
                eventRows(end+1,:)={trial,F.flyId(k),c,p,d.t(p),localSpeed,keep(j),numel(find(window))}; %#ok<AGROW>
            end
        end
        pooled{f,c}=cat(2,pooled{f,c},x(:,keep,:));counts(f,c)=counts(f,c)+nnz(keep);
    end
end
means=nan(numel(rel),nf,3,8);
for f=1:nf,for c=1:8,if counts(f,c)>0,means(:,f,:,c)=mean(pooled{f,c},2);end,end,end
% Confirm original four selections, their event counts, and their fly means
% match the existing bar-only output, including absent-event NaNs.
assert(isequal(counts(:,1:4),P.counts));
assert(isequal(isnan(means(:,:,:,1:4)),isnan(P.flyMeans)));
delta=abs(means(:,:,:,1:4)-P.flyMeans);assert(all(delta(isfinite(delta))<1e-10));
assert(isequal(counts(:,[1:4 7:8]),originalCounts(:,[1:4 7:8])));
assert(all(counts(:,5:6)<=originalCounts(:,5:6),'all'));
normalized=means;
for f=1:nf,for v=1:3,normalized(:,f,v,:)=(means(:,f,v,:)-N.center(f,v))/N.scale(f,v);end,end
delta=abs(normalized(:,:,:,1:4)-N.flyMeans);assert(all(delta(isfinite(delta))<1e-10));
labels={{'VES041 peaks'},{'All slowdowns'},{'Slowdowns without turns'},{'Slowdowns with turns'}, ...
    {'Low-speed rotation peaks','Forward drop <1 mm/s'}, ...
    {'Low-speed rotation peaks','Forward drop >=1 mm/s'}, ...
    {'Speed-ups','Without rotation (<15 deg)'},{'Speed-ups','With rotation (>=15 deg)'}};
titleText='VES041 bar-only trials | rotation peaks restricted to slow forward movement';
subtitleText={sprintf('%d recordings | %d estimated flies (%d inferred) | Gray: fly means; black: equal-fly mean',height(F),nf,nnz(contains(lower(flies),'inferred'))), ...
    'Only columns 5/6 restricted: median |forward| <2 mm/s over +/-0.25 s | Other selections and full-recording z-score references unchanged'};
newStem=fullfile(destination,'ves041_bar_all_selections_rotation_low_fwd_lt2');
for mode=1:2
    X=cell(1,8);for c=1:8,if mode==1,X{c}=means(:,:,:,c);else,X{c}=normalized(:,:,:,c);end,end
    target=newStem;if mode==2,target=[newStem '_normalized'];end
    pizza_ves041_render_eight_columns(X,rel,titleText,subtitleText,target,mode==2,flies,sum(counts,1),labels);
end
rotationEvents=cell2table(eventRows,'VariableNames',{'trial','flyId','column','sampleIndex', ...
    'eventTimeSec','medianAbsoluteForwardMmS','retained','gateSamples'});
writetable(rotationEvents,fullfile(destination,'rotation_event_filter.csv'));
summaryRows=cell(0,4);
for c=1:8,summaryRows(end+1,:)={c,strjoin(string(labels{c}),' | '),sum(originalCounts(:,c)),sum(counts(:,c))};end %#ok<AGROW>
summary=cell2table(summaryRows,'VariableNames',{'column','selection','originalEvents','retainedEvents'});
writetable(summary,fullfile(destination,'selection_counts.csv'));
weightRows=cell(0,5);contributingFlies=sum(counts>0,1);
for c=1:8,for f=1:nf
    weight=0;if contributingFlies(c)>0,weight=double(counts(f,c)>0)/contributingFlies(c);end
    weightRows(end+1,:)={flies(f),c,counts(f,c),counts(f,c)>0,weight}; %#ok<AGROW>
end,end
weights=cell2table(weightRows,'VariableNames',{'flyId','column','events','contributes','weight'});
writetable(weights,fullfile(destination,'fly_weights.csv'));
methods=struct('cohort','Existing bar-only mapping; no added dark trials','cutoffMmS',2, ...
    'gateWindowSec',[-.25 .25],'gateStatistic','Median absolute 0.25-s median-smoothed forward velocity', ...
    'onlyFilteredColumns',[5 6],'unchangedColumns',[1:4 7:8], ...
    'normalization','Existing full-recording within-fly bar-only mean and sample SD for each metric', ...
    'averaging','Pool events within fly, then average contributing flies equally');
center=N.center;scale=N.scale;
save(fullfile(destination,'analysis.mat'),'F','flies','rel','counts','originalCounts','means', ...
    'normalized','center','scale','weights','methods','rotationEvents','summary','-v7.3');
fid=fopen(fullfile(destination,'methods.txt'),'w');assert(fid>=0);
fprintf(fid,'Bar-only crunch: %d recordings, %d estimated flies.\n',height(F),nf);
fprintf(fid,'Only rotation-peak columns 5/6 use median absolute forward velocity <2 mm/s over -0.25:+0.25 s. Forward is first median-smoothed over 0.25 s, separately in valid contiguous blocks. This is a median gate, not an all-frames gate; lateral speed is not restricted.\n');
fprintf(fid,'Rotation detector, full -2:+3 s valid windows, and forward-drop <1 versus >=1 mm/s split are unchanged. Columns 1-4 and 7-8 are unchanged saved selections. No fluorescence-based gate is added.\n');
fprintf(fid,'Events pooled within each fly, then contributing fly means weighted equally. Normalized metrics retain the original full-recording bar-only within-fly mean/sample-SD references, not low-speed-event references. Three fly identities are inferred in the original map. Source files and older plots are unchanged.\n');
fclose(fid);
result=struct('destination',destination,'summary',summary,'contributingFlies',contributingFlies);disp(result);disp(summary);
end
