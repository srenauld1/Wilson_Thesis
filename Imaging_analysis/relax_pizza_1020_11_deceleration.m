function result=relax_pizza_1020_11_deceleration()
% Remove minimum deceleration magnitude only from the current two-selection plot.
root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_heat_analysis/barrier/20261020_11_2';
previous=fullfile(root,'presentation_control_vs_heat_20261102');out=fullfile(previous,'no_deceleration_cutoff');if ~isfolder(out),mkdir(out);end
S=load(fullfile(previous,'comparison_data.mat'));B=load(fullfile(root,'analysis.mat'));A=load(fullfile(root,'control_and_heat_bins','binned_analysis.mat'),'P');d=B.d;P=A.P;t=d.t;fs=1/median(diff(t));valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
cut=[true;diff(t)>2/fs];beg=find(valid&(cut|[true;~valid(1:end-1)]));last=find(valid&([~valid(2:end);true]|[cut(2:end);true]));
f=nan(size(t));acc=f;
for k=1:numel(beg),ix=beg(k):last(k);f(ix)=smoothdata(d.fwd(ix),'movmedian',max(3,round(.25*fs)));if numel(ix)>2,acc(ix)=gradient(f(ix),t(ix));end,end
% Negative runs locate slowing episodes; there is no magnitude cutoff.
candidate=valid&acc<0;change=diff([false;candidate;false]);starts=find(change==1);ends=find(change==-1)-1;
idx=zeros(0,1);angles=idx;before=idx;drop=idx;lastAccepted=-inf;
for k=1:numel(starts)
 [~,j]=min(acc(starts(k):ends(k)));i=starts(k)+j-1;t0=t(i);
 if t0-lastAccepted<2,continue,end
 lo=find(t<=t0-2,1,'last');hi=find(t>=t0+3,1,'first');if isempty(lo)||isempty(hi),continue,end
 win=lo:hi;if ~all(valid(win))||any(diff(t(win))>2/fs),continue,end
 pre=t>=t0-.75&t<t0-.25;post=t>=t0+.25&t<=t0+.75;a=median(f(pre));b=a-median(f(post));
 if a<1||b<1,continue,end
 idx(end+1,1)=i;before(end+1,1)=a;drop(end+1,1)=b;lastAccepted=t0;
end
angles=pizza_ves041_turn_angle(d,idx,[0 .5]);withTurn=angles>=15;
conditions=strings(numel(idx),1);conditions(:)="Excluded";reason=conditions;reason(:)="";voltage=nan(numel(idx),1);
for k=1:numel(idx)
 i=idx(k);lo=find(t<=t(i)-2,1,'last');hi=find(t>=t(i)+3,1,'first');win=lo:hi;core=t>=t(i)-.5&t<=t(i)+.5;voltage(k)=median(P.heat(core));
 if ~withTurn(k),reason(k)="rotation below 15 degrees";
 elseif ~all(P.active(win))||any(~isfinite(P.control(win)))||any(~isfinite(P.trial(win))),reason(k)="unplaced barrier or missing protocol";
 elseif numel(unique(P.trial(win)))~=1||numel(unique(P.control(win)))~=1||numel(unique(P.geometry(win)))~=1,reason(k)="window crosses trial/control/geometry boundary";
 elseif P.control(i)==1,conditions(k)="Control";
 elseif P.control(i)~=0,reason(k)="unknown control flag";
 elseif ~all(isfinite(P.heat(core))),reason(k)="missing measured heat";
 elseif voltage(k)<-.05,reason(k)="unexpected negative measured command";
 else,conditions(k)="Heat";end
end
newEvents=table(idx,t(idx),acc(idx),before,drop,angles,P.trial(idx),P.control(idx),voltage,conditions,reason, ...
 'VariableNames',{'sampleIndex','eventTime','decelerationMmS2','forwardBefore','forwardDrop','absoluteTurnDeg','barrierTrial','controlTrial','medianHeatV','condition','exclusionReason'});
oldCounts=S.counts;vars={d.dff,d.fwd,abs(d.yaw)};
for s=1:2
 name="Control";if s==2,name="Heat";end
 ix=idx(conditions==name);X=nan(numel(S.rel),numel(ix),3);
 for k=1:numel(ix),for v=1:3,X(:,k,v)=interp1(t,vars{v},t(ix(k))+S.rel);end,end
 S.D{s,2}=X;S.counts(s,2)=numel(ix);S.barrierCounts(s,2)=numel(unique(P.trial(ix)));S.mu(:,:,s,2)=squeeze(mean(X,2));S.sem(:,:,s,2)=squeeze(std(X,0,2))/sqrt(numel(ix));
 assert(all(isfinite(X(:))));
 q=conditions==name;assert(all(before(q)>=1&drop(q)>=1&angles(q)>=15));assert(all(P.control(ix)==(s==1)));
end
% Verify the fluorescence-peak arrays, means, SEM and counts are unchanged.
original=load(fullfile(previous,'comparison_data.mat'));
for s=1:2,assert(isequaln(S.D{s,1},original.D{s,1}));end
assert(isequaln(S.mu(:,:,:,1),original.mu(:,:,:,1))&&isequaln(S.sem(:,:,:,1),original.sem(:,:,:,1))&&isequal(S.counts(:,1),original.counts(:,1)));
S.decelerationCutoffRemoved=true;S.methods.slowdownDetector=rmfield(S.methods.slowdownDetector,'maximumAccelerationMmS2');S.methods.slowdownDetector.decelerationMagnitudeCutoff='none';S.methods.slowdownDetector.alignment='strongest deceleration within each negative-acceleration run';
S.methods.change='Remove the -2 mm/s^2 magnitude cutoff; retain speed-before>=1, speed-drop>=1, turn>=15 deg in 0:+0.5 s, 2-second spacing and prior full-window protocol rules';
for r=1:height(S.T),s=find(["Control","Heat"]==string(S.T.condition(r)));c=S.T.displayColumn(r);S.T.events(r)=S.counts(s,c);S.T.barrierTrials(r)=S.barrierCounts(s,c);end
writetable(S.T,fullfile(out,'selection_counts.csv'));
S.methods.candidateRuns='Contiguous negative-acceleration runs replace the old <=-2 runs; merging runs can change alignment and event membership';
S.newEvents=newEvents;S.oldCounts=oldCounts;save(fullfile(out,'comparison_data.mat'),'-struct','S','-v7.3');writetable(newEvents,fullfile(out,'slowdown_events.csv'));
plot_pizza_1020_11_side_by_side(out);
fid=fopen(fullfile(out,'methods.json'),'w');fprintf(fid,'%s',jsonencode(S.methods,PrettyPrint=true));fclose(fid);
result=struct('output',out,'oldSlowdownCounts',oldCounts(:,2)','newSlowdownCounts',S.counts(:,2)','acceptedWithAbsDecelBelow2',nnz(conditions~="Excluded"&acc(idx)>-2));disp(result);
end
