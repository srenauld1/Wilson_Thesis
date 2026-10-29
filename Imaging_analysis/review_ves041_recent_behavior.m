function T=review_ves041_recent_behavior()
% Behavior-only review; fluorescence and stimulus responses do not select trials.
base=fileparts(mfilename('fullpath'));out=fullfile(imaging_analysis_artifact_root(),'reports','ves041_behavior_review_20261028');if ~isfolder(out),mkdir(out);end
I=readtable(fullfile(imaging_analysis_artifact_root(),'reports','ves041_trial_inventory_20261028.csv'),'TextType','string','Delimiter',',');
keep=true(height(I),1);
for j=1:height(I)
 if ismissing(I.analysis_folder(j)),keep(j)=isfile(fullfile(I.folder(j),'smr_analysis','fly_data.mat'));end
end
I=I(keep,:);
rows=cell(height(I),20);D=cell(height(I),1);
for k=1:height(I)
 folder=I.analysis_folder(k);if ismissing(folder)||strlength(folder)==0,folder=fullfile(I.folder(k),'smr_analysis');end
 file=fullfile(folder,'fly_data.mat');
 try
 S=load(file,'fly_data');a=S.fly_data;q=a.dq(1);t=double(q.t(:));f=double(q.bvf(:));
 if isfield(q,'bvy_deg'),r=double(q.bvy_deg(:));else,r=rad2deg(double(q.bvy(:)));end
 dt=median(diff(t));valid=isfinite(t)&isfinite(f)&isfinite(r);maskAvailable=false;
 if isfield(q,'bergii')&&isfield(q.bergii,'behaviorVolumeValid'),valid=valid&logical(q.bergii.behaviorVolumeValid(:));maskAvailable=true;end
 f(~valid)=nan;r(~valid)=nan;n=max(3,round(.25/dt));f=smoothdata(f,'movmedian',n,'omitnan');r=smoothdata(r,'movmedian',n,'omitnan');f(~valid)=nan;r(~valid)=nan;
 walk=valid&f>=1;turn=valid&abs(r)>=30;still=valid&abs(f)<.5&abs(r)<10;
 edges=t(1):30:t(end);active=nan(max(1,numel(edges)-1),1);
 for j=1:numel(edges)-1,z=t>=edges(j)&t<edges(j+1)&valid;if any(z),active(j)=mean(walk(z));end,end
 active=active(isfinite(active));longStill=longest(still,dt);longInvalid=longest(~valid,dt);
 bouts=diff([false;walk;false]);bs=find(bouts==1);be=find(bouts==-1)-1;dur=(be-bs+1)*dt;nb=nnz(dur>=1);
 if isfield(q,'heat'),h=double(q.heat(:));else,h=zeros(size(t));end
 if isfield(a,'protocol'),kind=string(a.protocol.kind);else,kind=I.context(k);end
 if isfield(q,'bh'),heading=double(q.bh(:));else,heading=cumsum(fillmissing(r,'constant',0))*dt*pi/180;end
 x=cumsum(fillmissing(f,'constant',0).*sin(heading))*dt;y=cumsum(fillmissing(f,'constant',0).*cos(heading))*dt;
 % Paths here are reconstructed behavioral paths, not new stimulus coordinates.
 if isfield(a,'protocol')&&isfield(a.protocol,'pos_x_mm'),x=a.protocol.pos_x_mm(:);y=a.protocol.pos_y_mm(:);end
 rows(k,:)={I.trial(k),I.status(k),kind,t(end)-t(1),100*mean(valid),maskAvailable,100*mean(walk(valid)),median(f(valid)),prctile(f(valid),95),100*mean(turn(valid)),100*mean(r(valid)>30),100*mean(r(valid)<-30),longStill,longInvalid,nb,100*mean(active>=.1),100*mean(h(valid)>.05),sum(max(f(valid),0))*dt,~ismissing(I.analysis_folder(k))&&strlength(I.analysis_folder(k))>0,file};
 D{k}=struct('t',t-t(1),'f',f,'r',r,'h',h,'x',x,'y',y,'valid',valid);
 fprintf('%s: valid %.1f%% walk %.1f%% active bins %.1f%% longest still %.1fs\n',I.trial(k),rows{k,5},rows{k,7},rows{k,16},longStill);
 catch ME,error('Review failed for %s: %s',I.trial(k),ME.message);end
end
T=cell2table(rows,'VariableNames',{'trial','priorStatus','context','durationSec','validPct','trackingMaskAvailable','walkingPct','medianForwardMmS','p95ForwardMmS','turningPct','positiveTurnPct','negativeTurnPct','longestStillSec','longestInvalidSec','walkingBoutsAtLeast1Sec','active30SecBinsPct','heatPct','forwardDistanceMm','recent','sourceFile'});
writetable(T,fullfile(out,'behavior_metrics.csv'));
ix=find(T.recent);
for page=1:ceil(numel(ix)/8)
 ids=ix((page-1)*8+1:min(page*8,numel(ix)));fig=figure('Visible','off','Color','w','Position',[10 10 1800 2200]);tl=tiledlayout(fig,numel(ids),4,'TileSpacing','compact','Padding','compact');
 for k=ids'
 d=D{k};label=erase(T.trial(k),'_ves041');
 ax=nexttile(tl);plot(d.t,d.f,'k');ylim([-5 max(10,min(20,T.p95ForwardMmS(k)*1.4))]);ylabel('Forward mm/s');title(sprintf('%s | walk %.0f%%',label,T.walkingPct(k)),'Interpreter','none');xlim([0 d.t(end)]);
 ax=nexttile(tl);plot(d.t,d.r,'Color',[.2 .3 .7]);ylim([-300 300]);ylabel('Yaw deg/s');title(sprintf('turn %.0f%% | still max %.0fs',T.turningPct(k),T.longestStillSec(k)));xlim([0 d.t(end)]);
 ax=nexttile(tl);plot(d.x,d.y,'k');axis equal;title('XY trajectory');xlabel('mm');ylabel('mm');
 ax=nexttile(tl);plot(d.t,d.h,'r');hold on;plot(d.t,~d.valid,'Color',[.4 .4 .4]);ylim([-.1 max(1.1,max(d.h)+.2)]);title(sprintf('Heat V (red); invalid (gray) | %.1f%% valid',T.validPct(k)));xlim([0 d.t(end)]);
 end
 exportgraphics(fig,fullfile(out,sprintf('behavior_review_%02d.png',page)),'Resolution',130);close(fig);
end
save(fullfile(out,'behavior_review_data.mat'),'D','T','-v7.3');
end
function n=longest(mask,dt)
z=diff([false;mask;false]);n=max([0;(find(z==-1)-find(z==1))*dt]);
end
