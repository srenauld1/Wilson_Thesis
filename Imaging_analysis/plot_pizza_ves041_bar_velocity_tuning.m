function result=plot_pizza_ves041_bar_velocity_tuning(root)
% Simple position curves separated by signed visual-bar velocity, not fly yaw.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
out=fullfile(root,'bar_tuning','bar_velocity_split');if ~isfolder(out),mkdir(out);end
M=readtable(fullfile(root,'manifest.csv'),'TextType','string');
names=["20261019_3_2","20261020_7_1","20261024_2_4_ves041_barjumps"];
M=M(ismember(M.trial,names)&M.status=="complete"&M.category=="bar",:);
assert(height(M)==3,'Expected the three accepted October bar recordings');
edges=-180:30:180;centers=edges(1:end-1)+15;R=struct([]);rows=cell(0,7);
for k=1:height(M)
 S=load(M.source(k),'fly_data');a=S.fly_data;q=a.dq(1);
 D=load(fullfile(root,'analysis_data',M.trial(k)+"_events.mat"),'d');d=D.d;
 t=double(q.t(:));vh=double(q.vh(:));assert(isequal(t,d.t));
 assert(abs(rad2deg(vh(find(isfinite(vh),1))))<1,'Initial angle is not near zero');
 angle=rad2deg(atan2(sin(vh),cos(vh))); % Preserve recorded zero, no recentering.
 fs=1/median(diff(t));width=max(1,round(.2*fs));
 velocity=movmedian(rad2deg(gradient(unwrap(vh),t)),width,'omitnan');
 assert(isfield(a,'jump_detected')&&numel(a.jump_detected)==numel(t),'Missing synchronized jump labels');
 jump=logical(a.jump_detected(:));away=true(size(t));
 for j=find(jump)',away(abs(t-t(j))<=.25)=false;end
 gaps=find(diff(t)>2*median(diff(t)));
 for j=gaps',away(t>=t(j)-.25&t<=t(j+1)+.25)=false;end
 valid=d.valid(:)&isfinite(vh)&isfinite(velocity)&isfinite(d.dff)&away;
 weights=diff([t(1);(t(1:end-1)+t(2:end))/2;t(end)]);
 weights=min(weights,2*median(diff(t)));
 masks=[valid,valid&velocity>1,valid&velocity< -1];
 b=discretize(angle,edges);mu=nan(12,3);occupancy=zeros(12,3);
 for c=1:3
  for p=1:12
   keep=masks(:,c)&b==p;occupancy(p,c)=sum(weights(keep));
   if occupancy(p,c)>=2,mu(p,c)=sum(d.dff(keep).*weights(keep))/occupancy(p,c);end
   rows(end+1,:)={M.trial(k),edges(p),edges(p+1),c,mu(p,c),occupancy(p,c),nnz(keep)}; %#ok<AGROW>
  end
 end
 r=struct('trial',M.trial(k),'source',M.source(k),'t',t,'barAngleDeg',angle, ...
  'barVelocityDegS',velocity,'valid',valid,'meanDff',mu,'occupancySec',occupancy,'jumpTimes',t(jump));
 if isempty(R),R=r;else,R(end+1)=r;end %#ok<AGROW>
 fig=figure('Visible','off','Color','w','Position',[20 20 1100 480]);tl=tiledlayout(fig,1,2,'TileSpacing','compact');
 draw_curve(nexttile(tl),centers,mu);draw_occupancy(nexttile(tl),centers,occupancy);
 title(tl,strrep(M.trial(k),'_','\_'));subtitle(tl,'Recorded 0 deg reference | 30 deg bins | Bar velocity sign, not fly turning | Original ROI dF/F');
 saveplots(fig,fullfile(out,M.trial(k)+"_position_by_bar_velocity"));
 save(fullfile(out,M.trial(k)+"_position_by_bar_velocity.mat"),'r','edges','centers');
 fprintf('%s: %d jumps; %.1f s positive, %.1f s negative\n',M.trial(k),nnz(jump),sum(occupancy(:,2)),sum(occupancy(:,3)));
end
fig=figure('Visible','off','Color','w','Position',[20 20 1500 800]);tl=tiledlayout(fig,2,3,'TileSpacing','compact');
for k=1:3,ax=nexttile(tl);draw_curve(ax,centers,R(k).meanDff);title(ax,strrep(R(k).trial,'_','\_'));end
for k=1:3,draw_occupancy(nexttile(tl),centers,R(k).occupancySec);end
title(tl,'VES041 dF/F versus bar position, separated by bar velocity');
subtitle(tl,'October trials: recorded 0 deg reference | Blue: increasing angle; orange: decreasing angle | Each column is one trial');
saveplots(fig,fullfile(out,'ves041_latest_position_by_bar_velocity'));
T=cell2table(rows,'VariableNames',{'trial','binStartDeg','binEndDeg','condition','meanDff','occupancySec','samples'});
T.condition=string(T.condition);T.condition(T.condition=="1")="all";T.condition(T.condition=="2")="barVelocityPositive";T.condition(T.condition=="3")="barVelocityNegative";
writetable(T,fullfile(out,'position_by_bar_velocity_bins.csv'));
methods=struct('trials',names,'position','Recorded vh wrapped to [-180,180] degrees; no per-trial offset subtraction', ...
 'velocity','Time derivative of unwrapped vh, then 0.2 s moving median; deg/s; independent of fly yaw', ...
 'sign','Positive = increasing recorded angle; negative = decreasing recorded angle; no anatomical left/right interpretation', ...
 'deadband','Positive >+1 deg/s, negative <-1 deg/s; near-stationary samples retained only in gray all-sample curve', ...
 'jumpExclusion','Exclude +/-0.25 s around saved jump_detected labels, and around time gaps >2 median sample intervals', ...
 'averaging','Original ROI dF/F, same-time bar position/velocity, no lag or event baseline subtraction; midpoint time-weighted means', ...
 'bins','30 degree bins; require >=2 s occupancy per bin and condition; no curve smoothing', ...
 'interpretation','Exploratory association; fly movement is not matched between bar-velocity directions');
save(fullfile(out,'position_by_bar_velocity_summary.mat'),'R','T','methods','edges','centers');
fid=fopen(fullfile(out,'methods.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('folder',out,'trials',height(M));disp(result);
end
function draw_curve(ax,x,y)
hold(ax,'on');plot(ax,x,y(:,1),'-','Color',[.65 .65 .65],'LineWidth',1.2);
plot(ax,x,y(:,2),'-o','Color',[0 .447 .741],'LineWidth',1.8,'MarkerSize',4);
plot(ax,x,y(:,3),'-o','Color',[.85 .325 .098],'LineWidth',1.8,'MarkerSize',4);
legend(ax,{'All','Bar velocity +','Bar velocity -'},'Location','best');ylabel(ax,'Mean VES041 \DeltaF/F');format(ax);title(ax,'Activity within position bins');
end
function draw_occupancy(ax,x,y)
bar(ax,x,y(:,2:3),'grouped');colororder(ax,[0 .447 .741;.85 .325 .098]);ylabel(ax,'Time in bin (s)');format(ax);title(ax,'Occupancy by bar-velocity direction');
end
function format(ax)
xlim(ax,[-180 180]);xticks(ax,-180:90:180);xlabel(ax,'Bar position (deg; 0 = front)');grid(ax,'on');
end
function saveplots(fig,stem)
exportgraphics(fig,stem+".png",'Resolution',160);savefig_visible_on_open(fig,stem+".fig");close(fig);
end
