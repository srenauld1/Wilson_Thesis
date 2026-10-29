function result=plot_pizza_ves041_bar_tuning(root)
% Exploratory angular tuning. Start-relative zero is provisional, not calibrated.
if nargin<1,root='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
out=fullfile(root,'bar_tuning');if ~isfolder(out),mkdir(out);end
M=readtable(fullfile(root,'manifest.csv'),'TextType','string');
F=readtable(fullfile(root,'bar_fly_mapping.csv'),'TextType','string');
M=M(M.status=="complete" & M.category=="bar",:);
edges=-180:30:180;centers=edges(1:end-1)+15;lags=[0 .5 1];
R=struct([]);audit=cell(0,6);
for k=1:height(M)
 name=M.trial(k);source=M.source(k);fprintf('BAR TUNING %s\n',name);
 if ~endsWith(source,'.mat')
  audit(end+1,:)={name,"unavailable: no saved synchronized visual heading",source,NaN,NaN,NaN};continue
 end
 S=load(source,'fly_data');a=S.fly_data;q=a.dq(1);
 if ~isfield(q,'vh')||isempty(q.vh)
  audit(end+1,:)={name,"unavailable: vh missing",source,NaN,NaN,NaN};continue
 end
 D=load(fullfile(root,'analysis_data',name+"_events.mat"),'d');d=D.d;
 t=double(q.t(:));vh=double(q.vh(:));assert(numel(vh)==numel(t));
 assert(numel(t)==numel(d.t)&&max(abs(t-d.t))<1e-5,'Source time mismatch');
 valid=d.valid(:)&isfinite(vh)&isfinite(d.dff(:));
 % Midpoint sample durations handle nonuniform synchronized volume times.
 dt=diff([t(1);(t(1:end-1)+t(2:end))/2;t(end)]);
 dt=min(dt,2*median(diff(t))); % Never assign a long missing interval to a sample.
 start=find(isfinite(vh),1);offset=vh(start);
 angle=rad2deg(atan2(sin(vh-offset),cos(vh-offset)));
 recorded=rad2deg(atan2(sin(vh),cos(vh)));
 z=(d.dff-mean(d.dff(valid)))/std(d.dff(valid));
 jump=false(size(t));if isfield(a,'jump_detected')&&numel(a.jump_detected)==numel(t),jump=logical(a.jump_detected(:));end
 % Remove peri-jump samples only in the explicitly labeled companion curve.
 away=true(size(t));for j=find(jump)',away(abs(t-t(j))<=1)=false;end
 masks=[valid,valid&abs(d.yaw)<20,valid&abs(d.yaw)>=20,valid&away];
 raw=nan(12,4);occ=raw;standard=raw;
 for c=1:4
  [raw(:,c),occ(:,c)]=curve(angle,d.dff,masks(:,c),edges,dt);
  standard(:,c)=curve(angle,z,masks(:,c),edges,dt);
 end
 rec=curve(recorded,z,valid,edges,dt);
 half=nan(12,2);mid=(t(1)+t(end))/2;
 half(:,1)=curve(angle,z,valid&t<=mid,edges,dt);
 half(:,2)=curve(angle,z,valid&t>mid,edges,dt);
 lagCurves=nan(12,3);
 for j=1:3
  % Calcium at t paired with bar angle at t-lag; interpolate on circle.
  earlier=interp1(t,unwrap(vh),t-lags(j),'linear',NaN)-offset;
  ang=rad2deg(atan2(sin(earlier),cos(earlier)));
  lagCurves(:,j)=curve(ang,z,valid&isfinite(ang),edges,dt);
 end
 fwd=curve(angle,d.fwd,valid,edges,dt);yaw=curve(angle,abs(d.yaw),valid,edges,dt);
 r=struct('trial',name,'flyId',F.flyId(F.trial==name),'source',source, ...
  'initialRecordedDeg',rad2deg(offset),'initialSampleTime',t(start), ...
  'raw',raw,'z',standard,'recordedZ',rec,'occupancySec',occ,'halfZ',half, ...
  'lagZ',lagCurves,'forward',fwd,'yaw',yaw,'detectedJumps',nnz(jump));
 if isempty(R),R=r;else,R(end+1)=r;end %#ok<AGROW>
 audit(end+1,:)={name,"included: start-relative reference provisional",source,rad2deg(offset),t(start),sum(dt(valid))};
 fig=figure('Visible','off','Color','w','Position',[20 20 1550 900]);tl=tiledlayout(fig,2,3,'TileSpacing','compact');
 ax=nexttile(tl);plot(ax,t(t<=5),recorded(t<=5),'.-');xlabel(ax,'Time (s)');ylabel(ax,'Recorded angle (deg)');title(ax,'First 5 s: reference audit');grid(ax,'on');
 ax=nexttile(tl);bar(ax,centers,occ(:,1));xlabel(ax,'Start-relative bar angle (deg)');ylabel(ax,'Occupancy (s)');title(ax,'All valid samples');xlim(ax,[-180 180]);
 ax=nexttile(tl);plot(ax,centers,raw,'LineWidth',1.5);legend(ax,{'All','|yaw| <20 deg/s','|yaw| >=20 deg/s','Outside +/-1 s of detected jumps'},'Location','best');ylabel(ax,'Mean original \DeltaF/F');title(ax,'Activity / occupied angle bin');format_angle(ax);
 ax=nexttile(tl);yyaxis(ax,'left');plot(ax,centers,fwd,'LineWidth',1.5);ylabel(ax,'Forward (mm/s)');yyaxis(ax,'right');plot(ax,centers,yaw,'LineWidth',1.5);ylabel(ax,'|yaw| (deg/s)');title(ax,'Behavior varies with bar angle');format_angle(ax);
 ax=nexttile(tl);plot(ax,centers,half,'LineWidth',1.5);legend(ax,{'First half','Second half'},'Location','best');ylabel(ax,'Mean activity (trial z score)');title(ax,'Does the angular pattern repeat?');format_angle(ax);
 ax=nexttile(tl);plot(ax,centers,lagCurves,'LineWidth',1.5);legend(ax,{'0 s','0.5 s','1 s'},'Location','best');ylabel(ax,'Mean activity (trial z score)');title(ax,'Bar precedes activity by fixed lag');format_angle(ax);
 title(tl,strrep(name,'_','\_'));subtitle(tl,'Exploratory | 30-degree bins; bins with <2 s occupancy hidden | 0 = first saved bar angle, provisionally front');
 saveplots(fig,fullfile(out,name+"_bar_tuning"));
 save(fullfile(out,name+"_bar_tuning.mat"),'r','edges','centers','lags');
end
A=cell2table(audit,'VariableNames',{'trial','status','source','initialRecordedDeg','initialSampleTime','validSeconds'});writetable(A,fullfile(out,'bar_tuning_audit.csv'));
assert(~isempty(R),'No synchronized visual heading found');
flies=unique(string({R.flyId}));nf=numel(flies);flyZ=nan(12,nf,4);flyOcc=zeros(12,nf);
for f=1:nf
 ix=find(string({R.flyId})==flies(f));
 % Pool occupied samples within each fly using occupancy weights; equal flies after.
 for c=1:4
  num=zeros(12,1);den=num;
  for i=ix
   w=R(i).occupancySec(:,c);v=R(i).z(:,c);w(~isfinite(v))=0;v(~isfinite(v))=0;
   num=num+w.*v;den=den+w;
  end
  flyZ(:,f,c)=num./den;
 end
 for i=ix,flyOcc(:,f)=flyOcc(:,f)+R(i).occupancySec(:,1);end
end
fig=figure('Visible','off','Color','w','Position',[20 20 1550 1000]);tl=tiledlayout(fig,2,3,'TileSpacing','compact');
labels=string({R.trial});
ax=nexttile(tl);heat(ax,centers,cat(2,R.recordedZ)',labels);title(ax,'Recorded coordinates: trial z score');
ax=nexttile(tl);v=cat(3,R.z);heat(ax,centers,squeeze(v(:,1,:))',labels);title(ax,'Start-relative coordinates: trial z score');
ax=nexttile(tl);v=cat(3,R.occupancySec);heat(ax,centers,squeeze(v(:,1,:))',labels);title(ax,'Start-relative occupancy (s)');
for c=1:3
 ax=nexttile(tl);plot(ax,centers,flyZ(:,:,c),'Color',[.7 .7 .7]);hold(ax,'on');plot(ax,centers,mean(flyZ(:,:,c),2,'omitnan'),'k','LineWidth',2);format_angle(ax);ylabel(ax,'Activity (within-trial z score)');names={'All samples','Low rotation: |yaw| <20 deg/s','High rotation: |yaw| >=20 deg/s'};title(ax,names{c});
end
title(tl,sprintf('VES041 exploratory bar tuning | %d recordings, %d estimated flies',numel(R),nf));
subtitle(tl,'Start-relative zero is provisional; sign conventions may differ | gray: fly means; black: equal-fly mean | <2 s bins hidden');
saveplots(fig,fullfile(out,'ves041_bar_tuning_overview'));
methods=struct('reference','Recorded vh radians, and wrap(vh-first finite vh); start-relative zero is provisionally front, not physical calibration', ...
 'sign','Recorded vh sign retained; left/right labels intentionally omitted until era-specific sign is verified', ...
 'bins','30 degrees, minimum 2 s valid occupancy per trial/condition; midpoint time weights (capped at twice median interval); time-weighted activity means; no circular smoothing', ...
 'activity','Original source ROI dF/F; z score uses full valid trial, not tuning curve normalization', ...
 'pool','Within fly occupancy-weight trial z scores; then equal mean across flies with data per bin', ...
 'lag','Fixed 0, 0.5, 1 s; bar angle precedes calcium; no best lag selection', ...
 'caution','Closed-loop bar angle covaries with behavior; angular association alone is not evidence of a purely visual response; start reference is sensitive to motion before first saved sample', ...
 'exclusions','Two March recordings have recovered calcium/behavior traces but no saved synchronized vh; omitted pending DAQ recovery');
save(fullfile(out,'ves041_bar_tuning_summary.mat'),'R','A','flies','flyZ','flyOcc','edges','centers','lags','methods');
fid=fopen(fullfile(out,'methods.json'),'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(methods,PrettyPrint=true));fclose(fid);
result=struct('folder',out,'recordings',numel(R),'estimatedFlies',nf);disp(result);
end
function [v,occ]=curve(angle,y,valid,edges,dt)
b=discretize(angle,edges);v=nan(numel(edges)-1,1);occ=zeros(size(v));
for k=1:numel(v),m=valid&b==k&isfinite(y);occ(k)=sum(dt(m));if occ(k)>=2,v(k)=sum(y(m).*dt(m))/occ(k);end,end
end
function format_angle(ax)
xlim(ax,[-180 180]);xticks(ax,-180:90:180);xlabel(ax,'Start-relative bar angle (deg)');grid(ax,'on');
end
function heat(ax,x,z,labels)
h=imagesc(ax,x,1:size(z,1),z);h.AlphaData=isfinite(z);set(ax,'Color',[.85 .85 .85]);yticks(ax,1:numel(labels));yticklabels(ax,strrep(labels,'_','\_'));ax.YAxis.FontSize=7;colorbar(ax);xlabel(ax,'Bar angle (deg)');
end
function saveplots(fig,stem)
exportgraphics(fig,stem+".png",'Resolution',150);savefig_visible_on_open(fig,stem+".fig");close(fig);
end
