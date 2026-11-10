function render_acr_reference_scale(dataDir)
% Match the edited Control_inhibition_acr_timeseries.svg coordinate gains.
% Full peaks remain visible even when outside the labeled rotation range.
if nargin<1,dataDir=fileparts(mfilename('fullpath'));end
s=load(fullfile(dataDir,'ACR_Fly22C_pulse31_original_detector.mat'));
m=s.metadata;win=reshape(m.window_zoom,1,[])-1.5;
dt=median(diff(s.time));sigma=.010;radius=ceil(4*sigma/dt);offset=(-radius:radius)*dt;
kernel=exp(-.5*(offset/sigma).^2);kernel=kernel/sum(kernel);
fw=conv([repmat(s.forward(1),radius,1);s.forward;repmat(s.forward(end),radius,1)],kernel,'valid');
yaw=conv([repmat(s.yaw(1),radius,1);s.yaw;repmat(s.yaw(end),radius,1)],kernel,'valid');
take=s.time>=win(1)&s.time<=win(2);t=s.time(take)-win(1);
cw=334.574;ch=420.403;left=34.063;width=295.451;height=159;tops=[.4 226.9];
fig=figure('Visible','off','Color','white','Units','points','Position',[50 50 cw ch], ...
 'Name','ACR Fly 22C: shifted 1.5 s earlier, reference scaling','NumberTitle','off');
set(fig,'PaperUnits','points','PaperPosition',[0 0 cw ch],'PaperSize',[cw ch],'PaperPositionMode','manual');
blue=[0 77 217]/255;shade=[240 225 250]/255;
values={fw(take),yaw(take)};limits={[-15 25],[-300 300]};colors={[0 0 0],blue};
for j=1:2
 pos=[left/cw 1-(tops(j)+height)/ch width/cw height/ch];
 ax=axes(fig,'Units','normalized','Position',pos,'XLim',[0 10],'YLim',limits{j}, ...
  'FontName','Arial','FontSize',12.24,'Box','off','TickDir','out','LineWidth',.8, ...
  'TickLength',[5.615/width .025],'XColor','none','XTick',[],'YColor',colors{j});
 hold(ax,'on');yl=limits{j};ab=[m.light_start_s m.light_end_s]-win(1);
 patch(ax,[ab(1) ab(2) ab(2) ab(1)],[yl(1) yl(1) yl(2) yl(2)],shade,'EdgeColor','none','Tag','light_on');
 if j==1,set(ax,'YTick',[-10 0 10 20]);else
  set(ax,'YTick',[-300 -150 0 150 300],'YTickLabel',{'-300','','0','','300'});
  plot(ax,[0 10],[0 0],'Color',blue*.702+[1 1 1]*.298,'LineWidth',1.1,'Tag','rotation_zero_axis');
 end
 plot(ax,t,values{j},'Color',colors{j},'LineWidth',1.1,'Clipping','off','Tag','smoothed_velocity');
 set(ax,'Position',pos);
 assert(abs(ax.Position(3)*cw/diff(ax.XLim)-29.5451)<1e-8);
 if j==1,assert(abs(ax.Position(4)*ch/diff(ax.YLim)-3.975)<1e-8);
 else,assert(abs(ax.Position(4)*ch/diff(ax.YLim)-.265)<1e-8);end
end
annotation(fig,'line',[36.708 66.253]/cw,[1 1]*(1-402.822/ch),'Color','black','LineWidth',2.2);
annotation(fig,'textbox',[44/cw 0 30/cw 17/ch],'String','1 s','FontName','Arial', ...
 'FontSize',12.24,'LineStyle','none','Margin',0,'HorizontalAlignment','center');
fig.UserData=struct('window_recording_s',win,'window_relative_light_s',win-m.light_start_s, ...
 'gaussian_sigma_s',sigma,'time_points_per_s',29.5451,'forward_points_per_mm_s',3.975,'rotation_points_per_deg_s',.265);
dest=fullfile(dataDir,'ACR_Fly22C_shift1p5s_rot300_no_circle');
metadata=fig.UserData;time=s.time;forward_smoothed=fw;yaw_smoothed=yaw;
save([dest '.mat'],'metadata','time','forward_smoothed','yaw_smoothed','-v7');
print(fig,[dest '_matlab.svg'],'-dsvg','-painters');
set(fig,'Visible','on');savefig(fig,[dest '.fig']);close(fig);
fprintf('Verified reference gains and saved editable FIG: %s\n',dest);
end
