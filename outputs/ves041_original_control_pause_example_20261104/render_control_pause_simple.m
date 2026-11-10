function render_control_pause_simple(dataDir)
% Minimal kinematics with Gaussian smoothing (sigma 10 ms).
if nargin<1,dataDir=fullfile(pwd,'outputs','ves041_original_control_pause_example_20261104');end
s=load(fullfile(dataDir,'Fly30_pulse31_original_detector.mat'));
m=s.metadata;win=reshape(m.window_zoom,1,[]);take=s.time>=win(1)&s.time<=win(2);
t=s.time(take)-m.light_start_s;
% Smooth the full available 25 s segment before taking the displayed crop.
dt=median(diff(s.time));sigmaSeconds=.010;
radius=ceil(4*sigmaSeconds/dt);offset=(-radius:radius)*dt;
kernel=exp(-.5*(offset/sigmaSeconds).^2);kernel=kernel/sum(kernel);
forwardSmooth=conv([repmat(s.forward(1),radius,1);s.forward;repmat(s.forward(end),radius,1)],kernel,'valid');
yawSmooth=conv([repmat(s.yaw(1),radius,1);s.yaw;repmat(s.yaw(end),radius,1)],kernel,'valid');
assert(numel(forwardSmooth)==numel(s.forward) && numel(yawSmooth)==numel(s.yaw));
f=figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1000 850], ...
 'Name','Control Fly 30: simple forward and rotation','NumberTitle','off');
set(f,'PaperUnits','inches','PaperPosition',[0 0 10 8.5],'PaperSize',[10 8.5],'PaperPositionMode','manual');
blue=[0 .30 .85];black=[0 0 0];shade=[240 225 250]/255;
pos=[.15 .57 .78 .26;.15 .20 .78 .26];
values={forwardSmooth(take),yawSmooth(take)};colors={black,blue};lims={[-15 25],[-400 400]};
labels={'Forward (mm/s)','Rotation (deg/s)'};
for j=1:2
 ax=axes(f,'Units','normalized','Position',pos(j,:),'XLim',win-m.light_start_s, ...
  'YLim',lims{j},'FontName','Arial','FontSize',17,'Box','off','TickDir','out', ...
  'XColor','none','YColor',colors{j},'XTick',[],'LineWidth',.8);
 hold(ax,'on');yl=lims{j};duration=m.light_end_s-m.light_start_s;
 patch(ax,[0 duration duration 0],[yl(1) yl(1) yl(2) yl(2)],shade,'EdgeColor','none','Tag','light_on');
 if j==2
  yline(ax,0,'-','Color',[.55 .55 .55],'LineWidth',.7,'Tag','rotation_zero_axis');
 end
 plot(ax,t,values{j},'Color',colors{j},'LineWidth',1.1,'Tag','smoothed_velocity');
 ylabel(ax,labels{j},'FontName','Arial','FontSize',18,'Color',colors{j});
 if j==1,set(ax,'YTick',[-10 0 10 20]);else,set(ax,'YTick',[-400 -200 0 200 400]);end
 if j==2
  % One-second scale bar, separate from the trace.
  barX=[5.8 6.8];barY=-570;
  plot(ax,barX,[barY barY],'k-','LineWidth',2.2,'Clipping','off','Tag','time_scale_bar');
  text(ax,mean(barX),barY-40,'1 s','FontName','Arial','FontSize',17, ...
   'HorizontalAlignment','center','VerticalAlignment','top','Clipping','off');
 end
end
annotation(f,'textbox',[.15 .90 .78 .045],'String','Control Fly 30', ...
 'FontName','Arial','FontSize',21,'FontWeight','bold','LineStyle','none','Margin',0);
annotation(f,'textbox',[.15 .845 .3 .035],'String','Light on', ...
 'FontName','Arial','FontSize',16,'Color',[.5 .15 .7],'LineStyle','none','Margin',0);
figData=struct('gaussian_sigma_s',sigmaSeconds,'time_s',s.time,'forward_smoothed',forwardSmooth,'yaw_smoothed',yawSmooth,'original',s);
save(fullfile(dataDir,'Control_Fly30_simple_black_blue_smooth10ms.mat'),'-struct','figData');
dest=fullfile(dataDir,'Control_Fly30_simple_black_blue_smooth10ms');
print(f,[dest '.svg'],'-dsvg','-painters');
print(f,[dest '.png'],'-dpng','-r600');
set(f,'Visible','on');savefig(f,[dest '.fig']);close(f);
fprintf('Saved simple SVG, PNG and FIG: %s\n',dest);
end
