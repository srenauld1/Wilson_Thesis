function render_original_control_pause_example(dataDir)
% Reproduce a control light-on reorientation from the original spatial rule.
% The MAT file is portable. Raw velocity samples are plotted without smoothing.
if nargin<1,dataDir=fileparts(mfilename('fullpath'));end
matPath=fullfile(dataDir,'Fly30_pulse31_original_detector.mat');
if isfile(matPath)
    load(matPath,'metadata','time','forward','yaw');
else
    metadata=jsondecode(fileread(fullfile(dataDir,'selected_example.json')));
    raw=readmatrix(fullfile(dataDir,'Fly30_pulse31_raw_velocity.csv'));
    time=raw(:,1);forward=raw(:,2);yaw=raw(:,3);
    save(matPath,'metadata','time','forward','yaw','-v7');
end
m=metadata;
assert(m.verification.passed && m.rotation_deg>30 && m.displacement_mm<=2 && m.directionality>=.7);
assert(m.start_s>=m.light_start_s && m.end_s<=m.light_end_s);
light=[150 61 206]/255;shade=[240 225 250]/255;orange=[197 82 22]/255;
ink=[32 33 36]/255;gray=[104 107 112]/255;
for view=1:2
    if view==1,win=m.window_25s;tag='25s';else,win=m.window_zoom;tag='10s_zoom';end
    win=reshape(win,1,[]);
    fig=figure('Visible','off','Color','w','Units','pixels','Position',[100 100 1050 1050], ...
        'Name',['Control Fly 30: original-detector light-on event (' tag ')'],'NumberTitle','off');
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 10.5 10.5], ...
        'PaperSize',[10.5 10.5],'PaperPositionMode','manual');
    note(fig,[.075 .944 .89 .04],'Control Fly 30 | VES041_1 x GFP',23,'bold',ink);
    note(fig,[.075 .910 .89 .03],sprintf('Original detector | pulse 31 | reorientation at %.2f s',m.peak_yaw_time_s),15,'normal',gray);
    annotation(fig,'rectangle',[.078 .878 .02 .015],'Color',shade,'FaceColor',shade);
    note(fig,[.108 .871 .14 .027],'Light on',14,'normal',ink);
    annotation(fig,'line',[.26 .29],[.885 .885],'Color',orange,'LineWidth',2);
    note(fig,[.30 .871 .31 .027],'Detected 1 s window',14,'normal',ink);
    annotation(fig,'line',[.60 .63],[.885 .885],'Color',ink,'LineWidth',1);
    note(fig,[.64 .871 .22 .027],'Raw velocity',14,'normal',ink);
    note(fig,[.075 .828 .89 .03],'Stored raw velocity fields | 1 kHz sampling | no added smoothing',13,'normal',gray);
    positions=[.12 .53 .82 .25;.12 .16 .82 .25];
    signals={forward,yaw};ylims={[-15 25],[-400 400]};yticks={[-10 0 10 20],[-400 -200 0 200 400]};
    labels={'Forward velocity (mm/s)','Yaw velocity (deg/s)'};
    axesHandles=gobjects(1,2);use=time>=win(1)&time<=win(2);
    relativeTime=time-m.light_start_s;relativeWindow=win-m.light_start_s;
    for j=1:2
        axesHandles(j)=axes(fig,'Units','normalized','Position',positions(j,:), ...
            'FontName','Arial','FontSize',15,'Box','off','TickDir','out', ...
            'XLim',relativeWindow,'YLim',ylims{j},'YTick',yticks{j}, ...
            'XColor',gray,'YColor',gray,'LineWidth',.7,'Layer','top');
        ax=axesHandles(j);hold(ax,'on');
        ybounds=ylims{j};duration=m.light_end_s-m.light_start_s;
        patch(ax,[0 duration duration 0],[ybounds(1) ybounds(1) ybounds(2) ybounds(2)],shade,'EdgeColor','none');
        yline(ax,0,'Color',[.85 .86 .87],'LineWidth',.6);
        plot(ax,relativeTime(use),signals{j}(use),'Color',ink,'LineWidth',.9,'Tag','raw_velocity');
        event=time>=m.start_s&time<=m.end_s;
        plot(ax,relativeTime(event),signals{j}(event),'Color',orange,'LineWidth',1.9,'Tag','original_detection_window');
        xline(ax,m.peak_yaw_time_s-m.light_start_s,':','Color',orange,'LineWidth',.85);
        note(fig,[.12 positions(j,2)+.26 .82 .033],labels{j},17,'bold',ink);
        if j==1,set(ax,'XTick',[]);else
            if view==1,xt=-5:5:15;else,xt=-2:1:7;end
            set(ax,'XTick',xt);
            xlabel(ax,'Time from light onset (s)','FontName','Arial','FontSize',15);
        end
        assert(isequal(ax.XLim,relativeWindow));
        assert(min(signals{j}(use))>=ybounds(1)&&max(signals{j}(use))<=ybounds(2));
    end
    note(fig,[.075 .037 .89 .04],sprintf('Net turn %.1f deg | net displacement %.2f mm | directionality %.2f', ...
        m.rotation_deg,m.displacement_mm,m.directionality),14,'normal',ink);
    dest=fullfile(dataDir,['Control_Fly30_original_reorientation_' tag]);
    print(fig,[dest '.svg'],'-dsvg','-painters');
    print(fig,[dest '.png'],'-dpng','-r600');
    set(fig,'Visible','on');savefig(fig,[dest '.fig']);close(fig);
    fprintf('Saved original control example: %s\n',dest);
end
end

function note(fig,pos,str,fontSize,fontWeight,color)
annotation(fig,'textbox',pos,'String',str,'LineStyle','none','FontName','Arial', ...
    'FontSize',fontSize,'FontWeight',fontWeight,'Color',color,'Interpreter','none', ...
    'Margin',0,'FitBoxToText','off');
end
