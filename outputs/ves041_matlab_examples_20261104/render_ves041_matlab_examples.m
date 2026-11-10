function render_ves041_matlab_examples(dataDir)
% Recreate editable 25-second square figures from the bundled MAT data.
% Raw forward/yaw fields are sampled at 1 kHz, without added smoothing.
% Heading and event detections are retained from the original analysis.
if nargin<1, dataDir=fileparts(mfilename('fullpath')); end
loaded=load(fullfile(dataDir,'ves041_example_timeseries.mat'),'cases');
for k=1:numel(loaded.cases)
    c=loaded.cases(k);
    fig=figure('Visible','off','Color','w','Units','pixels', ...
        'Position',[100 100 1050 1050],'Name',c.name,'NumberTitle','off');
    set(fig,'PaperUnits','inches','PaperPosition',[0 0 10.5 10.5], ...
        'PaperSize',[10.5 10.5],'PaperPositionMode','manual');
    purple=[150 61 206]/255; light=[240 225 250]/255;
    orange=[197 82 22]/255; ink=[32 33 36]/255; gray=[104 107 112]/255;
    label=sprintf('Fly %d  |  %s',c.fly,c.genotype);
    if c.event_id>0
        subtitle=sprintf('Event #%d  |  25-second window  |  full light-on pulse',c.event_id);
    else
        subtitle=sprintf('Light-on bout %d  |  %.1f-%.1f s  |  25-second window', ...
            c.pulse,c.light_intervals(1,1),c.light_intervals(1,2));
    end
    note(fig,[.075 .946 .88 .035],label,23,'bold',ink);
    note(fig,[.075 .912 .88 .025],subtitle,15,'normal',gray);
    annotation(fig,'rectangle',[.078 .879 .02 .015],'Color',light,'FaceColor',light);
    note(fig,[.108 .872 .12 .027],'Light on',14,'normal',ink);
    annotation(fig,'line',[.26 .29],[.886 .886],'Color',orange,'LineWidth',2);
    note(fig,[.30 .872 .17 .027],'Reorientation',14,'normal',ink);
    annotation(fig,'line',[.48 .51],[.886 .886],'Color',ink,'LineWidth',1.2);
    note(fig,[.52 .872 .25 .027],'Recorded behavior',14,'normal',ink);
    note(fig,[.075 .835 .88 .027], ...
        'Raw velocity fields | sampled at 1 kHz | no additional smoothing',13,'normal',gray);
    axesPositions=[.12 .65 .82 .15; .12 .41 .82 .15; .12 .17 .82 .15];
    ax=gobjects(1,3);
    values={c.forward,c.yaw,c.heading}; times={c.time,c.time,c.heading_time};
    hl=[floor(min(c.heading)/30)*30-10,ceil(max(c.heading)/30)*30+10];
    limits={[-20 50],[-1000 1000],hl};
    ticks={[-20 0 20 40],[-1000 -500 0 500 1000],linspace(hl(1)+10,hl(2)-10,4)};
    labels={'Forward velocity (mm/s)','Yaw velocity (deg/s)', ...
        'Heading change (deg; original processed heading)'};
    assert(all(c.forward>=-20 & c.forward<=50),'Forward axis clips the signal.');
    assert(all(abs(c.yaw)<=1000),'Yaw axis clips the signal.');
    for j=1:3
        ax(j)=axes(fig,'Units','normalized','Position',axesPositions(j,:), ...
            'FontName','Arial','FontSize',14,'Box','off','TickDir','out', ...
            'XLim',c.window,'YLim',limits{j},'YTick',ticks{j}, ...
            'XColor',gray,'YColor',gray,'LineWidth',.7,'Layer','top');
        hold(ax(j),'on');
        for q=1:size(c.light_intervals,1)
            ab=c.light_intervals(q,:);yl=limits{j};
            patch(ax(j),[ab(1) ab(2) ab(2) ab(1)],[yl(1) yl(1) yl(2) yl(2)], ...
                light,'EdgeColor','none','HandleVisibility','off');
        end
        yline(ax(j),0,'Color',[.86 .87 .88],'LineWidth',.5,'HandleVisibility','off');
        plot(ax(j),times{j},values{j},'Color',ink,'LineWidth',.85);
        for q=1:size(c.event_windows,1)
            ev=c.event_windows(q,:);sel=times{j}>=ev(1)&times{j}<=ev(2);
            plot(ax(j),times{j}(sel),values{j}(sel),'Color',orange,'LineWidth',1.8);
            xline(ax(j),ev(3),':','Color',orange,'LineWidth',.8,'HandleVisibility','off');
            if j==2
                text(ax(j),ev(3),limits{j}(2)+100,sprintf('#%d',ev(4)), ...
                    'Color',orange,'FontName','Arial','FontSize',14, ...
                    'FontWeight','bold','HorizontalAlignment','center','Clipping','off');
            end
        end
        note(fig,[.12 axesPositions(j,2)+.158 .82 .027],labels{j},15,'bold',ink);
        if j<3
            set(ax(j),'XTick',[]);
        else
            set(ax(j),'XTick',ceil(c.window(1)/5)*5:5:floor(c.window(2)/5)*5);
            xlabel(ax(j),'Time from recording start (s)','FontName','Arial','FontSize',14);
        end
        assert(isequal(ax(j).XLim,c.window) && diff(ax(j).XLim)==25);
    end
    linkaxes(ax,'x');
    if isempty(c.event_metrics)
        footer='No saved reorientation detections in this window.';
    else
        m=c.event_metrics(1,:);
        footer=sprintf('Net turn %.1f deg | net displacement %.2f mm | directionality %.2f',m);
    end
    note(fig,[.075 .035 .88 .04],footer,14,'normal',ink);
    dest=fullfile(dataDir,c.name);
    print(fig,[dest '.svg'],'-dsvg','-painters');
    print(fig,[dest '.png'],'-dpng','-r600');
    % Keep the FIG visible when subsequently opened in MATLAB.
    set(fig,'Visible','on');
    savefig(fig,[dest '.fig']);
    close(fig);
    fprintf('Saved editable FIG, SVG and 6300-pixel PNG: %s\n',dest);
end
end

function note(fig,pos,str,fontSize,fontWeight,color)
annotation(fig,'textbox',pos,'String',str,'LineStyle','none', ...
    'FontName','Arial','FontSize',fontSize,'FontWeight',fontWeight, ...
    'Color',color,'Interpreter','none','Margin',0,'FitBoxToText','off');
end
