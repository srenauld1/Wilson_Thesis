function results=plot_ves041_reorientation_threshold_sensitivity(folder)
% Compare post-VES041 behavior with paired pre-peak and pseudo controls.
% Columns in input matrices are mean forward speed, absolute net heading
% displacement, and path straightness for a 1 s window. Post is 0 to +1 s,
% paired pre is -2 to -1 s, and pseudo controls are stimulus-free windows
% sampled throughout the trial and kept >2 s from a VES041 peak.

if nargin<1 || isempty(folder)
    folder='/Users/sophiarenauld/stacks/ves041_behavior_quantification/spontaneous_reorientation_phase_space';
end
Post=readmatrix(fullfile(folder,'pooled_peak_windows_fwd_heading.csv'));
Pre=readmatrix(fullfile(folder,'pooled_prepeak_windows_fwd_heading.csv'));
Pseudo=readmatrix(fullfile(folder,'pooled_control_windows_fwd_heading.csv'));
headingThresholds=[5 10 15 20 30 45];straightThresholds=[.95 .9 .8 .7 .5];

% Require low forward speed in both definitions. This asks whether the
% combined behavioral state broadens when reorientation is defined more
% permissively, without converting ordinary fast turns into positives.
[hPost,hPre,hPseudo]=heading_curves(Post,Pre,Pseudo,headingThresholds);
[sPost,sPre,sPseudo]=straight_curves(Post,Pre,Pseudo,straightThresholds);

fig=figure('Visible','off','Color','w','Position',[80 80 1350 820]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
plot_prevalence(nexttile(tl),headingThresholds,100*hPost,100*hPre,100*hPseudo, ...
    'Minimum |heading change| (deg)','Low-speed reorientation prevalence');
plot_enrichment(nexttile(tl),headingThresholds,hPost./hPre,hPost./hPseudo, ...
    'Minimum |heading change| (deg)','Heading-threshold enrichment');
plot_prevalence(nexttile(tl),straightThresholds,100*sPost,100*sPre,100*sPseudo, ...
    'Maximum path straightness','Low-speed nonlinear-path prevalence');
plot_enrichment(nexttile(tl),straightThresholds,sPost./sPre,sPost./sPseudo, ...
    'Maximum path straightness','Straightness-threshold enrichment');
sgtitle(tl,'VES041 behavioral association: paired pre-peak versus trial-wide occupancy controls');
savefig(fig,fullfile(folder,'reorientation_threshold_sensitivity.fig'));
exportgraphics(fig,fullfile(folder,'reorientation_threshold_sensitivity.png'),'Resolution',250);close(fig);

results.heading=table(headingThresholds',hPost',hPre',hPseudo',(hPost./hPre)', ...
    (hPost./hPseudo)','VariableNames',{'thresholdDeg','postFraction','preFraction', ...
    'pseudoFraction','postVsPre','postVsPseudo'});
results.straightness=table(straightThresholds',sPost',sPre',sPseudo',(sPost./sPre)', ...
    (sPost./sPseudo)','VariableNames',{'maximumStraightness','postFraction','preFraction', ...
    'pseudoFraction','postVsPre','postVsPseudo'});
writetable(results.heading,fullfile(folder,'heading_threshold_sensitivity.csv'));
writetable(results.straightness,fullfile(folder,'straightness_threshold_sensitivity.csv'));
disp(results.heading);disp(results.straightness);
end

function [a,b,c]=heading_curves(A,B,C,t)
a=arrayfun(@(x)mean(A(:,1)<=2&A(:,2)>=x,'omitnan'),t);
b=arrayfun(@(x)mean(B(:,1)<=2&B(:,2)>=x,'omitnan'),t);
c=arrayfun(@(x)mean(C(:,1)<=2&C(:,2)>=x,'omitnan'),t);
end

function [a,b,c]=straight_curves(A,B,C,t)
a=arrayfun(@(x)mean(A(:,1)<=2&A(:,3)<=x,'omitnan'),t);
b=arrayfun(@(x)mean(B(:,1)<=2&B(:,3)<=x,'omitnan'),t);
c=arrayfun(@(x)mean(C(:,1)<=2&C(:,3)<=x,'omitnan'),t);
end

function plot_prevalence(ax,x,a,b,c,xlab,ttl)
plot(ax,x,a,'-o','LineWidth',2.2,'DisplayName','0 to +1 s after peak');hold(ax,'on');
plot(ax,x,b,'-o','LineWidth',2.2,'DisplayName','-2 to -1 s before peak');
plot(ax,x,c,'-o','LineWidth',2.2,'DisplayName','Trial-wide pseudo windows');
xlabel(ax,xlab);ylabel(ax,'Windows meeting definition (%)');title(ax,ttl);grid(ax,'on');legend(ax,'Location','best');
end

function plot_enrichment(ax,x,pre,pseudo,xlab,ttl)
plot(ax,x,pre,'-o','LineWidth',2.2,'DisplayName','Post / paired pre');hold(ax,'on');
plot(ax,x,pseudo,'-o','LineWidth',2.2,'DisplayName','Post / pseudo occupancy');
h=yline(ax,1,'k--','No enrichment');h.HandleVisibility='off';
xlabel(ax,xlab);ylabel(ax,'Fold enrichment');title(ax,ttl);grid(ax,'on');legend(ax,'Location','best');
end
