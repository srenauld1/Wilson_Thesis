% spikerate_v_velocityHeat
% analysis function for generating a summary plot of spike rate versus
% directional velocity as heatmaps
%
% OUTPUTS:
% binnedData - binned spikerate and directional velocity data
%
% INPUTS:
% int_panelps - downsampled panel positions
% int_forward - downsampled forward velocities
% int_forward - downsampled forward velocities
% int_angular - downsampled angular velocities
% int_sideway - downsampled sideways velocities
% int_spikerate - downsampled spikerate
% optPlot - 0 no plot, 1 to plot
%
% ORIGINAL: 06/24/2022 - MC
% UPDATED:  03/14/2024 - MC added optional plotting, output variables
%

function binnedData = dff_v_velocityHeat(int_forward,int_angular,int_sideway,dff,optPlot)
%% settings
% set bin size
fs = 1;
as = 30;
ss = 0.25;

%% reshape and discretize the data set
% reshape each to single column
rsp_forward = reshape(int_forward,[],1);
rsp_angular = reshape(int_angular,[],1);
rsp_sideway = reshape(int_sideway,[],1);
rsp_dff = reshape(dff,[],1);

% discretize each variable
f_max = 10;
f_min = -10;
f_edge = f_min-fs/2:fs:f_max+fs/2; %bin edges
f_bins = f_min:fs:f_max; %bin labels (center)
disc_forward=discretize(rsp_forward,f_edge,f_bins);

a_max = 200;
a_edge = -a_max-as/2:as:a_max+as/2; %bin edges
a_bins = -a_max:as:a_max; %bin labels (center)
disc_angular=discretize(rsp_angular,a_edge,a_bins);

s_max = 2.5;
s_edge = -s_max-ss/2:ss:s_max+ss/2; %bin edges
s_bins = -s_max:ss:s_max; %bin labels (center)
disc_sideway=discretize(rsp_sideway,s_edge,s_bins);

rs = 0.05; %bin size
r_max = round(max(rsp_dff),3); %max
r_edge = -rs/2:rs:r_max+rs/2; %bin edges
r_bins = 0:rs:r_max; %bin labels (center)
disc_spikerate=discretize(rsp_dff,r_edge,r_bins);

% compile data and remove NaN values
trialData_pre = [disc_spikerate disc_forward disc_angular disc_sideway];
[rNan, ~] = find(isnan(trialData_pre));
trialData_pre(rNan,:)=[];

% convert to table
colNames = {'dff','Forward','Angular','Sideway'};
binnedData = array2table(trialData_pre,'VariableNames',colNames);

%% plot

if optPlot
    % initialize
    figure; set(gcf,'Position',[100 100 1000 500])
    tiledlayout(1,2,'TileSpacing','compact')
    
    % plot forward and angular
    nexttile
    h(1) = heatmap(binnedData,'Angular','Forward','ColorVariable','dff');
    h(1).YDisplayData = flipud(h(1).YDisplayData);
    h(1).CellLabelColor = 'None';
    h(1).GridVisible = 'off';
    colormap('spring');
    
    % plot forward and sideways
    nexttile
    h(2) = heatmap(binnedData,'Sideway','Forward','ColorVariable','dff');
    h(2).YDisplayData = flipud(h(2).YDisplayData);
    h(2).CellLabelColor = 'None';
    h(2).GridVisible = 'off';
    colormap('spring');
end

end

