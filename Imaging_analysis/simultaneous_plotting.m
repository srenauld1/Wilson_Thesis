close all
clear all

baseFolder = '/Users/sophiarenauld/stacks/';

% Prompt the user to select a file from the base folder
[fileName, pathName] = uigetfile('fly_data_*.mat', 'Select Processed Behavior + imaging Data', baseFolder);
[fileName2, pathName2] = uigetfile('*cmrg_.mat', 'Select tiff', baseFolder);
% Check if the user clicked "Cancel"
if isequal(fileName, 0)
    disp('User canceled the file selection.');
else
    % Load the selected file into the workspace
    fullFilePath = fullfile(pathName, fileName);
    fullFilePath2 = fullfile(pathName2, fileName2);
    load(fullFilePath);
    load(fullFilePath2);
    disp(['User selected "', fullFilePath, '".']);
   
end

%% simultaneous plotting
% mov= stack;
% mov1 = squeeze(stack(:,:,1,:));
% mov2 = squeeze(stack(:,:,2,:));
% mov3 = squeeze(stack(:,:,3,:));
% % Then use above script
% trace = fly_data.ts{1,1}(1, :);% 1xT or Tx1 fluorescence vector
% 
% T = size(mov1, 3);
% 
% fig = figure;
% 
% % Use tiledlayout for custom arrangement (2 rows, 3 columns)
% tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % Top row: Large trace plot, span all 3 columns
% ax1 = nexttile([1 3]);
% hTrace = plot(fly_data.daq.t, trace, 'b');
% hold on;
% hBar = plot([1 1], ylim, 'r', 'LineWidth', 2); % vertical bar
% hold off;
% title('Fluorescence over time');
% xlabel('Frame (time)');
% ylabel('Fluorescence');
% 
% % Bottom row: 3 movie frames side by side
% ax2 = nexttile(4);
% hImg1 = imagesc(mov1(:,:,1)); axis image off; colormap(ax2, gray);
% title('Plane 1');
% 
% ax3 = nexttile(5);
% hImg2 = imagesc(mov2(:,:,1)); axis image off; colormap(ax3, gray);
% title('Plane 2');
% 
% ax4 = nexttile(6);
% hImg3 = imagesc(mov3(:,:,1)); axis image off; colormap(ax4, gray);
% title('Plane 3');
% 
% % Set color scale for all planes if you want
% caxis(ax2, [min(mov1(:)) max(mov1(:))]);
% caxis(ax3, [min(mov2(:)) max(mov2(:))]);
% caxis(ax4, [min(mov3(:)) max(mov3(:))]);
% 
% % Animation Loop
% for t = 1:T
%     set(hBar, 'XData', [t t]); % move vertical bar
%     set(hImg1, 'CData', mov1(:,:,t));
%     set(hImg2, 'CData', mov2(:,:,t));
%     set(hImg3, 'CData', mov3(:,:,t));
%     drawnow;
%     pause(0.03);
% end
% 
% %% now with averages
% ave_frames = 2;  % for example
% 
% [height, width, num_planes, num_frames] = size(mov);
% num_avged_frames = floor(num_frames / ave_frames);
% mov_avg = zeros(height, width, num_planes, num_avged_frames);
% 
% % Average each block:
% for k = 1:num_avged_frames
%     idx_start = (k-1)*ave_frames + 1;
%     idx_end = k*ave_frames;
%     mov_avg(:,:,:,k) = mean(mov(:,:,:,idx_start:idx_end), 4);
% end
% 
% % Average trace for completeness (optional, usually average like movie)
% if length(trace) == num_frames
%     trace_avg = zeros(1, num_avged_frames);
%     for k = 1:num_avged_frames
%         idx_start = (k-1)*ave_frames + 1;
%         idx_end = k*ave_frames;
%         trace_avg(k) = mean(trace(idx_start:idx_end));
%     end
% else
%     trace_avg = trace;
% end
% 
% T = num_avged_frames;
% 
% fig = figure;
% n_planes = num_planes;
% 
% % Use tiledlayout: 2 rows, max(n_planes,3) cols just in case you get many planes
% tl = tiledlayout(2, n_planes, 'TileSpacing', 'compact', 'Padding', 'compact');
% 
% % Top: Trace over all columns (span)
% ax1 = nexttile([1 n_planes]);
% hTrace = plot(trace_avg, 'b');
% hold on;
% hBar = plot([1 1], ylim, 'r', 'LineWidth', 2);
% hold off;
% title('Fluorescence over time');
% xlabel('Averaged Frame');
% ylabel('Fluorescence');
% 
% % Bottom: All planes, one in each column
% hImgs = gobjects(n_planes,1);
% cmin = min(mov_avg(:));
% cmax = max(mov_avg(:));
% for p = 1:n_planes
%     ax = nexttile(n_planes + p); % row 2, column p
%     hImgs(p) = imagesc(mov_avg(:,:,p,1)); 
%     axis image off; colormap(ax, gray); 
%     caxis([cmin cmax]);
%     title(['Plane ' num2str(p)]);
% end
% 
% % Animation loop
% for t = 1:T
%     set(hBar, 'XData', [t t]);
%     for p = 1:n_planes
%         set(hImgs(p), 'CData', mov_avg(:,:,p,t));
%     end
%     drawnow;
%     pause(0.03);
% end

%% variables
trace1 = fly_data.ts{1,1}(1, :); % Fluorescence vector [1 x T]
trace2 = fly_data.daq.bfv;       % Behavioral velocity [1 x T]
trace3 = fly_data.daq.t;         % Time vector [1 x T] (seconds)
mov = stack;                     % Movie [X, Y, n_planes, n_frames]

%% average variables
ave_frames = 2;  % frames per block average
[x, y, n_planes, n_frames] = size(mov);
n_blocks = floor(n_frames / ave_frames);

% Block average movie and traces
dFoF_avg = zeros(x, y, n_planes, n_blocks);
trace1_avg = zeros(n_blocks, 1);
trace2_avg = zeros(n_blocks, 1);
trace3_avg = zeros(n_blocks, 1); % mean time per block

for k = 1:n_blocks
    idx_start = (k-1)*ave_frames + 1;
    idx_end = k * ave_frames;
    idx = idx_start:idx_end;

    dFoF_avg(:, :, :, k) = mean(mov(:, :, :, idx), 4);
    trace1_avg(k) = mean(trace1(idx));
    trace2_avg(k) = mean(trace2(idx));
    trace3_avg(k) = mean(trace3(idx));
end

% Set image display contrast limits (robust)
valid_vals = dFoF_avg(isfinite(dFoF_avg));
lo = prctile(valid_vals, 1);
hi = prctile(valid_vals, 99);

% --- Plotting ---
fig = figure;
tl = tiledlayout(3, n_planes, 'TileSpacing', 'compact', 'Padding', 'compact');

% TOP ROW: Fluorescence, time on X, one trace (block-averaged)
ax1 = nexttile(1, [1 n_planes]);
hTrace1 = plot(trace3_avg, trace1_avg, 'b');
hold on;
hBar1 = plot([trace3_avg(1) trace3_avg(1)], ylim, 'r', 'LineWidth', 2);
hold off;
title('Trace 1 (Block-Averaged)');
xlabel('Time (s)');
ylabel('Fluorescence');

% MIDDLE ROW: Behavior, time on X
ax2 = nexttile(n_planes+1, [1 n_planes]);
hTrace2 = plot(trace3_avg, trace2_avg, 'g');
hold on;
hBar2 = plot([trace3_avg(1) trace3_avg(1)], ylim, 'r', 'LineWidth', 2);
hold off;
title('Trace 2 (Block-Averaged)');
xlabel('Time (s)');
ylabel('Behavior');

% BOTTOM ROW: Movies, one plane per tile
hImgs = gobjects(n_planes,1);
for p = 1:n_planes
    ax = nexttile(2*n_planes + p);
    hImgs(p) = imagesc(dFoF_avg(:,:,p,1));
    axis image off; colormap(ax, gray);
    caxis(ax, [lo hi]);
    title(['Plane ' num2str(p)]);
end
colorbar(ax);

% --- Animation loop: vertical line moves in time! ---
for t = 1:n_blocks
    set(hBar1, 'XData', [trace3_avg(t) trace3_avg(t)], 'YData', ylim(ax1));
    set(hBar2, 'XData', [trace3_avg(t) trace3_avg(t)], 'YData', ylim(ax2));
    for p = 1:n_planes
        set(hImgs(p), 'CData', dFoF_avg(:,:,p,t));
    end
    drawnow;
    pause(0.03);
end