function [daq, triggerIdx, rho, Meno_chunks, not_Meno_chunks,ts_rm] = SegmentMenovsNotMeno_2p(daq, savepath, window, minVel,highThres,lowThres)
% Performs index by index assignment into menotaxing category based on
% calculated rho value assigned to that index. Rho is calculated over a
% window of set length centered on the index including only data points
% when the fly was moving. 

% INPUTS

    % ts: trial's fictrac data
    % window: window to calculate idx rho over -start with 1s (ask Janki)
    % minVel: minimum total velocity value, below = not moving - in mm/s 3 (forward + sideways + rotational) was less than 0.67 rad s−1
    % highThres: threshold when crossed with a + slope indicates start of
    % meno bought ρ = 0.88 elena paper
    % lowThres: threshold when crossed with a - slope indicates end of meno
    % bought maybe as low as p=0.7

% Here we excluded datapoints where the fly’s cumulative speed (forward + sideways + rotational) was less than 0.67 rad s−1

% OUTPUTS

% triggerIdx: category of each index, boolean
% rho: rho values of each index, vector
% Meno_chunks: cell containing indices categorized as 'menotaxing', cell
% not_Meno_chunks: cell containing indices categorized as not 'menotaxing', cell
% bData_interp: downsampled behaviour data, used to assign correct category
% to 1000Hz sampled data in Group_Meno_notMeno_ephys.m

Meno_chunks = {}; 
not_Meno_chunks = {};

total_mov_mm = abs(daq.bfv) + abs(daq.bsv) + abs(daq.byv);        
no0vel_idx = find(total_mov_mm >= minVel); 


count = 1; 

% find jump idxes 
daq = compute_absolute_circular_diff_2p(daq);
daq = detect_local_peaks_2p(daq);
jump_idx = daq.jump_detected;

% calculate & assign rho value to each datapoint
mean_headingVectors = [];

all_idx =  no0vel_idx;
if isempty(jump_idx)
    noJump_idx = all_idx;
else
    noJump_idx = all_idx(~ismember(all_idx, jump_idx));
end

keep_idx = intersect(noJump_idx, no0vel_idx);
newSampRate = length(daq.t)/max(daq.t);

ts_rm = daq;
% Get the field names of the structure
fieldNames = fieldnames(daq);

% Loop through each field to check for emptiness, character arrays, and nested structs
for i = 1:length(fieldNames)
    field = fieldNames{i};  % Current field name
    fieldValue = daq.(field);  % Get the value of the current field
    
    % Check if the field is empty or if it's a character array or a nested struct
    if isempty(fieldValue) || ...
       (ischar(fieldValue)) || ... % For empty character arrays
       isstruct(fieldValue) || ... % For empty nested structs
       iscell(fieldValue)
        ts_rm = rmfield(ts_rm, field);  % Remove the field if the conditions are met
    else
        if ~(contains(field, 'supp'))
            ts_rm.(field) = daq.(field)(:,keep_idx);
        end
    end

end

window = window * newSampRate; 
if window > length(keep_idx)
    window = length(keep_idx); 
end
window = round(window);


total_points = length(ts_rm.vy);
mean_headingVectors = nan(2, total_points);
idx_windows = cell(total_points, 1);
count = 1;

for i = 1:total_points
    % Calculate window bounds
    half_window = floor(window/2);
    window_start = max(1, i - half_window);
    window_end = min(total_points, i + half_window);
    
    % Get angles for this window
    angles_flyFor = ts_rm.vy(window_start:window_end);
    
    % Calculate mean vector
    x = cos(angles_flyFor);
    y = sin(angles_flyFor);
    idx_windows{count} = window_start:window_end;
    mean_headingVectors(1,count) = sum(x)/length(x);
    mean_headingVectors(2,count) = sum(y)/length(y);
    count = count + 1;
end
rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 

% assign idx 0 or 1 values based on assigned rho values & thresholds
triggerIdxRho = schmittTrigger(rho',highThres,lowThres);

triggerIdx = zeros(size(triggerIdxRho)); 
triggerIdx(triggerIdxRho == 1) = 1; 


% connect meno chunks if non meno is shorter than x seconds (optional)
pastMenoIdx = 1;
pastNotMenoIdx = 1;
for i = 1:length(triggerIdx)
iValue = triggerIdx(i);
if i > 1
    if iValue == 1
        timeFromMeno = i - pastMenoIdx; 
        if timeFromMeno < 0.5 * newSampRate % last transition from menotaxis occured less than x seconds before
            triggerIdx(pastMenoIdx:i) = 1; 
        end
        pastMenoIdx = i; 
    elseif iValue == 0
        timeFromNotMeno = i - pastNotMenoIdx;
        if timeFromNotMeno < 0.5 * newSampRate
            triggerIdx(pastNotMenoIdx:i) = 0; 
        end
        pastNotMenoIdx = i; 
    end
end
end

% Group meno & not meno data chunks in arrays

chunk = []; 
Mcount = 1; 
nMcount = 1; 
for idx = 1:length(triggerIdx)
    if isempty(chunk)
        if ~isnan(triggerIdx(idx))
            chunk(1) = idx;
            chunkValue = triggerIdx(idx);
        end
    else
        if (triggerIdx(idx) == chunkValue || isnan(triggerIdx(idx))) && idx ~= length(triggerIdx)
            chunk = [chunk, idx];
        else
            if chunkValue == 1
                Meno_chunks{Mcount} = chunk;
                Mcount = Mcount + 1; 
                chunk = [];
            else
                not_Meno_chunks{nMcount} = chunk;
                nMcount = nMcount + 1; 
                chunk = [];
            end
        end
    end
end

% Create 2D path trajectory meno = red not meno = black

notNan_idx = find(~isnan(ts_rm.vy) & ~isnan(ts_rm.bfv) & ~isnan(ts_rm.bsv));
yawAngPos = rad2deg(ts_rm.vy(notNan_idx));
fwdAngVel = ts_rm.bfv;
slideAngVel = ts_rm.bsv(notNan_idx);
triggerIdx = triggerIdx(notNan_idx); 
transition = zeros(size(ts_rm.bsv));

% N = keep_idx';
% V = find(transitions == 1);
% 
% A = repmat(N,[1 length(V)]);
% [~,closestIndex] = min(abs(A-V'));
% 
% transition(closestIndex) = 1; 
% 
% transition = transition(notNan_idx); 
 
% conversion factor between degrees and mm
circum = 9 * pi; % circumference of ball, in mm
mmPerDeg = circum / 360; % mm per degree of ball

% position incorporating heading - as if fly were walking on x-y plane,
%  x-y coordinates at each time point
% start with fly at (0,0) and facing 0 deg
zeroedYawAngPos = yawAngPos - yawAngPos(1); 

% movement in x (in degrees) at each time point
xChangePos = (fwdAngVel ./ newSampRate) .* sind(zeroedYawAngPos) + ...
    (slideAngVel ./ newSampRate) .* sind(zeroedYawAngPos + 90); 

% x position in mm (i.e. x-coordinate of fly's position at each time 
%  point), starts at 0
xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
minX = min(xPos);
maxX = max(xPos);

% movement in y (in degrees) at each time point
yChangePos = (fwdAngVel ./ newSampRate) .* cosd(zeroedYawAngPos) + ...
    (slideAngVel ./ newSampRate) .* cosd(zeroedYawAngPos + 90);

% y position in mm (i.e. y-coordinate of fly's position at each timepoint), starts at 0
yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
minY = min(yPos);
maxY = max(yPos);

% summary plots
time = ts_rm.t;
nMenoTime = time(triggerIdx == 0);
MenoTime = time(triggerIdx ==1); 
timeStart = time(1); 
timeEnd = max(time); 
nRho = rho(triggerIdx == 0); 

mAngle = -ts_rm.vy(triggerIdx ==1);
mRho = rho(triggerIdx == 1); 

figure(77);clf;
set(gcf,'color','w','renderer','painters')
h(1) =  subplot(2,1,1);
hold on
a = plot(daq.t, daq.vy,'k');
b = plot(MenoTime(MenoTime > timeStart & MenoTime < timeEnd),mAngle(MenoTime > timeStart & MenoTime < timeEnd),'r');
try
    b.XData(abs(diff(b.XData)) > 20) = nan;
catch
end
ylabel('HD')
ylim([0,360])
xlim([timeStart,timeEnd])
box off
h(2) = subplot(2,1,2);
hold on
c = plot(nMenoTime(nMenoTime > timeStart & nMenoTime < timeEnd),nRho(nMenoTime > timeStart & nMenoTime < timeEnd),'k');
try
    % remove lines connecting cue angle transitions from 0 to 180 & jumps
    c.XData(abs(diff(c.XData)) > 10) = nan;
catch
end
d = plot(MenoTime(MenoTime > timeStart & MenoTime < timeEnd),mRho(MenoTime > timeStart & MenoTime < timeEnd),'r');
try
    d.XData(abs(diff(d.XData)) > 10) = nan;
catch
end
ylabel('rho')
box off
xlim([timeStart,timeEnd])
linkaxes(h,'x')
save_plot_with_title_as_filename('rho', 'heading', savepath);

% plot the 2D path trajectory
figure(88);clf;
patch('XData',ts_rm.px(triggerIdx == 0),'YData',ts_rm.py(triggerIdx == 0),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
hold on
patch('XData',ts_rm.px(triggerIdx == 1),'YData',ts_rm.py(triggerIdx == 1),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
patch('XData',ts_rm.px(1),'YData',ts_rm.py(1),'EdgeColor','g','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',20);
patch('XData',ts_rm.px(transition == 1),'YData',ts_rm.py(transition == 1),'EdgeColor','b','FaceColor','none','LineStyle','none','Marker','o', 'MarkerSize',7);
set(gcf,'color','w');
xlabel('mm')
xlim([min(min(ts_rm.px),min(ts_rm.py)),max(max(ts_rm.px),max(ts_rm.py))])
ylim([min(min(ts_rm.px),min(ts_rm.py)),max(max(ts_rm.px),max(ts_rm.py))])
% Add the jump detected points
patch('XData', daq.px(daq.jump_detected == 1), ...
      'YData', daq.py(daq.jump_detected == 1), ...
      'EdgeColor', 'bl', ...
      'FaceColor', 'none', ...
      'LineStyle', 'none', ...
      'Marker', '.', ...
      'MarkerSize', 10);
title(['window: ',num2str(window/newSampRate),' highThres: ', num2str(highThres),' lowThres: ', num2str(lowThres)], 'Interpreter', 'none')
save_plot_with_title_as_filename('menotaxis', 'path', savepath);

% % plot the 2D path trajectory - ball
% figure(89);clf;
% patch('XData',ts_rm.px_ball(triggerIdx == 0),'YData',ts_rm.py_ball(triggerIdx == 0),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
% hold on
% patch('XData',ts_rm.px_ball(triggerIdx == 1),'YData',ts_rm.py_ball(triggerIdx == 1),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
% patch('XData',ts_rm.px_ball(1),'YData',ts_rm.py_ball(1),'EdgeColor','g','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',20);
% patch('XData',ts_rm.px_ball(transition == 1),'YData',ts_rm.py_ball(transition == 1),'EdgeColor','b','FaceColor','none','LineStyle','none','Marker','o', 'MarkerSize',7);
% set(gcf,'color','w');
% xlabel('mm')
% xlim([min(min(ts_rm.px_ball),min(ts_rm.py_ball)),max(max(ts_rm.px_ball),max(ts_rm.py_ball))])
% ylim([min(min(ts_rm.px_ball),min(ts_rm.py_ball)),max(max(ts_rm.px_ball),max(ts_rm.py_ball))])
% % Add the jump detected points
% patch('XData', daq.px(daq.jump_detected == 1), ...
%       'YData', daq.py(daq.jump_detected == 1), ...
%       'EdgeColor', 'bl', ...
%       'FaceColor', 'none', ...
%       'LineStyle', 'none', ...
%       'Marker', '.', ...
%       'MarkerSize', 10);
% title(['window: ',num2str(window/newSampRate),' highThres: ', num2str(highThres),' lowThres: ', num2str(lowThres)], 'Interpreter', 'none')
% save_plot_with_title_as_filename('menotaxis_ball', 'path', savepath);

end