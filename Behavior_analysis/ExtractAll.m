function [triggerIdx, rho, Saccade_turns, no_Saccade, bData_interp] = ExtractAll(exptData, window, minVelforward, minVelyaw, Thres)
% Performs index by index selection to pick out saccadic turns using a
% calculated rho value assigned to that index. Rho is calculated over a
% window of set length centered on the index including only data points
% when the fly was moving. 

% saccadic turns are defined as short periods of high rotational velocity
% in one direction while the fly is moving

% INPUTS

    % exptData: trial's fictrac data
    % window: window to calculate idx rho over starting with 3 s to pick out
    % short bursts of turns
    % minVelforward: minimum forward velocity value, below = not moving -
    % in mm/s 1
    % minVelyaw: minimum  rotational velocity value, below = not turning -
    % 50 deg/s (its inverted in the code???)
    % Thres: threshold indicates rho value at which a turn occurs 0.3
    


% OUTPUTS

% triggerIdx: category of each index, boolean
% rho: rho values of each index, vector
% Saccade_turns: cell containing indices categorized as a saccade turn
% not_Saccade: cell containing indices categorized as not a saccade turn
% bData_interp: downsampled behaviour data, used to assign correct category
% to 1000Hz sampled data in Group_Meno_notMeno_ephys.m

Saccade_turns = {}; 
not_Saccade = {};

% ephys data had to be downsampled from 1000Hz to 60Hz due to speed & memory constraints
newSampRate = 60; 


%downsample time and opto and g4
exptData.t = downsample(exptData.t, 30);
exptData.optoStim = downsample(exptData.optoStim, 30);
exptData.g4displayXPos = downsample(exptData.g4displayXPos, 30);

% turn heading from upwrapped to wrapped
heading = exptData.headingPosition;
heading_wrapped = mod(heading, 360);
exptData.headingPosition = heading_wrapped;

% add in jumps
% find jump idxes 
exptData = compute_absolute_circular_diff(exptData);
exptData = detect_local_peaks(exptData);

% Define new sampling rate
newSampRate = 60; 


% Get field names of exptData for interpolation
fn = fieldnames(exptData);


bData_interp = struct();  % Initialize bData_interp as an empty structure
vq = linspace(min(exptData.t), max(exptData.t), length(exptData.t) / 30);  % New time vector

bData_interp.optoStim = interp1(exptData.t, exptData.optoStim, exptData.t, 'nearest');
for field = 1:numel(fn)
    if isempty(regexp(fn{field}, 'disp'))
        bData_interp.(fn{field}) = interp1(exptData.t, exptData.(fn{field}), exptData.t, 'nearest');
    end
end


if isfield(exptData, 't') && ~isempty(exptData.t)
    disp(['Length of exptData.t: ', num2str(length(exptData.t))]);  % Debugging: Print length of t field
else
    error('exptData.t is missing or empty!');
end


for field = 1:numel(fn)
    if isempty(regexp(fn{field}, 'disp')) || strcmp(fn{field}, 'g4displayXPos')
        bData_interp.(fn{field}) = interp1(exptData.t,exptData.(fn{field}),vq,'nearest');
    end
end  


% threshold for movement and high rotational speed
no0vel_idx = find(abs(bData_interp.forwardVelocity) >= minVelforward & abs(bData_interp.angularVelocity) >= minVelyaw); 


count = 1; 

% calculate & assign rho value to each datapoint
mean_headingVectors = [];


keep_idx = no0vel_idx;

% add in buffer around jump
% Find indices where ones are located
bData_interp = compute_absolute_circular_diff(bData_interp);
bData_interp = detect_local_peaks(bData_interp);
jump_idx = bData_interp.jump_detected;
one_indices = find(jump_idx == 1);

% Define the padding size on each side
buffer_before_and_after_jump = 2;
padding = buffer_before_and_after_jump*newSampRate;
block_size = 2 * padding + 1;

% Expand each 1 into a block of ones
for i = 1:length(one_indices)
    idx = one_indices(i);
    start_idx = max(1, idx - padding);
    end_idx = min(43200, idx + padding);
    jump_idx(start_idx:end_idx) = 1;
end

if isempty(jump_idx)
    noJump_idx = no0vel_idx;
else
    noJump_idx = no0vel_idx(~ismember(no0vel_idx, jump_idx));
end

keep_idx = noJump_idx;

% remove 0 velocity & jump idxes from data to evaluate
fn = fieldnames(bData_interp);
% if istable(bData_interp)
%     fn = fn(1:end - 3); 
% end
for field = 1:numel(fn)
        bData_interp.(fn{field}) = bData_interp.(fn{field})(keep_idx);
end 

window = window * newSampRate; 
if window > length(keep_idx)
    window = length(keep_idx); 
end
disp(bData_interp);  % Print entire structure to see if it's populated correctly



% slide window over remaining data idx by idx to calculate rho
for i = 1 - window/2:1:length(bData_interp.headingPosition) - window/2
    idx = i:i + window; 

    if idx(end) > length(bData_interp.headingPosition) && idx(1) < 1 
        idx = 1:1:length(bData_interp.headingPosition);
    elseif idx(end) > length(bData_interp.headingPosition)
        idx = idx(1):1:length(bData_interp.headingPosition);
    elseif idx(1) < 1 
        idx = 1:idx(end); 
    end

    angles_flyFor = bData_interp.headingPosition(idx);
    if ~isempty(angles_flyFor)
        x = cosd(angles_flyFor); 
        y = sind(angles_flyFor); 
        idx_windows{count,1} = idx;
        mean_headingVectors(1,count)= sum(x)/length(x); 
        mean_headingVectors(2,count)= sum(y)/length(y);
        count = count + 1;
    else
        mean_headingVectors(1,count)= nan; 
        mean_headingVectors(2,count)= nan; 
        idx_windows{count,1} = idx;
        slope(count) = nan;
        count = count + 1;
    end
end

rho = sqrt(mean_headingVectors(1,:).^2 + mean_headingVectors(2,:).^2); 

%calculate thres as a function of the rho in the experiment (unless too
%high)
Thres=prctile(rho, 5);

if Thres >0.5
       Thres = 0.5;
end

% assign idx 0 or 1 values based on assigned rho values & thresholds
triggerIdxRho = schmittTrigger(rho',Thres, Thres);

triggerIdx = zeros(size(triggerIdxRho)); 
triggerIdx(triggerIdxRho == 0) = 1; 



% Group saccade & not saccade data chunks in arrays

chunk = []; 
Saccade_count = 1; 
no_saccade_count = 1; 
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
                Saccade_turns{Saccade_count} = chunk;
                Saccade_count = Saccade_count + 1; 
                chunk = [];
            else
                not_Saccade{no_saccade_count} = chunk;
                no_saccade_count = no_saccade_count + 1; 
                chunk = [];
            end
        end
    end
end

if ~isfield(bData_interp, 'optoStim')
    error('optoStim field is missing in bData_interp.');
end
disp(['Unique optoStim values: ', num2str(unique(bData_interp.optoStim))]);


% Create 2D path trajectory saccade = red not saccade = black

notNan_idx = find(~isnan(bData_interp.headingPosition) & ~isnan(bData_interp.forwardVelocity) & ~isnan(bData_interp.sidewaysVelocity));
yawAngPos = bData_interp.headingPosition(notNan_idx);
fwdAngVel = bData_interp.forwardVelocity(notNan_idx);
slideAngVel = bData_interp.sidewaysVelocity(notNan_idx);
triggerIdx = triggerIdx(notNan_idx); 
transition = zeros(size(bData_interp.sidewaysVelocity));


% Ensure xPos and yPos are computed before plotting
newSampRate = 60; 
circum = 9 * pi; % circumference of ball, in mm
mmPerDeg = circum / 360; % mm per degree of ball

% Extract relevant movement data
yawAngPos = bData_interp.headingPosition;
fwdAngVel = bData_interp.forwardVelocity;
slideAngVel = bData_interp.sidewaysVelocity;

% Adjust heading position
zeroedYawAngPos = yawAngPos - yawAngPos(1); 

% Compute x and y position changes
xChangePos = (fwdAngVel ./ newSampRate) .* sind(zeroedYawAngPos) + ...
    (slideAngVel ./ newSampRate) .* sind(zeroedYawAngPos + 90);  
xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;

yChangePos = (fwdAngVel ./ newSampRate) .* cosd(zeroedYawAngPos) + ...
    (slideAngVel ./ newSampRate) .* cosd(zeroedYawAngPos + 90);
yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;

% x position in mm (i.e. x-coordinate of fly's position at each time 
%  point), starts at 0
xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
minX = min(xPos);
maxX = max(xPos);

% y position in mm (i.e. y-coordinate of fly's position at each timepoint), starts at 0
yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
minY = min(yPos);
maxY = max(yPos);

% summary plots
time = bData_interp.t;
nMenoTime = time(triggerIdx == 0);
MenoTime = time(triggerIdx ==1); 
timeStart = time(1); 
timeEnd = max(time); 
nRho = rho(triggerIdx == 0); 

mAngle = -bData_interp.headingPosition(triggerIdx ==1);
mRho = rho(triggerIdx == 1); 

figure(77);clf;
set(gcf,'color','w','renderer','painters')
h(1) =  subplot(2,1,1);
hold on
a = plot(exptData.t, exptData.headingPosition,'k');
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


% plot the 2D path trajectory
figure(88);clf;patch('XData',xPos(triggerIdx == 0),'YData',yPos(triggerIdx == 0),'EdgeColor','k','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
hold on
patch('XData',xPos(triggerIdx == 1),'YData',yPos(triggerIdx == 1),'EdgeColor','r','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',3);
patch('XData',xPos(1),'YData',yPos(1),'EdgeColor','g','FaceColor','none','LineStyle','none','Marker','.', 'MarkerSize',20);
patch('XData',xPos(transition == 1),'YData',yPos(transition == 1),'EdgeColor','b','FaceColor','none','LineStyle','none','Marker','o', 'MarkerSize',7);
set(gcf,'color','w');
xlabel('mm')
xlim([min(minX,minY),max(maxX, maxY)])
ylim([min(minX,minY),max(maxX, maxY)])

title(['window: ',num2str(window/newSampRate),' Thres: ', num2str(Thres)], 'Interpreter', 'none')


%% Identify optoStim ON and OFF indices
opto_on_idx = find(bData_interp.optoStim ~= 0);
opto_off_idx = find(bData_interp.optoStim == 0);


figure(99); clf; hold on;
title("2D Path Trajectory with OptoStim, Saccades, & Bar Jumps")

% Define colors
opto_on_color = [0 1 0]; % Green for Opto ON
opto_off_color = [0 0 0]; % Black for Opto OFF
saccade_color = [1 0 0]; % Red for Saccades
both_color = [0.5, 0, 0.5]; % Purple for both opto and saccades

for i = 1:length(xPos)-1
    if bData_interp.optoStim(i) ~= 0 && triggerIdx(i) == 1
        this_color = both_color;  % Purple for both opto and saccade
    elseif triggerIdx(i) == 1
        this_color = saccade_color;  % Red for saccade
    elseif bData_interp.optoStim(i) ~= 0
        this_color = opto_on_color;  % Green for opto ON
    else
        this_color = opto_off_color;  % Black for opto OFF
    end

    plot(xPos(i:i+1), yPos(i:i+1), 'Color', this_color, 'LineWidth', 1.5);
end

% Mark start position in green
plot(xPos(1), yPos(1), 'g.', 'MarkerSize', 20);
% Find indices where jumps occur in the interpolated data
jump_indices_interp = find(bData_interp.jump_detected == 1);

% Get corresponding x and y positions
jump_x = xPos(jump_indices_interp);
jump_y = yPos(jump_indices_interp);

% Overlay jump points on figure 99
plot(jump_x, jump_y, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); % Blue circles for jumps

% Update legend to include jumps
legend({'Both OptoStim & Saccade (Purple)', 'Saccade (Red)', 'OptoStim ON (Green)', 'No Saccade & Opto OFF (Black)', 'Start Position (Green)', 'Jumps (Blue)'}, 'Location', 'best');


xlabel('X Position (mm)');
ylabel('Y Position (mm)');
set(gcf,'color','w');
axis equal;

end 