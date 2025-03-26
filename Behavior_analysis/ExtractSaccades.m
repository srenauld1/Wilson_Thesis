function [triggerIdx, rho, Saccade_turns, no_Saccade, bData_interp] = ExtractSaccades(exptData, window, minVelforward, minVelyaw, Thres)
% Performs index by index selection to pick out saccadic turns using a
% calculated rho value assigned to that index. Rho is calculated over a
% window of set length centered on the index including only data points
% when the fly was moving. 

% saccadic turns are defined as short periods of high rotational velocity
% in one direction while the fly is moving

% INPUTS

    % exptData: trial's fictrac data
    % window: window to calculate idx rho over starting with 0.5s to pick out
    % short bursts of turns
    % minVelforward: minimum forward velocity value, below = not moving -
    % in mm/s 3
    % minVelyaw: minimum rotational velocity value, below = not turning -
    % 10 deg/s (its inverted in the code???)
    % Thres: threshold indicates rho value at which a turn occurs 0.5
    


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

%downsample time and optoopto and g4
exptData.t = downsample(exptData.t, 30);
exptData.optoStim = downsample(exptData.optoStim, 30);

% turn heading from upwrapped to wrapped
heading = exptData.headingPosition;
heading_wrapped = mod(heading, 360);
exptData.headingPosition = heading_wrapped;

vq = exptData.t(1):1/newSampRate:exptData.t(end);
fn = fieldnames(exptData);
% if istable(exptData)
%     fn = fn(1:end - 3); 
% end
for field = 1:numel(fn)
    if isempty(regexp(fn{field},'disp'))
        bData_interp.(fn{field}) = interp1(exptData.t,exptData.(fn{field}),vq,'linear');
    end
end  

        % threshold for movement and high rotational speed
no0vel_idx = find(abs(bData_interp.forwardVelocity) >= minVelforward & abs(bData_interp.angularVelocity) >= minVelyaw); 


count = 1; 

% find jump idxes 
% [jump_array, transitions, ~] = detect_jumps_ephys(bData_interp.frY, 5,5, 60);
% jump_idx = [];
% for jump = 1:size(jump_array,1)
%     jump_idx = [jump_idx , jump_array(jump,2):jump_array(jump,3)];
% end


% calculate & assign rho value to each datapoint
mean_headingVectors = [];

% all_idx =  no0vel_idx;
% if isempty(jump_idx)
%     noJump_idx = all_idx;
% else
%     noJump_idx = all_idx(~ismember(all_idx, jump_idx));
% end

keep_idx = no0vel_idx;

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
%hight)
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

% Create 2D path trajectory saccade = red not saccade = black

notNan_idx = find(~isnan(bData_interp.headingPosition) & ~isnan(bData_interp.forwardVelocity) & ~isnan(bData_interp.sidewaysVelocity));
yawAngPos = bData_interp.headingPosition(notNan_idx);
fwdAngVel = bData_interp.forwardVelocity(notNan_idx);
slideAngVel = bData_interp.sidewaysVelocity(notNan_idx);
triggerIdx = triggerIdx(notNan_idx); 
transition = zeros(size(bData_interp.sidewaysVelocity));

% N = keep_idx';
% V = find(transitions == 1);
% 
% A = repmat(N,[1 length(V)]);
% [~,closestIndex] = min(abs(A-V'));
% 
% transition(closestIndex) = 1; 
% 
% transition = transition(notNan_idx); 

newSampRate = 60; 
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


% plot the 2D path trajectory colored by forward velocity 
cmap=jet(100);
fwd_norm = round((exptData.forwardVelocity-min(exptData.forwardVelocity))/(max(exptData.forwardVelocity)-min(exptData.forwardVelocity))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(xPos)-1
    plot(xPos(i:i+1), yPos(i:i+1), 'Color', cmap(fwd_norm(i), :))
end
colorbar
plot(xPos(1), yPos(1), 'r.', 'MarkerSize', 20)
title("Forward Velocity colored")


% plot the 2D path trajectory colored by rotational speed 
cmap=jet(100);
rot_norm = round((exptData.angularSpeed-min(exptData.angularSpeed))/(max(exptData.angularSpeed)-min(exptData.angularSpeed))*(length(cmap)-1)) +1;
figure
hold on
for i = 1:length(xPos)-1
    plot(xPos(i:i+1), yPos(i:i+1), 'Color', cmap(rot_norm(i), :))
end
colorbar
plot(xPos(1), yPos(1), 'r.', 'MarkerSize', 20)
title("Rotational Speed colored")

% plot the 2D path trajectory colored by rho 
cmap=jet(100);
rho_norm = round(rho*100) +1;
figure
hold on
for i = 1:length(xPos)-1
    plot(xPos(i:i+1), yPos(i:i+1), 'Color', cmap(rho_norm(i), :))
end
colorbar
plot(xPos(1), yPos(1), 'r.', 'MarkerSize', 20)
title("Rho Speed colored")
end