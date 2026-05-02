function code = heat_barrier_nav_short_range_control
% linearTrack_v3   Code for the ViRMEn experiment linearTrack_v3.
%   code = linearTrack_v3   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)
%% TODO:
% what info needs to be stored for each experiment? System of mouse id and
% session ID is perhaps not needed?

global fictrac_string
fictrac_string = '';

% set whether in debug mode:
% currently these flags are not used
vr.debugMode = 0;
vr.imaging = 0; % what does this do?
vr.drawText = 0;
daqreset;

vr = commonInit(vr);

% names of variables (fields in vr) that will be saved on each iteration,
% followed by the number of elements in each of those.
% Each variable will be flattened and saved in a
% binary file in the order specified in saveVar.
vr.saveOnIter = {...
    % UNIVERSAL FIELDS
    'rawMovement',3;... % empty in this case
    'position',4;...
    'velocity',4;...
    'dp',4;...
    'iN',1;...
    'tN',1;...
    'isITI',1;...
    'isVisible',1;...
    'dt',1;...
    'sync1',1;... % sync 1 just whether the iteration is odd or even
    'sync2',1; ... % sync 2 changes randomly every second
    'collision',1;...
    'currentWorld',1;...
    'iterations',1;...
    'last_iter_position',4;...
    'keyPressed',1;...
    'timeElapsed',1;...
    % FICTRAC VARIABLES
    'ft_count',1;...
    'ft_dr_cam',3;...
    'ft_err',1;...
    'ft_dr_lab', 3;...
    'ft_r_cam', 3;...
    'ft_r_lab', 3;...
    'ft_posx', 1;...
    'ft_posy', 1;...
    'ft_heading', 1;...
    'ft_step_dir', 1;...
    'ft_step_mag', 1;...
    'ft_inty', 1;...
    'ft_intx', 1;...
    'ft_ts', 1;...
    'ft_seq', 1;...
    'ft_dt_ms', 1;...
    'ft_alt_ts', 1; ...
    'ft_iN',1; ...
    ... % custom fields
    'speed',1;...
    'heat_level',1; ...
    'output_vec',5;
    'bar_jump',1;
    'speed_forward',1;
    'speed_smoothed',1;
    'heat_barrier_time',1;
    'travel_integration_time',1;
    'heat_barrier_radius',1;
    'trial_start_position',4;
    'heat_barrier_active',1;
    'heat_barrier_distance',1;
    'time_in_trial',1;
    'heat_barier_point',2;
    'heat_barrier_normal',2;
    'control_trial_no_heat',1;
    'control_trial_virtual_heat_level',1;
    };

vr = initSaveOnIterFields(vr);
vr.isVisible = 1;
    
% initialize vr.trial
vr.trial(1:5000,1) = struct(...
    ...% the standard fields 
    'N',0,...
    'start',now(),...
    'trial_start_tic',tic,...
    'type',1);


%% initialize condition structure with default values
% condition is used to reset these fields in vr at the start and when a new
% trial is initialized
vr.tt(1:1000) = struct(...
    'ind',0,...
    'tt_name','base',...
    'trial_duration_seconds',60,...
    'condition_heading',0,...
    'condition_speed',0,...
    'condition_speed_type','linear_lookup_forward_smoothed',...
    'condition_heading_type','step_90deg',...
    'target_heading_deg',0,...
    'luminance_gain',1,...
    'min_heat_level',0,...
    'max_heat_level',eval(vr.exper.variables.high_heat),...
    'heat_level_list',[eval(vr.exper.variables.high_heat)/2 0 0],...
    'speed_list',[0 2 1000000],...
    'bar_jump_at_start',0,...
    'bar_jump_range_deg',[],...
    'bar_jump_list_deg',[-90 90],...
    'bar_jump_repeat_interval_seconds',0,...
    'open_loop_heat',0,...
    'open_loop_heat_trace',[],...
    'world',vr.currentWorld);

i = 1;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'bar_jump_60_sec_trial';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
vr.tt(i).max_heat_level = 0;
% vr.tt(i).trial_duration_seconds = 120; 

i = 2;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'condition_speed_bar_jump_60_sec_trial';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
% vr.tt(i).trial_duration_seconds = 120; 
vr.tt(i).condition_speed = 1;
vr.tt(i).condition_speed_type = 'linear_lookup_forward_smoothed';
vr.tt(i).min_heat_level = 0;
vr.tt(i).max_heat_level = eval(vr.exper.variables.high_heat);

i = 3;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'condition_speed_and_heading0_bar_jump_60_sec_trial';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
% vr.tt(i).trial_duration_seconds = 120; 
vr.tt(i).condition_speed = 1;
vr.tt(i).min_heat_level = 0;
vr.tt(i).max_heat_level = eval(vr.exper.variables.high_heat);
vr.tt(i).condition_heading = 1;
vr.tt(i).target_heading_deg = 0;

i = 4;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'condition_speed_and_heading180_bar_jump_60_sec_trial';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
% vr.tt(i).trial_duration_seconds = 120; 
vr.tt(i).condition_speed = 1;
vr.tt(i).min_heat_level = 0;
vr.tt(i).max_heat_level = eval(vr.exper.variables.high_heat);
vr.tt(i).condition_heading = 1;
vr.tt(i).target_heading_deg = 180;

i = 5;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'condition_speed_and_heading270_bar_jump_60_sec_trial';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
% vr.tt(i).trial_duration_seconds = 120; 
vr.tt(i).condition_speed = 1;
vr.tt(i).min_heat_level = 0;
vr.tt(i).max_heat_level = eval(vr.exper.variables.high_heat);
vr.tt(i).condition_heading = 1;
vr.tt(i).target_heading_deg = 270;

i = 6;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'condition_speed_and_heading90_bar_jump_60_sec_trial';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
% vr.tt(i).trial_duration_seconds = 120; 
vr.tt(i).condition_speed = 1;
vr.tt(i).min_heat_level = 0;
vr.tt(i).max_heat_level = eval(vr.exper.variables.high_heat);
vr.tt(i).condition_heading = 1;
vr.tt(i).target_heading_deg = 90;

i = 7;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'low_heat_level_bar_jump_60_sec_trial';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
% vr.tt(i).trial_duration_seconds = 120; 
vr.tt(i).condition_speed = 0;
vr.tt(i).min_heat_level = 0.1;
vr.tt(i).max_heat_level = eval(vr.exper.variables.high_heat);
vr.tt(i).condition_heading = 0;
vr.tt(i).target_heading_deg = 0;

% another conditon where heat does not come back on for a few seconds after
% stopping
i = 8;
vr.tt(i).ind = i;
vr.tt(i).tt_name = 'condition_speed_bar_jump_60_sec_trial_heat_smoothed';
vr.tt(i).bar_jump_at_start = 1;
vr.tt(i).bar_jump_list_deg = [-90 90];
vr.tt(i).bar_jump_repeat_interval_seconds = 0;
% vr.tt(i).trial_duration_seconds = 120; 
vr.tt(i).condition_speed = 1;
vr.tt(i).condition_speed_type = 'linear_lookup_forward_smoothed';
vr.tt(i).min_heat_level = 0;
vr.tt(i).max_heat_level = eval(vr.exper.variables.high_heat);

% remove all structures not used
vr.tt([vr.tt.ind]==0) = [];


%% trial_list

vr.trial_type_list = [1 1 5 5 5 5 5 1 7 1 7 1 6 6 6 6 6 1 7 1 7];
%vr.trial_type_list = [1 1 1 1 1 1 1 1 1 1 5 1 5 1 5 1 5 1 5 1 7 1 7 1 6 1 6 1 6 1 6 1 6 1 7 1 7];
%vr.trial_type_list = [1 7 1 8 8 8 8 8 8 8 8 8 8 1 7 1 7 1];
%vr.trial_type_list = [1 1 1 1 1 1 1 1 1 1 8 8 8 8 8 1 7 1 5 5 8 5 5 8 5 5 8 5 5 8 1 7 1 6 6 8 6 6 8 6 6 8 6 6 8 1 7 1];
% vr.trial_type_list = [1 7 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 1 1 7 1 1 1 7];

% vr.trial_type_list = [repmat(1,1,5) repmat(3,1,5) repmat(1,1,5) repmat(4,1,5) repmat(1,1,5) repmat(5,1,5) repmat(1,1,5)];
% vr.trial_type_list = [repmat(2,1,10) repmat(1,1,10)];

%[1 1 1 1 2 2 2 2 1 1 1 1 3 3 3 3 1 1 1 1 4 4 4 4 1 1 1 1 1 1];

vr.trial_start_tic = tic;
vr.tN = 0;
vr.current_trial_type = vr.trial_type_list(1);

%% Initialize fields in vr

vr.iN = 0;
vr.heat_level = 0;
vr.speed = 0;
vr.bar_jump = 0;
vr.speed_forward = 0;

vr.draw_distance = 10000;
vr.high_heat = eval(vr.exper.variables.high_heat);
vr.heat_levels = [0.3 0 0];
vr.speed_levels = [0 10 1000];

vr.target_heading = 0;
vr.condition_heading = 0;
vr.condition_speed = 0;
vr.speed_smoothed_dt = zeros(1,60*10);
vr.speed_smoothed_trace = zeros(1,60*10).*NaN;
vr.speed_smoothed = 0;
vr.smoothing_window_seconds = 2;

vr.output_vec = [0 0 0 0 0];
vr.last_output_vec = vr.output_vec;
vr.daqout.write(vr.output_vec);
vr.currentWorld = eval(vr.exper.variables.startingWorld);

% remove red and green from the world
vr = removeRedAndGreen(vr);
% dim the world a bit
vr.luminance_gain = 1;
vr.worlds{vr.currentWorld}.surface.colors = vr.worlds{vr.currentWorld}.surface.colors.*vr.luminance_gain;

vr.position(4) = 0;
vr.lastPosition = vr.position;

n_samp = 600;
figure; hold on;
for k = 1:4
    subplot(4,1,k);
    vr.sp(k) = plot(zeros(n_samp,1), zeros(n_samp,1),'k-');
    xlabel('time');
end
drawnow();

figure; 
vr.f2d = plot(zeros(n_samp,1), zeros(n_samp,1));
xlabel('x')
ylabel('y')

% Add heat barrier visualization
hold on
% Create line for heat barrier (will be updated with position/orientation)
vr.heat_barrier_line = plot([-100 100], [0 0], 'r-', 'LineWidth', 2);
% Create patch for heat gradient visualization (extending forward from barrier line)
x = [-100 100 100 -100];
y = [0 0 50 50];
vr.heat_patch = patch('XData', x, 'YData', y, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Create arrow to show integrated travel vector from start position
vr.travel_vector_arrow = quiver(0, 0, 1, 0, 'b-', 'LineWidth', 2, 'MaxHeadSize', 1);
hold off



vr.trial_start_tic = tic;
vr.trial_start_position = vr.position;
vr.last_bar_jump = tic;
vr.heat_barrier_active = 0;
vr.distance_to_barrier = 100000;
vr.time_in_trial = toc(vr.trial_start_tic);
vr.tN = 1;
vr.travel_integration_time = 5;
vr.heat_barrier_time = 20;
vr.heat_barrier_radius = 80;
vr.heat_barrier_distance_list = [0 vr.heat_barrier_radius 1000000000];
vr.heat_barrier_distance_heat = [eval(vr.exper.variables.high_heat) 0 0];
vr.heat_barrier_point = [0 0];
vr.heat_barrier_normal = [0 0];
vr.heat_barrier_set = 0;
vr.control_trial_no_heat = 0;
vr.control_trial_virtual_heat_level = 0;
vr = saveVr(vr);






% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

vr = commonRuntime(vr,'iterStart');

vr.bar_jump = 0;
% if toc(vr.last_bar_jump)>120
%     vr.last_bar_jump = tic;
%     deg_offset = randsample([-180 -90 90],1);
%     vr.position(4) = vr.position(4)+deg2rad(deg_offset);
%     vr.bar_jump = 1;
% end

vr.time_in_trial = toc(vr.trial_start_tic);

if vr.heat_barrier_active == 0
    if toc(vr.trial_start_tic) > vr.travel_integration_time
        % set up the heat barrier
        % position is a 1x4 vector x,y,z,theta
        travel_vector = vr.position(1:2) - vr.trial_start_position(1:2);
        
        % normalize travel vector
        travel_vector = travel_vector / norm(travel_vector);
        
        % calculate point 20 units ahead along travel direction
        heat_barrier_point = vr.position(1:2) + vr.heat_barrier_radius * travel_vector;
        
        % get perpendicular vector (rotate 90 degrees counterclockwise)
        %perpendicular_vector = [-travel_vector(2), travel_vector(1)];
        
        % store line parameters for heat barrier (point and normal vector)
        vr.heat_barrier_point = heat_barrier_point;
        vr.heat_barrier_normal = travel_vector;
        
        % compute signed distance from fly to heat barrier line
        % positive distance means fly is behind the line
        fly_to_barrier = vr.heat_barrier_point - vr.position(1:2);
        distance_to_barrier = dot(fly_to_barrier, vr.heat_barrier_normal);
        
        % activate heat barrier
        vr.heat_barrier_active = 1;
        vr.distance_to_barrier = distance_to_barrier;
        vr.heat_barrier_set = 1;


        % print that a heat barrier is set 
        disp('heat_barrier_set');
        disp(vr.heat_barrier_point);

    end
else
    % heat barrier is active

        % compute signed distance from fly to heat barrier line
        % positive distance means fly is behind the line
        if (abs(vr.distance_to_barrier)>(vr.heat_barrier_radius + 5))
            disp('trial aborted due to distance from heat barrier')
        end

        if (toc(vr.trial_start_tic) > vr.travel_integration_time + vr.heat_barrier_time) || (abs(vr.distance_to_barrier)>(vr.heat_barrier_radius + 5))
            % make the heat barrier disappear and start a new trial
            vr.heat_barrier_active = 0;
            vr.tN = vr.tN+1;
            vr.control_trial_no_heat = mod(vr.tN,2)==0;
            vr.trial_start_tic = tic;
            vr.trial_start_position = vr.position; 
        end
end

if vr.heat_barrier_set
    fly_to_barrier = vr.heat_barrier_point - vr.position(1:2);
    distance_to_barrier = dot(fly_to_barrier, vr.heat_barrier_normal);
    
    % activate heat barrier
    vr.distance_to_barrier = distance_to_barrier;
end


if vr.heat_barrier_active
    vr.heat_level = interp1(vr.heat_barrier_distance_list, vr.heat_barrier_distance_heat, abs(vr.distance_to_barrier));
else
    vr.heat_level = 0;
end

vr.control_trial_virtual_heat_level = vr.heat_level;
if vr.control_trial_no_heat
    vr.heat_level = 0;
end




% write to the daq
vr.output_vec = [vr.sync1 vr.sync2 0 0 vr.heat_level*5];
if max(abs(vr.output_vec-vr.last_output_vec))>0
    % if anything in output vec has changed, then write to daq
    vr.daqout.write(vr.output_vec);
    vr.last_output_vec = vr.output_vec;
end

if mod(vr.iterations,10)==0
% update the plot
y3 = vr.distance_to_barrier;
y3(y3>100000) = NaN;

yvar =  [wrapTo360(rad2deg(vr.position(4))), vr.target_heading, y3, vr.heat_level];
X = vr.sp(1).XData;
X = circshift(X,-1,2);
X(end) = vr.timeElapsed;

for k = 1:4
Y = vr.sp(k).YData;
Y = circshift(Y,-1,2);
Y(end) = yvar(k);
set(vr.sp(k),'YData',Y,'Xdata',X);
end

% update the 2d plot
X = vr.f2d.XData;
X = circshift(X,-1,2);
X(end) = vr.position(1);

Y = vr.f2d.YData;
Y = circshift(Y,-1,2);
Y(end) = vr.position(2);
set(vr.f2d,'YData',Y,'Xdata',X);
axis equal


% Update heat barrier visualization
if vr.heat_barrier_active
    % Calculate endpoints for barrier line (perpendicular to travel direction)
    barrier_length = 100; % Length of visualized barrier line
    perp_vector = [-vr.heat_barrier_normal(2), vr.heat_barrier_normal(1)];
    line_start = vr.heat_barrier_point + barrier_length * perp_vector;
    line_end = vr.heat_barrier_point - barrier_length * perp_vector;
    
    % Update barrier line
    set(vr.heat_barrier_line, 'XData', [line_start(1) line_end(1)], ...
                             'YData', [line_start(2) line_end(2)], ...
                             'Visible', 'on');
    
    % Update heat patch - extending forward from barrier line in direction of normal
    patch_depth = -vr.heat_barrier_radius; % Visual depth of heat zone
    patch_points_x = [line_start(1), line_end(1), ...
                     line_end(1) + patch_depth * vr.heat_barrier_normal(1), ...
                     line_start(1) + patch_depth * vr.heat_barrier_normal(1)];
    patch_points_y = [line_start(2), line_end(2), ...
                     line_end(2) + patch_depth * vr.heat_barrier_normal(2), ...
                     line_start(2) + patch_depth * vr.heat_barrier_normal(2)];
    
    set(vr.heat_patch, 'XData', patch_points_x, ...
                       'YData', patch_points_y, ...
                       'Visible', 'on');
    
    % Set color based on control trial status
    if vr.control_trial_no_heat
        set(vr.heat_barrier_line, 'Color', [0.7 0.7 0.7]); % Grey for control
        set(vr.heat_patch, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.1);
    else
        set(vr.heat_barrier_line, 'Color', 'r'); % Red for active
        set(vr.heat_patch, 'FaceColor', 'r', 'FaceAlpha', 0.2);
    end
    
    % Update travel vector visualization - from trial start position
    travel_vec = vr.heat_barrier_normal;
    if norm(travel_vec) > 0
        travel_vec = travel_vec / norm(travel_vec) * 20; % Normalize and scale to fixed length
        set(vr.travel_vector_arrow, 'XData', vr.trial_start_position(1), ...
                                   'YData', vr.trial_start_position(2), ...
                                   'UData', travel_vec(1), ...
                                   'VData', travel_vec(2), ...
                                   'Visible', 'on');
    end
else
    % Hide visualizations when barrier isn't set
    set(vr.heat_barrier_line, 'Visible', 'off');
    set(vr.heat_patch, 'Visible', 'off');
    set(vr.travel_vector_arrow, 'Visible', 'off');
end

% ylim([min(Y) max(Y)]);
% xlim([min(X) max(X)]);
% axis equal;

end



vr = commonRuntime(vr,'iterEnd');

% check to see if the experiment/session is over
if ((now-vr.sessionStart)*24*60)>eval(vr.exper.variables.experiment_max_duration_minutes)
    vr.experimentEnded = 1;
end




% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
vr = commonTermination(vr);