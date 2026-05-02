function code = heat_barrier_smr_test
% heat_barrier_smr — heat barrier only version
%
% Runs ONLY three trial types:
%   Type 1 : heat_barrier         — heat scales with distance to barrier
%   Type 2 : heat_barrier_control — no heat output; virtual heat recorded
%   Type 3 : accommodation        — closed-loop bar, no heat, no barrier
%
% Closed-loop bar, no bar jumps.
% Hardware/save structure from heat_basic_conditioning_wz
% (virmenInit/Run/End, vr.data arrays).
%
% Key parameters (adjust as needed):
%   travel_integration_time : seconds fly moves freely before barrier placed
%   heat_barrier_time       : seconds barrier stays active before trial ends
%   heat_barrier_radius     : distance (arena units) of the safe zone
%   n_repeats               : number of [heat, control] pairs in trial list

code.initialization = @initializationCodeFun;
code.runtime        = @runtimeCodeFun;
code.termination    = @terminationCodeFun;


% =========================================================================
%  INITIALIZATION
% =========================================================================
function vr = initializationCodeFun(vr)

vr = virmenInit(vr);

vr.currentWorld = 1;
vr.heat_on      = 1;

% ── Data arrays ───────────────────────────────────────────────────────────
% Pre-allocate to actual experiment duration to avoid array extension
n_pts = ceil(vr.trialDuration * 61);
vr.data.barJump               = nan(1, n_pts); % always 0 — kept for pipeline compatibility
vr.data.heat                  = nan(1, n_pts);
vr.data.currTt                = nan(1, n_pts);
vr.data.heat_barrier_distance = ones(1, n_pts) * 100000; % 100000 when barrier not active
vr.data.control_trial         = nan(1, n_pts);

% ── Live figure ───────────────────────────────────────────────────────────
vr.liveFig = figure('Position',[1200 200 700 800]);
subplot(3,2,1:4);
vr.h1 = plot(1:n_pts, nan(1,n_pts), 'k.');   % no heat
hold on;
vr.h2 = plot(1:n_pts, nan(1,n_pts), 'r.');   % heat on

% Barrier overlay — hidden until barrier is placed
vr.heat_barrier_line = plot([-100 100],[0 0],'r-', ...
    'LineWidth',2,'Visible','off');
vr.heat_patch = patch( ...
    'XData',[-100 100 100 -100],'YData',[0 0 50 50], ...
    'FaceColor','r','FaceAlpha',0.2,'EdgeColor','none','Visible','off');
vr.travel_vector_arrow = quiver(0,0,1,0,'b-', ...
    'LineWidth',2,'MaxHeadSize',1,'Visible','off');
axis equal;
xlabel('Floor Width');
ylabel('Floor Length');
hold off;

subplot(3,2,6);
vr.plot_direct = polarplot([0,nan],[0,1],'r');
ax = gca;
ax.ThetaZeroLocation = 'top';
title('Fly Heading');

% ── Infinite-cue setup (move_with_fly / move_with_y) ─────────────────────
vr.move_with_fly_vert_ind = {};
vr.move_with_fly_vert     = {};
vr.move_with_y_vert_ind   = {};
vr.move_with_y_vert       = {};
for k = 1:length(vr.worlds)
    fnames    = fieldnames(vr.worlds{k}.objects.indices);
    vert_ind1 = [];
    vert_ind2 = [];
    for fi = 1:length(fnames)
        if contains(fnames{fi},'move_with_fly')
            vert_ind1 = [vert_ind1, ...
                vr.worlds{k}.objects.vertices( ...
                    vr.worlds{k}.objects.indices.(fnames{fi}),1) : ...
                vr.worlds{k}.objects.vertices( ...
                    vr.worlds{k}.objects.indices.(fnames{fi}),2)];
        end
        if contains(fnames{fi},'move_with_y')
            vert_ind2 = [vert_ind2, ...
                vr.worlds{k}.objects.vertices( ...
                    vr.worlds{k}.objects.indices.(fnames{fi}),1) : ...
                vr.worlds{k}.objects.vertices( ...
                    vr.worlds{k}.objects.indices.(fnames{fi}),2)];
        end
    end
    vr.move_with_fly_vert_ind{k} = vert_ind1;
    vr.move_with_fly_vert{k}     = vr.worlds{k}.surface.vertices(:,vert_ind1);
    vr.move_with_y_vert_ind{k}   = vert_ind2;
    vr.move_with_y_vert{k}       = vr.worlds{k}.surface.vertices(:,vert_ind2);
end

% ── Trial type default struct ─────────────────────────────────────────────
vr.tt(1:1000) = struct( ...
    'ind',                        0, ...
    'tt_name',                    'base', ...
    'trial_duration_seconds',     60, ...
    'max_heat_level',             0, ...
    'world',                      vr.currentWorld, ...
    'travel_integration_time',    5, ...
    'heat_barrier_time',          20, ...
    'heat_barrier_radius',        80, ...
    'heat_barrier_distance_list', [0 80 1e9], ...
    'heat_barrier_distance_heat', [0 0 0]);

% ── Type 1: heat barrier — heat active ───────────────────────────────────
i = 1;
vr.tt(i).ind                        = i;
vr.tt(i).tt_name                    = 'heat_barrier';
vr.tt(i).trial_duration_seconds     = 60;   % hard cap; barrier logic ends trial earlier
vr.tt(i).travel_integration_time    = 5;    % sec fly moves freely before barrier placed
vr.tt(i).heat_barrier_time          = 20;   % sec barrier stays active
vr.tt(i).heat_barrier_radius        = 80;   % distance of safe zone ahead of barrier
vr.tt(i).max_heat_level             = eval(vr.exper.variables.high_heat);
vr.tt(i).heat_barrier_distance_list = [0,  80,  1e9];
vr.tt(i).heat_barrier_distance_heat = [eval(vr.exper.variables.high_heat), 0, 0];

% ── Type 2: heat barrier control — identical geometry, heat suppressed ────
% heat_barrier_distance_heat uses the SAME profile as type 1 so that
% control_trial_virtual_heat_level records what the fly would have felt.
% Heat output is zeroed in runtime via control_trial_no_heat flag.
i = 2;
vr.tt(i).ind                        = i;
vr.tt(i).tt_name                    = 'heat_barrier_control';
vr.tt(i).trial_duration_seconds     = 60;
vr.tt(i).travel_integration_time    = 5;
vr.tt(i).heat_barrier_time          = 20;
vr.tt(i).heat_barrier_radius        = 80;
vr.tt(i).max_heat_level             = 0;    % no actual heat
vr.tt(i).heat_barrier_distance_list = [0,  80,  1e9];
vr.tt(i).heat_barrier_distance_heat = [eval(vr.exper.variables.high_heat), 0, 0];

% ── Type 3: accommodation — closed-loop bar, no heat, no barrier ──────────
i = 3;
vr.tt(i).ind                        = i;
vr.tt(i).tt_name                    = 'accommodation';
vr.tt(i).trial_duration_seconds     = 120;
vr.tt(i).max_heat_level             = 0;    % no heat ever
vr.tt(i).travel_integration_time    = 999;  % effectively infinite — barrier never placed
vr.tt(i).heat_barrier_time          = 0;
vr.tt(i).heat_barrier_radius        = 0;
vr.tt(i).heat_barrier_distance_list = [0 1 1e9];
vr.tt(i).heat_barrier_distance_heat = [0 0 0];

vr.tt([vr.tt.ind]==0) = [];

% ── Trial list ────────────────────────────────────────────────────────────
n_repeats = 20;
n_accom   = 3;
if vr.accom == 1
    % Accommodation only — closed-loop bar, no heat, no barrier
    vr.trial_type_list = repmat(3, 1, 10);
else
    % Accommodation block then alternating heat barrier / control pairs
    vr.trial_type_list = [repmat(3, 1, n_accom), repmat([1 2], 1, n_repeats)];
end

vr.trial_start_tic    = tic;
vr.tN                 = 0;
vr.current_trial_type = vr.trial_type_list(1);

% ── vr field initialisation ───────────────────────────────────────────────
vr.heat_on                          = 1;
vr.heat_output                      = 0;
vr.heat_level                       = 0;
vr.speed_forward                    = 0;
vr.speed_smoothed                   = 0;
vr.smoothing_window_seconds         = 2;
vr.high_heat                        = eval(vr.exper.variables.high_heat);

% Required by virmenRun to keep bar in closed-loop mode
vr.jumpStart                        = 0;
vr.fixedBarAngle                    = 0;
vr.barJumpTic                       = tic;
vr.position(4)                      = 0;

% Barrier state
vr.heat_barrier_active              = 0;
vr.heat_barrier_set                 = 0;
vr.heat_barrier_point               = [0 0];
vr.heat_barrier_normal              = [0 0];
vr.distance_to_barrier              = 100000;
vr.trial_start_position             = vr.position;
vr.control_trial_no_heat            = 0;
vr.control_trial_virtual_heat_level = 0;

disp(['Trial duration cap: ' num2str(vr.trialDuration) ' s']);
if vr.accom == 1
    disp('Accommodation mode — no barrier trials');
else
    disp(['Total trials: ' num2str(numel(vr.trial_type_list)) ...
          '  (' num2str(n_accom) ' accom + ' ...
          num2str(n_repeats) ' heat + ' num2str(n_repeats) ' control)']);
end


% =========================================================================
%  RUNTIME
% =========================================================================
function vr = runtimeCodeFun(vr)

% Keep bar in closed-loop mode — must be before virmenRun
vr.jumpStart = 0;
vr = virmenRun(vr);

% ── Infinite cue position update ─────────────────────────────────────────
if ~isempty(vr.move_with_fly_vert_ind{vr.currentWorld})
    vr.worlds{vr.currentWorld}.surface.vertices(1, ...
        vr.move_with_fly_vert_ind{vr.currentWorld}) = ...
        vr.move_with_fly_vert{vr.currentWorld}(1,:) + vr.position(1);
    vr.worlds{vr.currentWorld}.surface.vertices(2, ...
        vr.move_with_fly_vert_ind{vr.currentWorld}) = ...
        vr.move_with_fly_vert{vr.currentWorld}(2,:) + vr.position(2);
end
if ~isempty(vr.move_with_y_vert_ind{vr.currentWorld})
    vr.worlds{vr.currentWorld}.surface.vertices(2, ...
        vr.move_with_y_vert_ind{vr.currentWorld}) = ...
        vr.move_with_y_vert{vr.currentWorld}(2,:) + vr.position(2);
end

% ── Forward speed (monitoring only) ──────────────────────────────────────
speed_window = floor(vr.smoothing_window_seconds / vr.dt);
if size(vr.data.locomotion,1) > speed_window + 80
    vr.speed_forward = diff(vr.data.locomotion(end-speed_window:end,5)) / vr.dt;
else
    vr.speed_forward = 10;  % large default avoids spurious heat at startup
end
vr.speed_forward(vr.speed_forward < 0) = 0;
vr.speed_smoothed = mean(vr.speed_forward,'omitnan');

newTrial = 0;

% Skip first 80 iterations while fictrac stabilises
if vr.iterations == 80
    newTrial = 1;
    vr.tN    = 0;
end

% Trial end: hard duration cap
if toc(vr.trial_start_tic) > vr.tt(vr.current_trial_type).trial_duration_seconds
    newTrial = 1;
end

% Trial end: barrier expired or fly escaped boundary
if vr.heat_barrier_active
    time_since_barrier = toc(vr.trial_start_tic) - ...
        vr.tt(vr.current_trial_type).travel_integration_time;
    barrier_expired = time_since_barrier > ...
        vr.tt(vr.current_trial_type).heat_barrier_time;
    barrier_escaped = abs(vr.distance_to_barrier) > ...
        (vr.tt(vr.current_trial_type).heat_barrier_radius + 5);

    if barrier_escaped
        disp('Trial ended: fly escaped heat barrier radius');
    end
    if barrier_expired || barrier_escaped
        vr.heat_barrier_active = 0;
        vr.heat_barrier_set    = 0;
        newTrial               = 1;
    end
end

% ── New trial setup ───────────────────────────────────────────────────────
if newTrial == 1
    vr.tN = vr.tN + 1;
    if vr.tN > numel(vr.trial_type_list)
        vr.experimentEnded = 1;
        return
    end
    vr.current_trial_type = vr.trial_type_list(vr.tN);
    vr.trial_start_tic    = tic;
    disp('--- New trial ---');
    disp(struct2table(vr.tt(vr.current_trial_type),'asarray',1));

    % Reset all barrier state for the new trial
    vr.heat_barrier_active  = 0;
    vr.heat_barrier_set     = 0;
    vr.heat_barrier_point   = [0 0];
    vr.heat_barrier_normal  = [0 0];
    vr.distance_to_barrier  = 100000;
    vr.trial_start_position = vr.position;

    % Type 1 = heat active
    % Type 2 = control (no heat output, virtual heat still recorded)
    % Type 3 = accommodation (no heat, no barrier)
    vr.control_trial_no_heat = (vr.current_trial_type == 2) || ...
                               (vr.current_trial_type == 3);
    disp(['Control/accom trial (no heat output): ' ...
          num2str(vr.control_trial_no_heat)]);
end

% =========================================================================
%  HEAT LEVEL CALCULATION
% =========================================================================

if ~vr.heat_barrier_active

    % Phase 1: integration window — fly moves freely, no heat
    if toc(vr.trial_start_tic) > vr.tt(vr.current_trial_type).travel_integration_time

        % Direction from trial start → current position
        travel_vector = vr.position(1:2) - vr.trial_start_position(1:2);
        if norm(travel_vector) > 0
            travel_vector = travel_vector / norm(travel_vector);
        else
            travel_vector = [1 0];  % fallback if fly has not moved
        end

        % Place barrier ahead along travel direction
        vr.heat_barrier_point  = vr.position(1:2) + ...
            vr.tt(vr.current_trial_type).heat_barrier_radius * travel_vector;
        vr.heat_barrier_normal = travel_vector;
        vr.heat_barrier_active = 1;
        vr.heat_barrier_set    = 1;

        % FIX: compute initial distance immediately on placement iteration
        % so escape check on the NEXT iteration sees a valid value (~radius)
        % rather than the stale 100000 which caused instant trial termination
        fly_to_barrier         = vr.heat_barrier_point - vr.position(1:2);
        vr.distance_to_barrier = dot(fly_to_barrier, vr.heat_barrier_normal);

        disp('heat_barrier_set at:');
        disp(vr.heat_barrier_point);
    end
    vr.heat_level = 0;

else

    % Phase 2: barrier active — update distance and compute heat
    fly_to_barrier         = vr.heat_barrier_point - vr.position(1:2);
    vr.distance_to_barrier = dot(fly_to_barrier, vr.heat_barrier_normal);

    vr.heat_level = interp1( ...
        vr.tt(vr.current_trial_type).heat_barrier_distance_list, ...
        vr.tt(vr.current_trial_type).heat_barrier_distance_heat, ...
        abs(vr.distance_to_barrier));

    % Guard against NaN or out-of-range interp1 output
    if isnan(vr.heat_level) || vr.heat_level < 0
        vr.heat_level = 0;
    end
    vr.heat_level = min(vr.heat_level, vr.high_heat);  % hard ceiling

end

% Record what heat would have been, then zero for control/accommodation
vr.control_trial_virtual_heat_level = vr.heat_level;
if vr.control_trial_no_heat
    vr.heat_level = 0;
end

% Write to DAQ
vr.heat_output = vr.heat_level * 5;

% ── Per-iteration data save ───────────────────────────────────────────────
vr.data.barJump(vr.iterations)  = 0;               % always 0 — no bar jumps
vr.data.heat(vr.iterations)     = vr.heat_level;
vr.data.currTt(vr.iterations)   = vr.current_trial_type;
vr.data.control_trial(vr.iterations) = vr.control_trial_no_heat;

% Save  10,000 for distance when barrier is not active (pre-barrier + accommodation)
% This keeps data clean for downstream analysis
if vr.heat_barrier_active
    vr.data.heat_barrier_distance(vr.iterations) = vr.distance_to_barrier;
else
    vr.data.heat_barrier_distance(vr.iterations) = vr.distance_to_barrier; % should be 10,000
end

% ── Live plot ─────────────────────────────────────────────────────────────
if mod(vr.iterations,20) == 0

    % Trajectory dot: red = heat on, black = off/control/accommodation
    if abs(vr.heat_output) > 0.2
        vr.h2.XData(vr.iterations) = vr.position(1);
        vr.h2.YData(vr.iterations) = vr.position(2);
    else
        vr.h1.XData(vr.iterations) = vr.position(1);
        vr.h1.YData(vr.iterations) = vr.position(2);
    end

    % Barrier overlay
    if vr.heat_barrier_active
        barrier_length = 100;
        perp_vec   = [-vr.heat_barrier_normal(2), vr.heat_barrier_normal(1)];
        line_start = vr.heat_barrier_point + barrier_length * perp_vec;
        line_end   = vr.heat_barrier_point - barrier_length * perp_vec;

        set(vr.heat_barrier_line, ...
            'XData', [line_start(1) line_end(1)], ...
            'YData', [line_start(2) line_end(2)], ...
            'Visible','on');

        patch_depth = -vr.tt(vr.current_trial_type).heat_barrier_radius;
        patch_x = [line_start(1), line_end(1), ...
                    line_end(1)   + patch_depth * vr.heat_barrier_normal(1), ...
                    line_start(1) + patch_depth * vr.heat_barrier_normal(1)];
        patch_y = [line_start(2), line_end(2), ...
                    line_end(2)   + patch_depth * vr.heat_barrier_normal(2), ...
                    line_start(2) + patch_depth * vr.heat_barrier_normal(2)];
        set(vr.heat_patch,'XData',patch_x,'YData',patch_y,'Visible','on');

        % Grey = control/accommodation, red = heat active
        if vr.control_trial_no_heat
            set(vr.heat_barrier_line,'Color',[0.7 0.7 0.7]);
            set(vr.heat_patch,'FaceColor',[0.7 0.7 0.7],'FaceAlpha',0.1);
        else
            set(vr.heat_barrier_line,'Color','r');
            set(vr.heat_patch,'FaceColor','r','FaceAlpha',0.2);
        end

        % Travel direction arrow from trial start position
        t_vec = vr.heat_barrier_normal * 20;
        set(vr.travel_vector_arrow, ...
            'XData', vr.trial_start_position(1), ...
            'YData', vr.trial_start_position(2), ...
            'UData', t_vec(1), 'VData', t_vec(2), ...
            'Visible','on');
    else
        set(vr.heat_barrier_line,   'Visible','off');
        set(vr.heat_patch,          'Visible','off');
        set(vr.travel_vector_arrow, 'Visible','off');
    end

elseif mod(vr.iterations,20) == 5
    vr.plot_direct.ThetaData(2) = vr.position(4);
end

% Terminate when experiment duration cap is reached
if vr.timeElapsed >= vr.trialDuration
    vr.experimentEnded = true;
end


% =========================================================================
%  TERMINATION
% =========================================================================
function vr = terminationCodeFun(vr)

% Use currTt as reference — it is written every iteration
% and never intentionally holds a non-NaN sentinel value
last_iter = find(~isnan(vr.data.currTt), 1, 'last');

if isempty(last_iter)
    warning('heat_barrier_smr: no data written — experiment may have ended before iteration 80');
    last_iter = 0;
end

% Trim all arrays to identical length
% Removes trailing pre-allocated NaNs only
vr.data.barJump               = vr.data.barJump(1:last_iter);
vr.data.heat                  = vr.data.heat(1:last_iter);
vr.data.currTt                = vr.data.currTt(1:last_iter);
vr.data.heat_barrier_distance = vr.data.heat_barrier_distance(1:last_iter);
vr.data.control_trial         = vr.data.control_trial(1:last_iter);

virmenEnd(vr);