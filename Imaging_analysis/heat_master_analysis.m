%% ANALYZE HEAT BARRIER DATA
% April 28, 2026

%% SETUP

baseFolder = '/Users/sophiarenauld/stacks/';  % change if needed
savepath   = '/Users/sophiarenauld/stacks/20260427-2_ves041_fly2_heat';

%% LOAD DATA

% --- FicTrac data ---
ftFile = uigetfile_full('*_ficTracData.mat', 'Select FicTrac Data', baseFolder);
if isempty(ftFile), return; end
load(ftFile, 'ftData');   % adjust if file contains different var names

% --- Virmen data ---
vmFile = uigetfile_full('*virmenData.mat', 'Select Virmen Data', baseFolder);
if isempty(vmFile), return; end
load(vmFile);   % adjust if needed

% --- Imaging data ---
imFile = uigetfile_full('*_or_z0z1_s_.mat', 'Select IMAGING Data', baseFolder);
if isempty(imFile), return; end
imaging_data = load(imFile);   % keep as struct

%% COMMON LABELS

% Time / signals
time_dff = daqVolumeTime.daqVolumeTime{1,1};          % imaging time
dff      = imaging_data.roi.ts;

% time_kin from vmData.Seconds (duration or cell)
t_raw = vmData.Seconds{1,1};  % often a duration vector
time_kin = seconds(t_raw);    % convert to double seconds
time_kin = time_kin(:);       % column
fwd_vel = smooth(vmData.vel_for{1,1});                         % forward velocity
rot_vel = smooth(vmData.vel_yaw{1,1});                         % yaw / rotational velocity

% Assume column vectors
time_dff = time_dff(:);
time_kin = time_kin(:);

tmax_kin = max(time_kin);

mask = time_dff <= tmax_kin;

time_dff_clipped = time_dff(mask);
dff_clipped      = dff(mask);   % if dff is same length as time_dff


%% QUICK CHECK PLOTS (OPTIONAL)

figure(1);
plot(time_kin, fwd_vel);
xlabel('Time (s)');
ylabel('Forward velocity (FicTrac)');
title('FicTrac forward velocity (smoothed)');

figure(2);
plot(ftData.pos_for{1,1}, ftData.pos_side{1,1});
xlabel('Forward position');
ylabel('Side position');
title('FicTrac trajectory');

figure; hold on;
plot(time_dff_clipped, dff_clipped*5, 'b');
plot(ftData.seconds{1,1}, ftData.fwVel{1,1}, 'k');
plot(ftData.seconds{1,1}, ftData.yawVel{1,1}, 'r');
xlabel('Time (s)');
ylabel('Signal');
legend({'dF/F (×5)','fwVel','yawVel'});
title('Imaging vs FicTrac signals');
grid on; hold off;

%% FLAT PATH: HEAT ON/OFF IN VM SPACE

figure(5);

x = vmData.pos_x{1,1};
y = vmData.pos_y{1,1};

% Heat may be stored as cell in table — extract to numeric vector
heat = vmData.heat;
if iscell(heat)
    heat = cell2mat(heat);
end
heat = heat(:);

heat_on  = heat > 0;
heat_off = ~heat_on;

hold on;
h_off = plot(x(heat_off), y(heat_off), 'k.', 'MarkerSize', 2);
h_on  = plot(x(heat_on),  y(heat_on),  'r.', 'MarkerSize', 2);
plot(x(1), y(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');  % start point

legend([h_off, h_on], {'Heat off', 'Heat on'}, 'Location', 'best');
xlabel('X position');
ylabel('Y position');
title('Flat path — heat on (red) vs off (black)');
axis equal; grid on;
hold off;
save_plot_with_title_as_filename('heat_colored', 'flatpath', savepath);


%% FIND CONTIGUOUS HEAT-ON SEGMENTS ON KINEMATIC TIME BASE

h = heat_on(:);

diff_array = diff([0; h; 0]);          % pad with zeros at ends
run_starts = find(diff_array ==  1);   % indices where heat turns on
run_ends   = find(diff_array == -1) - 1;

if numel(run_starts) ~= numel(run_ends)
    error('Mismatch between heat-on run starts and ends.');
end

% Optional: minimum run length (in samples) to keep
remove_if_shorter_than = 1;
for i = 1:numel(run_starts)
    run_length = run_ends(i) - run_starts(i) + 1;
    if run_length < remove_if_shorter_than
        h(run_starts(i):run_ends(i)) = 0;
    end
end

% Recompute runs after pruning short ones
diff_array = diff([0; h; 0]);
run_starts = find(diff_array ==  1);
run_ends   = find(diff_array == -1) - 1;

% Convert index runs to time intervals on kinematic time base
t_heat_start = time_kin(run_starts);
t_heat_end   = time_kin(run_ends);

%% PLOT: dF/F, FWD, ROT WITH HEAT-ON SHADING

figure;

% ---------- Subplot 1: dF/F ----------
subplot(3,1,1); hold on;

t       = time_dff_clipped;
sig_dff = dff_clipped;

yl = [min(sig_dff) max(sig_dff)];
if yl(1) == yl(2)
    yl = yl(1) + [-1 1];   % avoid zero-height y-limit
end

for i = 1:numel(t_heat_start)
    t1 = t_heat_start(i);
    t2 = t_heat_end(i);

    % Intersect with dF/F time range
    tt1 = max(t1, t(1));
    tt2 = min(t2, t(end));
    if tt2 > tt1
        patch([tt1 tt2 tt2 tt1], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 .8 .8], ...   % light red-ish
              'EdgeColor','none','FaceAlpha',0.35);
    end
end

plot(t, sig_dff, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('dF/F');
title('dF/F (red = heat on)');
grid on; hold off;

% ---------- Subplot 2: Forward Velocity ----------
subplot(3,1,2); hold on;

tk      = time_kin;
sig_fwd = fwd_vel;

yl = [min(sig_fwd) max(sig_fwd)];
if yl(1) == yl(2)
    yl = yl(1) + [-1 1];
end

for i = 1:numel(t_heat_start)
    t1 = t_heat_start(i);
    t2 = t_heat_end(i);

    tt1 = max(t1, tk(1));
    tt2 = min(t2, tk(end));
    if tt2 > tt1
        patch([tt1 tt2 tt2 tt1], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 .8 .8], ...
              'EdgeColor','none','FaceAlpha',0.35);
    end
end

plot(tk, sig_fwd, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
title('Forward Velocity');
grid on;
ylim([-5 20]);   % adjust as needed
hold off;

% ---------- Subplot 3: Rotational Velocity ----------
subplot(3,1,3); hold on;

sig_rot = rot_vel;
yl = [min(sig_rot) max(sig_rot)];
if yl(1) == yl(2)
    yl = yl(1) + [-1 1];
end

for i = 1:numel(t_heat_start)
    t1 = t_heat_start(i);
    t2 = t_heat_end(i);

    tt1 = max(t1, tk(1));
    tt2 = min(t2, tk(end));
    if tt2 > tt1
        patch([tt1 tt2 tt2 tt1], ...
              [yl(1) yl(1) yl(2) yl(2)], [1 .8 .8], ...
              'EdgeColor','none','FaceAlpha',0.35);
    end
end

plot(tk, sig_rot, 'r', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Rotational Velocity (deg/s)');
title('Rotational Velocity');
grid on; hold off;

% --- Link x-axes of all subplots and add title ---
linkaxes(findall(gcf, 'Type', 'axes'), 'x');
sgtitle('dF/F, Forward, and Rotational Velocity (shaded = heat on)');


% --- Save figure ---
if ~isempty(savepath)
    save_plot_with_title_as_filename('fwd_rot', 'dff_heat_shaded', savepath);
end

%% 
% --- Virmen position & time ---
x_vm = vmData.pos_x{1,1}(:);
y_vm = vmData.pos_y{1,1}(:);

t_vm = seconds(vmData.Seconds{1,1});  % duration -> double (sec)
t_vm = t_vm(:);

% --- Imaging dF/F and time (already clipped to overlap) ---
t_dff = time_dff_clipped(:);          % same length as dff_clipped
dff_vec = dff_clipped(:);

% --- Interpolate dF/F onto Virmen time base ---
% Use NaN outside overlap to avoid extrapolation artifacts
dff_on_vm = interp1(t_dff, dff_vec, t_vm, 'linear', NaN);

% Optional: remove NaNs so scatter only shows valid points
valid = ~isnan(dff_on_vm);

x_plot   = x_vm(valid);
y_plot   = y_vm(valid);
c_plot   = dff_on_vm(valid);

% --- Scatter plot: pos_x, pos_y colored by dF/F ---
figure;
scatter(x_plot, y_plot, 10, c_plot, 'filled', 'MarkerEdgeColor', 'none');
axis equal;
colormap(jet);      % or hot, parula, etc.
colorbar;
xlabel('X position');
ylabel('Y position');
title('Virmen trajectory colored by dF/F');
grid on;

save_plot_with_title_as_filename('dff_colored', 'flatpath', savepath);
%% HELPER FUNCTION: safe uigetfile + full path
function fullFilePath = uigetfile_full(pattern, prompt, baseFolder)
    [fileName, pathName] = uigetfile(pattern, prompt, baseFolder);
    if isequal(fileName, 0)
        disp(['User canceled: ' prompt]);
        fullFilePath = [];
    else
        fullFilePath = fullfile(pathName, fileName);
        disp(['Loaded: ' fullFilePath]);
    end
end