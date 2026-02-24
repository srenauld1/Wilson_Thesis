function a2p_data = openloop_bar_bidirectional_process(a2p_data, savepath, box, split)
%% calculate normalized dff if exists 2 rois
a=size(a2p_data.roi.ts);
if a(1)>1
    dff = a2p_data.roi.ts(1,:)-a2p_data.roi.ts(2,:) ;
else
    dff = a2p_data.roi.ts(1,:);
end
%% determine the kinematic variables
time = a2p_data.dq(1).t;
time_supp = a2p_data.dq(2).t;
if isfield(a2p_data.dq, 'vho')
    visual_heading = a2p_data.dq(1).vho;
else
    visual_heading = a2p_data.dq(1).vh;
end
if isfield(a2p_data.dq, 'vvyo')
    visual_pattern_velocity = a2p_data.dq(1).vvyo;
else
    visual_pattern_velocity = a2p_data.dq(1).vvy;
end
ball_forward_velocity = a2p_data.dq(1).bvf;
ball_forward_velocity_supp = a2p_data.dq(2).bvf;
ball_yaw_velocity = a2p_data.dq(1).bvy_deg;
ball_yaw_velocity_supp = a2p_data.dq(2).bvy_deg;
ball_side_velocity = a2p_data.dq(1).bvs;
ball_side_velocity_supp = a2p_data.dq(2).bvs;
totalspeed_supp = a2p_data.dq(2).totalspeed;
%% --- LABEL SWITCH ---
if split
    label_cw  = 'Back-to-Front';
    label_ccw = 'Front-to-Back';
else
    % changed bc bright bar and dark bar were opposite 12.1.2025
    label_cw  = 'Counterclockwise';
    label_ccw = 'Clockwise';
end

%% Plot the dff and pattern
    % Plot the dff and fwd
    figure;
    yyaxis left;
    plot(time, visual_heading, '-k');
    ylabel('pattern');  % Label for the left y-axis
    
    
    % Plot the second time series on the right y-axis
    yyaxis right;
    plot(time, dff, '-b');  % Plot second time series in red
    ylabel('dff');  % Label for the right y-axis
    
    % Add title and x-axis label
    xlabel('Time');
    title('Pattern and DFF');
    
    % Optional: Add legend
    legend('Pattern', 'dff');
    save_plot_with_title_as_filename('Pattern', 'dff', savepath);
    %% variables
    % frequency
    freq = time(end)/length(time);
    % Define the number of timepoints to include before and after motion
    len_sec = 2; % length in seconds of time to buffer with
    stimlen = 1; % how long the pattern moves for
    pre_motion_points = round(len_sec/freq);
    post_motion_points = round(len_sec/freq);
    remove_if_shorter_than = round(0.7/freq);

   
    %value for total speed (moving or no)
    still_cutoff = 1;
    
    %% using velocity to determine moving clockwise or counter
    % for 3 Hz

    
    positive_slope = zeros(size(visual_heading));
    positive_slope(visual_heading >= -3.1 & visual_pattern_velocity >= 0.02) = 1;
    % take care of less than 2pi case
    diff_array = diff([0, positive_slope, 0]);
    run_starts = find(diff_array == 1);
    run_ends = find(diff_array == -1) - 1;
    
    for i = 1:length(run_starts)
        run_length = run_ends(i) - run_starts(i) + 1;
        if run_length < remove_if_shorter_than
            positive_slope(run_starts(i):run_ends(i)) = 0;
        end
    end

    left_b2f = zeros(size(positive_slope));
    left_b2f(positive_slope == 1 & visual_heading <=0) = 1;
    right_f2b = zeros(size(positive_slope));
    right_f2b(positive_slope == 1 & visual_heading >=0) = 1;


    % neg slope
    negative_slope = zeros(size(visual_heading));
    negative_slope(visual_heading >= -3.1 & visual_pattern_velocity <= 0) = 1;
    if max(visual_pattern_velocity)<2.5
         negative_slope(visual_heading >= 0 & visual_heading <= 4) = 1;
    end
    diff_array = diff([0, negative_slope, 0]);
    run_starts = find(diff_array == 1);
    run_ends = find(diff_array == -1) - 1;
    
    for i = 1:length(run_starts)
        run_length = run_ends(i) - run_starts(i) + 1;
        if run_length < remove_if_shorter_than
            negative_slope(run_starts(i):run_ends(i)) = 0;
        end
    end


    % now ensure there are no tiny gaps
    vec = positive_slope; % Your original 1D row vector (1 x N)
    vec_padded = [1, vec, 1]; % Pad with ones at both ends to allow flanked checking at edges
    out = vec;
    
    z_idx = find(vec_padded == 0); % Indices of zeros (in the padded array)
    
    for i = 1:length(z_idx)
        idx = z_idx(i);
        % Check: if previous and next are 1s (in the padded array)
        if vec_padded(idx-1)==1 && vec_padded(idx+1)==1
            % It's a singleton, copy to original array
            out(idx-1)=1;
        elseif vec_padded(idx-1)==1 && vec_padded(idx+2)==1 && vec_padded(idx+1)==0
            % It's doublet (two zeros), check both are zero and both side are 1s
            out(idx-1:idx)=1;
        end
    end
    positive_slope = out;
    vec = negative_slope; % Your original 1D row vector (1 x N)
    vec_padded = [1, vec, 1]; % Pad with ones at both ends to allow flanked checking at edges
    out = vec;
    
    z_idx = find(vec_padded == 0); % Indices of zeros (in the padded array)
    
    for i = 1:length(z_idx)
        idx = z_idx(i);
        % Check: if previous and next are 1s (in the padded array)
        if vec_padded(idx-1)==1 && vec_padded(idx+1)==1
            % It's a singleton, copy to original array
            out(idx-1)=1;
        elseif vec_padded(idx-1)==1 && vec_padded(idx+2)==1 && vec_padded(idx+1)==0
            % It's doublet (two zeros), check both are zero and both side are 1s
            out(idx-1:idx)=1;
        end
    end
    negative_slope = out;

    % section out left and right b2f and f2b
    left_f2b = zeros(size(negative_slope));
    left_f2b(negative_slope == 1 & visual_heading <=0) = 1;
    right_b2f = zeros(size(negative_slope));
    right_b2f(negative_slope == 1 & visual_heading >=0) = 1;

    left_b2f = zeros(size(positive_slope));
    left_b2f(positive_slope == 1 & visual_heading <=0) = 1;
    right_f2b = zeros(size(positive_slope));
    right_f2b(positive_slope == 1 & visual_heading >=0) = 1;



    function plot_grouped_motion(time, dff, visual_heading, mask, pre_motion_points, post_motion_points, stimlen, len_sec, color, title_str)
    % Chop dff into segments via mask, add padding, and plot
    mask_padded = [0, mask(:)', 0];
    dmask = diff(mask_padded);
    seg_start = find(dmask==1);
    seg_end   = find(dmask==-1)-1;

    dff_segs = {};
    for s = 1:length(seg_start)
        idxs = (max(1,seg_start(s)-pre_motion_points)) : (min(seg_end(s)+post_motion_points, length(time)));
        dff_segs{end+1} = dff(idxs);
    end
    if isempty(dff_segs)
        title([title_str ' (No Segments)']);
        return
    end
    minlen = min(cellfun(@length, dff_segs));
    for s = 1:length(dff_segs)
        dff_segs{s} = dff_segs{s}(1:minlen); % Crop to minimum length
    end

    dff_mat = cell2mat(dff_segs'); % size [num_segments x minlen]
    ave_dff = mean(dff_mat,1);
    tvec = linspace(-len_sec, stimlen+len_sec, minlen);

    % Plot all single trials (light) + mean (bold)
    hold on
    for s = 1:size(dff_mat,1)
        plot(tvec, dff_mat(s,:), 'Color', [color 0.3], 'LineWidth', 0.5);
    end
    plot(tvec, ave_dff, 'Color', color, 'LineWidth', 2.5);

    xline(0, 'k--');
    xline(stimlen, 'k--');
    xlabel('Time (s)');
    ylabel('dF/F');
    title(title_str)
end

    figure;
    subplot(2,2,1);
    plot_grouped_motion(time, dff, visual_heading, left_f2b, pre_motion_points, post_motion_points, stimlen, len_sec, [0 0 1], 'Left Eye F2B');
    
    subplot(2,2,2);
    plot_grouped_motion(time, dff, visual_heading, left_b2f, pre_motion_points, post_motion_points, stimlen, len_sec, [0 0.8 0], 'Left Eye B2F');
    
    subplot(2,2,3);
    plot_grouped_motion(time, dff, visual_heading, right_f2b, pre_motion_points, post_motion_points, stimlen, len_sec, [1 0 0], 'Right Eye F2B');
    
    subplot(2,2,4);
    plot_grouped_motion(time, dff, visual_heading, right_b2f, pre_motion_points, post_motion_points, stimlen, len_sec, [1 0.7 0], 'Right Eye B2F');
    sgtitle('Grouped dFF by Presentation: Left/Right, F2B/B2F');
    
    save_plot_with_title_as_filename('bbar', 'f2b_b2f', savepath)


end

