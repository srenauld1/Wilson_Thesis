%% Assume you already have trialSummary (82x1 struct) in the workspace

% -------------------------------------------------------------------------
% 1. Fly IDs in the same order as rows of trialSummary
%    (these must match exactly what you listed)
% -------------------------------------------------------------------------
flyIDs = [ ...
     2  ,  3  ,  4  ,  5  ,  7  ,  8  , 12  , 14  , 15  , 16 , ...
    17  , 18  , 19  , 20  , 21  , 23  , 24  , 25  , 26  , 27 , ...
    29  , 31  , 34  , 35  , 39  , 40  , 44  , 46  , 50  , 52 , ...
    53  , 54  , 57  , 58  , 59  , 62  , 63  , 65  , 66  , 69 , ...
    70  , 72  , 73  , 74  , 75  , 77  , 79  , 80  , 81  , 84 , ...
    85  , 87  , 88  , 89  , 90  , 91  , 93  , 96  , 99  ,101 , ...
   100  ,102  ,103  ,105  ,107  ,108  ,110  ,111  ,113  ,114 , ...
   118  ,119  ,121  ,122  ,123  ,124  ,125  ,129  ,130  ,132 , ...
   133  ,137 ];

if numel(flyIDs) ~= numel(trialSummary)
    error('flyIDs length (%d) does not match trialSummary length (%d).', ...
          numel(flyIDs), numel(trialSummary));
end

% -------------------------------------------------------------------------
% 2. Define genotype fly-ID lists
% -------------------------------------------------------------------------
ves041_1_Chrim_IDs = [2, 7, 16, 19, 21, 26, 53, 58, 65, 70, 72, 74, 77, 79, 81, 85, 119, 130];
ves041_2_Chrim_IDs = [14, 18, 24, 29, 40, 44, 66, 73, 75, 80, 90, 100, 102, 105, 108, 110];
ves041_1_gfp_IDs   = [3, 4, 5, 8, 12, 15, 25, 27, 69, 99, 101, 103, 111, 114, 121, 124];
ves041_2_gfp_IDs   = [35, 52, 59, 63, 89, 93, 107, 113, 118, 122, 123, 125, 129, 132, 133, 137];
empty_Chrim_IDs    = [17, 20, 23, 31, 34, 39, 46, 50, 54, 57, 62, 84, 87, 88, 91, 96];

% -------------------------------------------------------------------------
% 3. Helper: subset trialSummary given a list of fly IDs
% -------------------------------------------------------------------------
subset_by_ids = @(id_list) trialSummary(ismember(flyIDs, id_list));

% -------------------------------------------------------------------------
% 4. Build genotype-specific structs (same fields as trialSummary)
% -------------------------------------------------------------------------
ves041_1_Chrim_trials = subset_by_ids(ves041_1_Chrim_IDs);
ves041_2_Chrim_trials = subset_by_ids(ves041_2_Chrim_IDs);
ves041_1_gfp_trials   = subset_by_ids(ves041_1_gfp_IDs);
ves041_2_gfp_trials   = subset_by_ids(ves041_2_gfp_IDs);
empty_Chrim_trials    = subset_by_ids(empty_Chrim_IDs);

% -------------------------------------------------------------------------
% 5. (Optional) check counts
% -------------------------------------------------------------------------
fprintf('ves041_1_Chrim: %d trials\n', numel(ves041_1_Chrim_trials));
fprintf('ves041_2_Chrim: %d trials\n', numel(ves041_2_Chrim_trials));
fprintf('ves041_1_gfp:   %d trials\n', numel(ves041_1_gfp_trials));
fprintf('ves041_2_gfp:   %d trials\n', numel(ves041_2_gfp_trials));
fprintf('empty_Chrim:    %d trials\n', numel(empty_Chrim_trials));

%% Choose one genotype struct
plot_genotype_opto(ves041_1_Chrim_trials, 'ves041_1_Chrim')
plot_genotype_opto(ves041_2_Chrim_trials, 'ves041_2_Chrim')
plot_genotype_opto(ves041_1_gfp_trials, 'ves041_1_gfp')
plot_genotype_opto(empty_Chrim_trials, 'empty_chrim')
plot_genotype_opto(ves041_2_gfp_trials, 'ves041_2_gfp')




function plot_genotype_opto(trialname, figurename)

    S = trialname;      % 18x1 struct with fields: id, mean_fwd_opto, mean_yaw_opto
    
    nFlies = numel(S);
    
    % Stack into matrices: [nFlies x nSamples]
    fwd_mat = vertcat(S.mean_fwd_opto);    % size: nFlies x nSamples
    yaw_mat = vertcat(S.mean_yaw_opto);    % size: nFlies x nSamples
    
    nSamples = size(fwd_mat, 2);
    t = 1:nSamples;                        % sample index; replace with time if you have dt
    
    mean_fwd = mean(fwd_mat, 1, 'omitnan');
    mean_yaw = mean(abs(yaw_mat), 1, 'omitnan');
    
    nSamples = size(fwd_mat, 2);   % or numel(S(1).mean_fwd_opto)
    
    t_sec = linspace(0, 3, nSamples);   % time from -1 to +3 s
    
    figure;
    
    %% Panel 1: Forward velocity during opto
    subplot(2,1,1); hold on;
    
    % All flies (thin, light lines)
    for i = 1:nFlies
        plot(t_sec, fwd_mat(i,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, ...
             'HandleVisibility','off');
    end
    
    % Across-fly mean (thick, dark line)
    plot(t_sec, mean_fwd, 'k', 'LineWidth', 2, 'DisplayName', 'Across-fly mean');
    xline(1, "LineWidth",2, "LineStyle","--")
    xline(2, "LineWidth",2, "LineStyle","--")
    
    xlabel('Samples (opto-aligned)');
    ylabel('Forward velocity (mean_fwd\_opto)');
    title('Opto-on forward velocity');
    grid on;
    legend('Location','best');
    
    hold off;
    
    %% Panel 2: Yaw velocity during opto
    subplot(2,1,2); hold on;
    
    for i = 1:nFlies
        plot(t_sec, abs(yaw_mat(i,:)), 'Color', [1 0.6 0.6], 'LineWidth', 0.8, ...
             'HandleVisibility','off');
    end
    
    plot(t_sec, mean_yaw, 'Color', [0.7 0 0], 'LineWidth', 2, 'DisplayName', 'Across-fly mean');
    xline(1,"LineWidth",2, "LineStyle","--")
    xline(2, "LineWidth",2, "LineStyle","--")
    ylim([0 100]);   % <-- lock rotation y-axis to [-100, 100]
    xlabel('Samples (opto-aligned)');
    ylabel('Yaw velocity (mean\_yaw\_opto)');
    title('Opto-on yaw velocity');
    grid on;
    legend('Location','best');
    
    hold off;
    save_plot_with_title_as_filename(figurename, 'kinetics', '/Users/sophiarenauld/Downloads')

end