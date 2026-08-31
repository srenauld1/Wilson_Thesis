function [resultsTable, flySummary] = summarize_nonopto_sinuosity(baseFolder, outCsvPath, winPts, minVel)
% SUMMARIZE_NONOPTO_SINUOSITY
%
%   [resultsTable, flySummary] = summarize_nonopto_sinuosity(baseFolder, outCsvPath, winPts, minVel)
%
% Recursively searches baseFolder for files matching:
%   *_processed_final.mat
%
% For each file it isolates the NON-OPTO (opto-off) portions of the flat
% path (exptData.optoStim == 0), chops them into fixed-length windows, and
% computes the ecological sinuosity (Benhamou, 2004) of each window using
% benhamou_sinuosity.m. It then reports a mean sinuosity per trial and a
% mean sinuosity per fly, and produces a bar chart with one value per fly.
%
% WHY FIXED-LENGTH WINDOWS:
%   Benhamou sinuosity is sensitive to the number of samples in a path
%   (more points -> more measured turning). To make the averaged epochs
%   comparable, every window used here contains exactly `winPts` samples;
%   opto-off runs shorter than one window are discarded, and leftover tails
%   shorter than winPts are dropped. winPts is locked to a single value for
%   the whole run so every epoch that goes into the averages is the same
%   length. - claude
%
% Inputs:
%   baseFolder : root directory to search (string or char)
%   outCsvPath : full path to output CSV; if empty/omitted, defaults to
%                'nonopto_sinuosity_summary.csv' in baseFolder
%   winPts     : (optional) number of samples per sinuosity window. If empty
%                or omitted, it is derived from the FIRST valid file as
%                round(5 s * sampling_rate) and then reused for every file.
%   minVel     : (optional) minimum mean speed (exptData.totSpeed units) for
%                a window to be counted; windows slower than this are skipped
%                so sinuosity is not dominated by a stationary fly.
%                Default = 3 (matches master_behavior_analysis). Ignored if
%                exptData.totSpeed is absent.
%
% Outputs:
%   resultsTable : per-TRIAL table with variables:
%                  - id            (trial id, e.g. 2026_03_04_fly01_cell01_trial01A)
%                  - fly           (fly id, e.g. 2026_03_04_fly01)
%                  - mean_sinuosity(mean over that trial's opto-off windows)
%                  - n_windows     (number of windows used)
%                  - win_pts       (samples per window)
%   flySummary   : per-FLY struct array with fields:
%                  - fly
%                  - mean_sinuosity (pooled over all of that fly's windows)
%                  - sem_sinuosity
%                  - n_windows
%                  - n_trials
%
% Written to match summarize_chrimson_trials.m conventions. - claude

    %% ---------------- input handling ----------------
    if nargin < 1 || isempty(baseFolder)
        error('You must provide baseFolder.');
    end
    if nargin < 2 || isempty(outCsvPath)
        outCsvPath = fullfile(baseFolder, 'nonopto_sinuosity_summary.csv');
    end
    if nargin < 3
        winPts = [];                 % derive from first valid file
    end
    if nargin < 4 || isempty(minVel)
        minVel = 3;
    end

    winSec_default = 5;              % window length (s) used only if winPts empty

    if ~isfolder(baseFolder)
        error('Base folder does not exist: %s', baseFolder);
    end

    %% ---------------- find files ----------------
    pattern  = fullfile(baseFolder, '**', '*_processed_final.mat');
    fileList = dir(pattern);

    if isempty(fileList)
        warning('No *_processed_final.mat files found under: %s', baseFolder);
        resultsTable = table(string.empty(0,1), string.empty(0,1), [], [], [], ...
            'VariableNames', {'id','fly','mean_sinuosity','n_windows','win_pts'});
        flySummary = struct('fly',strings(0,1), 'mean_sinuosity',[], ...
            'sem_sinuosity',[], 'n_windows',[], 'n_trials',[]);
        return;
    end

    nFiles = numel(fileList);
    fprintf('Found %d *_processed_final.mat files.\n', nFiles);

    % Per-trial accumulators
    trial_id        = strings(nFiles,1);
    trial_fly       = strings(nFiles,1);
    trial_meanSin   = nan(nFiles,1);
    trial_nWin      = zeros(nFiles,1);
    trial_winPts    = nan(nFiles,1);
    trial_sinWindows = cell(nFiles,1);   % per-window sinuosity values

    validCount = 0;

    %% ---------------- loop over files ----------------
    for k = 1:nFiles
        fpath = fullfile(fileList(k).folder, fileList(k).name);
        fname = fileList(k).name;

        % --- Extract trial id (part before ' pattern'), as in summarize_chrimson_trials
        token = regexp(fname, '^(.*)\s+pattern', 'tokens', 'once');
        if isempty(token)
            baseNameNoExt = regexprep(fname, '_processed_final\.mat$', '');
            parts = strsplit(baseNameNoExt, ' ');
            idStr = parts{1};
        else
            idStr = token{1};
        end

        % --- Derive fly id (everything up to and including 'flyNN')
        flyTok = regexp(idStr, '^(.*fly\d+)', 'tokens', 'once');
        if isempty(flyTok)
            flyStr = idStr;   % fallback: use full id
        else
            flyStr = flyTok{1};
        end

        % --- Load
        S = load(fpath);
        if ~isfield(S, 'exptData')
            warning('File %s does not contain exptData; skipping.', fpath);
            continue;
        end
        exptData = S.exptData;

        % --- Required fields
        if ~isfield(exptData,'x') || ~isfield(exptData,'y')
            warning('exptData.x/.y missing in %s; skipping.', fpath);
            continue;
        end
        x = exptData.x(:);
        y = exptData.y(:);
        N = numel(x);

        % --- Opto-off mask (opto-off == optoStim == 0)
        if isfield(exptData,'optoStim') && ~isempty(exptData.optoStim) ...
                && numel(exptData.optoStim) == N
            off_mask = exptData.optoStim(:) == 0;
        else
            warning(['exptData.optoStim missing/mismatched in %s; ' ...
                     'treating entire trial as opto-off.'], fpath);
            off_mask = true(N,1);
        end

        % --- Exclude bar-jump buffer windows if available
        if isfield(exptData,'omit_jump_blocks') && numel(exptData.omit_jump_blocks) == N
            off_mask = off_mask & (exptData.omit_jump_blocks(:) == 0);
        end

        % --- Exclude any samples with NaN position (keeps window lengths honest)
        off_mask = off_mask & ~(isnan(x) | isnan(y));

        % --- Determine window length (points), locked across all files
        if isempty(winPts)
            if ~isfield(exptData,'t') || numel(exptData.t) < 2
                warning(['Cannot derive winPts (no exptData.t) for %s and none ' ...
                         'was provided; skipping.'], fpath);
                continue;
            end
            t = exptData.t(:);
            rate = (numel(t)-1) / (t(end) - t(1));       % samples/s
            winPts = max(3, round(winSec_default * rate));
            fprintf('Locked window length to winPts = %d samples (~%.1f s).\n', ...
                    winPts, winSec_default);
        end

        % --- Speed per sample for the optional moving filter
        if isfield(exptData,'totSpeed') && numel(exptData.totSpeed) == N
            spd = exptData.totSpeed(:);
        else
            spd = [];   % moving filter disabled for this file
        end

        % --- Build fixed-length windows from contiguous opto-off runs
        [Xwin, Ywin] = build_windows(x, y, off_mask, spd, winPts, minVel);

        if isempty(Xwin)
            warning('No valid opto-off windows (winPts=%d) in %s.', winPts, fpath);
            % still record the trial with NaN so it appears in the table
            validCount = validCount + 1;
            trial_id(validCount)      = string(idStr);
            trial_fly(validCount)     = string(flyStr);
            trial_meanSin(validCount) = NaN;
            trial_nWin(validCount)    = 0;
            trial_winPts(validCount)  = winPts;
            trial_sinWindows{validCount} = [];
            continue;
        end

        % --- Sinuosity per window (columns = trajectories)
        sinVals = benhamou_sinuosity(Xwin, Ywin);   % 1 x nWindows
        sinVals = sinVals(isfinite(sinVals));

        validCount = validCount + 1;
        trial_id(validCount)        = string(idStr);
        trial_fly(validCount)       = string(flyStr);
        trial_meanSin(validCount)   = mean(sinVals, 'omitnan');
        trial_nWin(validCount)      = numel(sinVals);
        trial_winPts(validCount)    = winPts;
        trial_sinWindows{validCount} = sinVals(:);
    end

    % Trim to valid entries
    trial_id        = trial_id(1:validCount);
    trial_fly       = trial_fly(1:validCount);
    trial_meanSin   = trial_meanSin(1:validCount);
    trial_nWin      = trial_nWin(1:validCount);
    trial_winPts    = trial_winPts(1:validCount);
    trial_sinWindows = trial_sinWindows(1:validCount);

    if validCount == 0
        warning('No valid trials processed.');
        resultsTable = table(string.empty(0,1), string.empty(0,1), [], [], [], ...
            'VariableNames', {'id','fly','mean_sinuosity','n_windows','win_pts'});
        flySummary = struct('fly',strings(0,1), 'mean_sinuosity',[], ...
            'sem_sinuosity',[], 'n_windows',[], 'n_trials',[]);
        return;
    end

    %% ---------------- per-trial table + CSV ----------------
    resultsTable = table(trial_id, trial_fly, trial_meanSin, trial_nWin, trial_winPts, ...
        'VariableNames', {'id','fly','mean_sinuosity','n_windows','win_pts'});

    try
        writetable(resultsTable, outCsvPath);
        fprintf('Saved per-trial sinuosity CSV to: %s\n', outCsvPath);
    catch ME
        warning('Could not write CSV: %s\n  Reason: %s', outCsvPath, ME.message);
    end

    %% ---------------- per-fly summary (pooled over windows) ----------------
    flies = unique(trial_fly, 'stable');
    nFly  = numel(flies);
    flySummary = repmat(struct('fly',"", 'mean_sinuosity',nan, ...
        'sem_sinuosity',nan, 'n_windows',0, 'n_trials',0), nFly, 1);

    fly_mean = nan(nFly,1);
    fly_sem  = nan(nFly,1);
    for f = 1:nFly
        sel = (trial_fly == flies(f));
        pooled = vertcat(trial_sinWindows{sel});   % all windows for this fly
        pooled = pooled(isfinite(pooled));

        flySummary(f).fly           = flies(f);
        flySummary(f).n_trials      = sum(sel);
        flySummary(f).n_windows     = numel(pooled);
        if isempty(pooled)
            flySummary(f).mean_sinuosity = NaN;
            flySummary(f).sem_sinuosity  = NaN;
        else
            flySummary(f).mean_sinuosity = mean(pooled);
            flySummary(f).sem_sinuosity  = std(pooled) / sqrt(numel(pooled));
        end
        fly_mean(f) = flySummary(f).mean_sinuosity;
        fly_sem(f)  = flySummary(f).sem_sinuosity;
    end

    %% ---------------- bar chart: one value per fly ----------------
    figure('Name','Non-opto path sinuosity per fly','Color','w');
    hold on;
    bar(1:nFly, fly_mean, 'FaceColor', [0.4 0.4 0.4]);
    errorbar(1:nFly, fly_mean, fly_sem, 'k', 'LineStyle','none', 'LineWidth', 1);
    set(gca, 'XTick', 1:nFly, 'XTickLabel', cellstr(flies), 'XTickLabelRotation', 45, ...
        'TickLabelInterpreter','none');
    ylabel('Mean sinuosity (opto-off)');
    xlabel('Fly');
    title(sprintf('Non-opto flat-path sinuosity (winPts = %d samples)', trial_winPts(1)));
    grid on;
    hold off;

    % Save figure next to the CSV
    [csvDir, csvName] = fileparts(outCsvPath);
    figBase = fullfile(csvDir, [csvName '_perfly']);
    try
        savefig(gcf, [figBase '.fig']);
        saveas(gcf,  [figBase '.png']);
        fprintf('Saved per-fly chart to: %s.(fig/png)\n', figBase);
    catch ME
        warning('Could not save chart: %s', ME.message);
    end

    %% ---------------- save summary MAT ----------------
    matPath = fullfile(csvDir, [csvName '_flySummary.mat']);
    try
        save(matPath, 'flySummary', 'resultsTable', '-v7.3');
        fprintf('Saved flySummary MAT to: %s\n', matPath);
    catch ME
        warning('Could not write flySummary MAT: %s\n  Reason: %s', matPath, ME.message);
    end
end


% =====================================================================
% Helper: chop contiguous opto-off runs into fixed-length windows.
% Every returned window has exactly winPts samples so Benhamou sinuosity
% is comparable across epochs. Windows whose mean speed is below minVel
% (when speed is available) are dropped. - claude
% =====================================================================
function [Xwin, Ywin] = build_windows(x, y, off_mask, spd, winPts, minVel)
    Xcols = {};
    Ycols = {};

    % Find contiguous runs of opto-off samples
    d = diff([false; off_mask(:); false]);
    runStarts = find(d == 1);
    runEnds   = find(d == -1) - 1;

    for r = 1:numel(runStarts)
        s = runStarts(r);
        e = runEnds(r);
        runLen = e - s + 1;

        % Number of whole windows in this run; leftover tail is discarded
        nWin = floor(runLen / winPts);
        for w = 1:nWin
            i0 = s + (w-1)*winPts;
            i1 = i0 + winPts - 1;
            idx = i0:i1;

            % Optional moving filter: require mean speed >= minVel
            if ~isempty(spd) && mean(spd(idx)) < minVel
                continue;
            end

            Xcols{end+1} = x(idx); %#ok<AGROW>
            Ycols{end+1} = y(idx); %#ok<AGROW>
        end
    end

    if isempty(Xcols)
        Xwin = [];
        Ywin = [];
    else
        Xwin = cell2mat(Xcols);   % winPts x nWindows (each col a window)
        Ywin = cell2mat(Ycols);
    end
end
