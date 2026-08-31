 function [trialTable, flyTable] = analyze_opto_heading_metrics(files, minVel, walkThresh)
% ANALYZE_OPTO_HEADING_METRICS
%
%   [trialTable, flyTable] = analyze_opto_heading_metrics(files, minVel)
%
% Loads any number of *_processed_final.mat files and, for each opto pulse,
% computes three heading/path metrics in two 1 s windows:
%       pre    = 1 s immediately BEFORE opto onset
%       during = first 1 s of the opto stim
%
% Metrics (per window):
%   (b) rho        : heading consistency = mean resultant length of heading
%                    unit vectors. STOP EPOCHS REMOVED (totSpeed < minVel).
%                    rho = |mean(exp(i*theta))|, bounded [0,1]. Not sensitive
%                    to turning-direction bias.
%   (c) hchange    : heading change = TOTAL angular displacement over the
%                    window (deg) = sum of |sample-to-sample heading change|
%                    (unwrapped). NOT net: both zig-zag and circling
%                    accumulate. Uses ALL samples (no stop removal).
%   (d) straight   : path straightness = displacement / path length, bounded
%                    [0,1]. Uses ALL samples (no stop removal).
%
% Aggregation:
%   - TRIAL = one opto pulse (one opto_sliced entry). One row / one dot per
%     pulse, so each fly contributes many trial dots.
%   - FLY   = mean over ALL of that fly's pulses (pooled across recordings).
%
% Plots (each produced at BOTH the per-trial and per-fly level):
%   Set 1 : for each metric, scatter x = during vs y = pre (color = fly)
%   Set 2 : metric-vs-metric scatters (all 3 pairings) for pre and during
%   Set 3 : pre-opto wiggliness (rho_pre / hchange_pre / straight_pre) vs
%           the opto response (delta fwd, delta heading change, and a
%           combined responder score), color = GENOTYPE, per-genotype linear
%           fit, plus a console Pearson-r table to pick the best predictor.
%   Set 4 : JOINT outcome response space - x = delta fwd (up=faster),
%           y = delta total heading change (down=straighter), color =
%           genotype, genotype means as stars, "faster & straighter"
%           quadrant shaded. This is the plot of your actual outcome.
% All sets are drawn TWICE: once for all pulses, and once restricted to
% pulses whose pre-opto (1 s) mean forward velocity >= 1 (see fwdThresh).
%
% Inputs:
%   files  : (optional) how to choose files. May be:
%            - omitted/[]     : open a FOLDER picker at baseFolder, then
%                               recursively load all matches in the chosen folder
%            - a folder path  : recursively load all matches under it
%            - 'pick'         : open a multi-select FILE dialog (cmd/shift-click)
%            - a char file    : that single file
%            - a cell of paths: exactly those files
%   minVel : (optional) stop threshold for rho (totSpeed units). Default 3,
%            matching master_behavior_analysis. If exptData.totSpeed is
%            absent, no stop removal is applied and a warning is issued.
%   walkThresh : (optional) forward-velocity threshold (mm/s) defining a
%            "walking" sample. fwd_pre/fwd_during are averaged over walking
%            samples only (fwd > walkThresh), so stops/turns don't drag the
%            mean toward 0. Default 1. NaN for a window with no walking.
%
% Outputs:
%   trialTable : one row per opto pulse: id, fly (genotype folder), genotype,
%                is_control, pulse, fwd_pre, fwd_during, and the 6 heading/path
%                metric columns (rho/hchange/straight x pre/during).
%   flyTable   : one row per fly: fly, genotype, is_control, the 8 averaged
%                metric+fwd columns (mean over that fly's pulses), n_pulses,
%                n_trials[recordings]. Computed over ALL pulses (unfiltered).
%
% NOTES / assumptions (flag if wrong):
%   - Heading is exptData.headingPosition, assumed to be in DEGREES.
%   - Opto is "on" wherever exptData.optoStim_ds ~= 0.
%   - x, y, headingPosition, totSpeed, optoStim_ds, t_ds are all at the
%     behavior (FicTrac) sample rate and the same length (~800k here); the
%     raw optoStim/t are ~30x longer and are NOT used. The 1 s window is
%     round(rate) samples.
%   - Pulses without a full 1 s pre window, or whose opto period is shorter
%     than 1 s, are skipped (so "during" is genuinely all opto).
%
% Written to match summarize_chrimson_trials.m conventions. - claude

    %% ---------------- input handling ----------------
    if nargin < 2 || isempty(minVel)
        minVel = 3;
    end
    if nargin < 3 || isempty(walkThresh)
        walkThresh = 1;   % forward-velocity threshold (mm/s) defining "walking"
    end
    win_sec = 1;   % length of pre and during windows (s)

    % Base folder for the heading analysis. - claude
    baseFolder = '/Volumes/Neurobio/wilsonlab/Sophia/behavior/unblinded RubyACR';

    % ----- Resolve the list of files to process -----
    % `files` may be:
    %   omitted / []         -> open a FOLDER picker starting at baseFolder, then
    %                           recursively load all *_processed_final.mat in the
    %                           chosen subfolder
    %   a folder path        -> recursively load ALL *_processed_final.mat under it
    %   'pick'               -> open a multi-select FILE dialog (cmd/shift-click)
    %   a char file path     -> that single file
    %   a cell array of paths -> exactly those files
    if nargin < 1 || isempty(files)
        chosen = uigetdir(baseFolder, 'Select a folder to analyze');
        if isequal(chosen, 0)
            disp('User canceled folder selection.');
            trialTable = table(); flyTable = table();
            return;
        end
        files = scan_folder(chosen);
    elseif ischar(files) && strcmpi(files, 'pick')
        files = pick_files(baseFolder);
        if isempty(files)
            disp('User canceled file selection.');
            trialTable = table(); flyTable = table();
            return;
        end
    elseif ischar(files) && isfolder(files)
        files = scan_folder(files);
    elseif ischar(files)
        files = {files};                                 % single explicit file
    end
    % (cell array passed in is used as-is)

    if isempty(files)
        warning('No *_processed_final.mat files to process.');
        trialTable = table(); flyTable = table();
        return;
    end
    nFiles = numel(files);
    fprintf('Processing %d file(s).\n', nFiles);

    metricNames = {'rho','hchange','straight'};

    % Per-PULSE accumulators. A "trial" here is ONE opto pulse (one
    % opto_sliced entry), so each fly contributes many rows. - claude
    p_id    = strings(0,1);   % recording/file id per pulse
    p_fly   = strings(0,1);   % fly id (genotype folder) per pulse
    p_geno  = strings(0,1);   % genotype string per pulse
    p_ctrl  = false(0,1);     % is-control flag per pulse
    p_pulse = zeros(0,1);     % pulse index within its recording
    p_fwdpre = zeros(0,1);    % mean forward velocity in the 1 s pre window
    p_fwddur = zeros(0,1);    % mean forward velocity in the 1 s opto window
    P = struct();
    for m = 1:numel(metricNames)
        P.(metricNames{m}).pre    = zeros(0,1);
        P.(metricNames{m}).during = zeros(0,1);
    end
    nFilesUsed = 0;

    %% ---------------- loop over files ----------------
    for k = 1:nFiles
        fpath = files{k};
        [~, fname, ext] = fileparts(fpath);
        fname = [fname ext];

        if exist(fpath, 'file') ~= 2
            warning('File not found; skipping: %s', fpath);
            continue;
        end

        % --- recording id from the .mat filename (before ' pattern') ---
        token = regexp(fname, '^(.*)\s+pattern', 'tokens', 'once');
        if isempty(token)
            baseNameNoExt = regexprep(fname, '_processed_final\.mat$', '');
            parts = strsplit(baseNameNoExt, ' ');
            idStr = parts{1};
        else
            idStr = token{1};
        end

        % --- genotype + fly identity from the containing folder ---
        % Folder looks like: fly01_VT044519_x_rubyACR_9
        %   flyId=whole folder (unique animal); genotype=VT044519_x_rubyACR;
        %   effector=rubyACR; isControl = effector contains 'GFP'. - claude
        genoInfo = parse_genotype(fpath);
        if genoInfo.found
            flyStr = genoInfo.flyId;          % unique animal = genotype folder
        else
            flyTok = regexp(idStr, '^(.*fly\d+)', 'tokens', 'once');
            if isempty(flyTok), flyStr = idStr; else, flyStr = flyTok{1}; end
            warning('%s: no genotype folder found in path; genotype = unknown.', fpath);
        end
        genoStr = genoInfo.genotype;
        isCtrl  = genoInfo.isControl;

        % --- load ONLY the behavior-rate fields we use (fast) ---
        % These files are v7.3/HDF5 and hold nine ~24M-sample raw arrays we
        % never touch; a full load() pulls all of them (~13x slower). Read
        % just the needed datasets directly instead. - claude
        exptData = load_behavior_fields(fpath);
        if isempty(fieldnames(exptData))
            warning('No readable exptData behavior fields in %s; skipping.', fpath);
            continue;
        end

        % --- required fields ---
        % NOTE: exptData holds two sample rates. Raw DAQ fields (optoStim, t,
        % ...) are ~30x longer than the FicTrac behavior fields (x, y,
        % headingPosition, totSpeed). We work entirely on the behavior rate,
        % so use the DOWNSAMPLED opto/time (optoStim_ds, t_ds) that match x/y.
        % - claude
        req = {'x','y','headingPosition','optoStim_ds','t_ds'};
        missing = req(~isfield(exptData, req));
        if ~isempty(missing)
            warning('%s missing field(s): %s; skipping.', fpath, strjoin(missing, ', '));
            continue;
        end

        x       = exptData.x(:);
        y       = exptData.y(:);
        heading = deg2rad(exptData.headingPosition(:));   % assumed degrees -> radians
        opto    = exptData.optoStim_ds(:);
        t       = exptData.t_ds(:);
        N       = numel(x);

        if numel(y)~=N || numel(heading)~=N || numel(opto)~=N
            warning('%s: x/y/heading/optoStim length mismatch; skipping.', fpath);
            continue;
        end

        % speed for stop removal (rho only)
        if isfield(exptData,'totSpeed') && numel(exptData.totSpeed)==N
            spd = exptData.totSpeed(:);
        else
            spd = [];
            warning('%s: no totSpeed; rho computed without stop removal.', fpath);
        end

        % forward velocity (behavior rate) for the pre-opto filter
        if isfield(exptData,'forwardVelocity') && numel(exptData.forwardVelocity)==N
            fwd = exptData.forwardVelocity(:);
        else
            fwd = [];
            warning('%s: no forwardVelocity; pre-opto fwd-vel filter set to NaN.', fpath);
        end

        % --- window length in samples ---
        if numel(t) < 2 || (t(end)-t(1)) <= 0
            warning('%s: bad time vector; skipping.', fpath);
            continue;
        end
        rate   = (numel(t)-1) / (t(end)-t(1));   % samples/s
        winPts = max(3, round(win_sec * rate));

        % --- opto onsets/offsets ---
        is_on   = opto ~= 0;
        d       = diff([false; is_on; false]);
        onsets  = find(d == 1);
        offsets = find(d == -1) - 1;

        % --- per-pulse metrics; each usable pulse becomes one trial row ---
        nSkip = 0; nUsed = 0;
        for p = 1:numel(onsets)
            s = onsets(p);
            e = offsets(p);
            optoLen = e - s + 1;

            pre_idx = (s-winPts):(s-1);
            dur_idx = s:(s+winPts-1);

            % need full pre window, full during window, and opto >= 1 s
            if pre_idx(1) < 1 || dur_idx(end) > N || optoLen < winPts
                nSkip = nSkip + 1;
                continue;
            end

            [rp, hp, sp] = window_metrics(pre_idx, heading, x, y, spd, minVel);
            [rd, hd, sd] = window_metrics(dur_idx, heading, x, y, spd, minVel);

            % Forward velocity conditioned on WALKING samples only
            % (forwardVelocity > walkThresh). Averaging over the whole window
            % includes stops/turns and drags the mean toward 0; restricting to
            % walking samples gives the true walking-forward speed. NaN if the
            % fly never walked forward in that window. - claude
            if isempty(fwd)
                fpre = NaN; fdur = NaN;
            else
                wp = fwd(pre_idx); wp = wp(wp > walkThresh);
                wd = fwd(dur_idx); wd = wd(wd > walkThresh);
                if isempty(wp), fpre = NaN; else, fpre = mean(wp, 'omitnan'); end
                if isempty(wd), fdur = NaN; else, fdur = mean(wd, 'omitnan'); end
            end

            p_id(end+1,1)    = string(idStr);   %#ok<AGROW>
            p_fly(end+1,1)   = string(flyStr);  %#ok<AGROW>
            p_geno(end+1,1)  = string(genoStr); %#ok<AGROW>
            p_ctrl(end+1,1)  = isCtrl;          %#ok<AGROW>
            p_pulse(end+1,1) = p;               %#ok<AGROW>
            p_fwdpre(end+1,1) = fpre;           %#ok<AGROW>
            p_fwddur(end+1,1) = fdur;           %#ok<AGROW>
            P.rho.pre(end+1,1)         = rp;    %#ok<AGROW>
            P.rho.during(end+1,1)      = rd;    %#ok<AGROW>
            P.hchange.pre(end+1,1)     = hp;    %#ok<AGROW>
            P.hchange.during(end+1,1)  = hd;    %#ok<AGROW>
            P.straight.pre(end+1,1)    = sp;    %#ok<AGROW>
            P.straight.during(end+1,1) = sd;    %#ok<AGROW>
            nUsed = nUsed + 1;
        end

        if nUsed == 0
            warning('%s: no usable opto pulses (skipped %d).', fpath, nSkip);
        else
            nFilesUsed = nFilesUsed + 1;
            fprintf('  %s: %d pulses used (%d skipped).\n', idStr, nUsed, nSkip);
        end
    end

    if isempty(p_id)
        warning('No usable opto pulses found in any file.');
        trialTable = table(); flyTable = table();
        return;
    end

    %% ---------------- build trial table (one row per opto pulse) ----------------
    trialTable = table(p_id, p_fly, p_geno, p_ctrl, p_pulse, p_fwdpre, p_fwddur, ...
        P.rho.pre, P.rho.during, ...
        P.hchange.pre, P.hchange.during, ...
        P.straight.pre, P.straight.during, ...
        'VariableNames', {'id','fly','genotype','is_control','pulse', ...
        'fwd_pre','fwd_during','rho_pre','rho_during', ...
        'hchange_pre','hchange_during','straight_pre','straight_during'});

    %% ---------------- PLOTS ----------------
    metricLabels = struct('rho','\rho (heading consistency)', ...
                          'hchange','Heading change (deg)', ...
                          'straight','Path straightness');

    % --- Full set: all usable pulses (also returns the per-fly output table) ---
    flyTable = plot_level_sets(trialTable, metricNames, metricLabels, 'all pulses');

    % --- Filtered set: only pulses with pre-opto forward velocity >= threshold ---
    fwdThresh = 1;   % pre-opto forward velocity threshold (forwardVelocity units)
    keep = trialTable.fwd_pre >= fwdThresh;
    fprintf('Pre-opto fwd vel >= %g: keeping %d / %d pulses.\n', ...
            fwdThresh, sum(keep), height(trialTable));
    if any(keep)
        plot_level_sets(trialTable(keep,:), metricNames, metricLabels, ...
            sprintf('pre-opto fwd vel >= %g', fwdThresh));
    else
        warning('No pulses with pre-opto forward velocity >= %g; skipping filtered plots.', fwdThresh);
    end
end


% =====================================================================
% Build per-fly table + colors for a (sub)set of pulses and draw all four
% figures (Set 1 & Set 2, at trial and fly level). Returns the per-fly
% table for this subset. `tag` is appended to figure names/titles. - claude
% =====================================================================
function flyT = plot_level_sets(T, metricNames, metricLabels, tag)
    % Columns averaged per fly (includes fwd for the predictor->response plots)
    fcols = {'rho_pre','rho_during','hchange_pre','hchange_during', ...
             'straight_pre','straight_during','fwd_pre','fwd_during'};

    [uFlies, ~, flyIdxPulse] = unique(T.fly, 'stable');
    nFly = numel(uFlies);

    flyVals = nan(nFly, numel(fcols));
    flyNpul = zeros(nFly,1);
    flyNtri = zeros(nFly,1);
    flyGeno = strings(nFly,1);
    flyCtrl = false(nFly,1);
    for f = 1:nFly
        sel = (flyIdxPulse == f);
        flyNpul(f) = sum(sel);
        flyNtri(f) = numel(unique(T.id(sel)));   % # recordings for this fly
        gi = find(sel, 1);
        flyGeno(f) = T.genotype(gi);             % constant within a fly
        flyCtrl(f) = T.is_control(gi);
        for c = 1:numel(fcols)
            flyVals(f,c) = mean(T.(fcols{c})(sel), 'omitnan');
        end
    end
    flyT = table(uFlies, flyGeno, flyCtrl, ...
        flyVals(:,1), flyVals(:,2), flyVals(:,3), flyVals(:,4), ...
        flyVals(:,5), flyVals(:,6), flyVals(:,7), flyVals(:,8), ...
        flyNpul, flyNtri, ...
        'VariableNames', ['fly','genotype','is_control', fcols, {'n_pulses','n_trials'}]);

    % colors: one per fly
    if nFly <= 7, cmap = lines(nFly); else, cmap = hsv(nFly); end
    flyNames    = cellstr(uFlies);
    trialColors = flyIdxPulse;    % one entry per pulse
    flyColors   = (1:nFly).';

    % Set 1: during (x) vs pre (y), per metric
    plot_set1(T,    trialColors, cmap, flyNames, metricNames, metricLabels, ['trial | ' tag]);
    plot_set1(flyT, flyColors,   cmap, flyNames, metricNames, metricLabels, ['fly | '   tag]);

    % Set 2: metric-vs-metric, for pre and during
    plot_set2(T,    trialColors, cmap, flyNames, metricNames, metricLabels, ['trial | ' tag]);
    plot_set2(flyT, flyColors,   cmap, flyNames, metricNames, metricLabels, ['fly | '   tag]);

    % Set 3: pre-opto wiggliness (predictor) -> opto response, per genotype
    plot_predictor_response(T,    ['trial | ' tag]);
    plot_predictor_response(flyT, ['fly | '   tag]);

    % Set 4: joint outcome response space (delta fwd vs delta heading change)
    plot_response_space(T,    ['trial | ' tag]);
    plot_response_space(flyT, ['fly | '   tag]);
end


% =====================================================================
% Per-window metric computation. - claude
% rho removes stop epochs; hchange and straight use all samples.
% =====================================================================
function [rho, hchange_deg, straight] = window_metrics(idx, heading, x, y, spd, minVel)
    hr = heading(idx);
    xk = x(idx);
    yk = y(idx);

    % (b) rho — heading consistency, stops removed
    if isempty(spd)
        mv = true(numel(idx),1);
    else
        mv = spd(idx) >= minVel;
    end
    hr_mv = hr(mv);
    hr_mv = hr_mv(~isnan(hr_mv));
    if numel(hr_mv) >= 2
        rho = abs(mean(exp(1i*hr_mv)));   % mean resultant length in [0,1]
    else
        rho = NaN;
    end

    % (c) heading change — TOTAL angular displacement over the window (deg),
    % i.e. sum of the absolute sample-to-sample heading changes (not net).
    % unwrap first so wraps are handled, then accumulate |d theta|. - claude
    hu = unwrap(hr);
    hchange_deg = sum(abs(rad2deg(diff(hu))));

    % (d) path straightness — displacement / path length, all samples
    disp_ = hypot(xk(end)-xk(1), yk(end)-yk(1));
    plen  = sum(hypot(diff(xk), diff(yk)));
    if plen > 0
        straight = disp_ / plen;
    else
        straight = NaN;
    end
end


% =====================================================================
% Set 1: for each metric, scatter x = during, y = pre, with unity line.
% =====================================================================
function plot_set1(T, colorIdx, cmap, flyNames, metricNames, metricLabels, level)
    figure('Name', sprintf('Set1 during-vs-pre (%s-level)', level), 'Color','w');
    nM = numel(metricNames);
    for m = 1:nM
        ax = subplot(1, nM, m); hold(ax,'on');
        mn = metricNames{m};
        xv = T.([mn '_during']);
        yv = T.([mn '_pre']);
        scatter_by_fly(ax, xv, yv, colorIdx, cmap, flyNames);

        % unity line over common range
        lo = min([xv; yv], [], 'omitnan');
        hi = max([xv; yv], [], 'omitnan');
        if isfinite(lo) && isfinite(hi) && hi > lo
            plot(ax, [lo hi], [lo hi], 'k--', 'HandleVisibility','off');
            axis(ax, [lo hi lo hi]); axis(ax,'square');
        end
        xlabel(ax, ['During: ' metricLabels.(mn)]);
        ylabel(ax, ['Pre: '    metricLabels.(mn)]);
        title(ax, metricLabels.(mn));
        grid(ax,'on');
        if m == 1
            legend(ax, 'Location','bestoutside', 'Interpreter','none');
        end
    end
    sgtitle(sprintf('During vs pre-opto  (%s-level, color = fly)', level));
end


% =====================================================================
% Set 2: metric-vs-metric scatters (3 pairings) x {pre, during}.
% Layout: rows = windows (pre, during), cols = metric pairs.
% =====================================================================
function plot_set2(T, colorIdx, cmap, flyNames, metricNames, metricLabels, level)
    pairs = {[1 2],[1 3],[2 3]};   % rho-hchange, rho-straight, hchange-straight
    windows = {'pre','during'};
    figure('Name', sprintf('Set2 metric-vs-metric (%s-level)', level), 'Color','w');
    nP = numel(pairs);
    for w = 1:numel(windows)
        win = windows{w};
        for pp = 1:nP
            ax = subplot(numel(windows), nP, (w-1)*nP + pp); hold(ax,'on');
            a = metricNames{pairs{pp}(1)};
            b = metricNames{pairs{pp}(2)};
            xv = T.([a '_' win]);
            yv = T.([b '_' win]);
            scatter_by_fly(ax, xv, yv, colorIdx, cmap, flyNames);
            xlabel(ax, metricLabels.(a));
            ylabel(ax, metricLabels.(b));
            title(ax, sprintf('%s window', win));
            grid(ax,'on');
            if w == 1 && pp == 1
                legend(ax, 'Location','bestoutside', 'Interpreter','none');
            end
        end
    end
    sgtitle(sprintf('Metric vs metric  (%s-level, color = fly)', level));
end


% =====================================================================
% Scatter helper: one series per fly so the legend maps color -> fly.
% =====================================================================
function scatter_by_fly(ax, xv, yv, colorIdx, cmap, flyNames)
    nF = size(cmap,1);
    if isscalar(xv) || isscalar(colorIdx)
        % fly-level: colorIdx is (1:nF)', one point each
    end
    for f = 1:nF
        sel = (colorIdx == f);
        if any(sel)
            scatter(ax, xv(sel), yv(sel), 45, cmap(f,:), 'filled', ...
                'MarkerFaceAlpha', 0.7, 'DisplayName', flyNames{f});
        end
    end
end


% =====================================================================
% Load only the behavior-rate fields the analysis needs. For v7.3/HDF5
% files this reads each dataset individually (~13x faster than load(),
% since it skips the huge raw arrays). Falls back to full load() for
% older .mat formats. Missing fields are simply absent from the struct
% (so isfield() checks downstream still work). - claude
% =====================================================================
function E = load_behavior_fields(fpath)
    fields = {'x','y','headingPosition','optoStim_ds','t_ds', ...
              'totSpeed','forwardVelocity','sidewaysVelocity'};
    E = struct();
    isHDF5 = true;
    try
        h5info(fpath);            % errors quickly if not HDF5 (v7 .mat)
    catch
        isHDF5 = false;
    end
    if isHDF5
        for i = 1:numel(fields)
            try
                E.(fields{i}) = h5read(fpath, ['/exptData/' fields{i}]);
            catch
                % dataset absent -> leave field out
            end
        end
    else
        S = load(fpath);
        if isfield(S, 'exptData'), E = S.exptData; end
    end
end


% =====================================================================
% Recursively find all *_processed_final.mat under a folder. - claude
% =====================================================================
function files = scan_folder(folder)
    if ~isfolder(folder)
        warning('Folder not found: %s', folder);
        files = {};
        return;
    end
    d = dir(fullfile(folder, '**', '*_processed_final.mat'));
    files = arrayfun(@(e) fullfile(e.folder, e.name), d(:), 'UniformOutput', false);
    fprintf('Found %d *_processed_final.mat file(s) under %s\n', numel(files), folder);
end


% =====================================================================
% Multi-select dialog (cmd/shift-click to pick several). Returns a cell
% array of full paths, or {} if canceled. - claude
% =====================================================================
function files = pick_files(baseFolder)
    startSpec = fullfile(baseFolder, '*_processed_final.mat');
    [fn, fp] = uigetfile(startSpec, ...
        'Select processed behavior files (cmd/shift-click for multiple)', ...
        'MultiSelect', 'on');
    if isequal(fn, 0)
        files = {};
        return;
    end
    if ischar(fn)          % single file selected -> wrap
        fn = {fn};
    end
    files = cellfun(@(f) fullfile(fp, f), fn, 'UniformOutput', false);
end


% =====================================================================
% Parse genotype/fly identity from a genotype folder in the path. Handles
% the messy real-cohort naming, e.g.:
%   fly10_VT044519_x_rubyACR_10   -> driver VT044519, effector rubyACR (EXPERIMENTAL)
%   fly11_VT044519_x_GFP_11       -> GFP control
%   fly03_VT044519_x_UAS-GFP_3    -> GFP control (hyphen tolerated)
%   fly15_VT044519_GFP_15         -> GFP control (no "_x_")
%   fly08_emptysplit_x_rubyACR_8  -> effector control (no driver = no expression)
%   fly41_R31_rubyACR_41          -> driver R31, effector rubyACR (EXPERIMENTAL)
% Group logic: control if effector is GFP OR the driver is an empty split
% (no functional silencing); experimental only when a real driver drives
% rubyACR. - claude
% =====================================================================
function info = parse_genotype(fpath)
    info = struct('found',false, 'flyId',"unknown", 'genotype',"unknown", ...
                  'driver',"unknown", 'effector',"unknown", 'isControl',false);
    parts = strsplit(fpath, filesep);
    for i = numel(parts):-1:1          % search deepest folder up
        pc = parts{i};
        if isempty(regexp(pc, '^fly\d+_', 'once')), continue; end   % genotype folder
        info.found = true;
        info.flyId = string(pc);

        % effector: search the whole name case-insensitively
        if ~isempty(regexpi(pc, 'gfp', 'once'))
            info.effector = "GFP";
        elseif ~isempty(regexpi(pc, 'acr', 'once'))     % rubyACR / ACR
            info.effector = "rubyACR";
        else
            info.effector = "unknown";
        end

        % driver
        if ~isempty(regexpi(pc, 'empty', 'once'))
            info.driver = "emptysplit";
        elseif ~isempty(regexp(pc, 'VT044519', 'once'))
            info.driver = "VT044519";
        else
            tok = regexp(pc, '^fly\d+_([A-Za-z0-9]+)', 'tokens', 'once');
            if ~isempty(tok), info.driver = string(tok{1}); else, info.driver = "unknown"; end
        end

        info.genotype  = info.driver + ">" + info.effector;
        % control = GFP effector, or a driverless (empty-split) line
        info.isControl = (info.effector == "GFP") || (info.driver == "emptysplit");
        return;
    end
end


% =====================================================================
% Set 3: pre-opto wiggliness (predictor) -> opto response, colored by
% genotype with a per-genotype linear fit. Also prints a correlation
% table so the data can pick the best wiggliness predictor. - claude
%   predictors (pre-opto): rho_pre, hchange_pre, straight_pre
%   responses            : delta_fwd, delta_hchange, response_score
%                          (the JOINT outcome: faster + less turning)
% Works on either the pulse table or the per-fly table (same column names).
% =====================================================================
function plot_predictor_response(T, tag)
    predictors = {'rho_pre','hchange_pre','straight_pre'};
    predLabels = {'pre \rho (low = wiggly)','pre heading change (deg)','pre straightness (low = wiggly)'};

    % Outcome is the JOINT response: forward velocity UP and total heading
    % change DOWN. Responses focus on that: delta_fwd, delta_hchange, and a
    % combined responder score = z(delta_fwd) - z(delta_hchange) (high =
    % faster AND straighter, i.e. turning decreased). - claude
    dFwd  = T.fwd_during     - T.fwd_pre;
    dHead = T.hchange_during - T.hchange_pre;
    score = znan(dFwd) - znan(dHead);
    Rname   = {'delta_fwd','delta_hchange','response_score'};
    Rlabel  = {'\Delta fwd vel (up=faster)','\Delta heading change (down=straighter)','responder score (up=both)'};
    Rvals   = {dFwd, dHead, score};

    % genotype grouping + colors (control = gray, experimental = colored)
    [uGeno, ~, gidx] = unique(cellstr(T.genotype), 'stable');
    nG   = numel(uGeno);
    gcol = geno_colors(uGeno);

    nP = numel(predictors); nR = numel(Rname);
    figure('Name', sprintf('Set3 predictor->response (%s)', tag), 'Color','w');
    for pi = 1:nP
        xv = T.(predictors{pi});
        for ri = 1:nR
            ax = subplot(nP, nR, (pi-1)*nR + ri); hold(ax,'on');
            yv = Rvals{ri};
            for g = 1:nG
                sel = (gidx==g) & isfinite(xv) & isfinite(yv);
                if ~any(sel), continue; end
                scatter(ax, xv(sel), yv(sel), 28, gcol(g,:), 'filled', ...
                    'MarkerFaceAlpha', 0.55, 'DisplayName', uGeno{g});
                if sum(sel) >= 3
                    pc = polyfit(xv(sel), yv(sel), 1);
                    xr = linspace(min(xv(sel)), max(xv(sel)), 20);
                    plot(ax, xr, polyval(pc, xr), '-', 'Color', gcol(g,:), ...
                        'LineWidth', 2, 'HandleVisibility','off');
                end
            end
            grid(ax,'on');
            xlabel(ax, predLabels{pi}, 'Interpreter','tex');
            ylabel(ax, Rlabel{ri},    'Interpreter','tex');
            if pi==1 && ri==1
                legend(ax, 'Location','best', 'Interpreter','none');
            end
        end
    end
    sgtitle(sprintf('Pre-opto wiggliness -> opto response  (%s, color = genotype)', tag), ...
        'Interpreter','none');

    % ---- console correlation table: which predictor predicts best ----
    fprintf('\n[%s] Pearson r  (predictor -> response), by genotype:\n', tag);
    for g = 1:nG
        fprintf('  %s  (n=%d)\n', uGeno{g}, sum(gidx==g));
        for pi = 1:nP
            xv = T.(predictors{pi});
            line = sprintf('     %-13s', predictors{pi});
            for ri = 1:nR
                yv = Rvals{ri};
                sel = (gidx==g) & isfinite(xv) & isfinite(yv);
                if sum(sel) >= 3
                    r = corr(xv(sel), yv(sel));
                    line = [line sprintf('  %s r=% .2f', Rname{ri}, r)]; %#ok<AGROW>
                else
                    line = [line sprintf('  %s r=  NaN', Rname{ri})];    %#ok<AGROW>
                end
            end
            fprintf('%s\n', line);
        end
    end
end


% =====================================================================
% Set 4: the JOINT outcome as a 2-D response space -
%   x = delta forward velocity (during - pre); up  = faster
%   y = delta total heading change (during - pre); down = straighter
% Colored by genotype, with each genotype's mean (star). The "responder"
% quadrant (faster AND straighter: x>0, y<0) is shaded. Under the
% hypothesis, rubyACR responders fall into that quadrant while GFP sits at
% the origin. - claude
% =====================================================================
function plot_response_space(T, tag)
    dFwd  = T.fwd_during     - T.fwd_pre;
    dHead = T.hchange_during - T.hchange_pre;
    [uGeno, ~, gidx] = unique(cellstr(T.genotype), 'stable');
    nG   = numel(uGeno);
    gcol = geno_colors(uGeno);

    figure('Name', sprintf('Set4 response space dFwd vs dHeading (%s)', tag), 'Color','w');
    ax = gca; hold(ax,'on');
    for g = 1:nG
        sel = (gidx==g) & isfinite(dFwd) & isfinite(dHead);
        if ~any(sel), continue; end
        scatter(ax, dFwd(sel), dHead(sel), 30, gcol(g,:), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'DisplayName', uGeno{g});
    end
    xl = xlim(ax); yl = ylim(ax);

    % shade "responder" quadrant: faster (dFwd>0) AND straighter (dHead<0)
    patch(ax, [0 xl(2) xl(2) 0], [0 0 yl(1) yl(1)], [0.2 0.7 0.2], ...
        'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');
    xline(ax, 0, 'k-', 'HandleVisibility','off');
    yline(ax, 0, 'k-', 'HandleVisibility','off');

    % genotype means (star markers on top)
    for g = 1:nG
        sel = (gidx==g) & isfinite(dFwd) & isfinite(dHead);
        if ~any(sel), continue; end
        plot(ax, mean(dFwd(sel),'omitnan'), mean(dHead(sel),'omitnan'), 'p', ...
            'MarkerSize', 16, 'MarkerFaceColor', gcol(g,:), 'MarkerEdgeColor','k', ...
            'HandleVisibility','off');
    end

    xlim(ax, xl); ylim(ax, yl);
    text(ax, xl(2), yl(1), 'faster & straighter  ', 'Color',[0 0.5 0], ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom', 'FontWeight','bold');
    xlabel(ax, '\Delta forward velocity (during - pre)', 'Interpreter','tex');
    ylabel(ax, '\Delta total heading change (during - pre)', 'Interpreter','tex');
    title(ax, sprintf('Joint response space (%s, color = genotype, star = mean)', tag), ...
        'Interpreter','none');
    legend(ax, 'Location','best', 'Interpreter','none');
    grid(ax,'on'); hold(ax,'off');
end


% =====================================================================
% Genotype -> color: control (GFP) = gray, experimental = colored. - claude
% =====================================================================
function gcol = geno_colors(uGeno)
    nG = numel(uGeno);
    gcol = zeros(nG,3);
    palette = lines(max(nG,1)); ei = 1;
    for g = 1:nG
        % controls (GFP or empty-split) = gray shades; experimental = colored
        if ~isempty(regexpi(uGeno{g}, 'gfp', 'once'))
            gcol(g,:) = [0.5 0.5 0.5];
        elseif ~isempty(regexpi(uGeno{g}, 'empty', 'once'))
            gcol(g,:) = [0.7 0.7 0.7];
        else
            gcol(g,:) = palette(ei,:); ei = ei + 1;
        end
    end
end


% =====================================================================
% z-score ignoring NaNs; returns zeros if std is 0/non-finite. - claude
% =====================================================================
function z = znan(x)
    m = mean(x, 'omitnan');
    s = std(x, 'omitnan');
    if s == 0 || ~isfinite(s)
        z = zeros(size(x));
    else
        z = (x - m) / s;
    end
end
