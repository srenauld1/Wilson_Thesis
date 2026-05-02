function [resultsTable, trialSummary] = summarize_chrimson_trials(baseFolder, outCsvPath)
% SUMMARIZE_CHRIMSON_TRIALS
%
%   resultsTable = summarize_chrimson_trials(baseFolder, outCsvPath)
%
% Recursively searches baseFolder for files matching:
%   *_processed_final.mat
%
% For each file, extracts an ID of the form:
%   2026_03_04_fly01_cell01_trial01A
% from the filename (before ' pattern...').
%
% Loads the file (expects variable `exptData`), computes:
%   - mean_heading_disp : mean(exptData.opto_sliced_headingdisp  )
%   - mean_fwd_opto     : mean of all exptData.opto_sliced.forwardVelocity_raw
%   - mean_yaw_opto     : mean of all exptData.opto_sliced.angularVelocity_raw
%
% Outputs:
%   CSV: id, mean_heading_disp
%   MAT: trialSummary struct array with fields:
%        - id
%        - mean_fwd_opto
%        - mean_yaw_opto
%
% Inputs:
%   baseFolder : root directory to search (string or char)
%   outCsvPath : full path to output CSV file; if empty or omitted,
%                defaults to 'chrimson_summary.csv' in baseFolder
%
% Output:
%   resultsTable : table with variables:
%                  - id
%                  - value (mean_heading_disp)

    if nargin < 1 || isempty(baseFolder)
        error('You must provide baseFolder.');
    end
    if nargin < 2 || isempty(outCsvPath)
        outCsvPath = fullfile(baseFolder, 'chrimson_summary.csv');
    end

    if ~isfolder(baseFolder)
        error('Base folder does not exist: %s', baseFolder);
    end

    % ---------------------------------------------------------------------
    % 1. Find all *_processed_final.mat files recursively
    % ---------------------------------------------------------------------
    pattern  = fullfile(baseFolder, '**', '*_processed_final.mat');
    fileList = dir(pattern);

    if isempty(fileList)
        warning('No *_processed_final.mat files found under: %s', baseFolder);
        resultsTable = table(string.empty(0,1), [], 'VariableNames', {'id','value'});
        return;
    end

    nFiles = numel(fileList);
    ids    = strings(nFiles, 1);
    values = nan(nFiles, 1);  % mean_heading_disp

    % Preallocate trialSummary (only id, mean_fwd_opto, mean_yaw_opto)
    trialSummary(nFiles,1) = struct('id',"", ...
                                    'mean_fwd_opto',nan, ...
                                    'mean_yaw_opto',nan);

    fprintf('Found %d *_processed_final.mat files.\n', nFiles);

    % ---------------------------------------------------------------------
    % 2. Loop over files, extract ID, load, compute values
    % ---------------------------------------------------------------------
    validCount = 0;

    for k = 1:nFiles
        fpath = fullfile(fileList(k).folder, fileList(k).name);
        fname = fileList(k).name;

        % --- Extract ID from filename: part before ' pattern'
        % Example:
        %   '2026_03_04_fly01_cell01_trial01A pattern20 stim2 function1_processed_final.mat'
        % We want: '2026_03_04_fly01_cell01_trial01A'
        token = regexp(fname, '^(.*)\s+pattern', 'tokens', 'once');

        if isempty(token)
            % Fallback: strip '*_processed_final.mat' and take part before first space
            baseNameNoExt = regexprep(fname, '_processed_final\.mat$', '');
            parts = strsplit(baseNameNoExt, ' ');
            idStr = parts{1};
        else
            idStr = token{1};
        end

        % Load the .mat file (expect it to contain exptData)
        S = load(fpath);
        if ~isfield(S, 'exptData')
            warning('File %s does not contain exptData; skipping.', fpath);
            continue;
        end
        exptData = S.exptData;

        % ---------- Compute mean heading displacement ----------
        mean_heading_disp = NaN;
        if isfield(exptData, 'opto_sliced_headingdisp') ...
                && ~isempty(exptData.opto_sliced_headingdisp )
            mean_heading_disp = mean(exptData.opto_sliced_headingdisp  (:), 'omitnan');
        else
            warning('No exptData.opto_sliced_headingdisp   in %s; value set to NaN.', fpath);
        end

        % ---------- Compute mean forward & yaw velocity during opto ----------
        mean_fwd_opto = NaN;
        mean_yaw_opto = NaN;

        if isfield(exptData, 'opto_sliced') && ~isempty(exptData.opto_sliced)
            % forwardVelocity_raw
            if isfield(exptData.opto_sliced, 'forwardVelocity_raw')
                fwd_cells = {exptData.opto_sliced.forwardVelocity_raw}';  % 80x1 cell
                fwd_mat   = vertcat(fwd_cells{:});                       % 80 x 60000
                mean_fwd_opto = mean(fwd_mat, 1, 'omitnan');             % 1 x 60000
            end

            % angularVelocity_raw
            if isfield(exptData.opto_sliced, 'angularVelocity_raw')
                yaw_cells = {exptData.opto_sliced.angularVelocity_raw}';
                yaw_mat   = vertcat(yaw_cells{:});
                mean_yaw_opto = mean(yaw_mat, 1, 'omitnan');
            end
        else
            warning('exptData.opto_sliced missing or empty in %s; fwd/yaw means set to NaN.', fpath);
        end

        % ---------- Store results ----------
        validCount = validCount + 1;

        ids(validCount)    = string(idStr);
        values(validCount) = mean_heading_disp;

        trialSummary(validCount).id            = string(idStr);
        trialSummary(validCount).mean_fwd_opto = mean_fwd_opto;
        trialSummary(validCount).mean_yaw_opto = mean_yaw_opto;
    end

    % Trim to valid entries only
    ids    = ids(1:validCount);
    values = values(1:validCount);
    trialSummary = trialSummary(1:validCount);

    % ---------------------------------------------------------------------
    % 3. Build table and write CSV
    % ---------------------------------------------------------------------
    resultsTable = table(ids, values, 'VariableNames', {'id','value'});

    try
        writetable(resultsTable, outCsvPath);
        fprintf('Saved summary CSV to: %s\n', outCsvPath);
    catch ME
        warning('Could not write CSV: %s\n  Reason: %s', outCsvPath, ME.message);
    end

    % ---------------------------------------------------------------------
    % 4. Save trialSummary as MAT file next to CSV
    % ---------------------------------------------------------------------
    [csvDir, csvName] = fileparts(outCsvPath);
    matPath = fullfile(csvDir, [csvName '_trialSummary.mat']);

    try
        save(matPath, 'trialSummary');
        fprintf('Saved trialSummary MAT to: %s\n', matPath);
    catch ME
        warning('Could not write trialSummary MAT: %s\n  Reason: %s', matPath, ME.message);
    end
end