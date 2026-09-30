%% Rerun all Finder-tagged good-behavior VES041 trials
% Standard trials use master_imaging_analysis. Heat trials use
% heat_master_analysis. A timestamped text log records successes, missing
% processed inputs, and errors without stopping the remainder of the batch.

clearvars
baseFolder='/Users/sophiarenauld/stacks';
analysisFolder=fileparts(mfilename('fullpath'));
addpath(analysisFolder);
addpath(genpath('/Users/sophiarenauld/Documents/GitHub/scopa'));

entries=dir(fullfile(baseFolder,'*ves041*')); entries=entries([entries.isdir]);
logRows=cell(0,4);
for q=1:numel(entries)
    trialFolder=fullfile(entries(q).folder,entries(q).name);
    [tagStatus,tags]=system(sprintf('/usr/bin/mdls -raw -name kMDItemUserTags "%s"',trialFolder));
    if tagStatus~=0 || ~contains(tags,'good behavior'),continue,end

    isHeat=contains(lower(entries(q).name),'heat');
    if isHeat
        mode="heat_master_analysis";
        inputFiles=dir(fullfile(trialFolder,'*_or_z0z1_s_.mat'));
    else
        mode="master_imaging_analysis";
        inputFiles=dir(fullfile(trialFolder,'*_s_.mat'));
    end

    if numel(inputFiles)~=1
        logRows(end+1,:)={entries(q).name,mode,"missing input", ... %#ok<SAGROW>
            sprintf('Expected 1 processed input; found %d',numel(inputFiles))};
        fprintf('SKIP %s: processed input count=%d\n',entries(q).name,numel(inputFiles));
        continue
    end

    fprintf('\nRUN %s with %s\n',entries(q).name,mode);
    try
        if isHeat
            heatTrialFolder=trialFolder; %#ok<NASGU>
            run(fullfile(analysisFolder,'heat_master_analysis.m'));
            clear heatTrialFolder
        else
            fullFilePath=fullfile(inputFiles(1).folder,inputFiles(1).name); %#ok<NASGU>
            run(fullfile(analysisFolder,'master_imaging_analysis.m'));
            clear fullFilePath
        end
        close all force
        logRows(end+1,:)={entries(q).name,mode,"success",""}; %#ok<SAGROW>
    catch ME
        close all force
        logRows(end+1,:)={entries(q).name,mode,"error", ... %#ok<SAGROW>
            string(getReport(ME,'extended','hyperlinks','off'))};
        fprintf(2,'ERROR %s: %s\n',entries(q).name,ME.message);
    end
end

logTable=cell2table(logRows,'VariableNames',{'trial','analysis','status','detail'});
logTable.trial=string(logTable.trial); logTable.analysis=string(logTable.analysis);
logTable.status=string(logTable.status); logTable.detail=string(logTable.detail);
logPath=fullfile(baseFolder,'ves041_behavior_quantification','rerun_good_behavior_log.csv');
if ~isfolder(fileparts(logPath)),mkdir(fileparts(logPath));end
writetable(logTable,logPath);
disp(logTable)
