notes from behavior exp:
```
% runExptBattery.m
%
% Top level function for running 1 or more trials of MULTIPLE electrophysiology experiment.
%
% Sets up folder structure for organizing data (date/fly/cell),
%  prompts whether this is new cell or additional trial(s) for cell. For
%  new cell, prompts for experimental condition, asks to run pre-expt
%  routines (pipette resistance, seal test, etc.), and runs through
%  experimental trials until no more new trials. For same cell (in case
%  previous run of runEphysExpt stopped unexpectedly), can run additional
%  trials.
%
% Original: 01/30/2021 -HY
% Updated:  01/02/2023 - MC converted from single to battery
% Updated: 10/17/2024 - SMR adapted


function runExptBattery()
                    
% clean up



close all


% clear function collectData() because it has persistent variable
clear collectData

% initialize persistent variables for this function
persistent cellDirPath trialNum;

% load constant settings
[dataDir, exptFnDir, settings] = ephysSettings();

% Asks whether this is a new cell
newCell = input('[INPUT] New trial batch? (y/n): ', 's');


%% NEW CELL
if (~strcmpi(newCell, 'n'))
    %% set up new cell
    % clear functions and persistent and global variables for new cell
    clear collectData % has persistent variable whichInScan

    % reset values of persistent variables
    cellDirPath = [];
    trialNum = [];

    % set up folder structure by date for organizing data
    dateDirPath = fullfile(dataDir, datestr(now,'yyyy_mm_dd'));

    % get just name of directory
    [~, dateDirName, ~] = fileparts(dateDirPath);
    % make date folder if it doesn't already exist
    if ~exist(dateDirPath, 'dir')
        mkdir(dateDirPath)
    end
    cd(dateDirPath);

    % figure out which fly we're on in the date directory
    dateDirContents = dir(dateDirPath);
    currFlyDirs = dateDirContents(...
        contains({dateDirContents.name},'fly'));

    % asks whether this is a new fly
    newFly = input('\n[INPUT] New Fly? (y/n): ', 's');
    if (strcmpi(newFly, 'n')) % NOT A NEW FLY
        flyNum = length(currFlyDirs);

        if (flyNum == 0)
            flyNum = 1;
            flyDirName = sprintf('fly%02d',flyNum); % fly folder name
            mkdir(flyDirName);
        else
            flyDirName = sprintf('fly%02d',flyNum); % fly folder name
        end

        % go to fly directory
        cd(flyDirName);
        flyDirPath = pwd;
    else % NEW FLY (this is default)
        flyNum = length(currFlyDirs) + 1;
        flyDirName = sprintf('fly%02d',flyNum); % fly folder name

        % creates fly folder in date directory
        mkdir(flyDirName); % make fly folder

        % go to fly directory
        cd(flyDirName);
        flyDirPath = pwd;
    end

    % create new cell folder
    % existing cells for this fly
    flyDirContents = dir(flyDirPath);
    currCellDirs = flyDirContents(contains({flyDirContents.name}, 'cell'));
    currCellDirs(contains({currCellDirs.name}, 'processed'))=[];
    cellNum = length(currCellDirs) + 1;
    cellDirName = sprintf('cell%02d', cellNum); % cell folder name
    % create cell folder
    mkdir(cellDirName);

    % go to cell directory
    cd(cellDirName);
    cellDirPath = pwd;

    % generate basic experimental info struct
    exptInfo.dateDir = dateDirName;
    exptInfo.flyDir = flyDirName;
    exptInfo.cellDir = cellDirName;
    exptInfo.exptDate = datestr(now, 'yymmdd');
    exptInfo.exptStartTime = datestr(now, 'HH:MM:SS');

    % save settings, fly metadata, experimental info to file
    save('metaDat.mat', 'settings', 'exptInfo', '-v7.3');


    %% **ASKS TO RUN PRE-PRE-EXPT ACCLIMITIZATION PERIOD**
    % if using display, remind user to check if panels are working


    % Asks whether to include acclimatization block, then runs until user declines
    acclimateDur = input('[INPUT] Add acclimatization period? (min): ');
    acclimateN = 1; %counter
    acclimatefictrac = input('[INPUT] dark acclim (0) or bright bar acclim (1)? : ');
    if acclimatefictrac == 0
        acclimateFTMode = 0; % (0) OL with nothing, (1) CL with dark bar, (2) CL with bright bar
    else
        acclimateFTMode = 2;
    end   
    while 1
        if ~isempty(acclimateDur) && acclimateDur>0
            % connect display during first run
            if acclimateN==1
                disp('[NOTICE] Connecting to G4 display...');
                connectHost; %required
                Panel_com('all_on'); %turn on
                checkPanels = input('[INPUT] Check display is working...');
                Panel_com('all_off'); %turn off

                %make g4 log dir if needed
                g4log_dir = fullfile(cellDirPath,'g4logs');
                if ~exist(g4log_dir, 'dir')
                    mkdir(g4log_dir)
                end
            end
            % configure/start fictrac
            startFicTrac_smr(cellDirPath,0)

            % begin acclimatization block
            [rawData, inputParams, rawOutput] = runAcclimatization(settings,acclimateDur*60,acclimateFTMode);
            % save and process pre-pre-expt data !!
            inputParams.filename = ['preExptAcclimate_' num2str(acclimateN)];
            save([inputParams.filename '_raw.mat'], 'rawData', 'inputParams','-v7.3');
            [~] = postExptRoutine(inputParams, rawData, rawOutput, settings);

            % asks about running another acclimatization block
            runAcclimateAgain = input('\n[INPUT] Run another acclimatization? (y/n): ','s');
            acclimateN = acclimateN+1; %update counter
            if ~strcmpi(runAcclimateAgain, 'y')
                disp('Acclimatization was not run again.');
                % close fictrac prior to main expt
                system('Taskkill/IM cmd.exe');
                disp('Fictrac terminated.')
                close all
                Panel_com('all_off');
                disconnectHost
                break;
            end

        else
            disp('Acclimatization was not run.');
            break;
        end
    end
end

%% NOT A NEW CELL
cd(cellDirPath) % make sure we're in the cell directory

%% prompt user to select an experiment

exptTypeFileName = 'batteryStimG4PanelsFictrac_ves041.m';

% get name without .m
exptTypeName = extractBefore(exptTypeFileName, '.');
% convert selected experiment file into function handle
exptFn = str2func(exptTypeName);

inUseFictrac = contains(exptTypeName,'fictrac','IgnoreCase',true); %determine if fictrac used
inUsePanels = contains(exptTypeName,'g4','IgnoreCase',true); %determine if panels used
inUseStim = contains(exptTypeName,'stim','IgnoreCase',true); %determine if opto pulse used


%% select a battery

% select which battery to be run ahead of time
% SMR CHANGE
fprintf(['\nBATTERY OPTIONS\n',...
    '(01) RubyACR1: 1s and 5sec opto bouts cl bar\n'...
    '(02) CsChrimson: 1s Closed bar \n']);
batterySelect = input('[INPUT] Select experimental battery: ');

% then pull matching patterns and functions
% SMR CHANGE
% current options
% stim 1 = 0.25, 5s on 5s off
% stim 2 = 0.25, 7 off, 1 on, 7 off
% stim 3 = 0.25, 7.5 off, 0.5 on, 7.5 off
% stim 4 = 0.25, 15 off 10 on 15 off
% stim 5 = 0.25, 12.5 off 5 on 12.5 off
switch batterySelect
    case 1 %rubyacr
        pattSelect = [20]; % cl bar
        funcSelect = [1]; % this is likely where i will go open vs closed 1=closed, 0 = open
        stimSelect = [2]; % 1 trial 1s stim (stim=2, v=0.25) 1 trial 5s stim (stim=5, v=0.25)
        trialBreak = 5; %sec, break between trial
    case 2 %gain battery
        pattSelect = [20 20 20]; % bright bar
        funcSelect = [1 1 1]; % closed loop only 1=closed, 0 = open, 2 = closed loop w bar jumps
        stimSelect = [2 3 2]; % 1 trial 1s stim (stim=2, v=0.25)
        trialBreak = 5; %sec, break between trial
end

%% check relevant external dependencies
cd(cellDirPath)

% if fictrac, configure and start based on exeriment type (open/closed)
if inUseFictrac
    %determine if open or closed loop experiment. thhis is outdated i think
    %if funcSelect == 1
        ftMode = 0; %set closed-loop
    %else %default
        %ftMode = 0; %set open-loop
    %end
    startFicTrac_smr(cellDirPath,ftMode)
end

% if opto stim, reminder user to check ND filters
if inUseStim
    if sum(stimSelect) > 0
        ndCheck = input('[INPUT] Check LED is on...');
    end
end

connectHost; %required
Panel_com('all_on'); %turn on
checkPanels = input('[INPUT] Check display is working...');
Panel_com('all_off'); %turn off

%% run experiment
cd(cellDirPath)

% general settings:
batteryN = length(pattSelect);
trialPerExpt = input('\n[INPUT] How many trials? (#) ');
trialPerExpt = ceil(trialPerExpt/batteryN); %divide total trials by N battery variants
%trialDuration = input('\n[INPUT] Duration of each trial? (sec) ');
alpha = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'; %for denoting each battery variant

% set necessary filepaths
baseFilename = strrep([erase(cellDirPath,dataDir) '\trial'],'\','_'); %set filename
pc_file = fullfile(cellDirPath,'patchcheck.mat');
if exist(pc_file,'file')
    load('patchcheck.mat','checkVm','checkInputR','checkSR')
else
    checkInputR=[];
    checkVm=[];
    checkSR=[];
end

% begin experiment trials
while 1 %repeat expt until user breaks
    for te = 1:trialPerExpt

        % update trial number with each itteration
        if isempty(trialNum)
            trialNum = 1;
        else
            trialNum = trialNum + 1;
        end

        % RUN EXPERIMENT !!
        % randomize battery order for each trial
        subt = 1; %keeping track of sub trials
        for rp = randperm(length(pattSelect))
            pattN = pattSelect(rp); %pull this pattern
            funcN = funcSelect(rp); %pull this open or closed loop
            if funcN==1 % closed loop
                % SMR changed 7/28, 10 min trials trialDuration = 12000; % 20 min
                % SMR - from 6/10 to 7/28, stim was 10 min
                trialDuration = 12000; % 20 min
            elseif funcN==0 % open loop
                trialDuration = 1200; % 5 min?
            elseif funcN==2 % closed loop bar jumps
                trialDuration = 12000; % 20 min
            end
            stim = stimSelect(:,rp); %pull this stim
            varN = alpha(rp); %pull this variant letter
            disp(['Starting trial ' num2str(trialNum,'%02.f') ' subset ' num2str(subt) '/' num2str(batteryN) ' ' varN ' stim ' num2str(stim) '...'])
            [rawData, inputParams, rawOutput] = exptFn(settings,trialDuration,pattN,funcN,stim, cellDirPath);

            % SAVE RAW TRIAL DATA !!
            inputParams.filename = [baseFilename num2str(trialNum,'%02.f') varN ' pattern' num2str(pattN) ' stim' num2str(stim) ' function' num2str(funcN)];
            save([inputParams.filename '_raw.mat'], 'rawData', 'rawOutput', 'inputParams','-v7.3');
            fprintf('Raw data for trial %02d saved! \n', trialNum);

            % RUN POST TRIAL ROUTINE !!
            [~] = postExptRoutine(inputParams, rawData, rawOutput, settings);
            subt=subt+1; %update counter
            pause(trialBreak) %break between trial
        end
    end

    % prompt user to select
    addTrial = input('\n[INPUT] Run another trial batch? (y/n): ','s');
    if strcmpi(addTrial, 'n')
        disp('Experiment ended by user.');
        close all
        break;
    end

end

%terminate external programs
if inUsePanels
    Panel_com('stop_display');
    Panel_com('all_off');
    disconnectHost;
end
if inUseFictrac
    system('Taskkill/IM cmd.exe');
    disp('Fictrac terminated.')
end



%% INVALID INPUT, DON'T DO ANYTHING
disp('[ERROR] Invalid input. Ending routine.');
return;
end

```
