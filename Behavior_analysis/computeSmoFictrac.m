% computeSmoFictrac.m
%
% Function that takes in one or more pData files and computes a very
%  smoothed version of the FicTrac velocity parameters. To be used for 
%  extracting turning bouts, not to be used as faithful representation of
%  actual velocities.
% Select pData files through GUI. These must have fictrac, fictracProc,
%  and fictracParams structs
% Saves smoothed data back into same pData file
%
% INPUTS:
%   pDataPath - full path to folder containing pData files
%   sigmaPos - sigma for smoothing position data using smoothdata(), in
%       fictrac samples
%   sigmaVel - sigma for smoothing velocity data using smoothdata(), in
%       fictrac samples
% 
% OUTPUTS:
%   none, but saves fictracSmo output struct back into same pData file
%
% CREATED: 6/22/23
%
% UPDATED: 
%   6/22/23 - HHY
%                   
%
function computeSmoFictrac(pDataPath, sigmaPos, sigmaVel)

    % default values of sigmaPos and sigmaVel
    if (isempty(sigmaPos))
        sigmaPos = 500;
    end
    if (isempty(sigmaVel))
        sigmaVel = 500;
    end

    % prompt user to select pData files
    [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
        'Select pData files', pDataPath, 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has fictrac, fictracProc, fictracParams
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'fictrac')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracParams')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')))
            continue;
        end

        % load variables from pData
        load(pDataFullPath, 'fictrac', 'fictracParams', 'fictracProc');

        % smooth position
        angPosUn = unwrap(fictrac.yawAngPosWrap);
        angPosUnSmo = smoothdata(angPosUn, 'gaussian',sigmaPos);
        fwdPosSmo = smoothdata(fictrac.fwdCumPos, 'gaussian', sigmaPos);
        slidePosSmo = smoothdata(fictrac.slideCumPos, 'gaussian', sigmaPos);

        % compute velocity
        angVel = gradient(angPosUnSmo, fictrac.t);
        fwdVel = gradient(fwdPosSmo, fictrac.t);
        slideVel = gradient(slidePosSmo, fictrac.t);

        % smooth velocity
        angVelSmo = smoothdata(angVel, 'gaussian',sigmaVel);
        fwdVelSmo = smoothdata(fwdVel, 'gaussian', sigmaVel);
        slideVelSmo = smoothdata(slideVel, 'gaussian', sigmaVel);

        % get bias for yaw and slide velocities - slope of line fit to
        % position over the whole trial
        pAng = polyfit(fictrac.t, angPosUnSmo,1);
        angVelBias = pAng(1);
        pSlide = polyfit(fictrac.t, slidePosSmo, 1);
        slideVelBias = pSlide(1);

        
        % interpolate to same time basis as fictracProc
        fictracSmo.yawAngVel = interp1(fictrac.t, angVelSmo, ...
            fictracProc.t, 'spline');
        fictracSmo.fwdVel = interp1(fictrac.t, fwdVelSmo, ...
            fictracProc.t, 'spline');
        fictracSmo.slideVel = interp1(fictrac.t, slideVelSmo, ...
            fictracProc.t, 'spline');

        % add time info to fictracSmo output struct
        fictracSmo.t = fictracProc.t;

        % copy over dropInd from fictracProc
        fictracSmo.dropInd = fictracProc.dropInd;

        % add bias info to fictracSmo output struct
        fictracSmo.angVelBias = angVelBias;
        fictracSmo.slideVelBias = slideVelBias;

        % add sigma info to fictracSmo output struct
        fictracSmo.sigmaPos = sigmaPos;
        fictracSmo.sigmaVel = sigmaVel;

        % append this fictracSmo output struct to same pData file
        save(pDataFullPath, 'fictracSmo', '-append');

        fprintf('Saved fictracSmo for %s!\n', pDataName);

    end
end