% dsFiltFictrac.m
%
% Function that downsamples and smooths fictrac data and returns new fwd,
%  slide, and yaw postion and velocity, x and y position, total speed, and
%  yaw speed, as well as t and dropInd in fictrac struct.
%
% Adapted from function of same name in 2PAnalysisCode-Helen
%
% INPUT:
%   fictracParams - struct of fictrac parameters, contains:
%       dsf - downsampling factor
%       filtParams - struct of parameters for smoothing all fictrac
%               data
%           padLen - padding length for convolving with Gaussian kernel
%           sigmaPos - standard deviation for Gaussian kernel for
%               position
%           sigmaVel - standard deviation for Gaussian kernel for
%               velocity
%       degPerMM - degrees per millimeter, conversion factor
%       mmPerDeg - millimeters per degree, conversion factor
%   fictrac - struct of fictrac data, output of preprocessFicTrac(),
%       selectDroppedFicTrac(), writePData()
%
% OUTPUT:
%   fictracProc - struct of processed fictrac data, contains:
%       fwdVel
%       slideCumPos
%       slideVel
%       yawAngCumPos
%       yawAngPosWrap
%       yawAngVel
%       yawAngSpd
%       totAngSpd
%       totSpd
%       xPos
%       yPos
%       t
%       dropInd
%
% CREATED: 9/11/20 - HHY
%
% UPDATED:
%   9/11/20 - HHY
%   8/2/22 - HHY - commented out portion interpolates to 1000 Hz,
%       regardless of input params (for when fictrac from vid)
%

function fictracProc = dsFiltFictrac(fictracParams, fictrac)

% downsample to 1000 Hz, all position variables are cumulative
fwdPosDS = downsample(fictrac.fwdCumPos, fictracParams.dsf);
yawPosDS = downsample(fictrac.yawAngPos, fictracParams.dsf);
slidePosDS = downsample(fictrac.slideCumPos, fictracParams.dsf);

% downsample timing vector
tDS = downsample(fictrac.t, fictracParams.dsf);
ftSampRate = 1/median(diff(tDS));

% perform smoothing on position and extract smoothed velocity
[fwdPosSmo, fwdVelSmo] = computeSmoothedVelocity(fwdPosDS, ...
    fictracParams.filtParams.padLen, ...
    fictracParams.filtParams.sigmaPos, ...
    fictracParams.filtParams.sigmaVel, ftSampRate);
[yawAngPosSmo, yawAngVelSmo] = computeSmoothedVelocity(yawPosDS, ...
    fictracParams.filtParams.padLen, ...
    fictracParams.filtParams.sigmaPos, ...
    fictracParams.filtParams.sigmaVel, ftSampRate);
[slidePosSmo, slideVelSmo] = computeSmoothedVelocity(slidePosDS, ...
    fictracParams.filtParams.padLen, ...
    fictracParams.filtParams.sigmaPos, ...
    fictracParams.filtParams.sigmaVel, ftSampRate);

% compute wrapped yaw position
yawAngPosSmoWrap = wrapTo360(yawAngPosSmo);

% compute yaw speed
yawAngSpd = abs(yawAngVelSmo);

% compute total speed, in degrees per second
% convert forward and slide velocities from mm to deg per sec
fwdVelDeg = fwdVelSmo .* fictracParams.degPerMM;
slideVelDeg = slideVelSmo .* fictracParams.degPerMM;

totAngSpd = yawAngSpd + abs(fwdVelDeg) + abs(slideVelDeg);
totSpd = totAngSpd .* fictracParams.mmPerDeg;

% compute x and y positions, using smoothed positions and velocities
% start with fly at (0,0) and facing 0 deg
zeroedYawAngPos = yawAngPosSmo - yawAngPosSmo(1);

% movement in x (in degrees) at each time point
xChangePos = (fwdVelDeg ./ ftSampRate) .* sind(zeroedYawAngPos) + ...
    (slideVelDeg ./ ftSampRate) .* sind(zeroedYawAngPos + 90);

% x position in mm (i.e. x-coordinate of fly's position at each time
%  point), starts at 0
xPos = (cumsum(xChangePos) - xChangePos(1)) .* fictracParams.mmPerDeg;

% movement in y (in degrees) at each time point
yChangePos = (fwdVelDeg ./ ftSampRate) .* cosd(zeroedYawAngPos) + ...
    (slideVelDeg ./ ftSampRate) .* cosd(zeroedYawAngPos + 90);

% y position in mm (i.e. y-coordinate of fly's position at each time
%  point), starts at 0
yPos = (cumsum(yChangePos) - yChangePos(1)) .* fictracParams.mmPerDeg;


% save all data into fictrac struct
fictracProc.fwdCumPos = fwdPosSmo;
fictracProc.fwdVel = fwdVelSmo;
fictracProc.slideCumPos = slidePosSmo;
fictracProc.slideVel = slideVelSmo;
fictracProc.yawAngCumPos = yawAngPosSmo;
fictracProc.yawAngPosWrap = yawAngPosSmoWrap;
fictracProc.yawAngVel = yawAngVelSmo;
fictracProc.yawAngSpd = yawAngSpd;
fictracProc.totAngSpd = totAngSpd;
fictracProc.totSpd = totSpd;
fictracProc.xPos = xPos;
fictracProc.yPos = yPos;
fictracProc.t = tDS'; % change to row vector

end