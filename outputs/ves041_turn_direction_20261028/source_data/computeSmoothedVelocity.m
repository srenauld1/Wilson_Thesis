% computeSmoothedVelocity.m
%
% Helper function called by dsFiltFictrac.m that takes in downsampled
%  FicTrac position and smoothing parameters and returns smoothed position
%  and velocity. Performs Gaussian kernel smoothing on position, computes
%  velocity, and performs Gaussian kernel smoothing on velocity.
%
% INPUT:
%   pos - input position values
%   padLen - padding length for Gaussian kernel smoothing, shared for
%       smoothing on position and on velocity
%   sigmaPos - standard deviation of Gaussian kernel for position
%   sigmaVel - standard deviation of Gaussian kernel for velocity
%   sampRate - sampling rate of pos, in Hz
%
% OUTPUT:
%   smoPos - smoothed position values
%   smoVel - smoothed velocity values
%
% CREATED: 8/28/19 - HHY
%
% UPDATED: 8/28/19 - HHY
%   10/4/21 - HHY - matlab code for smoothing, not python
%

function [smoPos, smoVel] = computeSmoothedVelocity(pos, padLen, ...
    sigmaPos, sigmaVel, sampRate)
   
    smoPos = gaussSmooth(pos, padLen, sigmaPos);

    % compute velocity
    vel = gradient(smoPos) .* sampRate;
    
    smoVel = gaussSmooth(vel,padLen,sigmaVel);
end