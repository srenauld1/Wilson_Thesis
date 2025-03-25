% fetchTimeIdx
% This function converts a desired time point (tsize) to the closest index in a given 
% time array (ttime). It finds the index corresponding to the time that is nearest 
% to the specified time size.
%
% INPUTS:
%   ttime - Array of time points (in seconds)
%   tsize - Desired time point (in seconds) to find the closest index for
%
% OUTPUTS:
%   tidx  - Index of the closest time point in the time array
%
% CREATED: [Date] MC
%
function [tidx] = fetchTimeIdx(ttime,tsize)
%find nearest index
[~,tidx] = min(abs(ttime - tsize));
end