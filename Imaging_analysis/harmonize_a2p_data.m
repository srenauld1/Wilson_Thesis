function [a2p_data, compatibility] = harmonize_a2p_data(a2p_data)
%HARMONIZE_A2P_DATA Normalize old and new A2P files for shared analyses.
% Older files contain dq(1) aligned to imaging and dq(2) at the original
% behavioral sampling rate. Newer integrated files may contain only the
% imaging-aligned dq(1). In that case dq(1) is copied to dq(2), allowing
% existing analyses to run at the available sampling rate. Old two-stream
% files are left unchanged.

% Spatial trajectories and visual heading are deliberately not fabricated.
% Callers can use the returned flags to skip plots that require them.

% compatibility fields:
%   synthesizedSupplementalDq - dq(2) was copied from dq(1)
%   hasPosition               - px/py or pxb/pyb are populated
%   hasVisualHeading          - vh is populated

if ~isstruct(a2p_data) || ~isfield(a2p_data,'dq') || isempty(a2p_data.dq)
    error('harmonize_a2p_data:MissingDq','Input does not contain a nonempty dq structure.');
end

compatibility = struct();
compatibility.synthesizedSupplementalDq = numel(a2p_data.dq) < 2;
if compatibility.synthesizedSupplementalDq
    a2p_data.dq(2) = a2p_data.dq(1);
end

compatibility.hasPosition = ...
    has_nonempty_pair(a2p_data.dq(1),'px','py') || ...
    has_nonempty_pair(a2p_data.dq(1),'pxb','pyb');
compatibility.hasVisualHeading = ...
    isfield(a2p_data.dq(1),'vh') && ~isempty(a2p_data.dq(1).vh);

a2p_data.compatibility = compatibility;
end

function tf = has_nonempty_pair(s,xField,yField)
tf = isfield(s,xField) && isfield(s,yField) && ...
    ~isempty(s.(xField)) && ~isempty(s.(yField));
end
