function result=require_heat_exploration(a,savepath,recordingName)
% Mandatory heat-exploration stage for every cell type in the master pipeline.
% Do not catch plot failures or silently skip recordings lacking geometry.
result=struct('required',false,'complete',false,'outputFolder','');
isHeat=isfield(a,'protocol')&&isfield(a.protocol,'kind')&&startsWith(a.protocol.kind,'heat-');
if isfield(a.dq(1),'heat'),isHeat=isHeat||any(isfinite(a.dq(1).heat(:))&a.dq(1).heat(:)>.05);end
if ~isHeat,return,end
result.required=true;
assert(isfield(a,'protocol')&&isfield(a.protocol,'geometryAligned')&&a.protocol.geometryAligned, ...
    'imaging:RequiredHeatGeometry','Heat processing requires aligned spatial protocol geometry before the mandatory heat_exploration plots can be generated.');
result.outputFolder=fullfile(savepath,'heat_exploration');
plot_ves041_heat_exploration(a,result.outputFolder,recordingName);
result.files=verify_heat_exploration_outputs(result.outputFolder);
result.complete=true;
end
