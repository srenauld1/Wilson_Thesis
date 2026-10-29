function files=verify_heat_exploration_outputs(out)
% Required equivalent of the October 20 heat-exploration deliverables.
stems={'slowdown_heat_comparison','slowdown_transition_control','activity_vs_movement'};
paths=dir(fullfile(out,'xy_heat_paths_*.png'));
assert(~isempty(paths),'imaging:MissingHeatExploration','Missing XY heat-geometry figure in %s',out);
for k=1:numel(paths),[~,stem]=fileparts(paths(k).name);stems{end+1}=stem;end %#ok<AGROW>
files=strings(0,1);
for k=1:numel(stems)
    for ext=[".png" ".fig"]
        file=fullfile(out,string(stems{k})+ext);info=dir(file);
        assert(~isempty(info)&&info.bytes>0,'imaging:MissingHeatExploration','Missing or empty required heat-exploration output: %s',file);
        files(end+1,1)=file; %#ok<AGROW>
    end
end
for name=["slowdown_events.csv" "analysis_settings.json"]
    file=fullfile(out,name);info=dir(file);
    assert(~isempty(info)&&info.bytes>0,'imaging:MissingHeatExploration','Missing or empty required heat-exploration output: %s',file);
    files(end+1,1)=file; %#ok<AGROW>
end
end
