function run_blinded_5s_server_master(manifest_path,worker,workers)
% Run the existing master once per copied recording, with resumable status.
if nargin<2,worker=1;end
if nargin<3,workers=1;end
records = jsondecode(fileread(manifest_path));
root = '/Volumes/Neurobio/wilsonlab/Sophia/behavior/RubyACR Data/5s stim';
status_path = fullfile(root,sprintf('master_status_worker%d.json',worker));
status = struct('fly_id',{},'source',{},'status',{},'message',{});
set(groot,'defaultFigureVisible','off');
for i = 1:numel(records)
    if mod(i-1,workers)~=worker-1,continue;end
    r = records(i);
    item = struct('fly_id',r.cohort_ids{1},'source',r.destination,'status','pending','message','');
    % Copy phase writes recording_info only AFTER checksum verification.
    started=tic;
    while ~isfile(fullfile(fileparts(r.destination),'recording_info.json'))
        assert(toc(started)<18000,'Timed out awaiting verified source copy.');
        pause(2);
    end
    if ~strcmp(r.status,'ok')
        item.status = 'unassessable';item.message = 'Source processed MAT file truncated; copied unchanged, not reconstructed.';
    else
        try
            [parent,name] = fileparts(r.destination);
            parts = strsplit(name,' ');
            output = fullfile(parent,'trial01Bstim5');
            final = fullfile(output,[parts{1} '_processed_final.mat']);
            done = fullfile(output,'MASTER_COMPLETE.txt');
            if isfile(done) && isfile(final)
                item.status = 'complete';item.message = 'Previously verified completed output retained.';
            else
                fprintf('\nSTART %d/%d ID %s\n',i,numel(records),item.fly_id);
                run_one(r.destination);
                assert(isfile(final),'Master did not produce its final MAT output.');
                required={'forward_vs_rotation.png','rotation_stop_vs_rotation_start.png','fwd_stop_vs_fwd_start.png','flatpath_vs_opto.png','flatpath_aligned_2s_rawY_vs_opto.png','flatpath_aligned_2s_mirrored_vs_opto.png'};
                for j=1:numel(required),assert(isfile(fullfile(output,required{j})),'Missing master figure.');end
                variables = whos('-file',final);
                assert(all(ismember({'exptData','exptMeta'},{variables.name})),'Final MAT missing data or metadata.');
                fid=fopen(done,'w');assert(fid>=0);fprintf(fid,'Master completed; expected PNGs and final exptData/exptMeta verified.\n');fclose(fid);
                item.status = 'complete';item.message = 'Master outputs saved and verified.';
            end
        catch ME
            item.status = 'failed';item.message = getReport(ME,'extended','hyperlinks','off');
            fprintf(2,'FAILED ID %s: %s\n',item.fly_id,ME.message);
        end
    end
    status(end+1)=item; %#ok<AGROW>
    fid=fopen(status_path,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(status,PrettyPrint=true));fclose(fid);
    fprintf('FINISH ID %s: %s\n',item.fly_id,item.status);
    close all;
end
end

function run_one(path)
master_behavior_config.files={path};
master_behavior_analysis;
end
