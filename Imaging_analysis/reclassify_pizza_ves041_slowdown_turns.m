function reclassify_pizza_ves041_slowdown_turns(out)
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
M=readtable(fullfile(out,'manifest.csv'),'Delimiter',',','TextType','string');ix=find(M.status=="complete");
backup=fullfile(out,'analysis_data','previous_symmetric_turn_definition');if ~isfolder(backup),mkdir(backup);end
for k=ix'
 trial=M.trial(k);file=fullfile(out,'analysis_data',trial+"_events.mat");S=load(file);idx=S.groups{2};
 angle=pizza_ves041_turn_angle(S.d,idx,[0 .5]);turn=angle>=15;
 for suffix=["_events.mat","_events.csv","_events_with_rotation_peaks.mat"]
  f=fullfile(out,'analysis_data',trial+suffix);[~,name,ext]=fileparts(f);dest=fullfile(backup,string(name)+ext);if ~isfile(dest),copyfile(f,dest);end
 end
 for g=1:2,S.groups{g+2}=idx(turn==(g==2));S.X{g+2}=S.X{2}(:,turn==(g==2),:);end
 slow=S.E.eventType=="slowdown";assert(nnz(slow)==numel(idx));assert(isequal(S.E.sampleIndex(slow),idx));
 S.E.absoluteTurnDeg(slow)=angle;S.E.hasTurn(slow)=double(turn);S.settings.turnWindowSec=[0 .5];
 S.settings.turnDefinition='Accumulated absolute yaw over exact 0 to +0.5 s after maximum deceleration';
 save(file,'-struct','S','-v7.3');writetable(S.E,fullfile(out,'analysis_data',trial+"_events.csv"));
 f=fullfile(out,'analysis_data',trial+"_events_with_rotation_peaks.mat");R=load(f);
 for c=1:4,R.groups{c}=S.groups{c};R.X{c}=S.X{c};end
 R.slowdownTurnWindowSec=[0 .5];save(f,'-struct','R','-v7.3');
 assert(isequal(sort([S.groups{3};S.groups{4}]),idx));
 M.slowdownsWithoutTurn(k)=nnz(~turn);M.slowdownsWithTurn(k)=nnz(turn);
end
writetable(M,fullfile(out,'manifest.csv'));
end
