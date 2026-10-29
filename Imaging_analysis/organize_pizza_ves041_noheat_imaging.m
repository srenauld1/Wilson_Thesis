function T=organize_pizza_ves041_noheat_imaging(out)
% Reorganize existing event figures without changing their selections/content.
if nargin<1,out='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging';end
T=readtable(fullfile(out,'manifest.csv'),'Delimiter',',','TextType','string');
T.category=strings(height(T),1);T.pathColorLow=nan(height(T),1);T.pathColorHigh=nan(height(T),1);
for cat=["bar","gratings"],if ~isfolder(fullfile(out,cat)),mkdir(fullfile(out,cat));end,end
flyRows=cell(0,4);
for k=find(T.status=="complete")'
 name=char(T.trial(k));cat=pizza_ves041_trial_category(name,T.trialFolder(k));T.category(k)=cat;
 dest=string(out);if any(cat==["bar","gratings"]),dest=fullfile(out,cat);end
 newEvent=fullfile(dest,string(name)+"_peaks_slowdowns");newPath=fullfile(dest,string(name)+"_flat_path_dff");
 for stemKey=["eventFigureStem","pathFigureStem"]
  old=T.(stemKey)(k);if stemKey=="eventFigureStem",new=newEvent;else,new=newPath;end
  if old~=new
   for ext=[".png",".fig"],assert(isfile(old+ext),'Missing existing figure %s',old+ext);assert(~isfile(new+ext),'Destination exists: %s',new+ext);movefile(old+ext,new+ext);end
  end
 end
 T.eventFigureStem(k)=newEvent;T.pathFigureStem(k)=newPath;
 source=fullfile(out,'analysis_data',[name '_events.mat']);S=load(source,'d','meta');
 info=plot_pizza_flat_path_scatter(S.d,name,char(newPath),S.meta.pathSource);
 T.pathColorLow(k)=info.low;T.pathColorHigh(k)=info.high;
 if cat=="bar"
  [flyId,basis]=fly_identity(name,char(T.trialFolder(k)));flyRows(end+1,:)={name,flyId,basis,T.trialFolder(k)}; %#ok<AGROW>
 end
 fprintf('%s -> %s | dF/F colors %.3g to %.3g | %d points\n',name,cat,info.low,info.high,info.sampleCount);
end
writetable(T,fullfile(out,'manifest.csv'));
F=cell2table(flyRows,'VariableNames',{'trial','flyId','identityBasis','trialFolder'});writetable(F,fullfile(out,'bar_fly_mapping.csv'));
settingsPath=fullfile(out,'analysis_settings.json');settings=jsondecode(fileread(settingsPath));settings.flatPath=struct('rendering','scatter','markerSize',14,'colorPercentiles',[5 95],'drawOrder','ascending original dF/F; higher values last','samplesRemovedForColor',false);settings.folderLayout='bar and gratings in subfolders; dark and odor in parent';
fid=fopen(settingsPath,'w');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(settings,PrettyPrint=true));fclose(fid);
notes=sprintf('\nApproved organization and flat-path update\nBar figures are in bar/; grating figures are in gratings/. Dark and odor figures remain in the parent. The 4x4 figures and event selections were moved without modification. All flat paths are now filled scatterplots (14-point marker area), with color limits at each trial''s fifth and 95th dF/F percentiles. Values outside those limits saturate at the endpoint colors; samples and their numeric dF/F are retained. Higher dF/F points are drawn last so low-activity revisits do not cover peaks.\n\nBar dataset: %d recordings, %d estimated flies. bar_fly_mapping.csv records explicit versus inferred identities. Fly IDs combine date and within-date animal identity; unnumbered recordings retain inferred labels. This is an estimated biological count, not 13 independent flies.\n',height(F),numel(unique(string(F.flyId))));
fid=fopen(fullfile(out,'README.txt'),'a');assert(fid>=0);fprintf(fid,'%s',notes);fclose(fid);
fprintf('BAR: %d recordings, %d estimated flies\n',height(F),numel(unique(string(F.flyId))));
end

function [id,basis]=fly_identity(name,folder)
token=regexp(name,'^(\d{8})','tokens','once');assert(~isempty(token));date=string(token{1});
token=regexp(lower(folder),'fly(\d+)','tokens','once');
if ~isempty(token),id=date+"_fly"+string(token{1});basis="Date and explicit fly number in recording folder";return,end
if startsWith(name,'20261019_'),id="20261019_fly1";basis="User confirmed: all October 19 recordings from one fly";
elseif startsWith(name,'20261024_'),id="20261024_fly1";basis="User confirmed: all October 24 recordings from one fly";
elseif startsWith(name,'20260217-'),id="20260217_fly1_inferred";basis="Prior fly mapping: sequential clbbar1-3; later recordings explicitly label fly2";
elseif startsWith(name,'20260213-2_'),id="20260213_fly1_inferred";basis="Unnumbered first fly inferred from later same-date fly2 folders";
else,id=date+"_unlabeled_inferred";basis="One included bar recording on this date; animal ID absent, distinct-date animal assumed";end
end
