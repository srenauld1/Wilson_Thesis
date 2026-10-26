function report=run_ves041_oct19_oct20(refreshPrefixes)
% Run four representative contexts first; batch only after all four pass.
% Optional recording-name prefixes refresh selected trials and replace their
% report rows, preserving the other completed trials in the report.
oldFigureVisible=get(groot,'defaultFigureVisible');
restoreFigures=onCleanup(@()set(groot,'defaultFigureVisible',oldFigureVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
here=fileparts(mfilename('fullpath'));addpath(here);addpath(genpath('/Users/sophiarenauld/Documents/GitHub/scopa'));
out='/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis';if ~isfolder(out),mkdir(out);end
files=[dir('/Users/sophiarenauld/stacks/2026-10-19/**/*_or_z0z1_s_.mat');dir('/Users/sophiarenauld/stacks/2026-10-20/**/*_or_z0z1_s_.mat')];
paths=string(fullfile({files.folder},{files.name}))';
pilotPrefixes=["20261020_4_2_","20261020_7_1_","20261019_4_1_","20261020_6_1_"];
pilot=[];for p=pilotPrefixes,k=find(startsWith(string({files.name}),p),1);assert(~isempty(k),'Missing pilot %s',p);pilot(end+1)=k;end
order=[pilot,setdiff(1:numel(files),pilot,'stable')];rows=cell(0,6);
refresh=nargin>0&&~isempty(refreshPrefixes);
if refresh
    order=find(startsWith(string({files.name}),string(refreshPrefixes)));
    previous=fullfile(out,'run_report.csv');
    if isfile(previous)
        opts=delimitedTextImportOptions('NumVariables',6);
        opts.DataLines=[2 Inf];opts.Delimiter=',';
        opts.VariableNames={'input','phase','status','seconds','error','lastWarning'};
        opts.VariableTypes={'string','string','string','double','string','string'};
        old=readtable(previous,opts);
        old.error=string(old.error);old.error(ismissing(old.error)|old.error=="NaN")="";
        old.lastWarning=string(old.lastWarning);old.lastWarning(ismissing(old.lastWarning))="";
        old=old(~ismember(old.input,paths(order)),:);rows=table2cell(old);
    end
end
for n=1:numel(order)
    k=order(n);phase="batch";if n<=4,phase="pilot";end;if refresh,phase="refresh";end
    log=fullfile(out,[files(k).name '_analysis.log']);diary(log);lastwarn('');
    fprintf('\nSTART %s %s\n',phase,paths(k));tic;
    try
        run_one(char(paths(k)),here);status="success";err="";
    catch ME
        status="error";err=string(getReport(ME,'extended','hyperlinks','off'));fprintf('%s\n',err);
    end
    elapsed=toc;[warn,~]=lastwarn;fprintf('END %s %s (%.1f s)\n',paths(k),status,elapsed);diary off;close all force;
    rows(end+1,:)={paths(k),phase,status,elapsed,err,string(warn)}; %#ok<AGROW>
    report=cell2table(rows,'VariableNames',{'input','phase','status','seconds','error','lastWarning'});writetable(report,fullfile(out,'run_report.csv'));
    if ~refresh&&n==4&&any(report.status=="error"),error('ves041:PilotFailed','Pilot errors recorded; full batch was not started.');end
end
end
function run_one(input,here)
fullFilePath=input; %#ok<NASGU>
run(fullfile(here,'master_imaging_analysis.m'));
end
