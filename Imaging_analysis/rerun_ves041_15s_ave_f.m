function report=rerun_ves041_15s_ave_f()
% Reanalyze updated 15-second local-baseline dF/F without replacing the old
% fixed-baseline outputs. Plotting/peak settings are unchanged for comparison.
oldVisible=get(groot,'defaultFigureVisible');
cleanup=onCleanup(@()set(groot,'defaultFigureVisible',oldVisible)); %#ok<NASGU>
set(groot,'defaultFigureVisible','off');
here=fileparts(mfilename('fullpath'));addpath(here);addpath(genpath('/Users/sophiarenauld/Documents/GitHub/scopa'));
out='/Users/sophiarenauld/stacks/ves041_oct19_oct20_analysis/15s_ave_f';if ~isfolder(out),mkdir(out);end
files=[dir('/Users/sophiarenauld/stacks/2026-10-19/**/*_or_z0z1_s_.mat');dir('/Users/sophiarenauld/stacks/2026-10-20/**/*_or_z0z1_s_.mat')];
pilots=["20261020_4_2_","20261020_7_1_","20261019_4_1_","20261020_6_1_"];order=[];
for p=pilots,order(end+1)=find(startsWith(string({files.name}),p),1);end
order=[order,setdiff(1:numel(files),order,'stable')];rows=cell(0,5);
for n=1:numel(order)
    k=order(n);input=fullfile(files(k).folder,files(k).name);output=fullfile(files(k).folder,'15s_ave_f');
    diary(fullfile(out,[files(k).name '.log']));fprintf('START %d/%d %s\n',n,numel(files),input);tic;
    try
        a=load(input,'roi');assert(strcmp(a.roi.opt.nrm.nrmstr,'dff005015'),'Expected 15 s fifth-percentile baseline');clear a
        run_one(input,output,here);status="success";err="";
    catch ME
        status="error";err=string(getReport(ME,'extended','hyperlinks','off'));fprintf('%s\n',err);
    end
    seconds=toc;close all force;fprintf('END %s %.1f s\n',status,seconds);diary off;
    rows(end+1,:)={string(input),string(output),status,seconds,err}; %#ok<AGROW>
    report=cell2table(rows,'VariableNames',{'input','output','status','seconds','error'});writetable(report,fullfile(out,'run_report.csv'));
    if n==4&&any(report.status=="error"),error('ves041:LocalBaselinePilotFailed','Pilot error; batch halted.');end
end
compare_ves041_15s_ave_f(files,out);
end
function run_one(input,output,here)
fullFilePath=input;analysisOutputFolder=output;analysisOmitImageStack=true; %#ok<NASGU>
run(fullfile(here,'master_imaging_analysis.m'));
end
