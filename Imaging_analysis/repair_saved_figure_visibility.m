function report=repair_saved_figure_visibility(folder)
% Make saved FIGs open normally; never show a window during the repair.
% With no argument, repair only smr_analysis FIGs from October 19 and 20.
% Writes through a temporary file; source is replaced only after validation.
% Delegate desktop calls to a no-display process. The desktop's root defaults
% and existing figures are left alone.
oldDefault=get(groot,'defaultFigureVisible');
restoreDefault=onCleanup(@()set(groot,'defaultFigureVisible',oldDefault)); %#ok<NASGU>
if feature('ShowFigureWindows')
    csv=[tempname '.csv'];here=fileparts(mfilename('fullpath'));
    if nargin==0,call='repair_saved_figure_visibility()';
    else,call=sprintf('repair_saved_figure_visibility(''%s'')',strrep(folder,'''',''''''));end
    expr=sprintf('addpath(''%s'');r=%s;writetable(r,''%s'');',strrep(here,'''',''''''),call,strrep(csv,'''',''''''));
    cmd=[shellquote(fullfile(matlabroot,'bin','matlab')) ' -nodisplay -batch ' shellquote(expr)];
    [status,msg]=system(cmd);assert(status==0,'Headless figure repair failed: %s',msg);
    report=readtable(csv,'TextType','string');delete(csv);return
end
if nargin==0
    files=[dir('/Users/sophiarenauld/stacks/2026-10-19/**/smr_analysis/*.fig'); ...
        dir('/Users/sophiarenauld/stacks/2026-10-20/**/smr_analysis/*.fig')];
else
    files=dir(fullfile(folder,'*.fig'));
end
rows=cell(numel(files),3);
for k=1:numel(files)
    path=fullfile(files(k).folder,files(k).name);fig=[];tmp='';
    try
        % Inspect legacy visibility without loading the saved HG2 objects.
        saved=load(path,'hgS_070000','-mat');
        if isfield(saved,'hgS_070000') && ...
                all(arrayfun(@(v)isfield(v.properties,'Visible') && ...
                strcmp(v.properties.Visible,'on'),saved.hgS_070000))
            rows(k,:)={string(path),"already visible",""};continue
        end
        fig=openfig(path,'new','invisible');
        nAxes=numel(findall(fig,'Type','axes'));nLines=numel(findall(fig,'Type','line'));
        tmp=[tempname(files(k).folder) '.fig'];
        savefig_visible_on_open(fig,tmp);close(fig);fig=[];
        % No-display MATLAB forces loaded figures off, even when saved on.
        % Verify the stored setting, then round-trip the plot objects hidden.
        stored=load(tmp,'hgS_070000','-mat');
        assert(all(arrayfun(@(v)strcmp(v.properties.Visible,'on'),stored.hgS_070000)), ...
            'Saved figure is still hidden.');
        check=openfig(tmp,'new','invisible');fig=check;
        assert(numel(findall(check,'Type','axes'))==nAxes && ...
            numel(findall(check,'Type','line'))==nLines,'Plot objects changed during round trip.');
        close(check);fig=[];
        [ok,msg]=movefile(tmp,path,'f');assert(ok,'%s',msg);tmp='';
        rows(k,:)={string(path),"repaired",""};
    catch ME
        if ~isempty(fig),close(fig);end
        if ~isempty(tmp)&&isfile(tmp),delete(tmp);end
        rows(k,:)={string(path),"error",string(ME.message)};
        warning('imaging:FigureVisibility','%s: %s',path,ME.message);
    end
    if mod(k,10)==0,fprintf('Figure visibility: %d/%d checked\n',k,numel(files));end
end
report=cell2table(rows,'VariableNames',{'file','status','error'});
assert(~any(report.status=="error"),'Some figures could not be repaired; inspect the report.');
end
function q=shellquote(s)
quote=char(39);q=[quote strrep(s,quote,[quote '"' quote '"' quote]) quote];
end
