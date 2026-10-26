function savefig_visible_on_open(fig,filename)
% Save with Visible=on, without displaying a window during hidden analysis.
% Hidden figures require a -nodisplay MATLAB process; do not toggle visibility
% in a desktop process, where doing so would briefly create a window.
oldVisible=get(fig,'Visible');
assert(~feature('ShowFigureWindows') || strcmp(oldVisible,'on'), ...
    'imaging:HeadlessRequired','Hidden figure saving requires MATLAB -nodisplay.');
cleanup=onCleanup(@()restoreState(fig,oldVisible)); %#ok<NASGU>
set(fig,'Visible','on');
savefig(fig,filename);
end
function restoreState(fig,oldVisible)
if isgraphics(fig),set(fig,'Visible',oldVisible);end
end
