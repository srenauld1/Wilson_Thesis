function context=classify_ves041_peak_heat_context(on,control)
% Classify a complete, finite +/-0.5 s command window around an activity peak.
% Entry/exit use measured-command threshold crossings, not geometric zones.
on=logical(on(:));control=control(:);steps=diff(double(on));
if all(control==1)
    context="control";
elseif sum(steps~=0)>1
    context="multiple switches";
elseif any(steps==1)
    context="heat entry";
elseif any(steps==-1)
    context="heat exit";
elseif any(control==1)
    context="control transition";
elseif all(on)
    context="heat";
else
    context="no heat";
end
end
