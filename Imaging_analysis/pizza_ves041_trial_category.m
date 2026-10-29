function category=pizza_ves041_trial_category(trial,folder)
% Identify sensory context before considering generic bar text in names.
s=lower(string(trial)+" "+string(folder));
if contains(s,'odor'),category="odor";
elseif contains(s,'dark'),category="dark";
elseif contains(s,'grating'),category="gratings";
elseif contains(s,'bar'),category="bar";
else,error('pizza:UnknownContext','Unrecognized trial context: %s',trial);end
end
