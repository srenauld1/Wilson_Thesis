% schmittTrigger
% implements a Schmitt trigger and pulls indices accordingly. Once x exceeds
% the high threshold, the trigger will remain active until x falls below the low
% threshold.
%
% INPUTS:
% - xdata : Dataset to apply the trigger against (2D array).
% - highT : High threshold to cross in order to activate the trigger.
% - lowT  : Low threshold to cross in order to deactivate the trigger.
%
% OUTPUT:
% - triggerIdx : Binary matrix indicating the state of the trigger (1 = active, 
%                 0 = inactive) for each element in the dataset.
%
% CREATED: 08/05/2022 by MC
%
function triggerIdx = schmittTrigger(xdata,highT,lowT)
% initialize
[x,n] = size(xdata);
triggerIdx = zeros(x,n);

for e = 1:n
    highIdx = xdata(:,e)>=highT;
    lowIdx = xdata(:,e)>=lowT;

    triggerIdx(:,e) = highIdx; %initialize w/all high points
    trig=0; %set trigger off
    for i=1:x
        % turn trigger ON if above HIGH
        if highIdx(i)==1
            trig=1;
        end
        % if above LOW add point
        if lowIdx(i)==1
            if trig==1
                triggerIdx(i,e)=1;
            end
            % else, turn trigger OFF if below LOW
        else
            trig=0;
        end
    end

end