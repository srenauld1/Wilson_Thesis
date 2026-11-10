function [idx,angle,before,after,info]=pizza_ves041_slowdowns_no_accel_threshold(d)
% Select forward-velocity drops directly, without an acceleration threshold.
% Keep original 0.25-s median smoothing; before median [-.75,-.25), after
% median [.25,.75]; before>=1 mm/s and before-after>=1 mm/s. Full -3:+3 s
% windows must be valid and continuous. One candidate per contiguous run of
% qualifying velocity drops; time zero is maximum deceleration within that
% run ONLY for alignment, with NO required acceleration magnitude or sign.
% Accepted centers >=2 s apart. Rotation split retains accumulated absolute
% 0.25-s median-smoothed yaw >=15 degrees during [0,.5] s after the center.
% No fluorescence is used to choose centers (finite dF/F is needed for plots).
t=d.t(:);fwd=d.fwd(:);fs=1/median(diff(t));
valid=d.valid(:)&isfinite(d.dff(:))&isfinite(fwd)&isfinite(d.yaw(:));
gap=[true;diff(t)>2/fs];s=find(valid&(gap|[true;~valid(1:end-1)]));
e=find(valid&([~valid(2:end);true]|[gap(2:end);true]));
f=nan(size(t));acc=f;pre=f;post=f;eligible=false(size(t));
for k=1:numel(s)
    ix=s(k):e(k);f(ix)=smoothdata(fwd(ix),'movmedian',max(3,round(.25*fs)));
    if numel(ix)>2,acc(ix)=gradient(f(ix),t(ix));end
    % Include native bracketing samples at both window endpoints.
    centers=ix(t(ix)>=t(s(k))+3&t(ix)<=t(e(k))-3);
    a1=s(k);a2=s(k)-1;b1=s(k);b2=s(k)-1;
    for i=centers
        while a1<=e(k)&&t(a1)<t(i)-.75,a1=a1+1;end
        while a2<e(k)&&t(a2+1)<t(i)-.25,a2=a2+1;end
        while b1<=e(k)&&t(b1)<t(i)+.25,b1=b1+1;end
        while b2<e(k)&&t(b2+1)<=t(i)+.75,b2=b2+1;end
        pre(i)=median(f(a1:a2));post(i)=median(f(b1:b2));eligible(i)=true;
    end
end
candidate=eligible&pre>=1&pre-post>=1;
change=diff([false;candidate;false]);on=find(change==1);off=find(change==-1)-1;
idx=zeros(0,1);before=idx;after=idx;last=-inf;
for k=1:numel(on)
    [~,j]=min(acc(on(k):off(k)));i=on(k)+j-1;
    if t(i)-last<2,continue,end
    idx(end+1,1)=i;before(end+1,1)=pre(i);after(end+1,1)=post(i);last=t(i); %#ok<AGROW>
end
angle=pizza_ves041_turn_angle(d,idx,[0 .5]);
info=struct('accelerationThreshold','none; acceleration used only to align within velocity-drop runs', ...
    'minimumVelocityBeforeMmS',1,'minimumVelocityDropMmS',1,'smoothSec',.25, ...
    'preWindowSec',[-.75 -.25],'postWindowSec',[.25 .75],'eventWindowSec',[-3 3], ...
    'spacingSec',2,'rotationWindowSec',[0 .5],'rotationThresholdDeg',15);
end
