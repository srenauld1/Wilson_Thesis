function angle=pizza_ves041_turn_angle(d,idx,window)
% Integrate absolute 0.25-second median-smoothed yaw over exact time bounds.
t=d.t;fs=1/median(diff(t));valid=d.valid&isfinite(d.yaw);cut=[true;diff(t)>2/fs];
s=find(valid&(cut|[true;~valid(1:end-1)]));e=find(valid&([~valid(2:end);true]|[cut(2:end);true]));r=nan(size(t));
for k=1:numel(s),ix=s(k):e(k);r(ix)=smoothdata(d.yaw(ix),'movmedian',max(3,round(.25*fs)));end
angle=nan(numel(idx),1);
for k=1:numel(idx)
 bounds=t(idx(k))+window;inside=t>bounds(1)&t<bounds(2);tt=[bounds(1);t(inside);bounds(2)];
 yy=interp1(t,abs(r),tt,'linear',NaN);assert(all(isfinite(yy)),'Invalid rotation window');angle(k)=trapz(tt,yy);
end
end
