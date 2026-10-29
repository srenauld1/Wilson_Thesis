function [idx,angle,before,after,info]=pizza_ves041_speedups(d)
% Mirror the slowdown selection: positive acceleration and forward-speed gain.
t=d.t;fs=1/median(diff(t));valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);
cut=[true;diff(t)>2/fs];s=find(valid&(cut|[true;~valid(1:end-1)]));e=find(valid&([~valid(2:end);true]|[cut(2:end);true]));
f=nan(size(t));r=f;acc=f;n=max(3,round(.25*fs));
for k=1:numel(s)
 ix=s(k):e(k);f(ix)=smoothdata(d.fwd(ix),'movmedian',n);r(ix)=smoothdata(d.yaw(ix),'movmedian',n);
 if numel(ix)>2,acc(ix)=gradient(f(ix),t(ix));end
end
candidate=valid&acc>=2;change=diff([false;candidate;false]);s=find(change==1);e=find(change==-1)-1;
idx=zeros(0,1);angle=idx;before=idx;after=idx;last=-inf;
for k=1:numel(s)
 [~,j]=max(acc(s(k):e(k)));i=s(k)+j-1;t0=t(i);
 if t0-last<2,continue,end
 lo=find(t<=t0-2,1,'last');hi=find(t>=t0+3,1,'first');
 if isempty(lo)||isempty(hi)||~all(valid(lo:hi))||any(diff(t(lo:hi))>2/fs),continue,end
 pre=t>=t0-.75&t<t0-.25;post=t>=t0+.25&t<=t0+.75;core=t>=t0-.5&t<=t0+.5;
 a=median(f(pre));b=median(f(post));if b<1||b-a<1,continue,end
 idx(end+1,1)=i;angle(end+1,1)=0;before(end+1,1)=a;after(end+1,1)=b;last=t0; %#ok<AGROW>
end
angle=pizza_ves041_turn_angle(d,idx,[-.5 0]);
info=struct('accelerationMinimumMmS2',2,'minimumSpeedGainMmS',1,'minimumSpeedAfterMmS',1, ...
 'smoothWindowSec',.25,'minimumAcceptedSpacingSec',2,'preWindowSec',[-.75 -.25], ...
 'postWindowSec',[.25 .75],'turnWindowSec',[-.5 0],'turnAbsoluteYawDeg',15,'eventWindowSec',[-2 3]);
end
