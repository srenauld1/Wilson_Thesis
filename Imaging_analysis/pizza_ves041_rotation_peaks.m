function [idx,drop,before,after,info]=pizza_ves041_rotation_peaks(d)
t=d.t;fs=1/median(diff(t));valid=d.valid&isfinite(d.dff)&isfinite(d.fwd)&isfinite(d.yaw);cut=[true;diff(t)>2/fs];
s=find(valid&(cut|[true;~valid(1:end-1)]));e=find(valid&([~valid(2:end);true]|[cut(2:end);true]));
r=nan(size(t));f=r;for k=1:numel(s),ix=s(k):e(k);r(ix)=abs(smoothdata(d.yaw(ix),'gaussian',max(3,round(.12*fs))));f(ix)=smoothdata(d.fwd(ix),'movmedian',max(3,round(.25*fs)));end
baseline=movmedian(r,max(3,round(5*fs)),'omitmissing');res=r-baseline;res=res(isfinite(res));prom=max(20,2*1.4826*mad(res,1));
info=struct('prominence',prom,'minimumHeight',20,'minimumDistanceSec',1.5);candidates=zeros(0,1);
for k=1:numel(s)
 ix=(s(k):e(k))';if numel(ix)<3,continue,end
 [~,loc]=findpeaks(r(ix),'MinPeakHeight',20,'MinPeakProminence',prom,'MinPeakDistance',max(1,round(1.5*fs)));candidates=[candidates;ix(loc)]; %#ok<AGROW>
end
idx=zeros(0,1);drop=idx;before=idx;after=idx;
for i=candidates'
 lo=find(t<=t(i)-2,1,'last');hi=find(t>=t(i)+3,1,'first');
 if isempty(lo)||isempty(hi)||~all(valid(lo:hi))||any(diff(t(lo:hi))>2/fs),continue,end
 pre=t>=t(i)-.75&t<t(i)-.25;post=t>=t(i)+.25&t<=t(i)+.75;
 a=median(f(pre));b=median(f(post));idx(end+1,1)=i;before(end+1,1)=a;after(end+1,1)=b;drop(end+1,1)=a-b; %#ok<AGROW>
end
end
