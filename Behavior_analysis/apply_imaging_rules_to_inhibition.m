function apply_imaging_rules_to_inhibition()
% Native-rate reproduction of the pizza-talk imaging behavioral event rules.
addpath(fullfile(pwd,'Imaging_analysis'));
out=fullfile(pwd,'outputs','ves041_inhibition_imaging_rules_20261104');
if ~isfolder(out),mkdir(out);end
cache=fullfile(out,'recordings');if ~isfolder(cache),mkdir(cache);end
records=jsondecode(fileread(fullfile(pwd,'outputs','ves041_blinded_pivots_gt30_20261102','manifest.json')));
assert(numel(records)==28);
% Minimal positive, negative and validity tests use the identical detector.
z.t=(0:.01:20)';z.fwd=6-4./(1+exp(15*(10-z.t)));z.yaw=100*exp(-((z.t-10.25)/.18).^2);z.valid=true(size(z.t));
[ix,ang]=detect(z);assert(numel(ix)==1&&ang>=15);
z.yaw(:)=0;[ix,ang]=detect(z);assert(numel(ix)==1&&ang==0);
z.valid(z.t>9&z.t<11)=false;ix=detect(z);assert(isempty(ix));
status=struct('recording_key',{},'status',{},'message',{});
for k=1:numel(records)
 rec=records(k);dest=fullfile(cache,[rec.recording_key '.mat']);
 try
  if isfile(dest)
   fprintf('CACHED %d/%d %s\n',k,numel(records),rec.recording_key);
  else
   file=rec.source;assert(isfile(file),'Source missing: %s',file);
   fields=h5info(file,'/exptData');names={fields.Datasets.Name};
   d.fwd=vec(file,'forwardVelocity');d.yaw=vec(file,'angularVelocity');speed=vec(file,'totSpeed');n=numel(d.fwd);
   assert(numel(d.yaw)==n&&numel(speed)==n);
   if ismember('t_ds',names),d.t=vec(file,'t_ds');else,d.t=vec(file,'tDS');end
   assert(numel(d.t)==n);dt=median(diff(d.t));assert(dt>0&&all(diff(d.t)>0));
   if ismember('optoStim_ds',names),light=vec(file,'optoStim_ds')>0;
   else,sz=h5info(file,'/exptData/optoStim');ratio=prod(sz.Dataspace.Size)/n;assert(abs(ratio-round(ratio))<1e-10);light=vec(file,'optoStim',round(ratio))>0;end
   assert(numel(light)==n);
   d.valid=isfinite(d.t)&isfinite(d.fwd)&isfinite(d.yaw)&isfinite(speed);
   [idx,angle,details,fsm,ysm]=detect(d);
   % Direct comparison against the existing imaging rotation integration.
   take=1:min(5,numel(idx));reference=pizza_ves041_turn_angle(d,idx(take),[0 .5]);
   assert(max(abs(reference-angle(take)),[],'omitnan')<1e-6||isempty(take));
   edges=diff([false;light;false]);starts=find(edges==1);ends=find(edges==-1);assert(numel(starts)==40);
   phases=struct('pulse',{},'phase',{},'definition',{},'events',{},'seconds',{},'moving_onset',{},'moving_fraction',{},'mean_forward',{},'onset_1s_pivot',{},'all_slowdowns',{});
   for p=1:numel(starts)
    a=starts(p);b=ends(p);len=b-a;pre=a-len;
    assert(pre>=1&&~any(light(pre:a-1)));
    assert(all(d.valid(pre:b-1))&&all(diff(d.t(pre:b-1))<=2*dt),'Invalid phase exposure requires review.');
    moving=mean(speed(max(1,a-round(1/dt)):a-1)>=3)>=.5;
    for definition={'primary','within_phase_support'}
     pad=0;if strcmp(definition{1},'within_phase_support'),pad=.875+dt;end
     for phase={'pre','inhibition'}
      if strcmp(phase{1},'pre'),lo=d.t(pre);hi=d.t(a);else,lo=d.t(a);if b<=n,hi=d.t(b);else,hi=d.t(end)+dt;end,end
      lo=lo+pad;hi=hi-pad;assert(hi>lo);
      selected=d.t(idx)>=lo&d.t(idx)<hi;turn=selected&angle>=15;
      data=d.t>=lo&d.t<hi;
      row=struct('pulse',p,'phase',phase{1},'definition',definition{1},'events',nnz(turn),'seconds',hi-lo,'moving_onset',moving, ...
       'moving_fraction',mean(speed(data)>=3),'mean_forward',mean(d.fwd(data)),'onset_1s_pivot',any(d.t(idx(turn))<lo+1),'all_slowdowns',nnz(selected));
      phases(end+1)=row; %#ok<AGROW>
     end
    end
   end
   events=struct('sample_index',num2cell(idx),'time_s',num2cell(d.t(idx)), ...
    'rotation_deg',num2cell(angle),'is_reorientation',num2cell(angle>=15), ...
    'pre_forward_median',num2cell(details(:,1)),'post_forward_median',num2cell(details(:,2)), ...
    'forward_drop',num2cell(details(:,3)),'acceleration',num2cell(details(:,4)));
   % Save small example traces, retaining unmodified native detection values.
   examples=struct('relative_time',{},'forward',{},'smoothed_forward',{},'yaw',{},'rotation',{},'smoothed_yaw',{},'detection',{},'phase',{});
   for phase={'pre','inhibition'}
    times=d.t(idx);eligible=false(size(idx));
    for p=1:numel(starts)
     a=starts(p);b=ends(p);len=b-a;
     if strcmp(phase{1},'pre'),lo=d.t(a-len);hi=d.t(a);else,lo=d.t(a);hi=d.t(min(b,n));end
     eligible=eligible|(times>=lo&times<hi&angle>=15);
    end
    take=find(eligible,2,'first');
    for j=take(:)'
     ii=idx(j);q=lower(d.t,d.t(ii)-2):max(1,upper(d.t,d.t(ii)+3)-1);q=q(1:max(1,round(.02/dt)):end);
     examples(end+1)=struct('relative_time',d.t(q)-d.t(ii),'forward',d.fwd(q),'smoothed_forward',fsm(q), ...
       'yaw',d.yaw(q),'rotation',abs(d.yaw(q)),'smoothed_yaw',ysm(q),'detection',events(j),'phase',phase{1}); %#ok<AGROW>
    end
   end
   metadata=struct('source',file,'sample_rate_hz',1/dt,'total_slowdowns',numel(idx),'total_reorientations',nnz(angle>=15), ...
    'finite_behavior_fraction',mean(d.valid),'native_samples',n,'direct_angle_validation_max_error',max([0;abs(reference-angle(1:numel(reference)))]));
   save(dest,'rec','phases','events','examples','metadata','-v7.3');
   payload=struct('record',rec,'presentations',phases,'events',events,'examples',examples,'metadata',metadata);
   fid=fopen(fullfile(cache,[rec.recording_key '.json']),'w');assert(fid>=0);fprintf(fid,'%s',jsonencode(payload));fclose(fid);
   fprintf('EXTRACTED %d/%d %s: %d slowdowns, %d with rotation\n',k,numel(records),rec.recording_key,numel(idx),nnz(angle>=15));
  end
  item=struct('recording_key',rec.recording_key,'status','complete','message','');
 catch ME
  item=struct('recording_key',rec.recording_key,'status','failed','message',getReport(ME,'extended','hyperlinks','off'));
  fprintf(2,'ERROR %s: %s\n',rec.recording_key,ME.message);
 end
 status(end+1)=item; %#ok<AGROW>
 fid=fopen(fullfile(out,'extraction_status.json'),'w');fprintf(fid,'%s',jsonencode(status,PrettyPrint=true));fclose(fid);
end
assert(all(strcmp({status.status},'complete')),'Some recordings failed; review extraction_status.json.');
fprintf('COMPLETE: %d native-rate recordings\n',numel(records));
end

function v=vec(file,field,stride)
if nargin<3,stride=1;end
path=['/exptData/' field];info=h5info(file,path);sz=info.Dataspace.Size;
if stride==1,v=double(h5read(file,path));else
 start=ones(size(sz));step=start;count=sz;[~,axis]=max(sz);step(axis)=stride;count(axis)=ceil(count(axis)/stride);
 v=double(h5read(file,path,start,count,step));end
v=v(:);
end

function [idx,angle,details,f,y]=detect(d)
t=d.t;fs=1/median(diff(t));n=max(3,round(.25*fs));valid=d.valid;
cut=[true;diff(t)>2/fs];a=find(valid&(cut|[true;~valid(1:end-1)]));b=find(valid&([~valid(2:end);true]|[cut(2:end);true]));
f=nan(size(t));y=f;acc=f;integral=f;
for k=1:numel(a)
 ix=a(k):b(k);f(ix)=smoothdata(d.fwd(ix),'movmedian',n);y(ix)=smoothdata(d.yaw(ix),'movmedian',n);
 if numel(ix)>2,acc(ix)=gradient(f(ix),t(ix));end
 integral(ix)=cumtrapz(t(ix),abs(y(ix)));
end
candidate=valid&acc<=-2;change=diff([false;candidate;false]);s=find(change==1);e=find(change==-1)-1;
bad=[0;cumsum(~valid)];gaps=[0;cumsum([false;diff(t)>2/fs])];
idx=zeros(0,1);angle=zeros(0,1);details=zeros(0,4);last=-inf;
for j=1:numel(s)
 [~,loc]=min(acc(s(j):e(j)));i=s(j)+loc-1;t0=t(i);
 if t0-last<2||t0-2<t(1)||t0+3>t(end),continue,end
 lo=lower(t,t0-2);hi=upper(t,t0+3)-1;q=upper(t,t0-2)-1;r=lower(t,t0+3);
 if q<1||r>numel(t)||bad(r+1)-bad(q)>0||gaps(hi+1)-gaps(lo+1)>0,continue,end
 pre=lower(t,t0-.75):(lower(t,t0-.25)-1);post=lower(t,t0+.25):(upper(t,t0+.75)-1);
 before=median(f(pre));after=median(f(post));drop=before-after;
 if before<1||drop<1,continue,end
 target=t0+.5;k=upper(t,target)-1;
 endpoint=abs(y(k));delta=target-t(k);
 if delta>0,endpoint=endpoint+(abs(y(k+1))-endpoint)*delta/(t(k+1)-t(k));end
 area=integral(k)+delta*(abs(y(k))+endpoint)/2-integral(i);
 idx(end+1,1)=i;angle(end+1,1)=area;details(end+1,:)=[before after drop acc(i)];last=t0; %#ok<AGROW>
end
end

function i=lower(t,x)
% First index with t(i)>=x, matching MATLAB logical-window boundaries.
lo=1;hi=numel(t)+1;while lo<hi,m=floor((lo+hi)/2);if t(m)<x,lo=m+1;else,hi=m;end,end;i=lo;
end
function i=upper(t,x)
% First index with t(i)>x.
lo=1;hi=numel(t)+1;while lo<hi,m=floor((lo+hi)/2);if t(m)<=x,lo=m+1;else,hi=m;end,end;i=lo;
end
