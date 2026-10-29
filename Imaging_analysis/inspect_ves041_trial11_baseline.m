function result=inspect_ves041_trial11_baseline()
% Read-only movie/ROI diagnostic. Original TIFFs, ROIs and analysis are unchanged.
% Sample both raw and registered four-plane volumes at 2 Hz. Compare fixed
% ROI, nearby annulus, bright non-ROI tissue, and dark non-ROI background.
% A calcium-only change should not be inferred from brightness alone: examine
% morphology, plane redistribution, motion and background alongside the trace.
set(groot,'defaultFigureVisible','off');
folder='/Users/sophiarenauld/stacks/2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo';
out=fullfile(imaging_analysis_artifact_root(),'diagnostics','ves041_20261020_11');
if ~isfolder(out),mkdir(out);end
a=load(fullfile(folder,'20261020_11_2_or_z0z1_s_.mat'),'roi','dq','t');
t=a.dq(1).t(:); mask=a.roi.mask; nt=numel(t);nz=size(mask,3);
idx=unique(arrayfun(@(s)find(abs(t-s)==min(abs(t-s)),1),0:.5:t(end)));
ts=t(idx); n=numel(idx); pre=ts>=375&ts<385;post=ts>=405&ts<415;
raw=Tiff(fullfile(folder,'20261020_11_2_o_.tif'),'r');reg=Tiff(fullfile(folder,'20261020_11_2_or_.tif'),'r');
cleanup=onCleanup(@()closefiles(raw,reg));
meanimg=single(a.roi.mnt); ann=false(size(mask));dark=ann;other=ann;
for z=1:nz
    m=mask(:,:,z);im=meanimg(:,:,z);safe=true(size(m));safe([1:12 end-11:end],:)=false;safe(:,[1:12 end-11:end])=false;
    ann(:,:,z)=imdilate(m,strel('disk',12,0))&~imdilate(m,strel('disk',3,0))&safe;
    outside=~imdilate(m,strel('disk',15,0))&safe;
    v=im(outside);dark(:,:,z)=outside&im<=prctile(v,25);other(:,:,z)=outside&im>=prctile(v,75);
end
% Dimensions: sample, plane, region, movie. Regions share fixed pixel masks.
F=nan(n,nz,4,2);centroid=nan(n,nz,2,2);shapeCorr=nan(n,nz,2);sharp=nan(n,nz,2);
means=zeros([size(mask) 2 2],'single');counts=zeros(2,1);snapTimes=[380 390 395 400 405 415];
snapIdx=arrayfun(@(s)find(abs(ts-s)==min(abs(ts-s)),1),snapTimes);snaps=zeros([size(mask) numel(snapIdx) 2],'single');
[xx,yy]=meshgrid(1:size(mask,2),1:size(mask,1));
for k=1:n
    for z=1:nz
        masks={mask(:,:,z),ann(:,:,z),other(:,:,z),dark(:,:,z)};
        for movie=1:2
            if movie==1,h=raw;else,h=reg;end
            h.setDirectory((idx(k)-1)*nz+z);im=single(h.read());
            for q=1:4,F(k,z,q,movie)=mean(im(masks{q}));end
            w=max(im-F(k,z,4,movie),0).*single(masks{1});
            centroid(k,z,:,movie)=[sum(w.*xx,'all') sum(w.*yy,'all')]/sum(w,'all');
            v=im(masks{1});ref=meanimg(:,:,z);shapeCorr(k,z,movie)=corr(double(v),double(ref(masks{1})));
            sm=imgaussfilt(im,1);[gx,gy]=gradient(sm);sharp(k,z,movie)=mean(gx(masks{1}).^2+gy(masks{1}).^2)/max(mean(v)^2,eps);
            for ep=1:2
                if (ep==1&&pre(k))||(ep==2&&post(k)),means(:,:,z,ep,movie)=means(:,:,z,ep,movie)+im;end
            end
            ss=find(snapIdx==k);if ~isempty(ss),snaps(:,:,z,ss,movie)=im;end
        end
    end
    if pre(k),counts(1)=counts(1)+1;end;if post(k),counts(2)=counts(2)+1;end
    if mod(k,200)==0,fprintf('Read %d/%d sampled raw and registered volumes\n',k,n);end
end
for ep=1:2,means(:,:,:,ep,:)=means(:,:,:,ep,:)/counts(ep);end
shifts=zeros(nt,2,nz);for z=1:nz,r=load(fullfile(folder,sprintf('20261020_11_2_or_reginfo_z%02d.mat',z)),'shifts');shifts(:,:,z)=r.shifts;end
names={'ROI','nearby annulus','bright outside ROI','dark background'};
for movie=1:2
    fprintf('\nMOVIE %d (1=unregistered,2=registered)\n',movie);
    for z=1:nz
        fprintf('plane %d ROI pixels %d\n',z,nnz(mask(:,:,z)));
        for q=1:4,x=F(:,z,q,movie);fprintf('%s pre %.3f post %.3f ratio %.3f\n',names{q},mean(x(pre)),mean(x(post)),mean(x(post))/mean(x(pre)));end
        fprintf('centroid change dx dy %s; shape corr pre/post %.3f %.3f\n',mat2str(squeeze(mean(centroid(post,z,:,movie),1)-mean(centroid(pre,z,:,movie),1))',3),mean(shapeCorr(pre,z,movie)),mean(shapeCorr(post,z,movie)));
    end
end
epochNames=["375-385 s","405-415 s","post minus pre"];movieNames=["raw","registered"];
fig=figure('Position',[50 50 1500 1000]);tl=tiledlayout(4,4,'TileSpacing','compact');
for z=1:nz
    before=means(:,:,z,1,2);after=means(:,:,z,2,2);lims=prctile([before(:);after(:)],[1 99.7]);
    for j=1:3
        ax=nexttile; if j==1,im=before;elseif j==2,im=after;else,im=after-before;end
        imagesc(ax,im);axis image off;hold on;contour(mask(:,:,z),[.5 .5],'g','LineWidth',.6);
        if j<3,clim(lims);colormap(ax,gray);else,v=max(abs(prctile(im(:),[1 99])));clim([-v v]);colormap(ax,turbo);end
        colorbar;title(sprintf('z%d %s',z,epochNames(j)));
    end
    nexttile;plot(ts,squeeze(F(:,z,1,2)),'b');hold on;plot(ts,squeeze(F(:,z,2,2)),'k');plot(ts,squeeze(F(:,z,3,2)),'Color',[.8 .4 0]);xlim([350 450]);xline(395,':');title(sprintf('z%d raw intensity',z));if z==1,legend(names(1:3));end
end
title(tl,'Registered movie: fixed ROI (green), common pre/post intensity scale');
exportgraphics(fig,fullfile(out,'spatial_pre_post.png'),'Resolution',130);savefig(fig,fullfile(out,'spatial_pre_post.fig'));close(fig);
fig=figure('Position',[50 50 1350 1000]);tl=tiledlayout(5,1,'TileSpacing','compact');
nexttile;plot(t,a.roi.ts,'b');xlim([350 450]);ylabel('saved dF/F');
nexttile;plot(ts,squeeze(mean(F(:,:,1,1),2)),'Color',[.6 .6 .6]);hold on;plot(ts,squeeze(mean(F(:,:,1,2),2)),'b');xlim([350 450]);ylabel('ROI intensity');legend('unregistered','registered');
nexttile;plot(ts,squeeze(mean(F(:,:,2:4,2),2)));xlim([350 450]);ylabel('control intensity');legend(names(2:4));
nexttile;plot(t,reshape(shifts,nt,[]));xlim([350 450]);ylabel('applied shift (px)');
nexttile;plot(t,a.dq(1).bvf);hold on;plot(t,a.dq(1).heat);xlim([350 450]);ylabel('behavior / command');xlabel('Time (s)');legend('forward mm/s','heat command V');
title(tl,'Trial 11 baseline shift diagnostic');exportgraphics(fig,fullfile(out,'temporal_diagnostic.png'),'Resolution',140);savefig(fig,fullfile(out,'temporal_diagnostic.fig'));close(fig);
fig=figure('Position',[50 50 1600 700]);tl=tiledlayout(2,6,'TileSpacing','compact');
for movie=1:2
    ims=squeeze(max(snaps(:,:,:,:,movie),[],3));lims=prctile(ims(:),[1 99.8]);
    for j=1:6,ax=nexttile;imagesc(ims(:,:,j),lims);colormap(ax,gray);axis image off;hold on;contour(any(mask,3),[.5 .5],'g');title(sprintf('%s %.1f s',movieNames(movie),ts(snapIdx(j))));end
end
title(tl,'Single-volume z projections; fixed intensity scale within each row');exportgraphics(fig,fullfile(out,'movie_contact_sheet.png'),'Resolution',140);close(fig);
result=struct('t',ts,'indices',idx,'F',F,'centroid',centroid,'shapeCorr',shapeCorr,'sharpness',sharp,'means',means,'shifts',shifts,'mask',mask,'annulus',ann,'other',other,'dark',dark);
save(fullfile(out,'diagnostic.mat'),'result','-v7.3');fprintf('Saved diagnostics: %s\n',out);
end
function closefiles(a,b)
a.close();b.close();
end
