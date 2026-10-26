function S=inspect_ves041_trial11_pixel_noise()
% Dense raw-pixel check around the abrupt step. Quantifies fixed dark-region
% intensity distributions independently of the calcium ROI. No data modified.
set(groot,'defaultFigureVisible','off');
out=fullfile(fileparts(mfilename('fullpath')),'diagnostics','ves041_20261020_11');
q=load(fullfile(out,'diagnostic.mat'),'result');d=q.result;
folder='/Users/sophiarenauld/stacks/2026-10-20/20261020_11_2_ves041_syt8m_fly1_heatbarrier_redo';
a=load(fullfile(folder,'20261020_11_2_or_z0z1_s_.mat'),'dq','roi');t=a.dq(1).t(:);
idx=find(t>=390&t<400);pre=t(idx)<394;post=t(idx)>=396;ts=t(idx);
h=Tiff(fullfile(folder,'20261020_11_2_o_.tif'),'r');c=onCleanup(@()h.close());
stats=nan(numel(idx),4,2,7);means=zeros([size(d.mask) 2],'single');counts=[sum(pre),sum(post)];
for k=1:numel(idx)
    for z=1:4
        h.setDirectory((idx(k)-1)*4+z);im=single(h.read());
        masks={d.mask(:,:,z),d.dark(:,:,z)};
        for r=1:2
            v=im(masks{r});p=prctile(v,[10 50 90 99]);
            stats(k,z,r,:)=[mean(v) std(v) p mean(v>0)];
        end
        if pre(k),means(:,:,z,1)=means(:,:,z,1)+im;elseif post(k),means(:,:,z,2)=means(:,:,z,2)+im;end
    end
end
for ep=1:2,means(:,:,:,ep)=means(:,:,:,ep)/counts(ep);end
names={'mean','SD','p10','median','p90','p99','positive fraction'};
for z=1:4
    fprintf('\nz%d dark raw pixels\n',z);
    for j=1:7,x=stats(:,z,2,j);fprintf('%s pre %.4f post %.4f\n',names{j},mean(x(pre)),mean(x(post)));end
end
fig=figure('Position',[50 50 1300 1000]);tl=tiledlayout(4,2,'TileSpacing','compact');
for j=1:7,nexttile;plot(ts,squeeze(stats(:,:,2,j)));xline(395,':');title(['Dark pixels: ' names{j}]);xlabel('s');end
nexttile;plot(t,a.roi.ts);xlim([390 400]);title('Saved dF/F');xlabel('s');
title(tl,'Unregistered movie: dark-background distributions around baseline step');exportgraphics(fig,fullfile(out,'raw_noise_distributions.png'),'Resolution',130);savefig(fig,fullfile(out,'raw_noise_distributions.fig'));close(fig);
fig=figure('Position',[50 50 1300 900]);tl=tiledlayout(4,3,'TileSpacing','compact');
for z=1:4
    ims=imgaussfilt(means(:,:,z,1),1);imp=imgaussfilt(means(:,:,z,2),1);lims=prctile([ims(:);imp(:)],[2 99.5]);
    for j=1:3,ax=nexttile;if j==1,im=ims;elseif j==2,im=imp;else,im=imp-ims;end
        imagesc(im);axis image off;hold on;contour(d.mask(:,:,z),[.5 .5],'g');colorbar;
        if j<3,clim(lims);colormap(ax,gray);else,v=prctile(abs(im(:)),99);clim([-v v]);colormap(ax,turbo);end
        title(sprintf('z%d %s',z,string(j)));
    end
end
title(tl,'Raw movie means (smoothed 1 px): 390-394 s, 396-400 s, difference');exportgraphics(fig,fullfile(out,'raw_spatial_step.png'),'Resolution',130);close(fig);
S=struct('t',ts,'stats',stats,'names',{names},'means',means);save(fullfile(out,'pixel_noise.mat'),'S');
end
