function plot_ves041_stack_maps(root,binFactor)
% Plot exact per-pixel mean fluorescence and maximum dF/F for every z slice.
% Run compute_ves041_stack_maps.py first. Numerical maps retain full values;
% only the display saturates at a shared recording-level 99.5th percentile.
if nargin<1,root='/Users/sophiarenauld/stacks/2026-10-20';end
if nargin<2,binFactor=1;end
mapFolder='pixel_fluorescence_maps';if binFactor==4,mapFolder='pixel_fluorescence_maps_binned4x4';end
folders=dir(fullfile(root,'*fly1*'));folders=folders([folders.isdir]);
for k=1:numel(folders)
    out=fullfile(folders(k).folder,folders(k).name,'smr_analysis',mapFolder);
    jsonfile=fullfile(out,'pixel_fluorescence_map_settings.json');
    if ~isfile(jsonfile),continue,end
    settings=jsondecode(fileread(jsonfile));file=fullfile(out,'pixel_fluorescence_maps.h5');
    % h5read reverses HDF5 dimensions: disk [y,x,z] -> MATLAB [z,x,y].
    meanFluorescence=permute(h5read(file,'/mean_fluorescence'),[3 2 1]);
    maxDff=permute(h5read(file,'/max_dff'),[3 2 1]);
    f0=permute(h5read(file,'/f0'),[3 2 1]);
    maxFluorescence=permute(h5read(file,'/max_fluorescence'),[3 2 1]);
    nz=settings.slices;
    assert(isequal(size(meanFluorescence),[settings.source_shape_tzxy(4)/binFactor,settings.source_shape_tzxy(3)/binFactor,nz]),'Unexpected map orientation');
    % Independently verify source orientation and exact statistics at pixels
    % spread across all slices, including registration padding.
    source=matfile(settings.source);
    for xy=[1 1;round(size(f0,2)/2) round(size(f0,1)/2);size(f0,2)-1 size(f0,1)-1]'
        x=xy(1);y=xy(2);
        iy=(y-1)*binFactor+(1:binFactor);ix=(x-1)*binFactor+(1:binFactor);
        raw=reshape(double(source.stack(iy,ix,:,:)),binFactor^2,nz,[]);
        samples=reshape(mean(raw,1),nz,[]);
        for z=1:nz
            v=samples(z,:);
            assert(abs(mean(v)-meanFluorescence(y,x,z))<1e-8,'Mean differs from source');
            assert(abs(prctile(v,5)-f0(y,x,z))<1e-8,'Baseline differs from fifth percentile');
            assert(max(v)==maxFluorescence(y,x,z),'Maximum differs from source');
            if f0(y,x,z)>0,assert(abs((max(v)-f0(y,x,z))/f0(y,x,z)-maxDff(y,x,z))<1e-8,'dF/F differs from source');end
        end
    end
    trial=folders(k).name;binLabel='';if binFactor>1,binLabel=' | 4x4 block-mean fluorescence';end
    data={meanFluorescence,maxDff};stems={'mean_fluorescence_slices','max_dff_slices'};
    labels={'Mean fluorescence (saved intensity units)','Maximum dF/F'};
    limits={settings.mean_display_limits,settings.max_dff_display_limits};
    for kind=1:2
        lim=limits{kind}(:)';if diff(lim)<=0,lim=[0 1];end
        fig=figure('Visible','off','Color','w','Position',[40 40 1500 950]);
        tl=tiledlayout(fig,ceil(nz/2),2,'TileSpacing','compact','Padding','loose');
        for z=1:nz
            ax=nexttile(tl);drawmap(ax,data{kind}(:,:,z),lim,kind,labels{kind},binFactor);title(ax,sprintf('Slice %d',z));
        end
        if kind==1,detail=sprintf('Mean over all %d volumes',settings.volumes);
        else,detail='Maximum over time of (F-F0)/F0 | F0: temporal 5th percentile of each displayed pixel/block | F0 <=0: black';end
        displayNote='Shared color scale across slices; values above the 99.5th percentile saturate for display';
        if kind==2 && ~any(isfinite(maxDff(:))),displayNote='No finite maximum dF/F: the fifth-percentile F0 is zero at every pixel';end
        title(tl,{[strrep(trial,'_',' ') binLabel],detail,displayNote});
        savefig(fig,fullfile(out,[stems{kind} '.fig']));
        exportgraphics(fig,fullfile(out,[stems{kind} '.png']),'Resolution',160);close(fig);
        for z=1:nz
            fig=figure('Visible','off','Color','w','Position',[40 40 1100 650]);ax=axes(fig);
            drawmap(ax,data{kind}(:,:,z),lim,kind,labels{kind},binFactor);
            title(ax,sprintf('%s | slice %d | %s',strrep(trial,'_',' '),z,[labels{kind} binLabel]));
            exportgraphics(fig,fullfile(out,sprintf('%s_z%02d.png',stems{kind},z)),'Resolution',160);close(fig);
        end
    end
    save(fullfile(out,'pixel_fluorescence_maps.mat'),'meanFluorescence','maxDff','f0','maxFluorescence','settings');
    fprintf('VERIFIED AND PLOTTED %s: %d slices, %d volumes\n',trial,nz,settings.volumes);
end
end
function drawmap(ax,v,lim,kind,label,binFactor)
cx=((1:size(v,2))-.5)*binFactor+.5;cy=((1:size(v,1))-.5)*binFactor+.5;
imagesc(ax,cx,cy,v,'AlphaData',isfinite(v));set(ax,'Color','k');axis(ax,'image');
if ~any(isfinite(v(:)))
    text(ax,.5,.5,'Undefined: F0 = 0 at every pixel','Units','normalized', ...
        'Color','w','HorizontalAlignment','center','FontSize',13);
end
xlabel(ax,'X (pixels)');ylabel(ax,'Y (pixels)');clim(ax,lim);
if kind==1,colormap(ax,turbo);else,colormap(ax,hot);end
if any(isfinite(v(:))),cb=colorbar(ax);cb.Label.String=label;end
end
