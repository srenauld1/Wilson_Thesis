function results=subtract_optomotor_random_background_roi(outputFolder)
% Add a reproducible random background ROI without changing any source ROI.
% ROI2: match ROI1 voxel count separately in each z plane; sample uniformly
% from outside all saved ROIs, a two-pixel ROI buffer, and a two-pixel image
% edge. Selection never uses fluorescence or stimulus responses.
% Both re-extracted signals use the original dff005000 convention, including
% its minimum subtraction. Corrected signal is SAVED ROI1 dF/F minus ROI2
% dF/F, not subtraction of raw fluorescence and not event-baseline subtraction.
% The independent normalizations are not a fitted bleedthrough coefficient;
% random non-ROI pixels may include other neurons/neuropil. Treat as a control.
% Existing events, ROI1, acquisition files, and earlier figures are untouched.
if nargin<1
 outputFolder='/Users/sophiarenauld/stacks/pizza_talk_analysis/data_crunches/random_background_roi';
end
assert(~feature('ShowFigureWindows'),'Run MATLAB -nodisplay.');
if ~isfolder(outputFolder),mkdir(outputFolder);end
source='/Users/sophiarenauld/stacks/pizza_talk_analysis/ves041_noheat_imaging/gratings/optomotor_exploration/optomotor_summary.mat';
S=load(source,'R');results=struct([]);
for k=1:numel(S.R)
 r=S.R(k);fprintf('BACKGROUND ROI %s\n',r.trial);
 a=load(char(r.source),'fly_data');roi=a.fly_data.roi;
 assert(strcmp(roi.opt.nrm.nrmstr,'dff005000')&&roi.opt.nrm.degdtr==0&&isempty(roi.opt.nrm.wavp));
 mask1=logical(roi.mask);allRoi=logical(roi.mask_allroi);
 assert(isequal(size(mask1),size(allRoi))&&roi.numroi==1);
 folder=fileparts(fileparts(char(r.source)));files=dir(fullfile(folder,'*_or_.tif'));assert(numel(files)==1);
 tiffPath=fullfile(files.folder,files.name);nt=numel(r.t);[ny,nx,nz]=size(mask1);
 rng(20261104+k,'twister');mask2=false(size(mask1));counts=zeros(nz,3);
 for z=1:nz
  excluded=conv2(double(allRoi(:,:,z)),ones(5),'same')>0;
  excluded([1:2,end-1:end],:)=true;excluded(:,[1:2,end-1:end])=true;
  candidates=find(~excluded);n=nnz(mask1(:,:,z));assert(numel(candidates)>=n,'Insufficient background pixels');
  chosen=candidates(randperm(numel(candidates),n));plane=false(ny,nx);plane(chosen)=true;mask2(:,:,z)=plane;
  counts(z,:)=[n,numel(candidates),nnz(plane)];
 end
 assert(~any(mask1&mask2,'all'));
 % Stream pages in z-fast, time-slow order; verify against saved ROI1 before
 % accepting either fluorescence trace. This catches wrong channel/order.
 tf=Tiff(tiffPath,'r');closer=onCleanup(@()tf.close());
 assert(tf.getTag('ImageLength')==ny&&tf.getTag('ImageWidth')==nx);
 raw=zeros(2,nt);meanImage=zeros(ny,nx,nz);n1=nnz(mask1);n2=nnz(mask2);
 for it=1:nt
  for z=1:nz
   im=double(tf.read());assert(isequal(size(im),[ny nx]));
   raw(1,it)=raw(1,it)+sum(im(mask1(:,:,z)))/n1;
   raw(2,it)=raw(2,it)+sum(im(mask2(:,:,z)))/n2;
   meanImage(:,:,z)=meanImage(:,:,z)+im/nt;
   if it<nt||z<nz
    assert(~tf.lastDirectory(),'TIFF shorter than volume clock');tf.nextDirectory();
   end
  end
  if mod(it,1000)==0,fprintf('  extracted %d/%d volumes\n',it,nt);end
 end
 assert(tf.lastDirectory(),'TIFF has unexpected extra frames');clear closer;
 shifted=raw-min(raw,[],2);f0=prctile(shifted,5,2);assert(all(f0>0));
 normalized=(shifted-f0)./f0;
 stored=double(roi.ts(1,:));assert(numel(stored)==nt);
 maxReconstructionError=max(abs(normalized(1,:)-stored));
 assert(maxReconstructionError<1e-3,'ROI1 extraction does not reproduce saved dF/F; do not subtract mismatched data.');
 assert(max(abs(double(r.d.dff(:))-stored(:)))<1e-6,'Saved analysis ROI1 differs from source ROI1');
 q=struct('trial',char(r.trial),'source',char(r.source),'tiff',tiffPath,'t',double(r.t(:)), ...
  'roi1_original',stored(:),'roi2_dff',normalized(2,:)','corrected_dff',stored(:)-normalized(2,:)', ...
  'roi1_mask',mask1,'roi2_mask',mask2,'raw_fluorescence',raw,'f0_after_min_subtraction',f0, ...
  'mean_image',meanImage,'pixel_counts',counts,'random_seed',20261104+k, ...
  'roi1_reconstruction_max_error',maxReconstructionError,'onsetSec',r.E.onsetSec, ...
  'turnSec',r.E.matchedTurnPeakSec,'onsetValid',r.E.onsetWindowValid, ...
  'turnValid',r.E.turnWindowValid,'valid',r.valid(:));
 q.methods=struct('normalization','dff005000: subtract each trace minimum, static 5th percentile F0, (F-F0)/F0', ...
  'subtraction','saved original ROI1 dF/F minus independently normalized ROI2 dF/F; no event-baseline subtraction', ...
  'selection','Uniform random pixels outside all ROIs plus 2-pixel buffer and 2-pixel image edge; match ROI1 pixel count in each plane', ...
  'warning','ROI2 may contain real neuropil activity. No fitted contamination coefficient; independent normalization may overcorrect.', ...
  'originals','ROI1 and all source files remain unchanged');
 if k==1,results=q;else,results(k)=q;end %#ok<AGROW>
 save(fullfile(outputFolder,[q.trial '_random_background_roi.mat']),'q','-v7');
 fprintf('  ROI1 reproduction max error %.6g; ROI2 pixels %d; ROI2 F0 %.4f\n',maxReconstructionError,n2,f0(2));
end
save(fullfile(outputFolder,'three_fly_random_background_roi.mat'),'results','-v7');
end
