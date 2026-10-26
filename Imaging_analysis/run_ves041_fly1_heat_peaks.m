function run_ves041_fly1_heat_peaks(root)
% Use the exact signals in the existing fwd_rot_vs_dff_separate figures.
if nargin<1,root='/Users/sophiarenauld/stacks/2026-10-20';end
files=[dir(fullfile(root,'*fly1*heatbarrier*','*_s_.mat'));dir(fullfile(root,'*fly1*stripes*','*_s_.mat'))];
for k=1:numel(files)
    a=load(fullfile(files(k).folder,files(k).name),'dq','roi');
    out=fullfile(files(k).folder,'smr_analysis');source=fullfile(out,'fwd_rot_vs_dff_separate.fig');
    fig=openfig(source,'invisible');axs=findall(fig,'Type','axes');
    assert(numel(axs)==3,'Expected three source axes');pos=cell2mat(get(axs,'Position'));[~,order]=sort(pos(:,2),'descend');
    vars=cell(3,1);
    for v=1:3
        lines=findall(axs(order(v)),'Type','line');[~,j]=max(arrayfun(@(z)numel(z.XData),lines));
        tt=lines(j).XData(:);vars{v}=lines(j).YData(:);
        assert(numel(tt)==numel(a.dq(1).t)&&max(abs(tt-a.dq(1).t(:)))<1e-6,'Source trace time differs from DAQ time');
    end
    close(fig);a.roi.ts(1,:)=vars{1}';a.dq(1).bvf=vars{2};a.dq(1).bvy_deg=vars{3};a.peakTraceSource=source;
    a=import_ves041_protocol(a,files(k).folder);
    plot_ves041_peak_heat_context(a,out);
end
end
