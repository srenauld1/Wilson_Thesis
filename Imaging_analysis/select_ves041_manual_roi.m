function a=select_ves041_manual_roi(a)
% Use full hand-drawn ROI traces, retaining the segmented alternatives.
% SCOPA stores both traces with identical normalization in pxfilter.
% Legacy files without pixel filtering already store hand-drawn roi.ts.
for k=1:numel(a.roi)
    r=a.roi(k);
    if isfield(r,'pxfilter') && ~isempty(r.pxfilter)
        assert(isfield(r.pxfilter,'manualts') && ~isempty(r.pxfilter.manualts), ...
            'ves041:MissingManualTrace','Pixel-filtered input lacks a saved hand-drawn trace.');
        assert(isequal(size(r.ts),size(r.pxfilter.manualts)), ...
            'ves041:ManualTraceShape','Manual and current traces have different dimensions.');
        a.roi(k).ts=r.pxfilter.manualts;
        if isfield(r.pxfilter,'manualmask'),a.roi(k).mask=r.pxfilter.manualmask;end
    end
end
a.analysisRoiSource='Full hand-drawn ROI (pxfilter.manualts when available)';
end
