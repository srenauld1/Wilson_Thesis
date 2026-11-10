function build_ves041_matlab_examples(dataDir)
% Package portable data and create editable MATLAB figures.
if nargin<1
    dataDir=fullfile(pwd,'outputs','ves041_matlab_examples_20261104');
end
cfg=jsondecode(fileread(fullfile(dataDir,'source_data','config.json')));
cases=struct([]);
for k=1:numel(cfg)
    c=cfg(k);
    velocity=readmatrix(fullfile(dataDir,'source_data',c.velocity_csv));
    heading=readmatrix(fullfile(dataDir,'source_data',c.heading_csv));
    row=struct('name',c.name,'fly',c.fly,'genotype',c.genotype, ...
        'event_id',c.event_id,'pulse',c.pulse,'window',[c.start_s c.end_s], ...
        'time',velocity(:,1),'forward',velocity(:,2),'yaw',velocity(:,3), ...
        'heading_time',heading(:,1),'heading',heading(:,2), ...
        'light_intervals',[[c.light.start_s]' [c.light.end_s]'], ...
        'event_windows',zeros(0,4),'event_metrics',zeros(0,3),'original_source',c.source, ...
        'velocity_processing','Stored raw velocity fields, sampled at 1 kHz, no additional smoothing', ...
        'heading_processing','Original processed heading, relative to recording start', ...
        'event_processing','Original primary >30-degree detections, unchanged');
    for q=1:numel(c.events)
        e=c.events(q);
        row.event_windows(q,:)=[e.start_s e.end_s e.marker_s e.event_id];
        row.event_metrics(q,:)=[e.rotation_deg e.displacement_mm e.directionality];
    end
    if k==1, cases=row; else, cases(k)=row; end %#ok<AGROW>
end
save(fullfile(dataDir,'ves041_example_timeseries.mat'),'cases','-v7');
copyfile(fullfile(fileparts(mfilename('fullpath')),'render_ves041_matlab_examples.m'),dataDir);
render_ves041_matlab_examples(dataDir);
end
