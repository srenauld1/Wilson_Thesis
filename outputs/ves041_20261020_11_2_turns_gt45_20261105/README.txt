Recording 20261020_11_2: strict >45 degree accumulated absolute yaw.

Two-page paper PDF, same axes and -2:+3 second alignment to strongest forward deceleration.
Page 1 preserves linked PDF rules: angle over 0:+0.5 s; heat off 0 events, heat on 8. Original >=15 degree counts were 20/19.
Page 2 preserves screenshot rules: angle over -0.5:+0.5 s and placed-control/heat context selection; 7 control, 12 heat events. Original >=30 degree counts were 15/18.
Only angle threshold changed within each definition. No additional display smoothing. Curves show mean +/- SEM over events from one fly. No eligible events is not zero activity.
Angle is accumulated absolute yaw, not net heading change.

Main SVG/FIG correspond to page 1; screenshot_window SVG/FIG correspond to page 2. Event tables and MAT arrays included.
Reproduce in MATLAB: addpath to Imaging_analysis, run plot_20261020_11_2_turns_gt45_paper. The Python merger needs pypdf.
