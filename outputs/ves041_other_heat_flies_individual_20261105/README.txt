Three-fly source cohort and individual plots

The screenshot combines two different source summaries:
1. Left pair: October 20 fly 1, recording 20261020_11_2. Original-criteria slowdown/reorientation averages, aligned to strongest forward deceleration. Conditions are placed control barriers (15 events) and heat (18 events). These are not calcium-peak-triggered averages.
2. Right pair: three-fly pooled heat-off/heat-on summary, from pooled_four_original_pizza_rules_20261102/baseline_subtracted_dff/excluding_april27_fly3/crunch_data.mat. Each event has its -2:-1 s neural baseline removed; a display offset of +1 was then added in the final SVG.

The three pooled flies are:
20260427_fly2: recording 20260427-2_ves041_fly2_heat; 6 heat-off, 3 heat-on turns.
20261019_fly1: recording 20261019_4_1; 17 heat-off, 7 heat-on turns.
20261020_fly1: recording 20261020_11_2; 20 heat-off, 19 heat-on turns.
Totals: 43 heat-off and 29 heat-on turns.
April 27 fly 3 (20260427-7_ves041_fly3_heat_20mm) was explicitly excluded in the saved three-fly source. It is not one of these three.

The other two flies are April 27 fly 2 and October 19 fly 1. Their *_heat_off_on_individual plots use the exact saved pooled-source event arrays, separately per fly; the groups were not redetected or filtered.
Heat off/on is the original split at the event sample. Heat off can include accommodation and zero-command portions of heat trials, so it is not a strict spontaneous/stimulus classification. Heat transitions elsewhere in the aligned window are retained in this source.
The single-example selection and the three-fly selection differ: single example uses >=30 degrees accumulated absolute median-smoothed yaw in +/-0.5 s and protocol control/heat contexts; pooled source uses >=15 degrees accumulated absolute median-smoothed yaw in 0:+0.5 s and heat at event time.

Additional *_calcium_peaks_and_turns_by_heat_state plots distinguish actual calcium-peak-triggered events (original analysis.mat group1) from slowdown+turn events (group4), each split by heat off/on for each fly. All original saved events are retained. These plots are included to resolve the calcium-peak versus turn-average labeling ambiguity.
April 27 fly 2 has 70/31 heat-off/on calcium peaks. October 19 fly 1 has 30/5 heat-off/on calcium peaks.
Neural traces in all new per-fly figures use per-event -2:-1 s baseline subtraction plus 1, matching the pooled screenshot. Kinematic units and event arrays are unchanged; no extra smoothing or signed-velocity averaging is added.
Shading in per-fly figures is descriptive SEM across events within that one fly; it is not between-fly SEM. The main plots share y-axis limits across both flies and conditions.

Why the example rotation mean is lower than individual event peaks:
The original arrays already contain absolute rotational speed before averaging. The traces are aligned to deceleration, so yaw maxima occur at different times across events.
Within +/-0.5 s, control event median peak is 109.37 deg/s, versus 63.44 deg/s for the maximum of the deceleration-aligned mean in that window. For heat, the corresponding values are 168.35 and 119.88 deg/s.
A separate diagnostic aligns these same events to their yaw maximum within +/-0.5 s. Mean peaks then become 106.16 (control) and 168.23 (heat) deg/s. This is an alignment diagnostic, not a replacement for the original analysis or new event selection.

Files include editable MATLAB FIGs, vector SVGs, PNGs, saved MAT arrays, event-count and rotation-audit CSVs, and plotting scripts. The original source data and figures are unchanged.
