# SMP544 inhibition and pivot-like turning

## Data and groups

- Source: `/Volumes/Neurobio/wilsonlab/Sophia/behavior/SMP544 unblinded data`
- Dates included: recordings from September 17 onward
- SMP544 experimental flies: global IDs 1, 3, 5, and 8
- SMP544 controls: global IDs 2, 10, and 12
- Each included fly contributed 80 optogenetic presentations: 320 experimental and 240 control presentations.
- The biologically matched reference is VES041 activation from the Chrimson cohort: 34 VES041-Chrimson flies and 47 readable GFP/empty-Chrimson controls, each with 80 presentations. One legacy control file, fly 39, was excluded from the presentation-level comparison because it is stored in an older non-HDF5 MATLAB format.
- VES041 activation is a contextual comparison rather than a contemporaneous control because the dates, effector, and stimulation protocol differ from the SMP544 experiment.

## Operational definition

A presentation was classified as pivot-like when the first second after light onset contained:

1. at least 45 degrees of absolute net heading change;
2. at most 2 mm of net displacement; and
3. directionality of at least 0.7, where directionality is absolute net heading change divided by total unwrapped heading change.

The same rule was applied to the matched one-second pre-light window. Results were also checked across 30, 45, and 60 degree turn thresholds and 1, 2, and 3 mm displacement thresholds.

The continuous analyses measured opto-minus-pre changes in absolute net heading and the 95th percentile of absolute angular velocity. Variability was quantified within each fly using the median absolute deviation of the pivot index and the standard deviation or median absolute deviation of opto-minus-pre net heading change.

SMP544 has four experimental and three control flies, so its results are presented descriptively without p-values, q-values, confidence intervals, or formal between-group tests. Fly-level values remain the unit of summary. Formal permutation tests are retained only for the larger VES041 activation cohort.

## Results

Pivot-like events were rare in the SMP544 groups. Experimental flies averaged 1.88% before light and 3.44% during the first second of light, a change of +1.56 percentage points. Controls averaged 3.33% before light and 2.92% during light, a change of -0.42 points. The descriptive difference in changes was +1.98 points.

The experimental-fly responses were heterogeneous:

- fly 1: 3.75% pre to 8.75% opto;
- fly 3: 2.50% pre to 1.25% opto;
- fly 5: 1.25% pre to 3.75% opto;
- fly 8: 0% pre to 0% opto.

All four SMP544 experimental flies had a positive mean opto-minus-pre change in absolute net heading. The experimental mean change was 8.87 degrees larger than the control mean change. The corresponding descriptive difference in the change in the 95th percentile of absolute angular velocity was +18.59 in the source units.

VES041 activation showed the same qualitative direction but a larger phenotype. VES041-Chrimson flies averaged 1.36% pivot-like presentations before light and 4.82% during light, a change of +3.46 percentage points. Its controls changed by +0.48 points. The activation-versus-control difference in changes was +2.98 points (permutation p = 0.00019, BH q = 0.00019). Descriptively, the SMP544 change was 1.89 points smaller than the VES041 activation change.

The continuous turning response was descriptively smaller under SMP544 inhibition than under VES041 activation. Mean opto-minus-pre absolute net heading was +10.26 degrees for SMP544 inhibition and +32.10 degrees for VES041 activation, a difference of -21.84 degrees. The mean change in the 95th percentile of absolute angular velocity was +21.87 versus +76.19, respectively, a difference of -54.32.

VES041 activation was also more variable by the pivot-index measure. Its mean within-fly median absolute deviation was 7.09, compared with 3.14 for SMP544 inhibition. The standard deviation of opto-minus-pre net heading was closer: 44.97 for VES041 activation and 40.89 for SMP544 inhibition. Across the nine threshold definitions, the descriptive SMP544 experimental-control differences were small and inconsistent.

## Interpretation

Some SMP544 inhibition trials visually resemble pivots, especially in flies 1 and 5, but similar high-scoring paths occur in controls. SMP544 inhibition and VES041 activation therefore point in the same qualitative direction, consistent with the proposed circuit relationship. The SMP544 effect is weaker and less consistently distinguishable from its controls, while VES041 activation produces a robust, larger, and generally more variable turning response. This comparison does not support describing SMP544 inhibition as more variable than VES041 activation.

## Outputs

- `pivot_rotation_vs_translation.png`: presentation-level rotation versus displacement
- `pivot_rate_pre_vs_opto.png`: paired fly-level pivot rates
- `pivot_response_variability.png`: within-fly variability comparisons
- `smp544_opto_aligned_kinematics.png`: fly-weighted angular and forward-velocity time courses
- `top_pivot_like_trajectories.png`: highest pivot-index SMP544 paths
- `fly_level_summary.csv`: fly-level estimates
- `ves041_activation_tests.csv`: fly-level permutation tests for the VES041 activation cohort only
- `pivot_threshold_sensitivity.csv`: threshold sensitivity analysis
- `presentation_metrics.csv`: one row per presentation
- `recording_manifest.csv`: included files and group assignments
