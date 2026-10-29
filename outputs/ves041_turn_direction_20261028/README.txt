VES041 csChrimson: turn direction per fly
Analysis date: 28 October 2026

QUESTION AND RESULT
Do VES041 activation flies turn more left or more right during light?
The pooled mean signed endpoint shift is -1.3559 degrees (leftward), with
a pointwise 95% fly-bootstrap interval of -7.9228 to +5.3596 degrees.
The two-sided fly-level sign-flip p-value is 0.695777. These data do not
establish a consistent population-wide left/right preference. This is not
an equivalence test and does not prove that the true directional bias is 0.

Line 1: 18 flies, mean -0.5595 degrees, 95% CI [-11.0830, +10.2605].
Line 2: 16 flies, mean -2.2518 degrees, 95% CI [-9.6678, +5.4300].
Observed fly means: 20 leftward and 14 rightward. These signs alone do not
establish statistically stable directional preferences for individual flies.

Of 2,720 light pulses, 1,338 have rightward net endpoints and 1,382 leftward
net endpoints. There are no exactly zero endpoints. The mean rightward
fraction is 49.1912%, with a 95% interval of 44.9632% to 53.4926%, obtained
by resampling flies. Pulses are not treated as independent population units.

WHAT THE FIGURES SHOW
VES041_turn_direction_per_fly.png/.svg/.pdf
  One colored dot per fly: its mean signed net heading change in the first
  1 s of light. Horizontal colored stems connect the point to 0 and are not
  error bars. Black diamonds: genotype means. Black error bars: pointwise
  95% percentile bootstrap confidence intervals across flies. Individual
  fly dots do not have uncertainty intervals. All 34 activation flies appear
  once, sorted by Fly ID within genotype.

VES041_rightward_pulse_fraction_per_fly.png/.svg/.pdf
  One colored dot per fly: percentage of its 80 pulses with positive net
  heading changes. The reference is 50%. There is no minimum turn-angle or
  walking-speed threshold. This classifies the direction of the net endpoint
  within a pulse, not the number of discrete left/right turn bouts within it.

METRIC AND COMPUTATION
For each pulse:
  signed turn = unwrapped heading at the end of the first light second
                minus unwrapped heading at light onset.
The extractor uses indices start through start + round(sample_rate) - 1.
At the saved sample rate of 666.6667 Hz, the two endpoints are about 0.999 s
apart. The initial heading is a coordinate reference, defined as 0 degrees.
The sign is retained. No absolute value or reflection is used for this
direction analysis. All 80 pulse estimates are averaged within each fly.
Population and genotype means give each fly equal weight.

The preceding 1 s is processed with the same signed endpoint definition.
Subtracting that fly's baseline mean from its light mean gives a pooled
shift of +1.5312 degrees, 95% CI [-4.9912, +8.6391], with an exploratory
two-sided sign-flip p of 0.676057. This sensitivity does not support a
consistent stimulation-induced directional shift.

STATISTICS
Seed 20261028. There are 100,000 bootstrap resamples of flies for each
pointwise 95% percentile interval. These are not simultaneous intervals.
The pooled VES041 interval resamples the 34 flies with equal fly weights.
The primary exploratory test uses 200,000 random sign flips of the 34 fly
means. It tests whether the mean signed light response is centered at 0,
assuming independent flies and a symmetric null distribution. The p-value
uses the plus-one Monte Carlo correction. Genotype intervals are descriptive
decompositions of the pooled cohort; no genotype-specific p-values are added.
An exact binomial sign-count sensitivity for 14 rightward versus 20 leftward
fly means gives p = 0.391528. The sign-count and baseline sensitivity tests
are unadjusted, exploratory checks, not additional confirmatory tests.

LEFT/RIGHT CONVENTION
Labels follow the existing lab analysis-code convention: positive = right,
negative = left. This is supported by the following copied source files:
  findYawVelPeaksFT.m: rightTurn uses positive angular-velocity peaks, while
  the left-turn branch multiplies angular velocity by -1 before peak finding.
  processExptData.m: headingPosition is yawAngCumPos and angularVelocity is
  yawAngVel, with no sign inversion between those assignments.
  dsFiltFictrac.m and computeSmoothedVelocity.m: yaw angular velocity is
  obtained from the gradient of smoothed yaw position, preserving its sign.
These labels follow the code convention; no new hardware calibration was run.

INPUTS AND SCOPE
The exact local pulse-export snapshot is in source_data/presentation_metrics.csv.
Original source:
  /Users/sophiarenauld/Documents/GitHub/Wilson_Thesis/Behavior_analysis/smp544_pivot_results/presentation_metrics.csv
Only dataset == ves041_activation is selected: 6,480 pulses from 81 flies,
80 pulses per fly. The main figure selects the 34 VES041 x csChrimson flies
(18 line 1, 16 line 2). GFP and empty-split control summaries remain in the
CSV outputs for context. Empty + csChrimson flies are negative controls and
are not included in the pooled VES041 activation estimate.
Genotype mapping comes from cschrim_behaviorcrunch.R, which records the
original blinded fly-ID assignments. The source CSV also contains SMP544
data; those records never enter this analysis.

There is no imputation, new fly exclusion, outlier removal, or movement
filter. No raw acquisition processing is rerun. Exact input hashes and
numerical results are in analysis_results.json. The saved summary includes
signed endpoint data, so a server connection is unnecessary.

REPRODUCING
  python3 analyze_turn_direction.py
  Rscript plot_turn_direction.R

Python dependencies: numpy, pandas.
R dependencies: ggplot2, dplyr, jsonlite, ragg, svglite.
Arial is used for figure text. R_session_info.txt records the environment.

fly_turn_direction.csv contains all 81 fly estimates, including controls.
group_turn_direction.csv contains genotype and pooled activation summaries.
analysis_results.json contains settings, results and source hashes.
