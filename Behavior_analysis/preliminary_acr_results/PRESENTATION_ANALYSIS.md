# Individual optogenetic-presentation analysis

## Bottom line

The presentation-level analysis reinforces the fly-level conclusion. The behavioral response depends strongly on whether the fly was moving when the light came on, but the response is not detectably different between ACR/rubyACR flies and GFP/empty controls in either cohort.

Among presentations that began while the fly was moving, experimental flies increased forward velocity and showed less total heading change during the light. Presentations that began while the fly was not moving also showed increased forward velocity, but total heading change increased as the fly initiated movement. Path straightness increased strongly in the not-moving category because pre-stimulus displacement/path length was near its floor; this should not be interpreted as an already-walking trajectory becoming straighter.

## Movement-state definition and sample sizes

“Moving at onset” means mean total speed was at least 3 during the final 250 ms before light onset. A 1-second-pre-window classification is included in the presentation table as a sensitivity variable.

| Cohort | Group | Not moving | Moving | Flies |
|---|---:|---:|---:|---:|
| Blinded round 1 | Control | 753 | 847 | 18 |
| Blinded round 1 | Experimental | 447 | 673 | 13 |
| Unblinded | Control | 749 | 611 | 16 |
| Unblinded | Experimental | 790 | 810 | 18 |

There were 5,680 usable presentations from 65 flies. One unblinded GFP recording remained unreadable and is listed in `skipped_recordings.csv`.

## Experimental-fly responses

Mixed-model estimated during-minus-pre changes are shown below. Positive heading-change values mean *less* total heading change. Positive forward-velocity values mean faster movement.

| Cohort | State at onset | Less heading change (deg) | Mean forward velocity | Path straightness |
|---|---|---:|---:|---:|
| Blinded round 1 | Moving | +9.78 | +2.95 | +0.069 |
| Blinded round 1 | Not moving | -12.54 | +2.60 | +0.316 |
| Unblinded | Moving | +8.56 | +1.87 | +0.031 |
| Unblinded | Not moving | -11.58 | +1.53 | +0.237 |

For moving presentations, reduced heading change and increased forward velocity were significant within experimental flies in both cohorts after BH correction. The same broad state-dependent pattern occurred in controls.

## Experimental versus control tests

- No experimental-control contrast was significant after BH correction for either movement state in either cohort.
- The smallest unadjusted experimental-control result was higher forward velocity in blinded-round-1 not-moving presentations: difference +0.91, p=0.072, BH q=0.359.
- Formal group-by-movement-state interactions were also null for every outcome: all BH q-values were at least 0.736.
- Thus, separating presentations by movement state clarifies the behavioral geometry but does not reveal a cell-type-specific group effect.

## Baseline wiggliness at the presentation level

Within experimental flies, moving presentations with more pre-opto heading change or lower pre-opto heading consistency tended to show greater straightening/heading stabilization and, for several predictor-outcome pairs, greater speed increases. However, controls showed similar relationships. After allowing each fly to have its own wiggliness slope and correcting within test families, no experimental-control wiggliness-slope difference had BH q<0.10 in either cohort or movement state.

This is important because models with only a random fly intercept produced overly precise presentation-level slopes. The final reported models include random intercepts and random slopes by fly.

## Statistical model

Response models were fit separately by cohort and outcome:

`response ~ experimental_group * movement_state + presentation_number + (1 + movement_state | fly)`

Baseline-wiggliness models were fit separately by cohort, movement state, predictor, and outcome:

`response ~ standardized_wiggliness * experimental_group + presentation_number + (1 + standardized_wiggliness | fly)`

P-values use Satterthwaite degrees of freedom. Benjamini-Hochberg q-values were calculated within cohort, movement state, and comparison family. Heading-consistency response models had singular random-slope fits, so their mixed-model uncertainty estimates deserve extra caution.

## Files

- `presentation_level_metrics.csv`: one row per opto presentation, including both movement classifications.
- `movement_state_counts.csv`: presentation and fly counts by cohort, group, and movement state.
- `presentation_mixed_model_response_tests.csv`: marginal responses, group contrasts, movement effects, and group-by-movement interactions.
- `presentation_baseline_wiggliness_models.csv`: experimental/control slopes and slope contrasts.
- `presentation_model_diagnostics.csv`: random-effect and singular-fit diagnostics.
- `presentation_level_changes_by_movement.png`: all presentations by cohort, group, state, and outcome.
- `presentation_sequence_by_movement.png`: fly-weighted response over presentation-number bins.
