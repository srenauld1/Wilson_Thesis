# Preliminary ACR/rubyACR silencing analysis

## Design and interpretation

- Unit of inference: fly. Pulses and repeated recordings were averaged within fly.
- Windows: 1 s immediately before opto onset and the first 1 s during opto.
- Experimental: names containing ACR/rubyACR, except `empty`; controls: names containing GFP or `empty`.
- Heading change is total absolute unwrapped angular displacement. Path straightness is displacement/path length. Heading consistency is circular mean resultant length (rho), excluding samples with total speed < 3.
- Forward walking speed is averaged over samples with forward velocity > 1; mean forward velocity over all samples is included as a sensitivity measure.
- All p-values are two-sided and exploratory. Benjamini-Hochberg q-values are provided within each cohort and test family.
- Positive hypothesis-oriented changes mean less heading change, more straightness/consistency, or faster walking.

## Sample sizes

- blinded_round1, control: 18 flies
- blinded_round1, experimental: 13 flies
- unblinded, control: 16 flies
- unblinded, experimental: 18 flies

Usable opto pulses: 5680 across 65 flies.

## Fly-level response summary

### blinded_round1

- less heading change (deg): experimental mean 0.661 (95% bootstrap CI -4.37 to 5.88; sign-flip p=0.811, BH q=0.811); control mean -0.381; experimental−control 1.04 (permutation p=0.768, BH q=0.778).
- more path straightness: experimental mean 0.167 (95% bootstrap CI 0.134 to 0.196; sign-flip p=0.00025, BH q=0.000417); control mean 0.159; experimental−control 0.00749 (permutation p=0.778, BH q=0.778).
- more heading consistency (rho): experimental mean 0.00904 (95% bootstrap CI 0.0027 to 0.0161; sign-flip p=0.021, BH q=0.0263); control mean 0.0124; experimental−control -0.00336 (permutation p=0.564, BH q=0.778).
- faster while walking: experimental mean 2.09 (95% bootstrap CI 1.22 to 3.02; sign-flip p=0.00025, BH q=0.000417); control mean 1.73; experimental−control 0.362 (permutation p=0.547, BH q=0.778).
- higher mean forward velocity: experimental mean 2.71 (95% bootstrap CI 1.98 to 3.53; sign-flip p=0.0002, BH q=0.000417); control mean 2.02; experimental−control 0.684 (permutation p=0.196, BH q=0.778).

### unblinded

- less heading change (deg): experimental mean -0.802 (95% bootstrap CI -5.5 to 4.01; sign-flip p=0.75, BH q=0.75); control mean -5.43; experimental−control 4.62 (permutation p=0.177, BH q=0.483).
- more path straightness: experimental mean 0.126 (95% bootstrap CI 0.0995 to 0.152; sign-flip p=5e-05, BH q=0.00025); control mean 0.133; experimental−control -0.00756 (permutation p=0.767, BH q=0.767).
- more heading consistency (rho): experimental mean 0.0097 (95% bootstrap CI 0.00331 to 0.0162; sign-flip p=0.0119, BH q=0.0149); control mean 0.00593; experimental−control 0.00377 (permutation p=0.488, BH q=0.61).
- faster while walking: experimental mean 1.4 (95% bootstrap CI 0.519 to 2.26; sign-flip p=0.00775, BH q=0.0129); control mean 0.9; experimental−control 0.498 (permutation p=0.394, BH q=0.61).
- higher mean forward velocity: experimental mean 1.73 (95% bootstrap CI 0.99 to 2.52; sign-flip p=0.0005, BH q=0.00125); control mean 1.07; experimental−control 0.659 (permutation p=0.193, BH q=0.483).

## Baseline wiggliness relationships in experimental flies

### blinded_round1

- high pre heading change → less heading change: Spearman rho=0.615, p=0.0286, BH q=0.0716, n=13.
- high pre heading change → more straightness: Spearman rho=-0.335, p=0.257, BH q=0.386, n=13.
- high pre heading change → more heading consistency: Spearman rho=-0.033, p=0.92, BH q=0.92, n=13.
- high pre heading change → faster while walking: Spearman rho=0.758, p=0.0037, BH q=0.0111, n=13.
- high pre heading change → higher mean forward velocity: Spearman rho=0.78, p=0.00315, BH q=0.0111, n=13.
- low pre straightness → less heading change: Spearman rho=-0.363, p=0.225, BH q=0.375, n=13.
- low pre straightness → more straightness: Spearman rho=0.808, p=0.00145, BH q=0.0111, n=13.
- low pre straightness → more heading consistency: Spearman rho=0.203, p=0.512, BH q=0.636, n=13.
- low pre straightness → faster while walking: Spearman rho=-0.181, p=0.554, BH q=0.636, n=13.
- low pre straightness → higher mean forward velocity: Spearman rho=-0.533, p=0.0631, BH q=0.135, n=13.
- low pre rho → less heading change: Spearman rho=0.484, p=0.0955, BH q=0.179, n=13.
- low pre rho → more straightness: Spearman rho=-0.209, p=0.491, BH q=0.636, n=13.
- low pre rho → more heading consistency: Spearman rho=0.165, p=0.593, BH q=0.636, n=13.
- low pre rho → faster while walking: Spearman rho=0.786, p=0.00205, BH q=0.0111, n=13.
- low pre rho → higher mean forward velocity: Spearman rho=0.775, p=0.0027, BH q=0.0111, n=13.

### unblinded

- high pre heading change → less heading change: Spearman rho=0.662, p=0.00315, BH q=0.0236, n=18.
- high pre heading change → more straightness: Spearman rho=0.259, p=0.3, BH q=0.346, n=18.
- high pre heading change → more heading consistency: Spearman rho=0.6, p=0.0105, BH q=0.0342, n=18.
- high pre heading change → faster while walking: Spearman rho=0.492, p=0.0402, BH q=0.0603, n=18.
- high pre heading change → higher mean forward velocity: Spearman rho=0.591, p=0.0114, BH q=0.0342, n=18.
- low pre straightness → less heading change: Spearman rho=-0.315, p=0.208, BH q=0.26, n=18.
- low pre straightness → more straightness: Spearman rho=0.0609, p=0.809, BH q=0.809, n=18.
- low pre straightness → more heading consistency: Spearman rho=-0.15, p=0.552, BH q=0.591, n=18.
- low pre straightness → faster while walking: Spearman rho=-0.406, p=0.0935, BH q=0.128, n=18.
- low pre straightness → higher mean forward velocity: Spearman rho=-0.519, p=0.0267, BH q=0.0572, n=18.
- low pre rho → less heading change: Spearman rho=0.633, p=0.0055, BH q=0.0275, n=18.
- low pre rho → more straightness: Spearman rho=0.494, p=0.0375, BH q=0.0603, n=18.
- low pre rho → more heading consistency: Spearman rho=0.785, p=0.00025, BH q=0.00375, n=18.
- low pre rho → faster while walking: Spearman rho=0.494, p=0.0391, BH q=0.0603, n=18.
- low pre rho → higher mean forward velocity: Spearman rho=0.546, p=0.0219, BH q=0.0549, n=18.

## Important caveats

- This is exploratory and uses several correlated outcomes/predictors; isolated p<0.05 results should not be treated as confirmatory.
- Relating a baseline metric to its own change score is vulnerable to regression to the mean. The baseline-response table therefore also includes an experimental-vs-control slope interaction (HC3 robust p-value) as a more informative check.
- Averaging across pulses gives each fly equal weight but does not model pulse-to-pulse habituation or recording date/batch.
- The walking-conditioned speed measure describes speed while moving forward; the all-sample forward-velocity sensitivity outcome also reflects stopping/turning.

## Output files

- `fly_level_metrics.csv`: analysis unit and derived responses.
- `pulse_level_metrics.csv`: per-pulse audit table.
- `fly_level_response_tests.csv`: directional changes and experimental-control permutation tests.
- `baseline_wiggliness_response_tests.csv`: Spearman associations and group-slope interactions.
- `joint_response_counts.csv`: flies that became both faster and less turning.
- PNG figures show fly-level change distributions and the joint response space.
