# Behavior analysis inventory and migration proposal

Status update November 10, 2026: at Sophia's request, 868 non-code files (4.499 GB) were moved from `Behavior_analysis` and 23 identified behavior output directories to `/Users/sophiarenauld/behavior_results`. The destination separates `results`, `inputs`, and `local_state`, preserving repository-relative paths. Every moved file passed full SHA-256 verification; 92 scripts, autosaves, and project files stayed in place and unchanged. Empty source directories were removed after their contents were moved. Imaging results and server data were not moved.

The [executed relocation plan](relocation_plan_20261110.json) and [external manifest](/Users/sophiarenauld/behavior_results/migration_20261110/manifest.json) map original paths to current paths. The [verification record](/Users/sophiarenauld/behavior_results/migration_20261110/verification.json) records completion and reversal guidance. Existing analysis scripts have not been patched and may fail on old paths or recreate old output folders until a later path update. No scientific analysis was rerun.

The inventory and proposal below describe the pre-move snapshot. Their original source links and proposed destinations are historical; use the executed plan for current file locations. Code consolidation and restructuring remain proposed, not implemented.

Start with the blinded round 1 VES041 5-second inhibition analysis. Preserve its existing inputs and results as reference checks, separate its shared functions from its run settings, and verify that the reorganized code reproduces those checks before touching other analyses.

## Inventory

The inventory covers every existing file under `Behavior_analysis`, including hidden and Git-ignored files. The 152 files occupy 1.079 GB using decimal units.

| File type | Count |
| --- | ---: |
| Source scripts | 72 |
| Generated results, accompanying reports, and figures | 41 |
| Processed recording or table inputs | 6 |
| Editor and runtime state, including an autosave | 32 |
| RStudio project configuration | 1 |

There are 69 scripts directly in the top-level folder. Across subfolders, the languages are MATLAB (34), Python (22), and R (16). The processed recording `2025_10_29_fly02_cell01_trial01B_processed_final.mat` alone is 869.7 MB. The SMP544 aligned-trace result table is another 182.7 MB. Neither belongs among source files.

Read the [file-by-file migration map](file_migration_map.md) for every proposed action and destination. The [machine-readable inventory](inventory.json) includes file sizes, hashes where available, Git tracking status, function definitions, explicit Python imports, candidate references, and source-line evidence for hard-coded paths. Eighteen existing repository output directories named in the behavior code are listed separately at directory level. This is not an inventory of everything on Stacks, the server, or the repository-wide `outputs` folder.

## Proposed boundaries

```text
Behavior_analysis/
    README.md
    pipelines/       Complete analysis workflows and thin entry points
    core/            Data loading, preprocessing, detectors and summaries
    plotting/        Cohort figures, example traces and export functions
    configs/         Named methods, reviewed cohort definitions and styles
    tests/           Synthetic checks and reference-result comparisons
    Behavior_analysis.Rproj

Data_management/behavior/       Server copy and batch-management utilities
legacy_analysis/behavior/      Older code retained after dependency review
```

These are proposed paths, not directories created by this inventory. Keep MATLAB, Python, and R initially. A language conversion is not needed to separate inputs, calculations, and plotting.

The file map preserves many existing script names inside the proposed folders to make the first migration traceable. Longer term, those scripts should be short entry points calling shared functions. The map is not a recommendation to keep each existing script as a separate numerical implementation.

`Behavior_analysis` should contain source, small hand-maintained configurations, tests, and documentation only. It should not contain recordings, derived measurement tables, example MAT bundles, saved workspaces, figures, caches, or generated reports. Hand-maintained fly/genotype definitions need a reviewed source and version; measured data must not be disguised as configuration.

## Proposed external locations

| Setting | Proposed local path | Contents |
| --- | --- | --- |
| `RESULTS_ROOT` | `/Users/sophiarenauld/stacks/behavior_analysis_results` | Analysis runs, tables, figures, caches and generated methods |
| `INPUT_ROOT` | `/Users/sophiarenauld/stacks/behavior_analysis_inputs` | Existing local input tables and processed-recording references where appropriate |
| `LOCAL_STATE_BACKUP` | `/Users/sophiarenauld/stacks/behavior_analysis_results/legacy_local_state` | Preserved editor history, runtime state and autosaves |

These locations have not been created or checked for availability, permissions, capacity, or backup. Approve the roots and local/server mapping before migration. Server recordings remain authoritative where applicable; this proposal does not require copying all recordings into a new local directory.

New analyses would save to `RESULTS_ROOT/<analysis_name>/<run_id>/`, with `figures/`, `tables/`, `cache/`, and `run_info.json`. Save the exact configuration, input manifest, inclusion and exclusion reasons, code version, units, relevant package versions, random seed when used, and failed-recording log with the run. Never overwrite a previous run silently. Keep each older result folder and its reports together under `legacy_snapshots/` without relabeling it as a new analysis.

Analysis entry points should require or validate an external results root before writing. Tests should reject a results destination inside the source tree. Output roots must also be separate from input recordings so that reprocessing cannot overwrite an input.

## Redundancies and scientific distinctions

The only byte-identical source pair found is `Behavior_Crunch/hfacr_behavior_crunch.R` and `Behavior_Crunch/ves041_behavior_crunch.R`, each 343 lines. Their matching SHA-256 values are recorded in the inventory. Choose a canonical table-summary implementation after checking callers and table provenance. Do not equate this older table analysis with the pulse-level preliminary analysis.

Flat-path exports repeat figure sizing, axis matching, legend placement, and SVG export. Their selected recordings and scaling policies should become recipes passed to a shared exporter. The trace renderers similarly repeat plotting with raw, 10 ms-smoothed, circled, and reference-scale displays. Preserve those choices explicitly; display smoothing must not silently change event detection.

The generic SVG flattener and the Chrimson-specific flattener overlap, but are not interchangeable yet. The specialized version labels individual flies, genotype means, and SEM bands; the generic version assigns generic object IDs. Preserve semantic IDs and independent editable objects when consolidating them.

| Scientific method | Distinction that must survive cleanup |
| --- | --- |
| Localized spatial turns | Net heading, net XY displacement and directionality within sliding windows; overlapping qualifying windows merge. Preserve each run's threshold and strict versus inclusive comparison. |
| Forward-speed-linked turns | Directional turn bouts with preceding slowdown or following speed-up. The primary recipe uses at least 30 degrees and 2 mm/s forward-speed changes in adjacent 0.5-second windows; it has no low-displacement requirement. |
| Peak-yaw forward gate | Angle/directionality candidates gated by forward velocity at peak yaw; the XY displacement condition is removed. Raw and processed-field variants remain distinct. |
| Imaging behavioral rule applied to inhibition | Slowdown-selected epochs, subsequent rotation, validity/context requirements and refractory ordering differ from the spatial rule. Preserve the reference validation against imaging code. |
| Early heading and wiggliness analyses | Matched windows, movement states, metrics and mixed models differ from the later blinded 5-second good-walker bout analysis. |
| SMP544 onset classification | First-second pre/on classification and descriptive summaries are not continuous merged-bout frequency. Preserve the no-statistical-tests policy for SMP544. VES041 activation reference tests are separate. |
| Activation scalar summaries and onset traces | Pre/on scalar values and time-resolved raw kinematic averages are separate products. Preserve presentation-within-fly averaging and equal-fly weighting where used. |

Genotypes, cohort membership, source fields, units, movement thresholds, time windows, missing-data rules, and sampling rates must not change as a side effect of reorganizing code. Do not pick a new default definition based on which version gives the clearest group difference.

## Dependencies to resolve before moving files

| Dependency | Required handling |
| --- | --- |
| Python scripts import `read_vector` from the walking audit and `aggregate`/`inference` from the spatial-turn analysis | Extract these functions into core modules with reference tests. Update imports together; moving scripts independently will break them. |
| Four example-selection/export scripts import `read` from the peak-forward scientific analysis | Put the HDF5 reader in core I/O. Example selection must not require importing an unrelated scientific method. |
| Example export/preparation imports `render_ves041_reorientations_clean` for loading/configuration | Separate the loader and example recipes from SVG drawing before consolidating renderers. |
| `apply_imaging_rules_to_inhibition.m` calls `Imaging_analysis/pizza_ves041_turn_angle.m` for validation | Keep that explicit reference dependency. Imaging source changes are outside this cleanup. |
| `save_plot_with_title_as_filename.m` is referenced by many imaging scripts | Keep compatibility at the old entry point or configure a shared MATLAB path before relocating its implementation. Do not treat it as behavior-only. |
| Scripts read prior dated outputs and caches | Replace these links with explicit input manifests/run references. Preserve cached field names, rates and units; a cache is not an independent source of cohort identity. |
| Activation code reads `.build_ves041_activation/reproducibility/source_data` | Recover and version the reviewed cohort mapping and input references. Do not discard that directory as a temporary build before resolving these dependencies. |
| Python uses `Path(__file__).resolve().parents[1]`; MATLAB/R use `pwd` or relative paths | Nesting files changes path resolution. Replace location assumptions before relocation and test from a different working directory. |

`master_behavior_analysis.m` currently combines file discovery, processing, plotting, and MAT saving. Its default starts at file 27 unless a configuration is supplied, and its discovery pattern is `*function1_pro.mat` despite comments describing processed-final files. Record those behaviors as migration risks; do not silently change them during a move.

`export_ves041_flatpaths_matched_scale.m` declares `export_ves041_flatpaths_common_axes()`. Align its function and filename during an approved refactor, then test the intended entry point. Several older R scripts also rely on workspace/package state; test them in clean sessions before calling them supported pipelines.

## Proposed migration order

1. Approve this map and the external roots. Preserve current dirty/untracked work with an agreed recoverable snapshot before any code refactoring. No commit, deletion, or file move is part of this inventory.
2. Freeze the VES041 blinded 5-second reference inputs, genotype/fly mapping, good-walker list, opto boundaries, per-presentation counts and exposure, and fly-level summaries. Keep all-fly and good-walker analyses separately identifiable.
3. Extract tested readers, opto-window segmentation, movement QC, localized-turn detection, and fly aggregation. Keep old entry points usable until their callers are updated. Preserve the original statistical procedure when inference is enabled.
4. Reproduce the existing spatial-turn analysis and figures using the same inputs and parameters. Compare event boundaries, counts, denominators, and fly-level values, not just the appearance of the final plot. Review numerical tolerances explicitly.
5. Redirect outputs through configuration. For any later approved relocation, copy and verify existing inputs/results first; retain the originals until paths and reproducibility checks pass. There is no deletion step authorized here.
6. Consolidate example plotting and exports, preserving recipes. Then migrate speed-linked and peak-forward methods, activation, SMP544, and early heading analyses in separate reviewable steps.
7. Review the ten older saccade/menotaxis/workspace scripts with their users before legacy relocation. Lack of a static caller does not prove that a function is unused outside this repository.

## Checks required before retirement or consolidation

Use synthetic cases for heading wraparound, signed rotations, directionality, overlapping-window merging, threshold equality, missing data and gaps, phase edges, duration/exposure, and movement classification. For real reference recordings, retain sample/event indices, source fields, units and sampling rates. Check that presentations are summarized within fly before calculating across-fly uncertainty where that is the established method.

Figure checks should cover unchanged data coordinates, equal XY scale, common bounds where requested, and independently editable SVG traces, means and SEM bands. The source-tree output guard should fail before any file is written. Failure on an unreadable recording should produce an explicit log, not a silent change in sample size.

The inventory is based on file metadata and static source inspection. It does not establish that every script runs in a clean environment, that candidate references are actual runtime calls, or that different implementations are scientifically equivalent. Files up to 2 MB have full hashes; larger files have size and modification time only. A later data migration needs full content verification.

The next approval should cover the VES041 reference checks and external output location. It should not yet approve deleting legacy code or changing any biological inclusion criterion.
