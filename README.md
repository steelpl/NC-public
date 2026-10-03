# NC-public

MATLAB code for **N-Tuple Collocation for Reference-Free Uncertainty Estimation from Multi-Source Hydrological Data: Beyond Triple Collocation and Three-Cornered Hat**, by Hung T. Pham and Seokhyeon Kim.

The repository includes scripts for reproducing the synthetic covariance experiments, TC–NC and TCH–NC comparisons, NC sensitivity tests, downstream weighted-averaging experiments, graphical abstract, and supplementary feasibility diagnostics. All input samples are generated internally. The repository contains code only; experiment outputs are generated locally.

## Requirements

- MATLAB R2025b was used for verification. Earlier releases have not been validated.
- The estimators, paired synthetic experiments, numerical checks, feasibility diagnostics, and graphical abstract use base MATLAB. No Optimization Toolbox is required.
- The four original experiment scripts require Statistics and Machine Learning Toolbox for their percentile/boxplot plotting functions.
- No external datasets, personal paths, or network drives are required.

## Start here

Open this folder in MATLAB and run the dependency check. A call with no arguments performs this check only:

```matlab
Run_All_Synthetic;
```

A reduced-replication software check explicitly starts calculations and writes a new timestamped folder under `results/`:

```matlab
options = struct('mode','run','nMC',2,'nSensitivity',2,'nLegacy',2, ...
    'runLegacy',true,'makePublicationOutputs',true);
check = Run_All_Synthetic([],options);
```

Reduced-replication outputs are software checks and must not be used as the reported Monte Carlo results. To reproduce the full synthetic design, explicitly request the production run:

```matlab
runInfo = Run_All_Synthetic([],struct('mode','run'));
saved = runInfo.outputRoot;
```

This uses 500 replicates for each original experiment and each of the 36 paired scenarios, plus 100 sensitivity realizations. Runtime depends on the computer; the runner records stage timings and progress. The workflow is sequential and does not require Parallel Computing Toolbox. Set `runLegacy=false` to run only the base-MATLAB paired/sensitivity suite. A new or empty output folder is required; existing runs are not overwritten or resumed.

The options `nMC`, `nSensitivity`, `nLegacy`, `runLegacy`, `makePublicationOutputs`, and `invisibleFigures` are centralized in `Run_All_Synthetic.m`. Changing replication counts or the scientific configuration creates a different experiment.

## Entry points and output figures

The integrated runner executes the following scripts in isolated stage folders and saves their workspaces and figures. They can also be run individually from a separate output working directory with the repository on the MATLAB path.

- `Sim_TC_vs_NC.m`: original covariance-level TC–NC comparison; `Sim_TC_vs_NC_results.png` (design corresponding to original Figure 1; current NC).
- `Sim_TCH_vs_NC.m`: original TCH–NC comparison and scaling test; `Sim_TCH_vs_NC_results.png` (design corresponding to original Figure 2; current NC).
- `Test_NC_Sensitivity.m`: original SNR, product-number, and sample-length tests; `Test_NC_Sensitivity_results.png` (design corresponding to original Figure 3; current NC).
- `Sim_WA_Comparison.m`: original covariance-based weighted-averaging experiment; `Sim_WA_Comparison.png` (design corresponding to original Figure 4; current NC).
- `run_revision_synthetic(nMC,outDir)`: the additional paired stationary/seasonal experiments with independent training/test samples, known-covariance and equal-weight controls, and raw/PSD diagnostics.
- `run_revision_identifiability(nMC,outDir)`: exact p=2/p=3 decomposition examples and six-initialization/five-floor sensitivity tests.
- `Build_Synthetic_Revision_Artifacts(mainDir,sensitivityDir,newOutputDir)`: generates tables and publication figures from completed CSV outputs, without fitting the models again.
- `Draw_JoH_Graphical_Abstract(outputDir)`: generates the graphical abstract in PNG, vector PDF, and editable MATLAB FIG. This entry point requires no simulations.
- `Figure_Supplementary_NC_PSD_Diagnostic.m`: optional original synthetic-pixel PSD diagnostic; `FigS_NC_PSD_Diagnostic.png/pdf`. Its same-sample validation proxy is a numerical diagnostic, distinct from the independent-test merging analysis.
- `Test_Revision_Core(outputDir)`: 16 numerical checks, including PSD feasibility, decomposition preservation, and weighted-averaging behavior.

All NC experiments call the single current `NC.m` implementation. The four covariance-level scripts retain the original experimental designs, but now use the current NC optimizer and its default residual PSD feasibility map. Their newly calculated values and figures can differ from the originally submitted Figures 1–4. For the historical pre-revision NC implementation and experiment sources, see source commit `1cb6c536a2f67724b4d23942228a6f100c00474e` (8 July 2026) in `steelpl/n-tuple-collocation`; the current file tree contains no separate submitted NC estimator. `WA_Submitted.m` remains the original ridge/backslash WA rule for the covariance-level WA comparison; additional paired experiments use current `WA.m`. TC, TCH and the covariance generator are unchanged from that source. The WA plotting script includes a display-only correction that preserves missing groups and method colors.

The paired experiments call the current `NC.m`/`WA.m`, add unit-gain independent-test comparisons, and retain every attempted realization. Their numerical values answer different experimental questions and should not be interchanged with the original figures. Both estimator implementations are needed to reproduce the original and additional experimental stages without silently changing the baseline.

The paired-experiment figure filenames in `09_publication_artifacts/` are:

- `Figure_R1_Variance_RMSE_alpha_0`: variance accuracy with signal–error orthogonality.
- `Figure_R1_Variance_RMSE_alpha_0p3`: deliberate signal–error correlation diagnostic.
- `Figure_R2_Heldout_Merging`: independent-test weighted-averaging MSE.
- `Figure_R3_Stress_Feasibility_and_PSD_Cost`: severe-correlation feasibility and PSD-map cost.
- `Figure_R4_Initialization_and_Eigenfloor` and `Figure_R4b_Initialization`: numerical sensitivity.
- `Figure_R5_Merging_Controls`: equal-weight and known-error-covariance controls.
- `Figure_R6_Paired_Seasonal_Change`: paired seasonal-minus-stationary differences.
- `Graphical_Abstract_NC`: conceptual graphical abstract.

`R1`–`R6` are filename identifiers, not final manuscript numbering. Supplementary CSVs include complete denominators, method contrasts, initialization/floor settings, and the p=2/p=3 nonuniqueness examples. Figure exports use MATLAB throughout: the paired figures use 600-dpi PNG, vector PDF and editable FIG; the graphical abstract uses 300-dpi PNG, vector PDF and FIG. `JoHFigureStyle.m` controls the common paired-figure fonts, palette and single legend.

To audit or regenerate figures from a completed integrated run:

```matlab
Run_All_Synthetic(saved,struct('mode','audit'));
p = Run_All_Synthetic(saved,struct('mode','postprocess'));
figureDir = p.publication.artifactDir;
```

To generate only the graphical abstract:

```matlab
Draw_JoH_Graphical_Abstract(fullfile(pwd,'results','graphical_abstract'));
```

## Reproducibility and interpretation

The original four scripts use seed 123. The paired suite uses base seed 20260917, and the sensitivity suite uses 20260918, each advanced by replicate index. Core numerical checks use 83955; the optional original PSD diagnostic uses 20260707. Configurations, source snapshots/hashes, all attempts, summaries, stage timings and the validation report are saved in each integrated run. Source snapshots are local outputs, not additional repository content.

The paired suite has 36 scenarios, 1095 training and 1095 independent test days, SNR 0.1/3/10 dB, base error correlations 0/0.4/0.99, and signal–error correlation parameters 0/0.3. Accuracy summaries are conditional on available estimates; failure rates retain attempted denominators. Pairwise contrasts use common valid realizations. Confidence intervals quantify Monte Carlo uncertainty across realizations.

`NC.m` uses an off-diagonal L1 fit with diagonal bounds, followed by a PSD feasibility map by default. This is a two-stage numerical procedure, not a globally optimal joint PSD-constrained fit or proof of unique physical covariance identification. Raw NC (`enforcePSD=false`) is a diagnostic. `WA.m` applies explicit SPD conditioning consistently to compared covariance estimates. Numerical feasibility does not imply universal accuracy superiority. The supplied TCH implementation is the repository penalty-optimization variant, and its first reconstructed covariance output is used by the experiment scripts.

## License

The original MIT license is retained in `LICENSE`.
