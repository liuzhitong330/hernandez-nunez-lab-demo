# From a spatial candidate to a falsifiable next experiment

An exploratory application project by Cathy Liu, inspired by the Hernandez-Nunez Lab at Stanford. Not affiliated with or endorsed by Stanford, the target lab, or the data authors. Developed with coding assistance.

[Live demo](https://liuzhitong330.github.io/hernandez-nunez-lab-demo/) · [Paper](https://www.nature.com/articles/s41467-026-76242-8) · [Versioned data](https://doi.org/10.25378/janelia.31967640.v1) · [Target lab research](https://www.hernandez-nunez-lab.com/research)

## Decision and scientific contribution

Which recorded cell segments are worth checking in a new experiment, and which apparently strong responses fail when the trial, time trend, or reference condition changes?

The interface links real acquisition-space centroids to measured calcium traces, reserved-trial predictions, source IDs, and a frozen discovery ranking. It draws no invented pathways. Exported candidates are **computational review records**, not biologically validated intervention targets.

Source: Chen, James, Ruetten et al. (2026), *Whole-brain, all-optical interrogation of neuronal dynamics underlying gut and vascular interoception in zebrafish*, Nature Communications 17, 9881, DOI [10.1038/s41467-026-76242-8](https://doi.org/10.1038/s41467-026-76242-8). Weiyu Chen deposited the associated Janelia Figshare dataset under CC BY 4.0. These are **external gut/vascular-interoception recordings, not Hernandez-Nunez Lab heart recordings**.

The paper already cross-validates nutrient-versus-UV cell selection. Its multimodal comparison uses full-data fits. This project's contribution is a bounded held-out-event audit with drift and reference-order sensitivity, not invention of the original method or replication of the complete study.

## Interpretation gate: position codes are not verified anatomy

The deposit lacks recording-specific position-to-organ mappings. The authors' example notebook supplies conventions, but these do not verify every recording's anatomical labels.

| Recording | Selected code | Primary reference | Additional comparison | Visual modulation |
|---|---:|---|---|---|
| Gut D-glucose, fish 3 | 2 | Code 3, earlier block | Code 2 versus code 1, later block | Present, open-loop |
| Blood D-glucose, fish 1 | 7 | Code 1 | None | Absent at imaging samples |
| Blood D-glucose, fish 2 | 4 | Code 1 | Other codes censored | Absent at imaging samples |

Dataset names identify the authors' experimental groups, not each position code's anatomy. Code 7 in one fish is not assumed equivalent to code 4 in another. Reference codes are not called verified shams, and nonreference codes are not pooled into one biological site.

Internal compatibility names such as `gut_event_times`, `nutrient_event_times`, `nutrient_error_reduction_pct`, and model `UG` mean the **selected position-code regressor**, not demonstrated nutrient specificity. `control_event_times` means chosen reference-code events, not an independently validated biological control.

## Data and files

A response-independent systematic subset is selected before observing fluorescence: 32 equally source-index-spaced blocks of 32 consecutive rows, or 1,024 ROIs per recording. This is not a random or representative population sample. Three ROIs fail the fixed camera-baseline denominator check.

| Recording | Source ROIs × frames | Retained sampled ROIs | Analysis frames | Reserved events: selected/reference |
|---|---|---:|---:|---|
| Gut fish 3 | 135,111 × 2,220 | 1,022 | 2,220 | 2 / 2 |
| Blood fish 1 | 203,231 × 1,870 | 1,023 | 623 | 2 / 2 |
| Blood fish 2 | 177,421 × 3,090 | 1,024 | 1,030 | 2 / 2 |

An ROI is an author-segmented cell footprint, not an independently confirmed neuron. There are no population-prevalence claims, cell-level p-values, or pooled gut-versus-blood effects.

- `index.html`, `styles.css`, `app.js`: static article and linked interactions.
- `data.js`: available recording summaries; the gut recording is first because visual modulation is present.
- `derived/*.json`: per-ROI metrics, frozen ranks, events, methods, and limitations.
- `derived/*-traces-*.json.gz`: on-demand measured/predicted trace shards. Numerical scoring uses the analysis grid, not the reduced display points.
- `provenance.json`: DOI, pinned code commit, manifest links, and public-subset SHA-256 values.
- `repro/*_subset.npz`: compact real inputs sufficient to recompute all model results. Larger raw fluorescence, baseline, and footprint arrays are omitted; derived traces and baseline QC remain.
- `repro/*_manifest.json`: source IDs, dimensions, selected rows, transformations, hashes, stimulus timing, and mapping limitations.
- `repro/model_audit.py`, `repro/run_model_audit.py`: models and three-recording runner, including the alternate gut reference.
- `repro/test_model_audit.py`: deterministic tests. Synthetic signals are test fixtures only and never displayed as experimental data.
- `repro/extract_subset.py`, `repro/range_access.py`: optional original-data extraction.

The earlier WARP marker-screen version is superseded and preserved in Git history only. It is not part of the current analysis or reproduction workflow.

## Reproduce from bundled inputs

Use Python 3.9–3.12. Analysis was run with Python 3.9.6 and NumPy 2.0.2; optional extraction also used SciPy 1.13.1, h5py 3.14.0, requests 2.32.5, and matplotlib 3.9.4. Requirements record those versions.

From repository root:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r repro/requirements.txt
git clone https://github.com/benJames1212/gutBrainPipeline.git gutBrainPipeline
git -C gutBrainPipeline checkout 4c62612ecd540ec786c43fe9ce67cc78259a6ec1
cd repro
python -m unittest -v test_model_audit.py
python run_model_audit.py
```

The pinned upstream clone belongs at **repository root `/gutBrainPipeline`**, alongside `repro/`, not inside it. The source-equivalence test and extractor expect this location. The numerical runner itself needs only NumPy and the NPZ inputs; cloning is needed to test the cited helper and re-extract source data.

Generated outputs beside the runner:

```text
fish1-analysis.json
fish2-analysis.json
gut3-analysis.json
gut3-reference1-analysis.json
```

The primary gut JSON also attaches alternate-reference results by source cell ID, with its original discovery ranks unchanged. Outputs contain fuller arrays than the website's compact summaries. To regenerate the browser summaries and compressed trace shards after the runner finishes, return to repository root and run `python repro/package_demo.py --output .`. Keep the pinned upstream clone present for its license file. Packaging records fresh hashes for the compact NPZs; ZIP-container timestamps can change their byte hashes without changing numerical arrays.

Nine tests check official lag-helper equivalence, intact-timeline lags, held-out outcome isolation, modality-boundary penalties, event censoring, binning, fallback labeling, unclassified-position exclusion, and a polynomial-only negative control. Passing tests does not validate a biological claim.

## Methods

### Signals and timing

Author-processed fluorescence is converted as `(cell_timeseries-cell_baseline)/(cell_baseline-100)`, following the source camera-background convention. ROIs with any baseline at or below 100 are excluded. This project does not rerun image registration or segmentation.

Imaging timestamps match recorded light-sheet falling edges; counts exactly match HDF5 frame counts. Motor input is the authors' Gaussian-smoothed sliding variance of both complete 6 kHz fictive-swim voltage channels, averaged and sampled at imaging times. It is voltage variance, not independently classified swim bouts. Low-amplitude fluctuations, particularly in blood fish 1, may be noise.

Gut fish 3 has 18 analog UV trains matching command counts and position codes. Recorded UV starts 1.410–1.503 seconds after its command; analysis uses actual recorded UV onset. Blood recordings do not cross the authors' analog detection threshold, so their times use documented digital UV-trial commands. These are **not independently validated optical-delivery times**. Neither channel provides a calibrated photon-dose measurement.

Visual velocity is binarized as `visual > 0`. Faster blood recordings are mean-binned to approximately 1 Hz; the slower gut recording is not binned for modeling. Binned visual values represent the fraction of source samples with positive visual velocity. Constant visual regressors are omitted. No visual modulation in blood is a limitation, not a negative visual-response finding.

### Models and protected split

Three smoothness-regularized models share an intercept and cubic session-time nuisance:

- `UVM`: all-UV timing + motor + available visual input.
- `UG`: all-UV timing + selected position-code timing.
- `FULL`: all-UV timing + selected position-code timing + motor + available visual input.

Fixed causal horizons are 55 seconds for both UV terms, 10 seconds for motor, and 5 seconds for visual input. Lags are constructed on the **intact timeline before masking**; deleted intervals are never spliced. Lag blocks use training-row RMS scaling. Time is centered/scaled on training rows. The penalty is `2I + 20 D.T @ D` within separate modality blocks; intercept and common time-polynomial terms are unpenalized. Degree and penalties are not tuned on held-out responses.

The final two eligible events per selected/reference code are reserved. Training excludes a 55-second buffer before and after each reserved response window. Scoring begins two seconds after onset, stopping at the next UV event, 55 seconds, or recording end. At least 20 scored seconds are required. Unclassified-code response intervals are excluded for 55 seconds; overlapping scored windows are not used. This prevents overlapping blood events, commonly only 40 seconds apart, from being counted as separate full 55-second responses.

Ranking uses **drift-adjusted training-only error reduction**. Entire reserved scoring windows, plus their buffers, are omitted from fitting; no held-out neural timepoints or alternate-reference test outcomes determine the rank. The protocol was revised after diagnostic review, including timing verification and drift-only tests. This is a retrospective analysis with protected numerical splits, not a prospectively registered experiment. Partial training windows can fit the fixed ridge model when full windows are scarce, but support is disclosed:

| Primary comparison | Complete discovery windows: selected/reference | Available discovery-response seconds: selected/reference |
|---|---|---|
| Gut fish 3, code 2 vs 3 | 2 / 3 | 117.303 / 163.106 |
| Blood fish 1, code 7 vs 1 | 1 / 3 | 59.583 / 114.400 |
| Blood fish 2, code 4 vs 1 | 0 / 1 | 48.858 / 59.583 |

Blood fish 2 has especially weak support: partial selected-code responses, not complete discovery trials. “Independent holdout” means test outcomes are protected from fitting/ranking, not independent animals, prospective deployment, or biological validation. A labeled leave-event-out exploratory fallback exists in the general module; all three packaged primary results use the reserved split.

### Metrics and sensitivity

Per reserved selected-position event, error reduction is `100 × (comparison MSE − full-model MSE) / comparison MSE`, averaged across events. Positive means lower held-out squared error; negative exposes failure. Percentages can be extremely negative when comparison MSE is small. The displayed correlation is also an arithmetic mean of per-trial Pearson correlations, not a pooled whole-recording correlation. It accompanies error because improving amplitude error does not ensure matching waveform shape.

`original_no_drift` repeats the same split without common time-polynomial terms. This sensitivity check was added after a blocked-design negative-control audit; it was not prospectively registered. Cubic adjustment does not eliminate every drift pattern or the structural confounding of site blocks with session time.

The gut analysis also contrasts selected code 2 with later reference code 1, separately from earlier code 3. These are different comparisons, not interchangeable controls. Primary ranks stay frozen. Inspecting multiple cells, budgets, or comparisons after seeing reserved outcomes remains exploratory.

## Checked examples and follow-up

Gut ROI **8725**, the top discovery-ranked candidate, reduces primary held-out error by 85.54%, but has **−24.77%** reduction against the later reference block. Its apparent promise is reference-dependent. ROI **15**, discovery rank 4, gives 41.48% and 9.20%, with full-model correlations 0.544 and 0.523. That consistency motivates further inspection, not a claim of chemical specificity or causal control.

A useful next experiment would verify position-to-organ labels, interleave selected/reference conditions, vary measured motor/visual conditions independently, and repeat the audit on new recordings. Molecular characterization, anatomical confirmation, reporter-specificity checks, and appropriate perturbation/off-target controls are still required before proposing a transgenic target.

## Optional original-data extraction

This is unnecessary for bundled-input model reproduction. It transfers approximately **2.88 GB**: complete ephys files totaling 2,794,152,400 bytes, about 83.14 MB of selected imaging ranges, and small headers/metadata. Full-band variance processing also needs additional memory and cache disk space.

After installing dependencies and cloning the pinned author repository:

```sh
cd repro
curl -L --fail --retry 3 https://ndownloader.figshare.com/files/63568599 -o blood_glucose_fish1_eP.26chFlt-v10
curl -L --fail --retry 3 https://ndownloader.figshare.com/files/63568617 -o blood_glucose_fish2_eP.26chFlt-v10
curl -L --fail --retry 3 https://ndownloader.figshare.com/files/63568746 -o glucose_fish3_eP.26chFlt-v10
python extract_subset.py blood_glucose_fish1
python extract_subset.py blood_glucose_fish2
python extract_subset.py glucose_fish3
python run_model_audit.py
```

The supplied `figshare-31967640.json` metadata must be beside the extractor. Ephys byte lengths and supplied MD5 values are checked. Imaging files are **not downloaded in full**: HTTP 206 requests retrieve selected fluorescence, baseline, and footprint arrays. Status, Content-Range, and byte length are checked, but whole imaging-file MD5 is not. Extracted full NPZs include extra arrays, so their hashes differ from compact public NPZs; consult the appropriate manifest field.

Extraction actually executes the authors' `load`, `pullCenters`, and `windowed_variance`. A compatibility alias supplies the older Gaussian-function import expected by unchanged upstream code. Lag construction is adapted and equivalence-tested against the pinned helper. The explicit solver, time-trend controls, protected masking, and sensitivity checks are documented extensions—not execution of the complete original pipeline.

## Limits and reuse

- Centroids are unweighted means of nonnegative footprint voxels, following the authors' helper. Physical voxel spacing and cross-fish atlas registration are not established. Do not label axes micrometres. Display z is expanded fourfold for legibility; geometry is not a metric anatomical reconstruction.
- The three recordings differ in route, unknown anatomical codes, schedules, and covariates. They are not equivalent biological replications, and have no established same-cell correspondence.
- Movement can mediate responses as well as confound them. Conditional predictability is not pure interoception, connectivity, organ control, or a stimulation prescription.
- Spatial location supplies no molecular identity. This artifact does not demonstrate Cathy's hands-on zebrafish transgenesis or spatial-transcriptomics expertise.
- Reusable inputs are imaging time, cell-by-time traces, source IDs/centroids, verified stimulus labels, and aligned motor/visual records. New datasets need acquisition metadata, preprocessing/QC review, and independent validation. Do not silently transfer these code conventions.

## Serve

From repository root:

```sh
python3 -m http.server 8000
```

Open `http://localhost:8000` in a modern browser. Use HTTP, not direct HTML-file opening, because JSON and gzip traces are fetched. There is no backend, account, analytics, camera access, or secret API key. Public contact details are intentionally included with Cathy's permission.

## Attribution and licensing

New code is MIT licensed. Adaptations retain Ben James's upstream MIT attribution (2025); the full notice is also in `repro/UPSTREAM_LICENSE`. Input subsets, derived summaries/traces, and provenance remain **CC BY 4.0**, attributed to Weiyu Chen and Chen, James, Ruetten et al. Changes include systematic subsetting, centroids, signal conversion, timing alignment, model fitting, sensitivity analysis, and presentation. See [LICENSE](LICENSE), [provenance](provenance.json), and recording manifests. No license grants institutional endorsement.
