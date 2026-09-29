# Molecular markers before circuit targeting

Personal exploratory application demo by Cathy Liu for Luis Hernandez-Nunez's Stanford lab. Not affiliated with or endorsed by either source lab.

Live: https://liuzhitong330.github.io/hernandez-nunez-lab-demo/

## Scientific decision

Choose chata co-markers for a first spatial-validation panel, trading anatomical concentration against coverage and cross-animal stability. No cardiac cell identity, causal result, transgenic driver fidelity or stimulation recommendation is inferred.

## Data and reproduction

Source: Marquez-Legorreta, Fleishman, Hesselink et al., WARP (2026), DOI 10.64898/2026.02.07.704095; publicly linked Figshare dataset https://figshare.com/s/d1d19b105c4f74865c32 (CC BY 4.0). This external reference is **not** the Hernandez-Nunez heart-brain dataset.

- `cells.npz`: derived granular data, one row per retained source cell, with binary calls for all 41 genes, fish, focal region and original row index; CC BY 4.0 with attribution to the WARP authors. Coordinate values are not redistributed, but their finiteness determines inclusion.
- `provenance.json`: original four download URLs, SHA256 checksums, inclusion counts and pinned reference mapping.
- `analyze.py`: independent preprocessing and exact aggregation. Python >=3.10 + NumPy. Run `python analyze.py` to rebuild `data.js` and `marker_counts.csv` from the bundled cell-level subset. Run `python analyze.py --download` to reacquire and hash-check the original arrays (~110 MB), normalize them and rebuild. Downloads fail closed on checksum changes.
- `marker_counts.csv`: complete co-marker × fish × focal-region counts and denominators.
- `app.js`: interaction, transparent ranking, leave-one-fish-out stress test and CSV export. No model fitting or invented data.

JaneliaSciComp/warp-atlas-viewer was fetched at commit `0f9451756b0a0bf3d051fea0f40789992a17e783` to interpret its documented gene order and focal-region mapping. That repository had no declared code license at inspection; its implementation is not copied. Gene labels and region meanings are scientific metadata; the analysis is independently implemented. The target paper links BrainBodyInteractionsLab, which had no public repositories at the September 29 check.

## Checked example

For the pooled medulla labels (1,2,3), chata+esrrb has in-area / all-pair-positive counts of 152/168, 115/139 and 87/106 for fish 59, 63 and 71. In-area chata-only denominators are 1243, 1129 and 387. Thus anatomical concentration ranges 82.1–90.5%, but coverage is only 10.2–22.5%. A more concentrated panel misses most chata+ cells. This is an anatomical reference result, not a cardiac-neuron specificity estimate.

## Ranking and limitations

Require each training fish to meet user-selected minimum in-area count and coverage. Rank eligible co-markers by their lowest training-fish concentration, then gene name. Exclude the chata-only baseline from the shortlist. For each leave-one-fish-out fold, select using the other two only and apply that fixed choice to the held-out fish. No information from the held-out fish enters that fold's selection. Report the held-out fractions and whether its floors pass; do not claim external validation.

Cells are nested within just three fish; no cell-level p-values or population confidence intervals. Author calls, segmentation, atlas registration, gene-positive preselection, the limited marker panel and developmental/study differences all constrain interpretation. Focal-region labels are coarse and do not isolate medial cardiac premotor cells. All cells outside the selected region, including unassigned ones, remain in the concentration denominator. An empty denominator yields “not estimable,” not zero. Activity-related cfos and potential glial cx43 are explicitly flagged.

Candidate panel requires independent spatial co-labeling, negative and single-marker controls, then physiological identification and reporter validation. No claim of Cathy's prior hands-on spatial transcriptomics or zebrafish transgenesis. Prepared with coding assistance.

## Serve

`python -m http.server 8000` then open http://localhost:8000. Static files only, no account, tracking, backend or API credentials.

New analysis/UI code: MIT license (see LICENSE). Derived data remain CC BY 4.0 as attributed above.
