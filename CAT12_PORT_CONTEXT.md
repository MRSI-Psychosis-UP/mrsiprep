# CAT12 Tissue Segmentation Port Context

Last updated: 2026-06-19

This file summarizes the working context for the MATLAB-free CAT12-like T1 tissue segmentation port in `mrsiprep`.

## Goal

Build a standalone Python implementation that approximates CAT12 T1 tissue segmentation into GM, WM, and CSF probability maps, without MATLAB runtime dependencies.

Preferred direction:

- Use Python-native preprocessing where practical.
- Use ANTs for registration/N4.
- Keep CSF outside the brain surface where CAT includes it; avoid HD-BET-style removal of outer CSF.
- Reuse low-level CAT12 C code where useful.
- Reimplement high-level MATLAB logic in Python.

## Important Paths

Repository:

```text
/home/flucchetti/Connectome/Dev/mrsiprep
```

CAT12 source:

```text
/home/flucchetti/Connectome/Dev/cat12
```

CAT12 standalone folder investigated earlier:

```text
/home/flucchetti/Connectome/Dev/cat-standalone-Linux/CAT_R2023b_MCR_Linux/cat_standalone/standalone
```

Sample T1:

```text
/home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_acq-mpragep3_run-01_T1w.nii.gz
```

CAT reference outputs for validation:

```text
/home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p0_T1w.nii.gz
/home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p1_T1w.nii.gz
/home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p2_T1w.nii.gz
/home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p3_T1w.nii.gz
```

## Committed Work

Recent commits:

```text
212e2fd ignre testing results
5cdc92e test(tissue): cover CAT12 AMAP and TPM scaffolding
89d474c feat(tissue): add standalone CAT12-like segmentation CLI
021df3e feat(tissue): add CAT12 TPM-seeded GMM initializer
73e3c55 feat(tissue): add Python wrapper for CAT12 AMAP core
e55da6e chore(cat12): vendor templates and AMAP sources
094fe94 feat(ants): add thread-aware registration facade
```

The generated result folders are intentionally not committed.

## Key Files Added Or Modified

ANTs support:

```text
mrsiprep/interfaces/ants.py
```

CAT12 data assets:

```text
mrsiprep/data/templates/cat12/TPM_Age11.5.nii.gz
mrsiprep/data/templates/cat12/T1.nii.gz
mrsiprep/data/templates/cat12/Template_T1.nii.gz
mrsiprep/data/templates/cat12/brainmask.nii.gz
mrsiprep/data/templates/cat12/brainmask_T1.nii.gz
mrsiprep/data/templates/cat12/cat.nii.gz
mrsiprep/data/templates/cat12/COPYING.CAT12
mrsiprep/data/templates/cat12/README.md
```

Vendored CAT12 AMAP C sources:

```text
mrsiprep/tissue/src/cat12_amap/
```

Python implementation:

```text
mrsiprep/tissue/cat12_amap.py
mrsiprep/tissue/tpm_gmm.py
mrsiprep/tissue/cat12_cleanup.py
mrsiprep/tissue/cat12_segment.py
```

Tests:

```text
tests/test_cat12_amap_scaffold.py
```

Package entry point:

```text
mrsiprep-cat12-tissue = "mrsiprep.tissue.cat12_segment:main"
```

## Current Pipeline

Main CLI:

```bash
mrsiprep-cat12-tissue T1.nii.gz OUTPUT_DIR --subject sub-... --session ses-... --use-amap
```

Core stages:

1. Register raw/native T1 to CAT template with ANTs.
2. Pull CAT TPM priors back into native T1 space.
3. N4 bias-correct the T1.
4. Fit six-class TPM-weighted GMM:
   - GM
   - WM
   - CSF
   - BONE
   - SOFT
   - BG
5. Use GM/WM/CSF maps to seed CAT AMAP.
6. Run vendored CAT12 AMAP C core through Python/ctypes.
7. Apply CAT-like cleanup and optional outer-CSF correction.

## Useful Commands

Default current pipeline on CHUVA009:

```bash
python -m mrsiprep.tissue.cat12_segment \
  /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_acq-mpragep3_run-01_T1w.nii.gz \
  temp_cat12_tpm_gmm_chuva009_catlike_preamap \
  --subject sub-CHUVA009 \
  --session ses-V5 \
  --use-amap \
  --overwrite \
  --ants-threads 22
```

CAT p-map seeded AMAP isolation test:

```bash
python -m mrsiprep.tissue.cat12_segment \
  /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_acq-mpragep3_run-01_T1w.nii.gz \
  temp_cat12_seeded_pmaps_chuva009_preamap \
  --subject sub-CHUVA009 \
  --session ses-V5 \
  --use-amap \
  --seed-p0 /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p0_T1w.nii.gz \
  --seed-p1 /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p1_T1w.nii.gz \
  --seed-p2 /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p2_T1w.nii.gz \
  --seed-p3 /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_desc-p3_T1w.nii.gz \
  --overwrite
```

Conservative seed-gating diagnostic:

```bash
python -m mrsiprep.tissue.cat12_segment \
  /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_acq-mpragep3_run-01_T1w.nii.gz \
  temp_cat12_tpm_gmm_chuva009_catlike_seedmin030 \
  --subject sub-CHUVA009 \
  --session ses-V5 \
  --use-amap \
  --amap-seed-min-probability 0.30 \
  --overwrite \
  --ants-threads 22
```

Six-class AMAP seed diagnostic:

```bash
python -m mrsiprep.tissue.cat12_segment \
  /home/flucchetti/Connectome/BIDS/Dummy-Project/sub-CHUVA009/ses-V5/anat/sub-CHUVA009_ses-V5_acq-mpragep3_run-01_T1w.nii.gz \
  temp_cat12_tpm_gmm_chuva009_catlike_sixclass_seed \
  --subject sub-CHUVA009 \
  --session ses-V5 \
  --use-amap \
  --use-six-class-amap-seed \
  --overwrite \
  --ants-threads 22
```

Run tests:

```bash
python -m unittest discover -s tests
```

Last passing result:

```text
Ran 22 tests: OK
```

## Validation Results

Default current pipeline:

```text
GM  MAE 0.01409  corr 0.93828  volume delta -28.17 ml
WM  MAE 0.00653  corr 0.96469  volume delta -20.23 ml
CSF MAE 0.01668  corr 0.74194  volume delta +147.89 ml
```

CAT p0/p1/p2/p3 seeded AMAP isolation:

```text
GM  MAE 0.00780  corr 0.97048  volume delta -37.65 ml
WM  MAE 0.00402  corr 0.98444  volume delta +22.12 ml
CSF MAE 0.00442  corr 0.93336  volume delta +12.28 ml
```

Conclusion from isolation:

- The vendored CAT AMAP C core and pre-AMAP prep can behave well.
- The main gap is upstream: the Python TPM/GMM initializer and generated AMAP seed labels are not yet CAT-like enough.

Conservative seed gate, `--amap-seed-min-probability 0.30`:

```text
CSF MAE 0.01374  corr 0.72021  volume delta +17.55 ml
```

This improves CSF volume and MAE but underestimates true CAT-CSF areas, so it is diagnostic/optional, not default.

Six-class CAT seed mapping:

```text
default CSF       MAE 0.01668  corr 0.74194  volume delta +147.89 ml
sixclass CSF      MAE 0.02510  corr 0.61287  volume delta +229.59 ml
```

Literal CAT six-class mapping is not safe yet because CAT maps bone to CSF-like label only inside a better CAT `Yb` support mask. Our current support mask is weaker.

## Important Technical Findings

CAT AMAP preparation in CAT source:

- `cat_main_amap.m` / `cat_main_amap1639.m` build `Yp0` from SPM/CAT class maps.
- Reordered max is effectively:

```text
Ycls order: [CSF, GM, WM, BONE, SOFT, BG, ...]
AMAP label mapping: [1, 2, 3, 1, 0, 0, ...]
```

- CAT then gates this through a brain/support mask `Yb`.
- CAT also prepares `Ymib`, applies CSF floor/noise, estimates auto-MRF, then calls `cat_amap`.

Current Python approximation:

- Implements CAT-like `Ymib/Yp0b`, CSF floor/noise, auto-MRF, cleanup, and AMAP call.
- Still does not fully replicate:
  - SPM unified segmentation.
  - CAT `Yb/Yb0`/gcut support construction.
  - CAT LAS/partvol details.
  - CAT nonlinear deformation/atlas cleanup behavior.

## Generated Output Folders

These are local validation/testing outputs and should remain uncommitted:

```text
temp_cat12_amap_chuva009/
temp_cat12_seeded_pmaps_chuva009_preamap/
temp_cat12_tpm_gmm_chuva009_*/
demo_fast_chuva009/
```

## Recommended Next Step

Next work should target CAT-like brain/support mask construction before more AMAP tuning.

Practical next implementation target:

1. Inspect/port more of CAT `Yb/Yb0` and gcut behavior.
2. Use the six-class TPM/GMM posteriors only after a better support mask exists.
3. Re-run the CAT-seeded isolation and default pipeline comparisons.

Do not spend much more time tuning AMAP cleanup until `Yb/Yb0` is closer to CAT.

