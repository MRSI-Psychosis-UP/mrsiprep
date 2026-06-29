# *MRSIPrep*: A Robust Preprocessing Pipeline for Whole-Brain MRSI Data

*MRSIPrep* is a preprocessing and derivative-generation pipeline for already
quantified whole-brain MRSI maps, run as a BIDS App via Docker.

[![Docker Pulls](https://img.shields.io/docker/pulls/fedlucchetti/mrsiprep)](https://hub.docker.com/r/fedlucchetti/mrsiprep)
[![Documentation Status](https://app.readthedocs.org/projects/mrsiprep/badge/?version=latest)](https://mrsiprep.readthedocs.io/en/latest/)
[![License: CHUV academic non-commercial](https://img.shields.io/badge/license-academic--non--commercial-blue)](https://github.com/MRSI-Psychosis-UP/mrsiprep/blob/main/LICENSE)

## About

`MRSIPrep` does not perform spectral fitting. It expects quantified MRSI maps,
quality maps, and T1w images as input. Its default light mode normalizes MRSI
maps and uses fast SynthSeg cortical parcellation for parcelwise anatomical
coverage and CRLB reporting. Full mode adds SynthSeg+FAST tissue maps,
PETPVC, and Chimera/MNI-atlas regional profile extraction. It is derived from
the implementation in `MRSI-Metabolic-Connectome` and preserves the CHUV
academic non-commercial research license.

## Pulling the published image

A prebuilt CPU image is published on Docker Hub:

```bash
docker pull fedlucchetti/mrsiprep:cpu
```

It bundles ANTs, FSL (FAST only), FreeSurfer (`recon-all`, `mri_synthseg`,
`mri_vol2vol`), PETPVC, and Chimera. It does not include a FreeSurfer license
file — mount your own and set `FS_LICENSE` as shown below. You will still
need to provide a BIDS dataset with already-quantified MRSI maps; see
"Minimal command" and "Docker / BIDS App usage" below for invocation.

## Minimal command

All invocations run through Docker — there is no supported host installation
of the `mrsiprep` CLI. Mount your BIDS dataset and derivatives directory and
run the published or locally built image:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --validate-only
```

Use `--validate-only` to check all selected subject/session inputs before
starting an expensive batch run. It reports invalid recordings and exits without
running SynthSeg, FAST, registration, parcellation, or PVC.

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --metabolites CrPCr GluGln GPCPCh NAANAAG Ins \
  --mode light \
  --output-spaces T1w MNI152NLin2009cAsym \
  --nthreads 16
```

`--nthreads` controls ANTs registration/transform-application thread count
per subject/session (default 16). `--mni-resolution` selects the MNI
template resolution used for both T1→MNI registration and final resampling:
`origres` (MRSI native resolution), `t1wres` (T1w resolution, default), or an
explicit `<N>mm`. CRLB, SNR, and FWHM(linewidth) maps for the configured
`--metabolites` are transformed into T1w/MNI space alongside the signal maps
whenever present; per-metabolite spike masks are only transformed if
`--transform-spikemask` is passed (the combined QC mask stays in MRSI
space).

### Verbosity, logging, and parallel subjects

`--verbose {0,1,2,3}` (default 1) controls console output:

- `0` — only the start/finish line and elapsed time for each subject/session.
- `1` — also announces each processing step as it starts (segmentation,
  registration, PVC, resampling, parcellation, connectivity, reports), without
  step-level detail.
- `2` — also shows step-level detail (info/success/warning/error messages).
- `3` — also lets ANTs, `recon-all`, and `mri_synthseg` print their own raw
  output instead of being captured.

Regardless of `--verbose`, full DEBUG-level detail is always written to a
timestamped log file under `<out>/mrsiprep/logs/`.

`--nproc` (default 1) processes that many subject/session recordings in
parallel, each one getting `--nthreads` threads. If `nproc * nthreads`
exceeds the machine's CPU count, `--nthreads` is automatically reduced (never
`nproc`) and a warning is shown in the preflight summary before processing
starts — e.g. on a 32-core machine, `--nproc 4 --nthreads 10` (40 threads)
is coerced down to `--nthreads 8` (32 threads).

Light mode uses `mri_synthseg --parc --robust --threads <nthreads>` for both
extraction and parcellation, registers MRSI to the extracted T1w, skips
tissue PVC, and does not invoke Chimera or `recon-all`. Use `--synthseg-mode
fast` for faster, lower-accuracy extraction; `--fast` and `--robust` are
never combined. The extraction always retains the whole brain (GM, WM,
ventricles, and inner/outer CSF) — only SynthSeg background (label 0) is
excluded from the brain mask.

Full processing is selected explicitly:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 --session-label V1 \
  --mode full \
  --tissue-backend synthseg-fast \
  --parcellation-mode chimera \
  --chimera-scheme LFMIHIFIS --chimera-scale 3
```

Use `--parcellation-mode mni --atlas chimera-LFMIHIFIS-3` for a bundled MNI
atlas. Full mode writes regional TSV output and a legacy-compatible
`desc-metprofiles_mrsi.npz` under `<out>/mrsi_parcel`.

## Command-line argument reference

`docker run --rm -v ... -v ... <image> bids_dir output_dir participant [options]`
— `bids_dir`/`output_dir` are the container-internal mount paths (e.g. `/data`
and `/out`), followed by the same options as below.

### Subjects, sessions, and acquisitions

| Argument | Default | Description |
| --- | --- | --- |
| `--participant-label` | (all) | One or more subject labels to process, e.g. `S001 S002`. |
| `--session-label` | (all) | One or more session labels to process, e.g. `V1 V2`. |
| `--participants` | none | Path to a TSV/CSV with `subject`/`session` columns for batch runs, as an alternative to `--participant-label`/`--session-label`. |
| `--b0` | `3.0` | Field strength: `3.0` or `7.0`. Selects the default metabolite list when `--metabolites` is not given. |
| `--metabolites` | (B0-dependent list) | Metabolite names to process, e.g. `CrPCr GluGln GPCPCh NAANAAG Ins`. |
| `--t1` | `desc-brain_T1w` | BIDS filename pattern used to locate the input T1w image. |

### Quality thresholds

| Argument | Default | Description |
| --- | --- | --- |
| `--quality-metrics` | `snr linewidth crlb` | Which quality maps to check/report. |
| `--snr-min` | (package default) | Minimum acceptable SNR. |
| `--linewidth-max` | (package default) | Maximum acceptable linewidth/FWHM. |
| `--crlb-max` | (package default) | Maximum acceptable per-metabolite CRLB. |

### Processing mode and tissue segmentation

| Argument | Choices / Default | Description |
| --- | --- | --- |
| `--mode` (alias `--processing-mode`) | `light`, `full` / `light` | `light`: SynthSeg extraction + parcellation only, no FAST/PVC/Chimera/`recon-all`. `full`: adds tissue probability maps, PVC, and Chimera/MNI-atlas parcellation. |
| `--tissue-backend` | `synthseg-fast`, `existing`, `none` / `synthseg-fast` | How GM/WM/CSF probability maps are produced in full mode. `synthseg-fast` segments with SynthSeg+FAST; `existing` reuses precomputed CAT12 maps from `derivatives/skullstrip`/`derivatives/cat12`; `none` disables tissue segmentation and PVC entirely. |
| `--synthseg-mode` | `fast`, `standard`, `robust` / `robust` | SynthSeg accuracy/speed trade-off; `fast` and `robust` are never combined. SynthSeg thread count is taken from `--nthreads`. |
| `--registration-t1-target` | `brain-csf`, `brain`, `raw` / `brain-csf` (full mode), `brain` (light mode) | Which T1w variant MRSI is registered to. |
| `--csf-pv-threshold` | `0.95` | CSF partial-volume threshold used when building registration masks. |
| `--fs-subjects-dir` | none | Existing FreeSurfer `SUBJECTS_DIR` to reuse/write into (used by Chimera parcellation). |

### Registration and normalization

| Argument | Choices / Default | Description |
| --- | --- | --- |
| `--registration-backend` | `ants` / `ants` | Registration engine (currently only ANTs). |
| `--normalization` | `simple`, `ants-syn`, `existing` / `simple` | T1→MNI normalization strategy; `existing` reuses a precomputed transform instead of registering. |
| `--output-spaces` | `T1w MNI152NLin2009cAsym` | Which space(s) to resample final MRSI maps into: `MRSI`, `T1w`, `MNI152NLin2009cAsym` (aliases `mrsi`, `t1`, `mni` accepted). |
| `--mni-resolution` | `origres`, `t1wres`, `<N>mm` / `t1wres` | MNI template resolution for both T1→MNI registration and final resampling. |
| `--ref-met` | `CrPCr` | Reference metabolite map used to build the MRSI registration target. |
| `--transform-spikemask` | off | Also transform per-metabolite spike masks into T1w/MNI space (the combined QC mask is never transformed). |
| `--transform` | `""` | Legacy output-transform override; prefer `--output-spaces`. |

### Filtering

| Argument | Default | Description |
| --- | --- | --- |
| `--no-filter` | off | Disable biharmonic spike-repair filtering of MRSI maps. |
| `--spikepc` | `99.0` | Percentile threshold used for spike detection. |
| `--no-pvc` | off | Disable partial-volume correction (full mode only; always disabled in light mode). |

### Parcellation and connectivity

| Argument | Choices / Default | Description |
| --- | --- | --- |
| `--parcellation-mode` | `synthseg`, `chimera`, `mni` / mode-dependent | Light mode only supports `synthseg`; full mode requires `chimera` or `mni`. |
| `--atlas` | `chimera-LFMIHIFIS-3` | Bundled MNI atlas name (used with `--parcellation-mode mni`). |
| `--custom-atlas` | none | Path to a custom atlas NIfTI. |
| `--custom-atlas-lut` | none | Path to the lookup table for `--custom-atlas`. |
| `--chimera-scheme` | `LFMIHIFIS` | Chimera parcellation scheme. |
| `--chimera-scale` | `3` | Chimera parcellation scale. |
| `--chimera-grow` | `2` | Chimera region-growing parameter. |
| `--write-connectivity` | off | Write a connectivity matrix from regional metabolite values. |
| `--connectivity-method` | `pearson`, `spearman`, `cosine`, `euclidean_distance` / `spearman` | Similarity/distance metric for the connectivity matrix. |
| `--connectivity-space` | `MRSI`, `T1w`, `MNI` / `MRSI` | Space the connectivity matrix is computed in. |
| `--connectivity-n-perturbations` | `50` | Number of CRLB-scaled noise perturbations per metabolite used to build the connectivity matrix. |
| `--connectivity-sigma-scale` | `2.0` | Scale factor applied to the CRLB-derived noise sigma when perturbing metabolite maps for connectivity. |
| `--regional-summary` | `mean`, `median`, `weighted_mean` / `mean` | How voxel values are summarized per parcel. |

### Performance, logging, and control flow

| Argument | Default | Description |
| --- | --- | --- |
| `--nthreads` | `16` | ANTs/ITK thread count per subject/session process. |
| `--nproc` | `1` | Number of subject/session recordings to process in parallel; combined with `--nthreads` this is capped at the host's CPU count (see [Verbosity, logging, and parallel subjects](#verbosity-logging-and-parallel-subjects)). |
| `--verbose`, `-v` | `0`-`3` / `1` | Console output detail level (see below). |
| `--work-dir` | `<output_dir>/work` | Scratch directory for intermediate files. |
| `--validate-only` | off | Check selected subject/session inputs and exit without processing. |
| `--check-external-libs` | off | Verify required external binaries are on `PATH`/installed and exit. |
| `--overwrite` | off | Force-rerun all steps, ignoring cached outputs. |
| `--overwrite-filt` | off | Force-rerun MRSI filtering only. |
| `--overwrite-seg` | off | Force-rerun tissue segmentation (SynthSeg brain extraction + dseg/probseg) only. |
| `--overwrite-pve` | off | Force-rerun tissue probability map generation only. |
| `--overwrite-t1-reg` | off | Force-rerun MRSI→T1w registration only. |
| `--overwrite-mni-reg` | off | Force-rerun T1w→MNI registration only. |
| `--overwrite-transform` | off | Force-rerun transform resampling only. |

## Docker / BIDS App usage

Run it like fMRIPrep, mounting your BIDS dataset, derivatives directory, and
a FreeSurfer license file:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/bids/derivatives:/out \
  -v /path/to/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --mode light \
  --nthreads 8
```

The PyQt import GUI is not included in the container and should be run on the
host. MRSIPrep writes native T1-space derivatives including
`space-T1w_label-GM/WM/CSF_probseg` GM/WM/CSF probability maps, the raw
SynthSeg parcellation (`desc-synthsegParcFast/Robust_dseg`), and the masked
T1w actually fed into FAST (`desc-synthsegFastInput_T1w`) — all under
`anat/`, alongside ANTs transforms under `transforms/`. Only the raw
pre-resampled `mri_synthseg` output and FAST's internal scratch files stay
in the configured work directory.

MRSIPrep outputs are grouped by processing space:

```text
<out>/mrsiprep/sub-*/ses-*/mrsi-orig/      native/imported-grid MRSI signal maps
<out>/mrsiprep/sub-*/ses-*/mrsi-orig-pvc/  PVC-corrected native-grid maps
<out>/mrsiprep/sub-*/ses-*/mrsi-t1w/       T1w-aligned MRSI maps
<out>/mrsiprep/sub-*/ses-*/mrsi-mni/       MNI-normalized MRSI maps
<out>/mrsiprep/sub-*/ses-*/tissue-mrsi/    MRSI-grid tissue probsegs and 4Dtissue
<out>/mrsiprep/sub-*/ses-*/qmasks/         QC, spike, and brain masks
<out>/mrsiprep/sub-*/ses-*/anat/           T1w tissue, SynthSeg, and registration files
<out>/mrsiprep/sub-*/ses-*/transforms/     ANTs MRSI→T1w and T1w→MNI transforms
<out>/chimera-atlases/sub-*/ses-*/anat/    Chimera atlas outputs
<out>/mrsiprep/sub-*/ses-*/connectomics/   matrices, nodes, and edges
<out>/mrsi_parcel/sub-*/ses-*/mrsi/        full-mode metabolite profile NPZ files
<out>/mrsiprep/logs/                       full-detail timestamped run logs (independent of --verbose)
```

## BIDS import utilities

```bash
mrsiprep-import /source/folder /path/to/bids --subject S001 --session V1
mrsiprep-skullstrip /path/to/bids --device cpu
mrsiprep-import-gui
```

The import helpers preserve the MRSI-Metabolic-Connectome derivative layout:
`derivatives/mrsi-orig`, `derivatives/cat12`, and `derivatives/skullstrip`.

## Attribution

Substantial implementation logic is cropped and refactored by Federico Lucchetti and Edgar Céléreau. The original
license is included in `LICENSE`.
