# Usage

All commands below run through Docker; there is no supported host
installation of the `mrsiprep` CLI.

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --mode light \
  --output-spaces T1w MNI152NLin2009cAsym \
  --synthseg-mode fast \
  --nthreads 8
```

Light mode registers the imported MRSI signal maps to a SynthSeg-extracted T1w image, writes
the requested T1w/MNI maps without filtering or PVC, and runs SynthSeg with
cortical parcellation. Its parcel QC table reports
the percentage of each anatomical T1w parcel covered by MRSI, parcelwise CRLB,
and valid-voxel fractions. Light mode does not run FAST, PETPVC,
Chimera, or `recon-all`.

SynthSeg-based brain extraction always retains the whole brain (GM, WM,
ventricles, and inner/outer CSF, including extra-ventricular CSF label 24) —
only SynthSeg background (label 0) is excluded from the brain mask, so FAST
sees the complete CSF compartment when estimating tissue probabilities.

Use full mode for tissue PVC and Chimera parcellation:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --mode full \
  --tissue-backend synthseg-fast \
  --parcellation-mode chimera \
  --chimera-scheme LFMIHIFIS \
  --chimera-scale 3
```

Full mode also writes a legacy-compatible parcel profile archive under
`<out>/mrsi_parcel/sub-*/ses-*/mrsi/*desc-metprofiles_mrsi.npz`.

Use a bundled MNI atlas instead of Chimera with:

```bash
--mode full --parcellation-mode mni --atlas chimera-LFMIHIFIS-3
```

Check all selected subject/session inputs without running preprocessing:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participants participants.tsv \
  --validate-only
```

This is intended for batch mode. It reports every invalid recording before any
expensive segmentation, registration, parcellation, or PVC step starts.

## Verbosity, logging, and parallel processing

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participants participants.tsv \
  --mode full \
  --verbose 1 \
  --nthreads 8 \
  --nproc 4
```

`--verbose` (`0`-`3`, default `1`):

- `0` — only the start/finish line and elapsed time per subject/session.
- `1` — also prints each processing step as it starts (tissue segmentation,
  anatomical prep, MRSI preprocessing, registration, tissue maps, PVC,
  resampling, parcellation, regional extraction, connectivity, reports), with
  no per-step detail.
- `2` — also prints step-level detail (info/success/warning/error messages).
- `3` — also lets ANTs, `recon-all`, and `mri_synthseg` print their own raw
  subprocess output instead of being captured.

A full-detail DEBUG log is always written to
`<out>/mrsiprep/logs/mrsiprep_<timestamp>.log`, independent of the console
`--verbose` level.

`--nproc` runs that many subject/session recordings concurrently; each one
gets `--nthreads` ANTs/ITK threads. MRSIPrep coerces `--nthreads` down (never
`--nproc`) if `nproc * nthreads` would exceed the host's CPU count, and shows
the resulting thread budget (or the coercion warning) in the preflight
summary before any recordings are processed.

Light mode requires `mri_synthseg` and ANTs. Full mode with the default
`synthseg-fast` tissue backend additionally requires FSL `fast`, PETPVC, and
(for Chimera parcellation) `recon-all` and a valid `FS_LICENSE`.

Use precomputed CAT12 tissue maps with:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --mode full \
  --tissue-backend existing
```

The `existing` backend requires:

- a skull-stripped T1w derivative in `derivatives/skullstrip`,
- the raw BIDS T1w acquisition,
- CAT12-style `p3` CSF probability map in `derivatives/cat12`.

If `p3` is missing, the recording fails. In batch processing, MRSIPrep logs the
error and continues with the next subject/session.

Disable tissue segmentation and PVC entirely with `--tissue-backend none`
(equivalent to also passing `--no-pvc`):

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/derivatives:/out \
  fedlucchetti/mrsiprep:cpu \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --mode full \
  --tissue-backend none
```

## Output layout

```text
<out>/mrsiprep/sub-*/ses-*/mrsi-orig/      native/imported-grid MRSI signal maps
<out>/mrsiprep/sub-*/ses-*/mrsi-orig-pvc/  PVC-corrected native-grid maps
<out>/mrsiprep/sub-*/ses-*/mrsi-t1w/       T1w-aligned MRSI maps
<out>/mrsiprep/sub-*/ses-*/mrsi-mni/       MNI-normalized MRSI maps
<out>/mrsiprep/sub-*/ses-*/tissue-mrsi/    MRSI-grid tissue probsegs and 4Dtissue
<out>/mrsiprep/sub-*/ses-*/qmasks/         QC, spike, and brain masks
<out>/mrsiprep/sub-*/ses-*/anat/           T1w tissue and registration files
<out>/chimera-atlases/sub-*/ses-*/anat/    Chimera atlas outputs
<out>/mrsiprep/sub-*/ses-*/parcellations/  regional metabolite tables and MNI atlases
<out>/mrsiprep/sub-*/ses-*/connectomics/   matrices, nodes, and edges
<out>/mrsi_parcel/sub-*/ses-*/mrsi/        full-mode metabolite profile NPZ files
<out>/mrsiprep/logs/                       full-detail timestamped run logs (independent of --verbose)
```
