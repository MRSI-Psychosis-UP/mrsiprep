# Usage

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --registration-t1-target brain-csf \
  --tissue-backend synthseg-fast
```

`brain-csf` is the default registration target. With the default SynthSeg+FAST
backend, MRSIPrep runs SynthSeg, builds a FAST mask from HD-BET brain plus
SynthSeg CSF/ventricles, runs FAST on the masked T1, derives native T1-space
`p1/p2/p3` maps, and uses `p3` to create the CSF-extended T1w registration
target.

Check all selected subject/session inputs without running preprocessing:

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participants participants.tsv \
  --validate-only
```

This is intended for batch mode. It reports every invalid recording before any
expensive segmentation, registration, parcellation, or PVC step starts.

SynthSeg+FAST requires `mri_synthseg`, `hd-bet`, and FSL `fast`.

FreeSurfer remains available with `--tissue-backend freesurfer` and requires
`recon-all`, `mri_vol2vol`, and a valid `FS_LICENSE`.

Use precomputed CAT12 tissue maps with:

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --tissue-backend existing
```

The `existing` backend requires:

- a skull-stripped T1w derivative in `derivatives/skullstrip`,
- the raw BIDS T1w acquisition,
- CAT12-style `p3` CSF probability map in `derivatives/cat12`.

If `p3` is missing, the recording fails. In batch processing, MRSIPrep logs the
error and continues with the next subject/session.

ANTs Atropos remains available for testing with:

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --tissue-backend ants-atropos
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
```
