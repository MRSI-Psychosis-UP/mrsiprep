# Usage

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --registration-t1-target brain-csf \
  --tissue-backend freesurfer
```

`brain-csf` is the default registration target. With the default FreeSurfer
backend, MRSIPrep runs or reuses `recon-all`, writes FreeSurfer intermediates to
`<output>/freesurfer`, derives native T1-space `p1/p2/p3` maps, and uses `p3` to
create the CSF-extended T1w registration target.

FreeSurfer requires `recon-all`, `mri_vol2vol`, and a valid `FS_LICENSE`.

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
- CAT12 `p3` CSF probability map in `derivatives/cat12`.

If `p3` is missing, the recording fails. In batch processing, MRSIPrep logs the
error and continues with the next subject/session.

ANTs Atropos remains available for testing with:

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --tissue-backend ants-atropos
```
