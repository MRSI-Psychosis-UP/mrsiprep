# Inputs

MRSIPrep reads the MRSI-Metabolic-Connectome derivative layout:

```text
derivatives/mrsi-orig/sub-<id>/ses-<id>/*_space-orig_met-<met>_desc-signal_mrsi.nii.gz
derivatives/mrsi-orig/sub-<id>/ses-<id>/*_desc-snr_mrsi.nii.gz
derivatives/mrsi-orig/sub-<id>/ses-<id>/*_desc-fwhm_mrsi.nii.gz
derivatives/mrsi-orig/sub-<id>/ses-<id>/*_met-<met>_desc-crlb_mrsi.nii.gz
derivatives/cat12/sub-<id>/ses-<id>/*_desc-p1_T1w.nii.gz
derivatives/cat12/sub-<id>/ses-<id>/*_desc-p2_T1w.nii.gz
derivatives/cat12/sub-<id>/ses-<id>/*_desc-p3_T1w.nii.gz
derivatives/skullstrip/sub-<id>/ses-<id>/*_desc-brain_T1w.nii.gz
```

CAT12 and skullstrip derivatives are required only for `--tissue-backend
existing`. The default `--tissue-backend freesurfer` requires raw BIDS T1w
images and a valid FreeSurfer installation/license, then writes the tissue and
skullstrip derivatives itself.
