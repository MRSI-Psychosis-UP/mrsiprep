CAT12 template assets vendored for the standalone tissue prototype.

Source:
`/home/flucchetti/Connectome/Dev/cat12/templates_MNI152NLin2009cAsym`

Files:
- `TPM_Age11.5.nii.gz`: six-class CAT12 TPM used as the first native
  prior source. Volumes are GM, WM, CSF, bone, non-brain soft tissue, and
  background.
- `T1.nii.gz`: 2 mm CAT12 MNI T1 template sharing the TPM orientation and
  physical space.
- `brainmask.nii.gz`: 2 mm mask matching `T1.nii.gz`.
- `Template_T1.nii.gz`: 1 mm CAT12 MNI T1 template, retained for future
  registration experiments.
- `brainmask_T1.nii.gz`: 1 mm mask matching `Template_T1.nii.gz`.

The files are copied from the CAT12 source tree and are distributed here
with the CAT12 license text in `COPYING.CAT12`.
