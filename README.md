# MRSIPrep

`MRSIPrep` is a preprocessing and derivative-generation pipeline for already
quantified whole-brain MRSI maps. It is derived from the implementation in
`MRSI-Metabolic-Connectome` and preserves the CHUV academic non-commercial
research license.

The package does not perform spectral fitting. It expects quantified MRSI maps,
quality maps, T1w images, skull-stripped T1w derivatives, and tissue probability
maps or a configured tissue-segmentation backend.

## Minimal command

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --metabolites CrPCr GluGln GPCPCh NAANAAG Ins \
  --registration-t1-target brain-csf \
  --tissue-backend existing \
  --parcellation-mode chimera \
  --chimera-scheme LFMIHIFIS \
  --chimera-scale 3
```

The default registration target is `brain-csf`, which adds CAT12 `p3` CSF voxels
to the skull-stripped T1w image before MRSI-to-T1 registration. If CAT12 `p3` is
missing, the current subject/session fails; in batch mode processing continues
with the next item.

## BIDS import utilities

```bash
mrsiprep-import /source/folder /path/to/bids --subject S001 --session V1
mrsiprep-skullstrip /path/to/bids --device cuda
mrsiprep-import-gui
```

The import helpers preserve the MRSI-Metabolic-Connectome derivative layout:
`derivatives/mrsi-orig`, `derivatives/cat12`, and `derivatives/skullstrip`.

## Attribution

Substantial implementation logic is cropped and refactored from
`MRSI-Metabolic-Connectome` by Federico Lucchetti and collaborators. The original
license is included in `LICENSE`.
