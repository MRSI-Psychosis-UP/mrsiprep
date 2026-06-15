# MRSIPrep

`MRSIPrep` is a preprocessing and derivative-generation pipeline for already
quantified whole-brain MRSI maps. It is derived from the implementation in
`MRSI-Metabolic-Connectome` and preserves the CHUV academic non-commercial
research license.

The package does not perform spectral fitting. It expects quantified MRSI maps,
quality maps, and T1w images. By default it runs FreeSurfer `recon-all` to
create the skull-stripped T1w image and GM/WM/CSF tissue maps needed downstream.

## Minimal command

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --metabolites CrPCr GluGln GPCPCh NAANAAG Ins \
  --registration-t1-target brain-csf \
  --tissue-backend freesurfer \
  --parcellation-mode chimera \
  --chimera-scheme LFMIHIFIS \
  --chimera-scale 3
```

The default registration target is `brain-csf`, which adds the `p3` CSF layer to
the skull-stripped T1w image before MRSI-to-T1 registration. With the default
FreeSurfer backend, `p1/p2/p3` are derived from `recon-all` outputs. With
`--tissue-backend existing`, precomputed CAT12 `p1/p2/p3` maps are required; if
they are missing, the current subject/session fails and batch processing
continues with the next item.

## Docker / BIDS App usage

Build the local image:

```bash
docker build -t mrsiprep:latest .
```

Run it like fMRIPrep:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/bids/derivatives:/out \
  -v /path/to/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  mrsiprep:latest \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --tissue-backend freesurfer \
  --registration-t1-target brain-csf
```

For the dummy dataset on this workstation:

```bash
docker run --rm \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project:/data:ro \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project/derivatives:/out \
  -v "${FS_LICENSE}:/opt/freesurfer/license.txt:ro" \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  mrsiprep:latest \
  /data /out participant \
  --participant-label CHUVUP013 \
  --session-label V1 \
  --tissue-backend freesurfer \
  --registration-t1-target brain-csf \
  --parcellation-mode mni \
  --verbose
```

FreeSurfer intermediates are kept under `/out/freesurfer`. MRSIPrep also writes
native T1-space derivatives from those intermediates, including
`desc-brain_T1w`, `desc-brain_mask`, `desc-p1_T1w`, `desc-p2_T1w`, and
`desc-p3_T1w`.

## BIDS import utilities

```bash
mrsiprep-import /source/folder /path/to/bids --subject S001 --session V1
mrsiprep-skullstrip /path/to/bids --device cuda
mrsiprep-import-gui
```

The import helpers preserve the MRSI-Metabolic-Connectome derivative layout:
`derivatives/mrsi-orig`, `derivatives/cat12`, and `derivatives/skullstrip`.

## Attribution

Substantial implementation logic is cropped and refactored by Federico Lucchetti and Edgar Céléreau. The original
license is included in `LICENSE`.
