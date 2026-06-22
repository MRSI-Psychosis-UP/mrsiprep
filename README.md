# MRSIPrep

`MRSIPrep` is a preprocessing and derivative-generation pipeline for already
quantified whole-brain MRSI maps. It is derived from the implementation in
`MRSI-Metabolic-Connectome` and preserves the CHUV academic non-commercial
research license.

The package does not perform spectral fitting. It expects quantified MRSI maps,
quality maps, and T1w images. By default it uses SynthSeg plus FSL FAST to
create GM/WM/CSF partial-volume tissue maps, with HD-BET and SynthSeg
CSF/ventricle labels defining the FAST working mask.

## Minimal command

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --validate-only
```

Use `--validate-only` to check all selected subject/session inputs before
starting an expensive batch run. It reports invalid recordings and exits without
running SynthSeg, HD-BET, FAST, registration, parcellation, or PVC.

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --metabolites CrPCr GluGln GPCPCh NAANAAG Ins \
  --registration-t1-target brain-csf \
  --tissue-backend synthseg-fast \
  --parcellation-mode chimera \
  --chimera-scheme LFMIHIFIS \
  --chimera-scale 3
```

The default registration target is `brain-csf`, which adds the `p3` CSF layer to
the skull-stripped T1w image before MRSI-to-T1 registration. With the default
SynthSeg+FAST backend, SynthSeg contributes CSF/ventricle labels, HD-BET
contributes the brain mask, and FAST estimates partial-volume `p1/p2/p3` maps.
With `--tissue-backend existing`, precomputed CAT12-style `p1/p2/p3` maps are
required; if they are missing, the current subject/session fails and batch
processing continues with the next item.

## Docker / BIDS App usage

If you already have a local `mrsiprep:latest` image available, run it like fMRIPrep:

```bash
docker run --rm \
  -v /path/to/bids:/data:ro \
  -v /path/to/bids/derivatives:/out \
  -v /usr/local/freesurfer:/usr/local/freesurfer:ro \
  -v /usr/local/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
  -e FREESURFER_HOME=/usr/local/freesurfer \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  mrsiprep:latest \
  /data /out participant \
  --participant-label S001 \
  --session-label V1 \
  --tissue-backend synthseg-fast \
  --registration-t1-target brain-csf
```

For the dummy dataset on this workstation:

```bash
docker run --rm \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project:/data:ro \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project/derivatives:/out \
  -v /usr/local/freesurfer:/usr/local/freesurfer:ro \
  -v /usr/local/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
  -e FREESURFER_HOME=/usr/local/freesurfer \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  mrsiprep:latest \
  /data /out participant \
  --participant-label CHUVUP013 \
  --session-label V1 \
  --tissue-backend synthseg-fast \
  --registration-t1-target brain-csf \
  --parcellation-mode mni \
  --verbose
```

To verify the external binaries from inside the container:

```bash
docker run --rm \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project:/data:ro \
  -v /home/flucchetti/Connectome/BIDS/Dummy-Project/derivatives:/out \
  -v /usr/local/freesurfer:/usr/local/freesurfer:ro \
  -v /usr/local/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
  -e FREESURFER_HOME=/usr/local/freesurfer \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  mrsiprep:latest \
  /data /out participant \
  --check-external-libs
```

For host-mounted FreeSurfer, the container image does not embed `/usr/local/freesurfer`; the license must be mounted at runtime. SynthSeg+FAST intermediates are kept under the configured work directory.
MRSIPrep writes native T1-space derivatives including GM/WM/CSF probsegs and
`desc-p1_T1w`, `desc-p2_T1w`, and `desc-p3_T1w`.

MRSIPrep outputs are grouped by processing space:

```text
<out>/mrsiprep/sub-*/ses-*/mrsi-orig/      native/imported-grid MRSI signal maps
<out>/mrsiprep/sub-*/ses-*/mrsi-orig-pvc/  PVC-corrected native-grid maps
<out>/mrsiprep/sub-*/ses-*/mrsi-t1w/       T1w-aligned MRSI maps
<out>/mrsiprep/sub-*/ses-*/mrsi-mni/       MNI-normalized MRSI maps
<out>/mrsiprep/sub-*/ses-*/tissue-mrsi/    MRSI-grid tissue probsegs and 4Dtissue
<out>/mrsiprep/sub-*/ses-*/qmasks/         QC, spike, and brain masks
<out>/mrsiprep/sub-*/ses-*/anat/           T1w tissue and registration files
<out>/chimera-atlases/sub-*/ses-*/anat/    Chimera atlas outputs
<out>/mrsiprep/sub-*/ses-*/connectomics/   matrices, nodes, and edges
```

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
