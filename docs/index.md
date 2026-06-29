# *MRSIPrep*: A Robust Preprocessing Pipeline for Whole-Brain MRSI Data

*MRSIPrep* is a preprocessing and derivative-generation pipeline for already
quantified whole-brain MRSI maps, run as a BIDS App via Docker.

[![Docker Pulls](https://img.shields.io/docker/pulls/fedlucchetti/mrsiprep)](https://hub.docker.com/r/fedlucchetti/mrsiprep)
[![Documentation Status](https://app.readthedocs.org/projects/mrsiprep/badge/?version=latest)](https://mrsiprep.readthedocs.io/en/latest/)
[![License: CHUV academic non-commercial](https://img.shields.io/badge/license-academic--non--commercial-blue)](https://github.com/MRSI-Psychosis-UP/mrsiprep/blob/main/LICENSE)

## About

`MRSIPrep` does not perform spectral fitting. It expects quantified MRSI
maps, quality maps, and T1w images as input, and produces BIDS-organized
derivatives: filtered/normalized MRSI maps, tissue probability maps,
anatomical parcellations, and regional metabolite tables. It is derived from
the implementation in `MRSI-Metabolic-Connectome` and preserves the CHUV
academic non-commercial research license.

## What it uses

- **[ANTs](http://stnava.github.io/ANTs/)** for MRSI↔T1w and T1w↔MNI registration.
- **[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)** (`mri_synthseg`, `recon-all`, `mri_vol2vol`) for brain
  extraction, cortical/subcortical parcellation, and surface reconstruction.
- **[FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) FAST** for tissue-class probability segmentation.
- **[PETPVC](https://github.com/UCL/PETPVC)** for partial-volume correction of MRSI maps.
- **[Chimera](https://github.com/connectomicslab/chimera)** for multi-atlas cortical/subcortical parcellation
  fusion.
- **[TemplateFlow](https://www.templateflow.org/)** for the bundled MNI152 reference templates and atlases.

## Pipelines

MRSIPrep runs in one of two modes, selected with `--mode`:

- **Light mode** — registers MRSI maps to a SynthSeg-extracted T1w image,
  resamples to the requested output spaces, and parcellates with SynthSeg
  cortical/subcortical labels. No tissue PVC, no Chimera, no `recon-all`.
  This is the fast default path for anatomical coverage and CRLB reporting.
- **Full mode** — adds SynthSeg+FAST tissue probability maps, partial-volume
  correction, and a choice of Chimera multi-atlas or bundled MNI-atlas
  parcellation, plus optional perturbation-based connectivity matrices.

Both modes share the same MRSI filtering, quality-masking, and T1w/MNI
normalization machinery; full mode is a superset of light mode's outputs.

## Principles

- **BIDS-native.** Inputs are read from and outputs are written to a BIDS
  dataset/derivatives layout, following BIDS-App conventions (`bids_dir
  output_dir participant [options]`).
- **Reproducible by default.** All external dependencies (ANTs, FreeSurfer,
  FSL, PETPVC, Chimera) are pinned inside a single Docker image; there is no
  supported host installation of the pipeline itself.
- **Idempotent steps.** Each processing step checks for existing outputs and
  skips recomputation unless an `--overwrite*` flag is passed, so batch runs
  can be safely resumed.
- **Quantification only.** MRSIPrep never performs spectral fitting — it
  starts from already-quantified MRSI maps and quality maps.

## License

MRSIPrep is distributed under the CHUV academic non-commercial research
license; see [LICENSE](https://github.com/MRSI-Psychosis-UP/mrsiprep/blob/main/LICENSE) for the full text.

## Acknowledgments

Substantial implementation logic is cropped and refactored by Federico
Lucchetti and Edgar Céléreau from `MRSI-Metabolic-Connectome`. MRSIPrep
builds on the work of the ANTs, FreeSurfer, FSL, PETPVC, Chimera, and
TemplateFlow projects.

## Publications

Code derived from this pipeline has been used in the following peer-reviewed
publications:

- Lucchetti, F., Céléreau, E., Steullet, P., Alemán-Gómez, Y., Hagmann, P.,
  Klauser, A., & Klauser, P. (2025). Constructing the human brain metabolic
  connectome with MR spectroscopic imaging reveals cerebral biochemical
  organization. *Nature Communications*, 16.
  [doi:10.1038/s41467-025-66124-w](https://doi.org/10.1038/s41467-025-66124-w)
- Céléreau, E., Lucchetti, F., Alemán-Gómez, Y., Dwir, D., Cleusix, M.,
  Ledoux, J.-B., Jenni, R., Conchon, C., Bach Cuadra, M., Schilliger, Z.,
  Solida, A., Armando, M., Plessen, K. J., Hagmann, P., Conus, P., Klauser,
  A., & Klauser, P. (2026). High-resolution whole-brain magnetic resonance
  spectroscopic imaging in youth at risk for psychosis. *Imaging
  Neuroscience*, 4.
  [doi:10.1162/imag.a.1276](https://doi.org/10.1162/imag.a.1276)

```{toctree}
:maxdepth: 2
:caption: 'Getting Started:'
:hidden:

installation
usage_basic
usage_normalization
usage_parcellation
```

```{toctree}
:maxdepth: 1
:hidden:

changelog
```
