# *MRSIPrep*: A Robust Preprocessing Pipeline for Whole-Brain MRSI Data

*MRSIPrep* is a preprocessing and derivative-generation pipeline for already
quantified whole-brain MRSI maps, run as a BIDS App via Docker.

[![Docker Pulls](https://img.shields.io/docker/pulls/fedlucchetti/mrsiprep)](https://hub.docker.com/r/fedlucchetti/mrsiprep)
[![Documentation Status](https://app.readthedocs.org/projects/mrsiprep/badge/?version=latest)](https://mrsiprep.readthedocs.io/en/latest/)
[![License: CHUV academic non-commercial](https://img.shields.io/badge/license-academic--non--commercial-blue)](LICENSE)

## About

`MRSIPrep` does not perform spectral fitting. It expects quantified MRSI maps,
quality maps, and T1w images as input. Its default light mode normalizes MRSI
maps and uses fast SynthSeg cortical parcellation for parcelwise anatomical
coverage and CRLB reporting. Full mode adds SynthSeg+FAST tissue maps,
PETPVC, and Chimera/MNI-atlas regional profile extraction.

**Full documentation, installation, and usage instructions are on
[Read the Docs](https://mrsiprep.readthedocs.io/en/latest/).**

## Use Cases

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

## License

MRSIPrep is distributed under the CHUV academic non-commercial research
license; see [LICENSE](LICENSE) for the full text.

## Attribution

Substantial implementation logic is cropped and refactored by Federico Lucchetti and Edgar Céléreau. The original
license is included in `LICENSE`.

## Acknowledgments

MRSIPrep builds on the work of the ANTs, FreeSurfer, FSL, PETPVC, Chimera,
and TemplateFlow projects.
