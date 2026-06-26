# Private Docker Build Workflow

MRSIPrep uses two runtime images:

```text
mrsiprep-deps:cpu  reduced private external/Python dependencies
mrsiprep:cpu       full application, including raw-to-Chimera support
```

The source dependency image containing complete upstream installations is kept
locally as `mrsiprep-deps:ubuntu22.04-cpu`. It is flattened and reduced into
`mrsiprep-deps:cpu`. Python source changes only require rebuilding the thin
`mrsiprep:cpu` application image.

## Build the full CPU image

After creating `mrsiprep-deps:ubuntu22.04-cpu` with the manual or private build
workflow below, build the reduced full image with:

```bash
docker/build_cpu_image.sh
```

This keeps SynthSeg, FAST, ANTs, PETPVC, Chimera, `recon-all`, and
`mri_vol2vol`. FSL is reduced to FAST and its shared libraries. FreeSurfer is
trimmed conservatively while retaining its core reconstruction executables,
libraries, models, registration atlases, and `fsaverage` subject.

## Recommended manual FSL/FreeSurfer workflow

Build an Ubuntu 22.04 CPU bootstrap containing ANTs, ANTsPy,
PETPVC, Chimera, and the Python dependencies:

```bash
docker/build_bootstrap.sh
```

Enter a persistent named container:

```bash
docker/enter_manual_deps.sh
```

Inside it, install FSL and FreeSurfer using the provided helper or your own
commands:

```bash
bash /root/manual_install_fsl_freesurfer.sh
exit
```

The helper uses the official FSL `getfsl.sh` installer and the FreeSurfer 8.2.0
Ubuntu 22 package. You can override either URL with `FSL_INSTALLER_URL` or
`FREESURFER_DEB_URL` inside the container.

Commit and verify the private dependency image:

```bash
docker/finalize_manual_deps.sh
```

Then build the thin MRSIPrep layer:

```bash
docker/update_mrsiprep_image.sh
```

This manual path is private and convenient, but the resulting dependency image
is less reproducible than the automated secret-based build below.

## 1. Configure private sources

```bash
cp docker/private-neurodeps.env.example docker/private-neurodeps.env
```

Populate the private file with download URLs:

```text
ANTS_URL=...
FSL_URL=...
FREESURFER_URL=...
PETPVC_URL=...
CHIMERA_URL=...
```

Alternatively set host archive paths such as:

```text
ANTS_ARCHIVE=/private/software/ants.tar.gz
FSL_ARCHIVE=/private/software/fsl.tar.gz
FREESURFER_ARCHIVE=/private/software/freesurfer.tar.gz
```

The configuration and archives are passed with BuildKit secrets. They are not
stored as Docker build arguments or copied into the final image.

## 2. Build private dependencies

```bash
docker/build_private_deps.sh
```

Configuration overrides:

```bash
DEPS_IMAGE=registry.private/mrsiprep-deps:2026-06 \
REQUIRE_PETPVC=1 \
REQUIRE_CHIMERA=1 \
docker/build_private_deps.sh
```

This builds `Dockerfile.deps` and runs the dependency verifier automatically.

## 3. Test an image

```bash
docker/test_container.sh mrsiprep-deps:cpu
docker/test_container.sh mrsiprep:cpu
```

The check verifies the winning tissue tools, ANTsPy/ANTs CLI, Python imports,
and optionally PETPVC and Chimera/FreeSurfer commands.

## 4. Update MRSIPrep only

After changing Python code:

```bash
docker/update_mrsiprep_image.sh
```

This builds `Dockerfile`, which starts from the existing dependency image and
only copies/reinstalls MRSIPrep. It does not rebuild FreeSurfer, FSL, ANTs,
PETPVC, or Chimera.

Override image names when needed:

```bash
DEPS_IMAGE=registry.private/mrsiprep-deps:2026-06 \
APP_IMAGE=registry.private/mrsiprep:cpu \
docker/update_mrsiprep_image.sh
```

If `pyproject.toml` dependency requirements change, rebuild the dependency
image. Ordinary MRSIPrep Python source changes only require the thin update.

## Runtime license

FreeSurfer is installed in the private image, but its license is mounted at
runtime:

```bash
docker run --rm \
  -v /path/license.txt:/opt/freesurfer/license.txt:ro \
  -e FS_LICENSE=/opt/freesurfer/license.txt \
  mrsiprep:cpu /data /out participant
```
