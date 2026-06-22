FROM python:3.11-slim

LABEL org.opencontainers.image.title="MRSIPrep"
LABEL org.opencontainers.image.description="BIDS App for preprocessing quantified whole-brain MRSI derivatives"
LABEL org.opencontainers.image.licenses="CHUV academic non-commercial research license"

ENV FREESURFER_HOME=/usr/local/freesurfer \
    SUBJECTS_DIR=/out/freesurfer \
    FS_LICENSE=/opt/freesurfer/license.txt \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg \
    PATH=/usr/local/freesurfer/bin:/usr/local/freesurfer/fsfast/bin:/usr/local/freesurfer/tktools:/usr/local/freesurfer/mni/bin:$PATH

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bc \
        ca-certificates \
        libgfortran5 \
        libgomp1 \
        libgl1 \
        libglib2.0-0 \
        libglu1-mesa \
        libquadmath0 \
        libx11-6 \
        libxext6 \
        libxmu6 \
        libxt6 \
        perl \
        tcsh \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/mrsiprep
COPY pyproject.toml README.md LICENSE ./
COPY mrsiprep ./mrsiprep

RUN python -m pip install --upgrade pip \
    && python -m pip install ".[ants]"

ENTRYPOINT ["mrsiprep"]
