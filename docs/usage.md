# Usage

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --registration-t1-target brain-csf \
  --tissue-backend existing
```

`brain-csf` is the default registration target. It requires:

- a skull-stripped T1w derivative in `derivatives/skullstrip`,
- the raw BIDS T1w acquisition,
- CAT12 `p3` CSF probability map in `derivatives/cat12`.

If `p3` is missing, the recording fails. In batch processing, MRSIPrep logs the
error and continues with the next subject/session.

Use ANTs Atropos tissue segmentation with:

```bash
mrsiprep /path/to/bids /path/to/derivatives participant \
  --participant-label S001 \
  --session-label V1 \
  --tissue-backend ants-atropos
```
