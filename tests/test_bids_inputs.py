import unittest

from mrsiprep.config.settings import MRSIPrepConfig
from mrsiprep.io.validators import validate_recording


class BIDSInputTests(unittest.TestCase):
    def test_brain_csf_requires_p3(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as td:
            import nibabel as nib
            import numpy as np

            tmp_path = Path(td)
            bids = tmp_path / "bids"
            anat = bids / "sub-S001" / "ses-V1" / "anat"
            mrsi = bids / "derivatives" / "mrsi-orig" / "sub-S001" / "ses-V1"
            skull = bids / "derivatives" / "skullstrip" / "sub-S001" / "ses-V1"
            anat.mkdir(parents=True)
            mrsi.mkdir(parents=True)
            skull.mkdir(parents=True)
            img = nib.Nifti1Image(np.ones((2, 2, 2), dtype="float32"), np.eye(4))
            nib.save(img, skull / "sub-S001_ses-V1_desc-brain_T1w.nii.gz")
            for met in ["CrPCr", "GluGln", "GPCPCh", "NAANAAG", "Ins"]:
                nib.save(img, mrsi / f"sub-S001_ses-V1_space-orig_met-{met}_desc-signal_mrsi.nii.gz")
            cfg = MRSIPrepConfig(bids, tmp_path / "derivatives", "participant", participant_label=["S001"], session_label=["V1"])
            with self.assertRaisesRegex(Exception, "p3"):
                validate_recording(cfg, "S001", "V1")
