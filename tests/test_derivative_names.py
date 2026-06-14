from pathlib import Path
import unittest

from mrsiprep.io.naming import anat_derivative, mrsi_derivative, parcellation_derivative


class DerivativeNameTests(unittest.TestCase):
    def test_derivative_names(self):
        root = Path("/out")
        self.assertTrue(str(mrsi_derivative(root, "sub-S001", "ses-V1", space="MRSI", met="CrPCr", desc="qcmask", suffix_override="mask")).endswith(
            "sub-S001/ses-V1/mrsi/sub-S001_ses-V1_space-MRSI_met-CrPCr_desc-qcmask_mask.nii.gz"
        ))
        self.assertTrue(str(anat_derivative(root, "S001", "V1", space="T1w", desc="brainCSF")).endswith(
            "sub-S001/ses-V1/anat/sub-S001_ses-V1_space-T1w_desc-brainCSF_T1w.nii.gz"
        ))
        self.assertTrue(str(parcellation_derivative(root, "S001", "V1", space="MRSI", atlas="chimera", scale="scale3", desc="regional", suffix_override="tsv")).endswith(
            "sub-S001/ses-V1/parcellations/sub-S001_ses-V1_space-MRSI_atlas-chimera_scale-scale3_desc-regional.tsv"
        ))
