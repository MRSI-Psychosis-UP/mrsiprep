from pathlib import Path
import unittest

from mrsiprep.io.naming import anat_derivative, chimera_derivative, mrsi_derivative, parcellation_derivative


class DerivativeNameTests(unittest.TestCase):
    def test_derivative_names(self):
        root = Path("/out")
        self.assertTrue(str(mrsi_derivative(root, "sub-S001", "ses-V1", space="MRSI", met="CrPCr", desc="qcmask", suffix_override="mask")).endswith(
            "sub-S001/ses-V1/qmasks/sub-S001_ses-V1_space-mrsi_met-CrPCr_desc-qcmask_mask.nii.gz"
        ))
        self.assertTrue(str(mrsi_derivative(root, "S001", "V1", space="MRSI", met="CrPCr", desc="preproc", suffix_override="mrsi")).endswith(
            "sub-S001/ses-V1/mrsi-orig/sub-S001_ses-V1_space-mrsi_met-CrPCr_desc-preproc_mrsi.nii.gz"
        ))
        self.assertTrue(str(mrsi_derivative(root, "S001", "V1", space="T1w", met="CrPCr", desc="preproc", suffix_override="mrsi")).endswith(
            "sub-S001/ses-V1/mrsi-t1w/sub-S001_ses-V1_space-T1w_met-CrPCr_desc-preproc_mrsi.nii.gz"
        ))
        self.assertTrue(str(mrsi_derivative(root, "S001", "V1", space="MNI152NLin2009cAsym", met="CrPCr", desc="preproc", suffix_override="mrsi")).endswith(
            "sub-S001/ses-V1/mrsi-mni/sub-S001_ses-V1_space-MNI152NLin2009cAsym_met-CrPCr_desc-preproc_mrsi.nii.gz"
        ))
        self.assertTrue(str(mrsi_derivative(root, "S001", "V1", space="MRSI", met="CrPCr", desc="pvc", suffix_override="mrsi")).endswith(
            "sub-S001/ses-V1/mrsi-orig-pvc/sub-S001_ses-V1_space-mrsi_met-CrPCr_desc-pvc_mrsi.nii.gz"
        ))
        self.assertTrue(str(mrsi_derivative(root, "S001", "V1", space="MRSI", label="GM", suffix_override="probseg")).endswith(
            "sub-S001/ses-V1/tissue-mrsi/sub-S001_ses-V1_space-mrsi_label-GM_probseg.nii.gz"
        ))
        self.assertTrue(str(mrsi_derivative(root, "S001", "V1", space="MRSI", desc="4Dtissue", suffix_override="mrsi")).endswith(
            "sub-S001/ses-V1/tissue-mrsi/sub-S001_ses-V1_space-mrsi_desc-4Dtissue_mrsi.nii.gz"
        ))
        self.assertTrue(str(anat_derivative(root, "S001", "V1", space="T1w", desc="brainCSF")).endswith(
            "sub-S001/ses-V1/anat/sub-S001_ses-V1_space-T1w_desc-brainCSF_T1w.nii.gz"
        ))
        self.assertTrue(str(parcellation_derivative(root, "S001", "V1", space="MRSI", atlas="chimera", scale="scale3", desc="regional", suffix_override="tsv")).endswith(
            "sub-S001/ses-V1/parcellations/sub-S001_ses-V1_space-mrsi_atlas-chimera_scale-scale3_desc-regional.tsv"
        ))
        self.assertTrue(str(chimera_derivative(root, "S001", "V1", space="MRSI", atlas="chimeraLFMIHIFIS", scale="scale3", suffix_override="dseg")).endswith(
            "chimera-atlases/sub-S001/ses-V1/anat/sub-S001_ses-V1_space-mrsi_atlas-chimeraLFMIHIFIS_scale-scale3_dseg.nii.gz"
        ))
