from pathlib import Path
import unittest

from mrsiprep.cli.parser import parse_args


class CLITests(unittest.TestCase):
    def test_cli_defaults_brain_csf(self):
        cfg = parse_args(["/tmp/bids", "/tmp/derivatives", "participant", "--participant-label", "sub-S001"])
        self.assertEqual(cfg.registration_t1_target, "brain-csf")
        self.assertEqual(cfg.participant_label, ["sub-S001"])
        self.assertEqual(cfg.tissue_backend, "synthseg-fast")
        self.assertEqual(cfg.derivative_dir, Path("/tmp/derivatives/mrsiprep"))

    def test_cli_ants_atropos_option(self):
        cfg = parse_args(["/tmp/bids", "/tmp/out", "participant", "--tissue-backend", "ants-atropos"])
        self.assertEqual(cfg.tissue_backend, "ants-atropos")

    def test_cli_synthseg_fast_option(self):
        cfg = parse_args(["/tmp/bids", "/tmp/out", "participant", "--tissue-backend", "synthseg-fast"])
        self.assertEqual(cfg.tissue_backend, "synthseg-fast")

    def test_cli_validate_only_option(self):
        cfg = parse_args(["/tmp/bids", "/tmp/out", "participant", "--validate-only"])
        self.assertTrue(cfg.validate_only)

    def test_cli_freesurfer_option(self):
        cfg = parse_args([
            "/tmp/bids",
            "/tmp/out",
            "participant",
            "--tissue-backend",
            "freesurfer",
            "--fs-subjects-dir",
            "/tmp/fs",
            "--overwrite-freesurfer",
        ])
        self.assertEqual(cfg.tissue_backend, "freesurfer")
        self.assertEqual(cfg.freesurfer_dir, Path("/tmp/fs"))
        self.assertTrue(cfg.overwrite_freesurfer)
