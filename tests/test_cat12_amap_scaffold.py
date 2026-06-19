from pathlib import Path
import tempfile
import unittest

import numpy as np

from mrsiprep.tissue.cat12_amap import AmapParameters, CAT12_AMAP_SOURCE_DIR, build_cat12_amap_library, prepare_amap_arrays, run_cat12_amap
from mrsiprep.tissue.cat12_cleanup import CATLabels, Cat12CleanupParameters, clean_gwc, correct_outer_csf_with_atlas
from mrsiprep.tissue.cat12_segment import Cat12PreAmapParameters, _initial_labels_from_probabilities, prepare_cat12_amap_input
from mrsiprep.tissue.tpm_gmm import Cat12TemplateAssets, TpmGmmParameters, _subject_to_template_transform_set, fit_tpm_weighted_gmm


class CAT12AMAPScaffoldTests(unittest.TestCase):
    def test_vendored_source_files_exist(self):
        expected = {
            "Amap.c",
            "Amap.h",
            "Kmeans.c",
            "MrfPrior.c",
            "Pve.c",
            "vollib.c",
            "cat_amap.c",
            "cat_amap_core.c",
            "cat_amap_core.h",
            "COPYING.CAT12",
        }
        self.assertTrue(CAT12_AMAP_SOURCE_DIR.exists())
        self.assertTrue(expected.issubset({path.name for path in Path(CAT12_AMAP_SOURCE_DIR).iterdir()}))

    def test_prepare_amap_arrays_uses_fortran_order(self):
        src = np.ones((3, 4, 5), dtype=np.float32)
        label = np.zeros((3, 4, 5), dtype=np.uint8)
        src_out, label_out = prepare_amap_arrays(src, label)
        self.assertEqual(src_out.dtype, np.float64)
        self.assertEqual(label_out.dtype, np.uint8)
        self.assertTrue(src_out.flags.f_contiguous)
        self.assertTrue(label_out.flags.f_contiguous)

    def test_initial_labels_match_cat_order(self):
        probs = np.zeros((2, 2, 2, 3), dtype=np.float32)
        probs[0, 0, 0, 0] = 0.9  # GM -> CAT AMAP label 2
        probs[0, 0, 1, 1] = 0.9  # WM -> CAT AMAP label 3
        probs[0, 1, 0, 2] = 0.9  # CSF -> CAT AMAP label 1
        labels = _initial_labels_from_probabilities(probs)
        self.assertEqual(labels[0, 0, 0], 2)
        self.assertEqual(labels[0, 0, 1], 3)
        self.assertEqual(labels[0, 1, 0], 1)
        self.assertEqual(labels[1, 1, 1], 0)

    def test_prepare_cat12_amap_input_lifts_low_csf_floor(self):
        src = np.zeros((12, 12, 8), dtype=np.float32)
        probs = np.zeros((*src.shape, 3), dtype=np.float32)
        probs[2:5, 2:10, 2:6, 2] = 0.95  # CSF
        probs[5:8, 2:10, 2:6, 0] = 0.95  # GM
        probs[8:11, 2:10, 2:6, 1] = 0.95  # WM
        src[probs[..., 2] > 0] = 0.01
        src[probs[..., 0] > 0] = 2.0 / 3.0
        src[probs[..., 1] > 0] = 1.0

        prepared = prepare_cat12_amap_input(
            src,
            probs,
            (1.0, 1.0, 1.0),
            AmapParameters(mrf_weight=0.0),
            Cat12PreAmapParameters(add_csf_noise=False, auto_mrf=False),
        )

        self.assertEqual(prepared.labels[3, 4, 4], 1)
        self.assertEqual(prepared.labels[6, 4, 4], 2)
        self.assertEqual(prepared.labels[9, 4, 4], 3)
        self.assertGreater(float(prepared.csf_floor[3, 4, 4]), 0.1)
        self.assertGreater(float(prepared.prepared[3, 4, 4]), float(src[3, 4, 4]))
        self.assertEqual(float(prepared.prepared[0, 0, 0]), 0.0)
        self.assertEqual(prepared.mrf_weight, 0.0)

    def test_prepare_cat12_amap_input_can_use_seed_labels(self):
        src = np.full((6, 6, 6), 0.5, dtype=np.float32)
        probs = np.zeros((*src.shape, 3), dtype=np.float32)
        probs[1:5, 1:5, 1:5, 0] = 0.9
        seed = np.zeros(src.shape, dtype=np.uint8)
        seed[1:5, 1:5, 1:5] = 3

        prepared = prepare_cat12_amap_input(
            src,
            probs,
            (1.0, 1.0, 1.0),
            pre_amap_parameters=Cat12PreAmapParameters(enabled=False),
            seed_labels=seed,
        )

        self.assertEqual(prepared.labels[2, 2, 2], 3)
        self.assertEqual(float(prepared.prepared[0, 0, 0]), 0.0)

    def test_prepare_cat12_amap_input_can_gate_weak_generated_seed_labels(self):
        src = np.full((5, 5, 5), 0.5, dtype=np.float32)
        probs = np.zeros((*src.shape, 3), dtype=np.float32)
        probs[1, 1, 1, 2] = 0.2  # weak CSF seed should be removed
        probs[2, 2, 2, 2] = 0.5  # confident CSF seed should remain
        probs[3, 3, 3, 0] = 0.6  # GM remains in CAT label order as 2

        prepared = prepare_cat12_amap_input(
            src,
            probs,
            (1.0, 1.0, 1.0),
            pre_amap_parameters=Cat12PreAmapParameters(
                enabled=False,
                seed_min_probability=0.3,
                csf_seed_min_probability=0.35,
            ),
        )

        self.assertEqual(prepared.labels[1, 1, 1], 0)
        self.assertEqual(prepared.labels[2, 2, 2], 1)
        self.assertEqual(prepared.labels[3, 3, 3], 2)

    def test_prepare_cat12_amap_input_can_seed_from_six_cat_classes(self):
        src = np.full((5, 5, 5), 0.5, dtype=np.float32)
        probs = np.zeros((*src.shape, 3), dtype=np.float32)
        probs[1, 1, 1, 2] = 0.8
        probs[2, 2, 2, 0] = 0.8
        probs[3, 3, 3, 1] = 0.8
        probs[1, 3, 1, 2] = 0.8
        probs[3, 1, 3, 2] = 0.8
        class_probs = np.zeros((*src.shape, 6), dtype=np.float32)
        class_probs[1, 1, 1, 5] = 0.95  # BG suppresses a three-class CSF seed
        class_probs[2, 2, 2, 0] = 0.95  # GM -> AMAP label 2
        class_probs[3, 3, 3, 1] = 0.95  # WM -> AMAP label 3
        class_probs[1, 3, 1, 3] = 0.95  # bone follows CAT's CSF-like label 1
        class_probs[3, 1, 3, 4] = 0.95  # soft tissue -> background label 0

        prepared = prepare_cat12_amap_input(
            src,
            probs,
            (1.0, 1.0, 1.0),
            pre_amap_parameters=Cat12PreAmapParameters(enabled=False, use_class_probabilities=True),
            class_probabilities=class_probs,
        )

        self.assertEqual(prepared.labels[1, 1, 1], 0)
        self.assertEqual(prepared.labels[2, 2, 2], 2)
        self.assertEqual(prepared.labels[3, 3, 3], 3)
        self.assertEqual(prepared.labels[1, 3, 1], 1)
        self.assertEqual(prepared.labels[3, 1, 3], 0)

    def test_prepare_cat12_amap_input_ignores_six_classes_by_default(self):
        src = np.full((5, 5, 5), 0.5, dtype=np.float32)
        probs = np.zeros((*src.shape, 3), dtype=np.float32)
        probs[1, 1, 1, 2] = 0.8
        class_probs = np.zeros((*src.shape, 6), dtype=np.float32)
        class_probs[1, 1, 1, 5] = 0.95

        prepared = prepare_cat12_amap_input(
            src,
            probs,
            (1.0, 1.0, 1.0),
            pre_amap_parameters=Cat12PreAmapParameters(enabled=False),
            class_probabilities=class_probs,
        )

        self.assertEqual(prepared.labels[1, 1, 1], 1)

    def test_ctypes_library_builds(self):
        lib = build_cat12_amap_library()
        self.assertTrue(lib.exists())

    def test_run_cat12_amap_small_volume(self):
        src = np.zeros((12, 12, 8), dtype=np.float64)
        label = np.zeros(src.shape, dtype=np.uint8)
        label[2:5, 2:10, 2:6] = 1
        label[5:8, 2:10, 2:6] = 2
        label[8:11, 2:10, 2:6] = 3
        src[label == 1] = 1.0 / 3.0
        src[label == 2] = 2.0 / 3.0
        src[label == 3] = 1.0
        result = run_cat12_amap(
            src,
            label,
            (1.0, 1.0, 1.0),
            AmapParameters(n_iters=1, sub=4, pve=0, mrf_weight=0.0, iters_icm=0),
        )
        self.assertEqual(result.probabilities.shape, (*src.shape, 3))
        self.assertEqual(result.means.shape, (3,))
        self.assertEqual(result.stds.shape, (3,))
        self.assertTrue(np.isfinite(result.probabilities).all())
        self.assertGreaterEqual(float(np.min(result.probabilities)), 0.0)
        self.assertLessEqual(float(np.max(result.probabilities)), 1.0)

    def test_cat12_template_assets_exist(self):
        assets = Cat12TemplateAssets().validated()
        self.assertTrue(assets.tpm.exists())
        self.assertTrue(assets.template.exists())
        self.assertTrue(assets.template_mask.exists())
        self.assertTrue(assets.atlas.exists())

    def test_tpm_weighted_gmm_uses_priors(self):
        image = np.zeros((6, 6, 6), dtype=np.float32)
        priors = np.zeros((*image.shape, 3), dtype=np.float32)
        image[:2, :, :] = 10.0
        image[2:4, :, :] = 50.0
        image[4:, :, :] = 90.0
        priors[:2, :, :, 2] = 1.0
        priors[2:4, :, :, 0] = 1.0
        priors[4:, :, :, 1] = 1.0
        posterior, stats = fit_tpm_weighted_gmm(
            image,
            priors,
            TpmGmmParameters(em_iterations=2, n4_bias_correct=False, support_threshold=0.01),
        )
        self.assertEqual(posterior.shape, (*image.shape, 3))
        self.assertGreater(float(np.mean(posterior[2:4, :, :, 0])), 0.9)
        self.assertGreater(float(np.mean(posterior[4:, :, :, 1])), 0.9)
        self.assertGreater(float(np.mean(posterior[:2, :, :, 2])), 0.9)
        self.assertLess(stats["CSF"]["mean"], stats["GM"]["mean"])
        self.assertLess(stats["GM"]["mean"], stats["WM"]["mean"])

    def test_tpm_weighted_gmm_supports_cat12_six_class_model(self):
        image = np.zeros((6, 6, 6), dtype=np.float32)
        priors = np.zeros((*image.shape, 6), dtype=np.float32)
        intensities = [50.0, 90.0, 10.0, 130.0, 170.0, 0.0]
        for idx, intensity in enumerate(intensities):
            image[idx, :, :] = intensity
            priors[idx, :, :, idx] = 1.0
        posterior, stats = fit_tpm_weighted_gmm(
            image,
            priors,
            TpmGmmParameters(
                em_iterations=2,
                n4_bias_correct=False,
                support_threshold=0.01,
                gaussians_per_class=(1, 1, 1, 1, 1, 1),
            ),
        )
        self.assertEqual(posterior.shape, (*image.shape, 6))
        for idx in range(6):
            self.assertGreater(float(np.mean(posterior[idx, :, :, idx])), 0.9)
        self.assertIn("BONE", stats)
        self.assertIn("SOFT", stats)
        self.assertIn("BG", stats)

    def test_subject_to_template_inverse_transform_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "sub-01_desc-T1w_to_cat12Template"
            affine = prefix.with_suffix(".affine.mat")
            warp = prefix.with_suffix(".syn.nii.gz")
            inverse_warp = prefix.with_suffix(".syn_inv.nii.gz")

            affine.write_text("affine\n")
            transforms = _subject_to_template_transform_set(prefix)
            self.assertEqual(transforms["forward"], [str(affine)])
            self.assertEqual(transforms["inverse"], [f"[{affine},1]"])

            warp.write_text("warp\n")
            transforms = _subject_to_template_transform_set(prefix)
            self.assertEqual(transforms, {"forward": [], "inverse": []})

            inverse_warp.write_text("inverse warp\n")
            transforms = _subject_to_template_transform_set(prefix)
            self.assertEqual(transforms["forward"], [str(warp), str(affine)])
            self.assertEqual(transforms["inverse"], [f"[{affine},1]", str(inverse_warp)])

    def test_cat12_clean_gwc_preserves_probability_contract(self):
        probs = np.zeros((24, 24, 16, 3), dtype=np.float32)
        probs[7:17, 7:17, 5:11, 1] = 1.0
        probs[6:18, 6:18, 4:12, 0] = np.maximum(probs[6:18, 6:18, 4:12, 0], 0.6)
        probs[5:19, 5:19, 3:13, 2] = np.maximum(probs[5:19, 5:19, 3:13, 2], 0.4)
        cleaned = clean_gwc(probs, (1.0, 1.0, 1.0), Cat12CleanupParameters(cleanup_strength=0.5))
        self.assertEqual(cleaned.shape, probs.shape)
        self.assertGreaterEqual(float(np.min(cleaned)), 0.0)
        self.assertLessEqual(float(np.max(cleaned)), 1.0)
        support = np.sum(cleaned, axis=-1) > 0
        self.assertTrue(np.allclose(np.sum(cleaned[support], axis=-1), 1.0, atol=1.0 / 255.0))

    def test_cat12_clean_gwc_removes_isolated_gmwm_island(self):
        probs = np.zeros((32, 32, 18, 3), dtype=np.float32)
        probs[10:22, 10:22, 6:12, 1] = 1.0
        probs[9:23, 9:23, 5:13, 0] = np.maximum(probs[9:23, 9:23, 5:13, 0], 0.7)
        probs[8:24, 8:24, 4:14, 2] = np.maximum(probs[8:24, 8:24, 4:14, 2], 0.3)
        probs[1:4, 1:4, 1:4, 0] = 1.0
        cleaned = clean_gwc(probs, (1.0, 1.0, 1.0), Cat12CleanupParameters(cleanup_strength=0.5))
        self.assertGreater(float(np.sum(cleaned[10:22, 10:22, 6:12])), 0.0)
        self.assertEqual(float(np.sum(cleaned[1:4, 1:4, 1:4])), 0.0)

    def test_cat12_atlas_cleanup_converts_boundary_gmwm_to_csf(self):
        probs = np.zeros((32, 32, 18, 3), dtype=np.float32)
        probs[8:24, 8:24, 5:13, 0] = 0.7
        probs[10:22, 10:22, 6:12, 1] = 1.0
        probs[7:25, 7:25, 4:14, 2] = np.maximum(probs[7:25, 7:25, 4:14, 2], 0.2)
        probs[4:8, 8:24, 5:13, 1] = 1.0
        intensity = np.zeros(probs.shape[:3], dtype=np.float32)
        intensity[probs[..., 2] > 0] = 1.0
        intensity[probs[..., 0] > 0] = 2.0
        intensity[probs[..., 1] > 0] = 3.0
        intensity[4:8, 8:24, 5:13] = 2.2
        atlas = np.zeros(probs.shape[:3], dtype=np.uint8)
        atlas[:, :16, :] = CATLabels.CT
        atlas[:, 16:, :] = CATLabels.CT + 1
        cleaned = correct_outer_csf_with_atlas(
            probs,
            intensity,
            atlas,
            (1.0, 1.0, 1.0),
            Cat12CleanupParameters(cleanup_strength=0.5),
        )
        self.assertEqual(cleaned.shape, probs.shape)
        self.assertGreater(float(np.mean(cleaned[4:8, 8:24, 5:13, 2])), 0.5)
        self.assertLess(float(np.mean(cleaned[4:8, 8:24, 5:13, 1])), 0.5)


if __name__ == "__main__":
    unittest.main()
