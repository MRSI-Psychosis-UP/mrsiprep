/*
 * Non-MEX CAT12 AMAP wrapper for MRSIPrep.
 *
 * This header defines a plain C API around CAT12's AMAP/PVE implementation.
 * It is intended for a Python extension layer and does not include MATLAB
 * headers.
 */

#ifndef MRSIPREP_CAT_AMAP_CORE_H
#define MRSIPREP_CAT_AMAP_CORE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int n_classes;
    int n_iters;
    int sub;
    int pve;
    int init_kmeans;
    double mrf_weight;
    int iters_icm;
    double bias_fwhm;
    int verbose;
} CatAmapOptions;

enum {
    CAT_AMAP_OK = 0,
    CAT_AMAP_ERR_NULL = 1,
    CAT_AMAP_ERR_DIMS = 2,
    CAT_AMAP_ERR_OPTIONS = 3,
    CAT_AMAP_ERR_ALLOC = 4
};

int cat_amap_output_classes(int n_classes, int pve);

const char *cat_amap_error_string(int code);

int cat_amap_run(
    const double *src0,
    const unsigned char *label0,
    const int dims[3],
    const double voxelsize[3],
    const CatAmapOptions *options,
    unsigned char *prob,
    double *means,
    double *stds,
    double *bias_corrected
);

#ifdef __cplusplus
}
#endif

#endif
