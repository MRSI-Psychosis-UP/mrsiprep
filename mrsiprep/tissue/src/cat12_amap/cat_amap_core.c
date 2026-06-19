/*
 * Non-MEX CAT12 AMAP wrapper for MRSIPrep.
 *
 * The call sequence mirrors CAT12's cat_amap.c MEX wrapper:
 * optional K-means initialization, AMAP, then PVE expansion.
 */

#include "cat_amap_core.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Amap.h"

int cat_amap_output_classes(int n_classes, int pve)
{
    if (n_classes <= 0) {
        return 0;
    }
    if (pve == 6) {
        return n_classes + 3;
    }
    if (pve == 5) {
        return n_classes + 2;
    }
    return n_classes;
}

const char *cat_amap_error_string(int code)
{
    switch (code) {
    case CAT_AMAP_OK:
        return "ok";
    case CAT_AMAP_ERR_NULL:
        return "null input pointer";
    case CAT_AMAP_ERR_DIMS:
        return "invalid image dimensions";
    case CAT_AMAP_ERR_OPTIONS:
        return "invalid AMAP options";
    case CAT_AMAP_ERR_ALLOC:
        return "memory allocation failed";
    default:
        return "unknown CAT AMAP error";
    }
}

static CatAmapOptions cat_amap_default_options(void)
{
    CatAmapOptions options;
    options.n_classes = 3;
    options.n_iters = 10;
    options.sub = 32;
    options.pve = 5;
    options.init_kmeans = 0;
    options.mrf_weight = 0.0;
    options.iters_icm = 0;
    options.bias_fwhm = 0.0;
    options.verbose = 0;
    return options;
}

int cat_amap_run(
    const double *src0,
    const unsigned char *label0,
    const int dims[3],
    const double voxelsize[3],
    const CatAmapOptions *options_in,
    unsigned char *prob,
    double *means,
    double *stds,
    double *bias_corrected
)
{
    CatAmapOptions options = options_in ? *options_in : cat_amap_default_options();
    int dims2[4];
    long nvox;
    int out_classes;
    double *src;
    unsigned char *label;
    double *mean;
    double *fmeans;
    double *fstds;
    double max_vol = -1e15;
    double offset;

    if (!src0 || !label0 || !dims || !voxelsize || !prob) {
        return CAT_AMAP_ERR_NULL;
    }
    if (dims[0] <= 0 || dims[1] <= 0 || dims[2] <= 0) {
        return CAT_AMAP_ERR_DIMS;
    }
    if (options.n_classes <= 0 || options.n_classes > MAX_NC || options.n_iters < 0 ||
        options.sub <= 0 || options.iters_icm < 0) {
        return CAT_AMAP_ERR_OPTIONS;
    }

    out_classes = cat_amap_output_classes(options.n_classes, options.pve);
    if (out_classes <= 0) {
        return CAT_AMAP_ERR_OPTIONS;
    }

    dims2[0] = dims[0];
    dims2[1] = dims[1];
    dims2[2] = dims[2];
    dims2[3] = out_classes;
    nvox = (long)dims[0] * (long)dims[1] * (long)dims[2];

    src = (double *)malloc(sizeof(double) * (size_t)nvox);
    label = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
    mean = (double *)calloc((size_t)options.n_classes + 3u, sizeof(double));
    fmeans = (double *)calloc((size_t)options.n_classes + 3u, sizeof(double));
    fstds = (double *)calloc((size_t)options.n_classes + 3u, sizeof(double));
    if (!src || !label || !mean || !fmeans || !fstds) {
        free(src);
        free(label);
        free(mean);
        free(fmeans);
        free(fstds);
        return CAT_AMAP_ERR_ALLOC;
    }

    memcpy(src, src0, sizeof(double) * (size_t)nvox);
    memcpy(label, label0, sizeof(unsigned char) * (size_t)nvox);
    memset(prob, 0, sizeof(unsigned char) * (size_t)nvox * (size_t)out_classes);

    for (long i = 0; i < nvox; i++) {
        if (src[i] > max_vol) {
            max_vol = src[i];
        }
    }
    offset = 0.2 * max_vol;
    for (long i = 0; i < nvox; i++) {
        if (label[i] > 0) {
            src[i] += offset;
        }
    }

    if (options.init_kmeans > 0) {
        unsigned char *mask = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)nvox);
        if (!mask) {
            free(src);
            free(label);
            free(mean);
            free(fmeans);
            free(fstds);
            return CAT_AMAP_ERR_ALLOC;
        }
        for (long i = 0; i < nvox; i++) {
            mask[i] = (src[i] > 0.0) ? 255 : 0;
        }
        (void)Kmeans(src, label, mask, 25, options.n_classes, (double *)voxelsize,
                     dims2, 0, 128, 0, KMEANS, options.bias_fwhm);
        (void)Kmeans(src, label, mask, 25, options.n_classes, (double *)voxelsize,
                     dims2, 0, 128, 0, NOPVE, options.bias_fwhm);
        free(mask);
    }

    Amap(src, label, prob, mean, options.n_classes, options.n_iters, options.sub,
         dims2, options.pve, options.mrf_weight, (double *)voxelsize,
         options.iters_icm, offset, options.bias_fwhm, options.verbose, fmeans, fstds);

    if (options.pve == 6) {
        Pve6(src, prob, label, mean, dims2);
    }
    if (options.pve == 5) {
        Pve5(src, prob, label, mean, dims2);
    }

    if (means) {
        for (int i = 0; i < options.n_classes; i++) {
            means[i] = fmeans[i];
        }
    }
    if (stds) {
        for (int i = 0; i < options.n_classes; i++) {
            stds[i] = fstds[i];
        }
    }
    if (bias_corrected) {
        for (long i = 0; i < nvox; i++) {
            bias_corrected[i] = src[i] - offset;
        }
    }

    free(src);
    free(label);
    free(mean);
    free(fmeans);
    free(fstds);
    return CAT_AMAP_OK;
}
