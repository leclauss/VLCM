#ifndef VLCM_PRECOMPUTE_H
#define VLCM_PRECOMPUTE_H

#include <cstddef>
#include <vector>

// Struct containing the precomputed statistics for an input time series
struct TimeSeriesInfo {
    TimeSeriesInfo(const std::vector<double> &time_series, double corr_threshold);

    ~TimeSeriesInfo();

    const double *__restrict data;
    const double *__restrict sums;
    const double *__restrict squaredSums;
    size_t len;
    double threshold;
};

// Struct containing the window specific precomputed statistics for an input time series
struct PrecomputedInfo {
    PrecomputedInfo(const TimeSeriesInfo &ts, int window, int windowMax);

    ~PrecomputedInfo();

    double *__restrict cov;
    const double *__restrict df;
    const double *__restrict dg;
    const double *__restrict norms;
    const double *__restrict sigSquaredInverse;
    const double *__restrict means;
    int w;
    int maxW;
    size_t n;
};

#endif //VLCM_PRECOMPUTE_H
