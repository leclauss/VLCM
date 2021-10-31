#include <vector>
#include <cmath>
#include <cstring>
#include "precompute.h"

constexpr const double FLATNESS_EPSILON = 1e-13;

TimeSeriesInfo::TimeSeriesInfo(const std::vector<double> &time_series, double corr_threshold) {
    threshold = corr_threshold;
    len = time_series.size();

    auto ts_a = new double[len];
    std::memcpy(ts_a, time_series.data(), sizeof(double) * len);
    data = ts_a;

    auto sums_a = new double[len + 1];
    sums_a[0] = 0;
    double sum = 0;
    auto squaredSums_a = new double[len + 1];
    squaredSums_a[0] = 0;
    double squaredSum = 0;
    for (size_t i = 0; i < len; i++) {
        sum += data[i];
        squaredSum += data[i] * data[i];
        sums_a[i + 1] = sum;
        squaredSums_a[i + 1] = squaredSum;
    }
    sums = sums_a;
    squaredSums = squaredSums_a;
}

TimeSeriesInfo::~TimeSeriesInfo() {
    delete[] data;
    delete[] sums;
    delete[] squaredSums;
}


// moving mean over sequence ts with window length w based on Ogita et. al, Accurate Sum and Dot Product
// A major source of rounding error is accumulated error in the mean values, so
// we use this to compensate. While the error bound is still ts function of the
// conditioning of ts very long dot product, we have observed ts reduction of 3 -
// 4 digits lost to numerical round-off when compared to older solutions.
void moving_mean(double *result, const TimeSeriesInfo &ts, int w) {
    double p = ts.data[0];
    double s = 0;
    for (int i = 1; i < w; ++i) {
        double x = p + ts.data[i];
        double z = x - p;
        s += (p - (x - z)) + (ts.data[i] - z);
        p = x;
    }
    result[0] = (p + s) / w;

    for (int i = w; i < ts.len; ++i) {
        double x = p - ts.data[i - w];
        double z = x - p;
        s += (p - (x - z)) - (ts.data[i - w] + z);
        p = x;

        x = p + ts.data[i];
        z = x - p;
        s += (p - (x - z)) + (ts.data[i] - z);
        p = x;
        result[i - w + 1] = (p + s) / w;
    }
}

// Computes the subsequence norms for ts time series, subsequence length w
// This is ts faster, less accurate method, but is good enough for nearly all cases.
void fast_sum_of_squared_difference(double *result, const TimeSeriesInfo &ts, const double *means, int w) {
    double sum = 0;
    for (int i = 0; i < w; ++i) {
        double val = ts.data[i] - means[0];
        sum += val * val;
    }
    result[0] = sum;

    for (int i = 1; i < ts.len - w + 1; ++i) {
        result[i] = result[i - 1] +
                    ((ts.data[i - 1] - means[i - 1]) + (ts.data[i + w - 1] - means[i])) *
                    (ts.data[i + w - 1] - ts.data[i - 1]);
    }
}

PrecomputedInfo::PrecomputedInfo(const TimeSeriesInfo &ts, int window, int windowMax) {
    w = window;
    maxW = windowMax;
    n = ts.len - window + 1;
    auto cov_a = new double[n];
    auto df_a = new double[n];
    auto dg_a = new double[n];
    auto norms_a = new double[n];
    auto sigSquaredInverse_a = new double[n];
    auto means_a = new double[n];

    // compute statistics
    moving_mean(means_a, ts, window);
    fast_sum_of_squared_difference(norms_a, ts, means_a, window);

    for (int i = 0; i < n; ++i) {
        // Check if the sum of differences from the mean is too small and this subsequence should be considered FLAT
        if (norms_a[i] <= FLATNESS_EPSILON) {
            sigSquaredInverse_a[i] = std::nan("NaN");
            norms_a[i] = std::nan("NaN");
        } else {
            sigSquaredInverse_a[i] = static_cast<double>(window) / norms_a[i];
            // Compute the inverse norm from the sum of squared differences
            norms_a[i] = static_cast<double>(1.0) / std::sqrt(norms_a[i]);
        }
    }

    for (int i = 0; i < n - 1; ++i) {
        df_a[i] = (ts.data[i + window] - ts.data[i]) / static_cast<double>(2);
        dg_a[i] = (ts.data[i + window] - means_a[i + 1]) + (ts.data[i] - means_a[i]);
    }

    // compute QT (cov)
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < window; ++j) {
            sum += (ts.data[i + j] - means_a[i]) * (ts.data[j] - means_a[0]);
        }
        cov_a[i] = sum;
    }

    double sqrt_w = std::sqrt(w);
    for (int i = 0; i < n; ++i) {
        means_a[i] *= sqrt_w;
    }

    df = df_a;
    dg = dg_a;
    norms = norms_a;
    cov = cov_a;
    sigSquaredInverse = sigSquaredInverse_a;
    means = means_a;
}

PrecomputedInfo::~PrecomputedInfo() {
    delete[] cov;
    delete[] df;
    delete[] dg;
    delete[] norms;
    delete[] sigSquaredInverse;
    delete[] means;
}
