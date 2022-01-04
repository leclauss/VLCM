#ifndef VLCM_UTILS_H
#define VLCM_UTILS_H

#include <fstream>
#include <vector>

#define TIME_INIT std::chrono::high_resolution_clock::time_point _TIME_start, _TIME_end;int _TIME_diff;
#define NOW std::chrono::high_resolution_clock::now()
#define START _TIME_start = NOW;
#define STOP _TIME_end = NOW;_TIME_diff = (std::chrono::duration_cast<std::chrono::microseconds>(_TIME_end - _TIME_start).count() / 1000);
#define TIME_MS _TIME_diff
#define END(msg) STOP std::cout << msg << " (" << TIME_MS << " ms)" << std::endl;

typedef std::vector<int> Clique;
typedef std::vector<double> TimeSeries;

// Reads input time series from file
void readFile(const std::string &filename, TimeSeries &v);

#endif //VLCM_UTILS_H
