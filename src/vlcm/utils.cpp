#include "utils.h"
#include "exception.h"

#include <cmath>
#include <iostream>


std::ifstream &read_value(std::ifstream &s, double &d, size_t count) {
    std::string line;
    double parsed;

    s >> line;
    if (line.empty()) {
        if (s.peek() != EOF) {
            std::cout << "Warning: ignoring empty line #" << count + 1 << " in input file" << std::endl;
        }
        d = NAN;
        return s;
    }

    try {
        parsed = std::stod(line);
    } catch (std::invalid_argument &e) {
        throw VLCMException("invalid argument: Could not parse line number "
                             + std::to_string(count + 1) + " from input file.");
    } catch (std::out_of_range &e) {
        throw VLCMException("out of range: Could not parse line number "
                             + std::to_string(count + 1) + " from input file.");
    }
    d = parsed;
    return s;
}

void readFile(const std::string &filename, std::vector<double> &v) {
    std::ifstream f(filename);
    if (f.fail()) {
        throw VLCMException("Unable to open " + filename + " for reading, please make sure it exists");
    }
    double num;
    while (read_value(f, num, v.size()) && f.peek() != EOF) {
        if (num != NAN) {
            v.push_back(num);
        }
    }
}
