#pragma once

class VarianceEstimator{
public:
    void Add(double val) {
        ++sampleCount;
        double delta = val - mean;
        mean += delta / sampleCount;
        double delta2 = val - mean;
        S += delta * delta2;
    }

    double Mean() const {
        return mean;
    }

    double Variance() const {
        return sampleCount > 1 ? S / (sampleCount - 1) : 0;
    }

    std::size_t Samples() const {
        return sampleCount;
    }

    double RelativeVariance() const {
        return mean == 0 ? 0 : 1.96 * std::sqrt(Variance()/sampleCount) / mean;
    }
private:
    double S = 0;
    double mean = 0;
    std::size_t sampleCount = 0;
};