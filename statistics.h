#include <cmath>
#include <random>
#include <vector>

// Helper Functions
double factorial(int n) {
    if (n <= 1) return 1;
    return n * factorial(n - 1);
}

// Random number generator
random_device rd;
mt19937 gen(rd());

// Discrete Distributions
int binomial(int trials, double p) {
    binomial_distribution<> binom(trials, p);
    return binom(gen);
}

int poisson(double lambda) {
    poisson_distribution<int> pois(lambda);
    return pois(gen);
}

int geometric(double p) {
    geometric_distribution<int> geom(p);
    return geom(gen);
}

// Continuous Distributions
double normal(double mean, double stddev) {
    normal_distribution<double> norm(mean, stddev);
    return norm(gen);
}

double exponential(double lambda) {
    exponential_distribution<double> expo(lambda);
    return expo(gen);
}

double uniform(double a, double b) {
    uniform_real_distribution<double> unif(a, b);
    return unif(gen);
}

// Descriptive Statistics
double mean(vector<double> data) {
    double sum = 0;
    for (double val : data) {
        sum += val;
    }
    return sum / data.size();
}

double variance(vector<double> data) {
    double avg = mean(data);
    double sum = 0;
    for (double val : data) {
        sum += pow(val - avg, 2);
    }
    return sum / data.size();
}

double standard_deviation(vector<double> data) {
    return sqrt(variance(data));
}
