#include <iostream>
#include <cmath>
using namespace std;

class Trigonometry {
public:
    // Primary Trigonometric Functions
    static double sin(double x) {
        return std::sin(x);
    }

    static double cos(double x) {
        return std::cos(x);
    }

    static double tan(double x) {
        return std::tan(x);
    }

    // Reciprocal Trigonometric Functions
    static double csc(double x) {
        return 1.0 / std::sin(x);
    }

    static double sec(double x) {
        return 1.0 / std::cos(x);
    }

    static double cot(double x) {
        return 1.0 / std::tan(x);
    }

    // Inverse Trigonometric Functions
    static double asin(double x) {
        return std::asin(x);
    }

    static double acos(double x) {
        return std::acos(x);
    }

    static double atan(double x) {
        return std::atan(x);
    }

    // Hyperbolic Trigonometric Functions
    static double sinh(double x) {
        return std::sinh(x);
    }

    static double cosh(double x) {
        return std::cosh(x);
    }

    static double tanh(double x) {
        return std::tanh(x);
    }

    // Degree and Radian Conversions
    static double deg_to_rad(double degrees) {
        return degrees * (M_PI / 180.0);
    }

    static double rad_to_deg(double radians) {
        return radians * (180.0 / M_PI);
    }
};
