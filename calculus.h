#include <iostream>
#include <cmath>

const int N_TERMS = 10;  // Number of terms in the Fourier series

// Numerical Integration: Trapezoidal Rule
double integral(double (*func)(double), double a, double b, int n) {
    double h = (b - a) / n;  // Step size
    double sum = 0.5 * (func(a) + func(b));  // Endpoints

    for (int i = 1; i < n; i++) {
        sum += func(a + i * h);  // Sum the inner points
    }
    
    return sum * h;
}

// Numerical Derivative: Using finite differences
double derivative(double (*func)(double), double x, double h = 1e-5) {
    return (func(x + h) - func(x - h)) / (2 * h);  // Central difference approximation
}

// Numerical Limit: Evaluate the function approaching x from both sides
double limit(double (*func)(double), double x, double h = 1e-5) {
    double left_limit = func(x - h);  // Approach from left
    double right_limit = func(x + h);  // Approach from right
    if (abs(left_limit - right_limit) < 1e-5) {
        return (left_limit + right_limit) / 2;  // Average the two
    } else {
        cout << "Limit does not exist!" << endl;
        return NAN;  // Not a number if limits don't match
    }
}

// Fourier series approximation: sum of sine and cosine terms
double fourier(double (*func)(double), double x, int terms = N_TERMS) {
    double a0 = 0;
    double sum = 0;
    double L = M_PI;  // Half-period for a 2π periodic function

    // Compute a0 term
    for (int k = 1; k <= terms; k++) {
        double ak = 0;
        double bk = 0;

        // Coefficients for each term
        for (double t = 0; t <= 2 * M_PI; t += 0.01) {
            ak += func(t) * cos(k * t);
            bk += func(t) * sin(k * t);
        }
        ak *= (1.0 / L);
        bk *= (1.0 / L);

        // Accumulate Fourier series terms
        sum += ak * cos(k * x) + bk * sin(k * x);
    }
    return a0 / 2 + sum;  // a0 term + sum of series
}

// Gamma function (recursive definition)
double gamma(int k) {
    if (k > 1) {
        return (k - 1) * gamma(k - 1);
    } else {
        return 1;
    }
}

// Beta function: B(x, y) = (Γ(x) * Γ(y)) / Γ(x + y)
double beta(double x, double y) {
    return (gamma(x) * gamma(y)) / gamma(x + y);
}

// Recursive Summation function
double summation(double (*func)(int), int start, int end) {
    if (start > end) {
        return 0;
    } else {
        return func(start) + summation(func, start + 1, end);
    }
}

// Recursive Infinite Product function
double infiniteProduct(double (*func)(int), int start, int end) {
    if (start > end) {
        return 1;
    } else {
        return func(start) * infiniteProduct(func, start + 1, end);
    }
}
