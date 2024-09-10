#include <iostream>
#include <cmath>
#include <complex>
using namespace std;

// Constants
const double h_bar = 1.054571817e-34;  // Reduced Planck constant (Joule-second)
const double m_e = 9.10938356e-31;     // Electron mass (kg)
const double pi = M_PI;

// Complex number type for wavefunctions
typedef complex<double> cdouble;

// Quantum Wave Function: ψ(x,t)
cdouble waveFunction(double (*psi_x)(double), double x, double t, double k, double omega) {
    return psi_x(x) * exp(cdouble(0, k * x - omega * t));  // ψ(x,t) = ψ(x) * exp(i(kx - ωt))
}

// Heisenberg Uncertainty Principle: Δx * Δp >= h_bar / 2
double uncertainty(double delta_x, double delta_p) {
    return delta_x * delta_p - h_bar / 2;
}

// Time-independent Schrödinger Equation: Hψ = Eψ
// H = - (h_bar^2 / 2m) * (d^2ψ / dx^2) + V(x)ψ
cdouble schrodinger(double (*psi_x)(double), double x, double V_x, double E) {
    double d2psi_dx2 = (psi_x(x + 1e-5) - 2 * psi_x(x) + psi_x(x - 1e-5)) / (1e-5 * 1e-5);
    return (-h_bar * h_bar / (2 * m_e) * d2psi_dx2 + V_x * psi_x(x)) - E * psi_x(x);
}

// Probability Amplitude: |ψ(x,t)|^2
double probabilityAmplitude(cdouble psi) {
    return norm(psi);  // |ψ(x,t)|^2 = ψ*ψ
}

// Pauli Matrices for spin-1/2 particles
cdouble pauliX[2][2] = {{0, 1}, {1, 0}};  // σx
cdouble pauliY[2][2] = {{0, cdouble(0, -1)}, {cdouble(0, 1), 0}};  // σy
cdouble pauliZ[2][2] = {{1, 0}, {0, -1}};  // σz

// Commutator: [A, B] = AB - BA
cdouble commutator(cdouble A[2][2], cdouble B[2][2]) {
    cdouble AB[2][2], BA[2][2], C[2][2];

    // Compute AB
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            AB[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j];
            BA[i][j] = B[i][0] * A[0][j] + B[i][1] * A[1][j];
            C[i][j] = AB[i][j] - BA[i][j];
        }
    }
    return C[0][0] + C[0][1] + C[1][0] + C[1][1];  // Returning sum for simplicity
}

// Quantum Harmonic Oscillator: E_n = (n + 1/2)ħω
double quantumHarmonicOscillator(double n, double omega) {
    return (n + 0.5) * h_bar * omega;
}

// Particle in a 1D box: ψ(x) = sqrt(2/L) * sin(nπx / L), E_n = (n^2 * h_bar^2 * π^2) / (2mL^2)
double particleInABox(double x, double L, int n) {
    if (x < 0 || x > L) return 0;
    return sqrt(2 / L) * sin(n * pi * x / L);
}

// Expectation value: ⟨ψ|A|ψ⟩
double expectationValue(double (*func)(double), double a, double b, int n) {
    double sum = 0;
    double h = (b - a) / n;

    for (int i = 0; i <= n; i++) {
        double x = a + i * h;
        sum += func(x) * x;
    }

    return sum * h;
}
