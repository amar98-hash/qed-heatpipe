#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>

// ============================================================
// Bare-minimum 3-level master-equation test
//
// System:
//
//      |3>
//       |
//      Ω23
//       |
//      |2>
//       |
//      Ω12
//       |
//      |1>
//
// Master equation:
//
// dρ/dt = -i [H,ρ]
//         + D[C21]ρ
//         + D[C32]ρ
//
// Here we use hbar = 1.
//
// This test intentionally has:
//   - no spatial domain
//   - no Gaussian beam
//   - no laser pulse envelope
//   - no Doppler distribution
//   - no ASL fields
//
// First we verify that the 3-level quantum dynamics work.
// ============================================================

using Complex = std::complex<double>;
using Matrix3 = Eigen::Matrix3cd;

constexpr Complex I{0.0, 1.0};

// ============================================================
// Physical parameters
//
// All frequencies/rates must use the SAME inverse-time units.
//
// For this test they are dimensionless simulation units.
// Later we will use SI angular frequencies [rad/s].
// ============================================================

struct Parameters
{
    // Laser couplings
    double omega12 = 1.0;
    double omega23 = 0.8;

    // Laser detunings
    double delta1 = 0.0;
    double delta2 = 0.0;

    // Population decay rates
    double gamma21 = 0.08;
    double gamma32 = 0.05;
};


// ============================================================
// Three-level rotating-frame Hamiltonian
//
// H / hbar =
//
// [ 0          Ω12/2              0             ]
// [ Ω12/2     -Δ1                Ω23/2          ]
// [ 0          Ω23/2            -(Δ1 + Δ2)     ]
//
// Since hbar = 1 numerically, this matrix has units of
// angular frequency.
// ============================================================

Matrix3 hamiltonian(const Parameters& p)
{
    Matrix3 H = Matrix3::Zero();
    H(0, 1) = p.omega12 / 2.0;
    H(1, 0) = p.omega12 / 2.0;
    H(1, 1) = -p.delta1;
    H(1, 2) = p.omega23 / 2.0;
    H(2, 1) = p.omega23 / 2.0;
    H(2, 2) = -(p.delta1 + p.delta2);
    return H;
}


// ============================================================
// Lindblad dissipator
//
// D[C] rho =
//      C rho C†
//    - 1/2 C†C rho
//    - 1/2 rho C†C
// ============================================================

Matrix3 dissipator(const Matrix3& rho,const Matrix3& C)
{
    const Matrix3 CdagC =C.adjoint() * C;
    returnC * rho * C.adjoint()- 0.5 * CdagC * rho- 0.5 * rho * CdagC;
}

// ============================================================
// Master-equation RHS
// ============================================================

Matrix3 master_rhs(const Matrix3& rho,const Parameters& p)
{
    const Matrix3 H =hamiltonian(p);
    // --------------------------------------------------------
    // Collapse operator:
    //
    // |2> -> |1>
    //
    // C21 = sqrt(gamma21) |1><2|
    // --------------------------------------------------------
    Matrix3 C21 = Matrix3::Zero();
    C21(0, 1) =std::sqrt(p.gamma21);

    // --------------------------------------------------------
    // Collapse operator:
    //
    // |3> -> |2>
    //
    // C32 = sqrt(gamma32) |2><3|
    // --------------------------------------------------------

    Matrix3 C32 = Matrix3::Zero();
    C32(1, 2) =std::sqrt(p.gamma32);

    // --------------------------------------------------------
    // Coherent Hamiltonian evolution
    // --------------------------------------------------------
    Matrix3 drho =-I * (H * rho - rho * H);
    // --------------------------------------------------------
    // Dissipation
    // --------------------------------------------------------
    drho += dissipator(rho, C21);
    drho += dissipator(rho, C32);

    return drho;
}


// ============================================================
// One RK4 step
//
// Temporary prototype integrator.
//
// Later replace this with:
//
//      C++QED Master
//      + ODE_EngineBoost
//
// ============================================================

void rk4_step(Matrix3& rho,double dt,const Parameters& p)
{
    const Matrix3 k1 =master_rhs(rho, p);
    const Matrix3 k2 =master_rhs(rho + 0.5 * dt * k1,p);
    const Matrix3 k3 = master_rhs(rho + 0.5 * dt * k2,p);
    const Matrix3 k4 =master_rhs(rho + dt * k3,p);
    rho +=(dt / 6.0)*(k1+ 2.0 * k2+ 2.0 * k3+ k4);
}


// ============================================================
// Diagnostics
// ============================================================

double hermiticity_error(const Matrix3& rho)
{
    return(rho - rho.adjoint()).norm();
}


// ============================================================
// main
// ============================================================
int main()
{
    Parameters p;

    // --------------------------------------------------------
    // Initial state
    //
    // All population in |1>
    //
    // rho(0) =
    //
    // [1 0 0]
    // [0 0 0]
    // [0 0 0]
    // --------------------------------------------------------

    Matrix3 rho =Matrix3::Zero();
    rho(0, 0) = 1.0;


    // --------------------------------------------------------
    // Time integration parameters
    // --------------------------------------------------------

    constexpr double dt =1.0e-3;
    constexpr double tMax =100.0;
    constexpr int outputEvery =100;

    // --------------------------------------------------------
    // Output
    // --------------------------------------------------------
    std::ofstream out( "results/three_level_master.dat");
    if (!out)
    {
        std::cerr<< "Could not open output file.\n";
        return 1;
    }

    out<< std::scientific << std::setprecision(12);
    out<< "# t "<< "rho11 "<< "rho22 "<< "rho33 "<< "trace "<< "hermiticity_error\n";

    // --------------------------------------------------------
    // Evolution
    // --------------------------------------------------------
    const int nSteps =static_cast<int>(std::round(tMax / dt));

    for (int step = 0;step <= nSteps;++step)
    {
        const double t =step * dt;
        if (step % outputEvery == 0)
        {
            const double rho11 =std::real(rho(0, 0));
            const double rho22 =std::real(rho(1, 1));
            const double rho33 =std::real(rho(2, 2));
            const double tr =std::real(rho.trace());
            out<< t << ' '<< rho11 << ' '<< rho22 << ' '<< rho33 << ' '<< tr << ' '<< hermiticity_error(rho)<< '\n';
        }
        if (step == nSteps)  break;
        rk4_step( rho,dt,p);
    }


    // --------------------------------------------------------
    // Final diagnostics
    // --------------------------------------------------------
    Eigen::SelfAdjointEigenSolver<Matrix3> solver(rho);
    std::cout << "\nDensity-matrix eigenvalues:\n"<< solver.eigenvalues()<< "\n";
    std::cout<< std::setprecision(12);
    std::cout<< "\nFinal density matrix:\n\n"<< rho<< "\n\n";
    std::cout<< "Trace = "<< rho.trace()<< '\n';
    std::cout<< "Hermiticity error = "<< hermiticity_error(rho)<< '\n';
    std::cout<< "\nOutput written to:\n"<< "three_level_master.dat\n";

    return 0;
}