#pragma once

#include <cmath>
#include <stdexcept>
#include <array>
#include <string>

#include "helper.hpp"

namespace laser
{

// ============================================================
// Physical constants
// ============================================================

inline constexpr double c0 =299792458.0; // m/s

inline constexpr double epsilon0 =8.8541878128e-12; // F/m

inline constexpr double pi =3.14159265358979323846;


// ============================================================
// LaserProfile
//
// Simplified pulsed Gaussian laser propagating along +z.
//
// Coordinates:
//      z       = propagation direction
//      r       = radial distance from beam axis
//      phi     = ignored for cylindrically symmetric TEM00 beam
//
// Spatial field envelope:
//
//      E(r) = E0 exp(-r^2 / w^2)
//
// Temporal field envelope:
//
//      f(t) = exp[-2 ln(2) (t-t0)^2 / tau^2]
//
// where tau is the INTENSITY FWHM.
//
// Therefore:
//
//      E_envelope(r,t)
//          = E0 exp(-r^2/w^2) f(t)
//
// For RWA calculations, use fieldEnvelope().
//
// ============================================================

struct EFieldMap2D
{
    helper::Matrix2D<double> Ex;
    helper::Matrix2D<double> Ey;
    helper::Matrix2D<double> Ez;

    std::size_t rows = 0;
    std::size_t cols = 0;
};


class LaserProfile
{
public:

    struct Parameters
    {
        // ----------------------------------------------------
        // Laser wavelength [m]
        // ----------------------------------------------------
        double wavelength = 600e-9;

        // ----------------------------------------------------
        // Peak intensity at beam center [W/m^2]
        // ----------------------------------------------------
        double peakIntensity = 1.0;

        // ----------------------------------------------------
        // Gaussian beam radius [m]
        //
        // Intensity:
        //      I(r) = I0 exp(-2 r^2 / w^2)
        //
        // Electric-field amplitude:
        //      E(r) = E0 exp(-r^2 / w^2)
        // ----------------------------------------------------
        double beamWaist = 1.0e-3;

        // ----------------------------------------------------
        // Intensity pulse FWHM [s]
        // ----------------------------------------------------
        double pulseFWHM = 10.0e-9;

        // ----------------------------------------------------
        // Time of pulse maximum [s]
        // ----------------------------------------------------
        double pulseCenter = 0.0;

        // ----------------------------------------------------
        // Carrier phase [rad]
        // ----------------------------------------------------
        double carrierPhase = 0.0;
    };


    explicit LaserProfile(const Parameters& p): p_(p)
    {
        validate();
    }


    // ========================================================
    // Fundamental laser quantities
    // ========================================================

    double wavelength() const
    {
        return p_.wavelength;
    }


    // Angular frequency:
    //
    //      omega = 2 pi c / lambda
    //
    double angularFrequency() const
    {
        return 2.0 * pi * c0/ p_.wavelength;
    }


    // Wave number:
    //
    //      k = 2 pi / lambda
    //
    double waveNumber() const
    {
        return 2.0 * pi/ p_.wavelength;
    }

    // ========================================================
    // Peak electric-field amplitude
    //
    // For a plane-wave-like optical field:
    //
    //      I = 1/2 c epsilon0 E0^2
    //
    // therefore
    //
    //      E0 = sqrt(2 I / c epsilon0)
    //
    // ========================================================

    double peakElectricField() const
    {
        return std::sqrt(2.0 * p_.peakIntensity/(c0 * epsilon0));
    }


    // ========================================================
    // Radial Gaussian field profile
    //
    //      exp(-r^2 / w^2)
    //
    // This is the FIELD amplitude profile.
    //
    // The corresponding intensity is:
    //
    //      I(r) = I0 exp(-2 r^2 / w^2)
    //
    // ========================================================

    double radialEnvelope(double r) const
    {
        const double q =r / p_.beamWaist;

        return std::exp(-q * q);
    }


    // ========================================================
    // Temporal Gaussian FIELD envelope
    //
    // tau = intensity FWHM
    //
    // Since intensity is proportional to E^2:
    //
    //      E(t)
    //      = E0 exp[
    //          -2 ln(2) (t-t0)^2 / tau^2
    //        ]
    //
    // and therefore
    //
    //      I(t)
    //      = I0 exp[
    //          -4 ln(2) (t-t0)^2 / tau^2
    //        ]
    //
    // ========================================================

    double temporalEnvelope(double t) const
    {
        const double q =(t - p_.pulseCenter)/p_.pulseFWHM;

        return std::exp(-2.0* std::log(2.0)* q * q);
    }


    // ========================================================
    // Slowly varying electric-field envelope
    //
    //      E_env(r,t)
    //
    // This is normally what you use when computing
    // the Rabi frequency in the rotating-wave approximation:
    //
    //      Omega(r,t) = d E_env(r,t) / hbar
    //
    // ========================================================

    double fieldEnvelope(double r,double t) const
    {
        return peakElectricField()*radialEnvelope(r)*temporalEnvelope(t);
    }


    // ========================================================
    // Same field envelope using Cartesian coordinates.
    //
    // Laser propagates along z:
    //
    //      r = sqrt(x^2 + y^2)
    //
    // z does not affect the envelope in this simplified model.
    // ========================================================

    double fieldEnvelopeXYZ(double x,double y,double /*z*/,double t) const
    {
        const double r =std::sqrt(x * x+y * y);

        return fieldEnvelope(r, t);
    }


    // ========================================================
    // Local intensity
    //
    //      I(r,t)
    //      = 1/2 c epsilon0 E_env^2
    //
    // ========================================================

    double intensity(double r,double t) const
    {
        const double E =fieldEnvelope(r, t);

        return 0.5* c0* epsilon0* E* E;
    }


    double intensityXYZ(double x,double y,double z,double t) const
    {
        const double E =
            fieldEnvelopeXYZ(x,y,z,t);

        return 0.5* c0* epsilon0* E* E;
    }


    // ========================================================
    // Full real optical electric field
    //
    //      E(r,z,t)
    //      =
    //      E_env(r,t)
    //      cos(kz - omega t + phi0)
    //
    // Usually DO NOT use this in the RWA master equation.
    //
    // It is included mainly for diagnostics / visualization.
    //
    // ========================================================

    double fullField(double r,double z,double t) const
    {
        const double phase =waveNumber() * z-angularFrequency() * t+p_.carrierPhase;

        return fieldEnvelope(r, t)*std::cos(phase);
    }


    double fullFieldXYZ(double x,double y,double z,double t) const
    {
        const double r =std::sqrt(x * x+y * y);
        return fullField(r,z,t);
    }


    // ========================================================
    // Linearly polarized electric-field map from measured
    // beam intensity profile.
    //
    // Assumptions:
    //
    // 1. Input file contains a NORMALIZED intensity map:
    //
    //        I_norm(x,y) = I(x,y) / I_max
    //
    //    typically:
    //
    //        0 <= I_norm <= 1
    //
    // 2. p_.peakIntensity is the actual peak intensity
    //    of the laser in W/m^2.
    //
    // 3. polVector gives the linear polarization direction:
    //
    //        {1,0,0} -> x polarized
    //        {0,1,0} -> y polarized
    //        {0,0,1} -> z polarized
    //
    //    The vector is automatically normalized.
    //
    // 4. Temporal pulse envelope is included.
    //
    // Result:
    //
    //      E(x,y,t)
    //          = epsilon_hat
    //            E0
    //            sqrt(I_norm(x,y))
    //            f(t)
    //
    // ========================================================
    
    EFieldMap2D linearlyPolarizedEField(
        const std::array<double, 3>& polVector,
        const std::string& intensityMapFile,
        double t) const
    {
        // ----------------------------------------------------
        // Read measured/normalized intensity map
        // ----------------------------------------------------
    
        const auto I_map =
            helper::readMatrix2D(
                intensityMapFile
            );
        
        
        // ----------------------------------------------------
        // Normalize polarization vector
        // ----------------------------------------------------
        
        const double polMagnitude =
            std::sqrt(
                polVector[0] * polVector[0]
                +
                polVector[1] * polVector[1]
                +
                polVector[2] * polVector[2]
            );
        
        if (polMagnitude <= 0.0)
        {
            throw std::invalid_argument(
                "Polarization vector cannot be zero."
            );
        }
    
    
        const double px =
            polVector[0] / polMagnitude;
    
        const double py =
            polVector[1] / polMagnitude;
    
        const double pz =
            polVector[2] / polMagnitude;
    
    
        // ----------------------------------------------------
        // Allocate vector-field maps
        // ----------------------------------------------------
    
        EFieldMap2D field;
    
        field.rows =
            I_map.rows;
    
        field.cols =
            I_map.cols;
    
    
        field.Ex.rows = field.rows;
        field.Ex.cols = field.cols;
    
        field.Ey.rows = field.rows;
        field.Ey.cols = field.cols;
    
        field.Ez.rows = field.rows;
        field.Ez.cols = field.cols;
    
    
        field.Ex.data.assign(
            field.rows,
            std::vector<double>(
                field.cols,
                0.0
            )
        );
    
        field.Ey.data.assign(
            field.rows,
            std::vector<double>(
                field.cols,
                0.0
            )
        );
    
        field.Ez.data.assign(
            field.rows,
            std::vector<double>(
                field.cols,
                0.0
            )
        );
    
    
        // ----------------------------------------------------
        // Peak electric-field amplitude
        //
        //      E0 = sqrt(2 Imax / c epsilon0)
        // ----------------------------------------------------
    
        const double E0 =
            peakElectricField();
    
    
        // ----------------------------------------------------
        // Temporal field envelope
        // ----------------------------------------------------
    
        const double temporal =
            temporalEnvelope(t);
    
    
        // ----------------------------------------------------
        // Convert intensity map -> electric-field map
        // ----------------------------------------------------
    
        for (std::size_t i = 0;
             i < field.rows;
             ++i)
        {
            for (std::size_t j = 0;
                 j < field.cols;
                 ++j)
            {
                const double normalizedIntensity =
                    I_map(i, j);
            
            
                if (normalizedIntensity < 0.0)
                {
                    throw std::runtime_error(
                        "Intensity map contains negative values."
                    );
                }
            
            
                // --------------------------------------------
                // Electric-field amplitude scales as sqrt(I)
                // --------------------------------------------
            
                const double E =
                    E0
                    *
                    std::sqrt(
                        normalizedIntensity
                    )
                    *
                    temporal;
                
                
                // --------------------------------------------
                // Apply polarization direction
                // --------------------------------------------
                
                field.Ex(i, j) =
                    px * E;
                
                field.Ey(i, j) =
                    py * E;
                
                field.Ez(i, j) =
                    pz * E;
            }
        }
    
    
        return field;
    }


    // ========================================================
    // Access to parameters
    // ========================================================

    const Parameters& parameters() const
    {
        return p_;
    }


private:

    Parameters p_;


    void validate() const
    {
        if (p_.wavelength <= 0.0)
            throw std::invalid_argument(
                "Laser wavelength must be > 0."
            );

        if (p_.peakIntensity < 0.0)
            throw std::invalid_argument(
                "Laser peak intensity must be >= 0."
            );

        if (p_.beamWaist <= 0.0)
            throw std::invalid_argument(
                "Laser beam waist must be > 0."
            );

        if (p_.pulseFWHM <= 0.0)
            throw std::invalid_argument(
                "Laser pulse FWHM must be > 0."
            );
    }
};

} // namespace laser