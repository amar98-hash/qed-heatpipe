#include "helper.hpp"
#include <iostream>
#include "beam_profile.hpp"

int main()
{
    const auto beam =
        helper::readMatrix2D(
            "data/beam_profile.dat"
        );

    helper::printMatrix2D(beam);

    laser::LaserProfile::Parameters p;

    p.wavelength = 600e-9;
    p.peakIntensity = 1.0e8;
    p.beamWaist = 1.0e-3;
    p.pulseFWHM = 10e-9;
    p.pulseCenter = 20e-9;

    laser::LaserProfile L1(p);

    std::array<double, 3> polarization =
    {
        0.3,  // Ex direction
        0.7,  // Ey direction
        0.0   // Ez direction
    };

    const auto Efield =
    L1.linearlyPolarizedEField(
        polarization,
        "data/beam_profile.dat",
        20e-9
    );

    helper::writeVectorField2D(
        "data/Efield_map.dat",
        Efield.Ex,
        Efield.Ey,
        Efield.Ez
    );

    return 0;
}