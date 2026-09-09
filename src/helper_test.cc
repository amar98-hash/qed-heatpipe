#include "helper.hpp"

int main()
{
    const auto beam =
        helper::readMatrix2D(
            "data/beam_profile.dat"
        );

    helper::printMatrix2D(beam);

    return 0;
}