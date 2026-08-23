#include "Composite.h"
#include "JaynesCummings.h"
#include "Mode.h"
#include "Qbit.h"
#include "QuantumJumpMonteCarlo.h"
#include "Pars.h"
#include "pcg_random.hpp"

using namespace cppqedutils;
using namespace quantumdata;
using namespace quantumtrajectory;
using namespace structure;

int main(int argc, const char* argv[])
{
    auto op = optionParser();

    trajectory::Pars<
        qjmc::Pars<
            qjmc::Algorithm::integrating,
            pcg64
        >
    > pt{op};

    qbit::Pars pq{op, "q"};
    mode::Pars pm{op, "m"};

    double g;

    add(
        op,
        parameters::_(
            "g",
            "Jaynes-Cummings coupling",
            1.0,
            g
        )
    );

    parse(op, argc, argv);

    StateVector<2> psi{
        qbit::init(pq) * mode::init(pm)
    };

    Broadcaster<2, 0> bQ{psi.extents};
    Broadcaster<2, 1> bM{psi.extents};

    auto [cutoffQ, haQ, liQ, evQ, descrQ]
        = qbit::make(pq);

    auto [cutoffM, haM, liM, evM, descrM]
        = mode::make(pm);

    auto ha = hana::make_tuple(
        bQ(haQ),
        bM(haM),
        jaynescummings::hamiltonian<0, 1>(
            pm.cutoff,
            g
        )
    );

    auto li = bQ(liQ);

    for (auto& l : bM(liM))
        li.push_back(std::move(l));

    quantumtrajectory::run<
        qjmc::Algorithm::integrating,
        ODE_EngineBoost
    >(
        ha,
        li,
        hana::make_tuple(
            bQ(evQ),
            bM(evM)
        ),
        json::object{
            {"system", "Jaynes-Cummings"},
            {"qubit", descrQ},
            {"mode", descrM},
            {"g", g}
        },
        std::move(psi),
        pt,
        trajectory::observerNoOp
    );

    return 0;
}