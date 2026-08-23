#include <vtkCleanPolyData.h>
#include <vtkCylinderSource.h>
#include <vtkFeatureEdges.h>
#include <vtkMassProperties.h>
#include <vtkPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

// ============================================================
// light-matter-heatpipe-tube
//
// Prototype domain model.
//
// Purpose:
//
//   1. Generate a closed cylindrical STL.
//   2. Load an arbitrary STL.
//   3. Gather basic cavity/domain geometry information.
//   4. Create Cartesian interior voxels.
//   5. Give every voxel storage for a future 3x3 density matrix.
//
// IMPORTANT:
//
//   This file does NOT evolve rho yet.
//   C++QED will be connected later.
//
// Coordinate convention:
//
//       cylinder axis = z
//
//       -L/2 <= z <= +L/2
//
//       x^2 + y^2 <= R^2
//
// Internal units:
//
//       meters
//
// ============================================================

namespace lmh
{

using Complex = std::complex<double>;


// ============================================================
// 3-level density matrix placeholder
// ============================================================
//
// Later this may be:
//   - replaced by a C++QED object,
//   - wrapped around a C++QED object,
//   - or retained as our fast fixed-size production backend.
//
// Storage is row-major:
//
// rho =
//
// [ 00 01 02 ]
// [ 10 11 12 ]
// [ 20 21 22 ]
//
// ============================================================

struct DensityMatrix3
{
    std::array<Complex, 9> data{};

    Complex& operator()(std::size_t row, std::size_t col)
    {
        return data[row * 3 + col];
    }

    const Complex& operator()(std::size_t row, std::size_t col) const
    {
        return data[row * 3 + col];
    }

    void setZero()
    {
        data.fill(Complex{0.0, 0.0});
    }

    // Prototype initial condition:
    //
    // |1><1|
    //
    // rho =
    //
    // [1 0 0]
    // [0 0 0]
    // [0 0 0]
    //
    void setInitialState1()
    {
        setZero();
        (*this)(0, 0) = Complex{1.0, 0.0};
    }

    Complex trace() const
    {
        return (*this)(0, 0)
             + (*this)(1, 1)
             + (*this)(2, 2);
    }
};


// ============================================================
// One discrete interaction-volume cell
// ============================================================

struct DomainCell
{
    std::uint64_t id = 0;

    // Structured-grid indices before the cylindrical/STL mask.
    std::size_t i = 0;
    std::size_t j = 0;
    std::size_t k = 0;

    // Position of voxel center [m].
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    // Nominal Cartesian voxel volume [m^3].
    double volume = 0.0;

    // --------------------------------------------------------
    // Future classical fields from ASL
    // --------------------------------------------------------

    double temperature = 0.0; // K

    double nAr  = 0.0; // m^-3
    double nLi  = 0.0; // m^-3
    double nLi2 = 0.0; // m^-3

    // --------------------------------------------------------
    // Future laser fields
    // --------------------------------------------------------

    double laser1Intensity = 0.0; // W/m^2
    double laser2Intensity = 0.0; // W/m^2

    // --------------------------------------------------------
    // Quantum state
    // --------------------------------------------------------

    DensityMatrix3 rho;
};


// ============================================================
// Geometry information obtained from STL
// ============================================================

struct CavityInfo
{
    // xmin xmax ymin ymax zmin zmax
    std::array<double, 6> bounds{};

    std::array<double, 3> center{};
    std::array<double, 3> dimensions{};

    double volume = 0.0;      // m^3
    double surfaceArea = 0.0; // m^2

    vtkIdType numberOfPoints = 0;
    vtkIdType numberOfTriangles = 0;

    bool closedAndManifold = false;

    // Number of boundary + non-manifold edge cells detected.
    vtkIdType problematicEdgeCount = 0;

    // Convenient values for our first z-axis cylindrical model.
    double radiusX = 0.0;
    double radiusY = 0.0;
    double approximateRadius = 0.0;
    double axialLengthZ = 0.0;
};


// ============================================================
// Result of voxelizing the STL interior
// ============================================================

struct VoxelDomain
{
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;

    std::vector<DomainCell> cells;

    double voxelVolume() const
    {
        return dx * dy * dz;
    }

    double totalVoxelizedVolume() const
    {
        return static_cast<double>(cells.size()) * voxelVolume();
    }
};


// ============================================================
// DomainModel
// ============================================================

class DomainModel
{
public:

    // ========================================================
    // Generate cylinder STL
    // ========================================================
    //
    // radius      [m]
    // length      [m]
    //
    // Cylinder is centered at (0,0,0) and points along z.
    //
    // VTK's vtkCylinderSource is oriented along y by default,
    // so we rotate it into the z direction.
    //
    // ========================================================

    static vtkSmartPointer<vtkPolyData>
    generateCylinderSTL(
        double radius,
        double length,
        const std::string& outputFile,
        int angularResolution = 256)
    {
        if (radius <= 0.0)
        {
            throw std::invalid_argument(
                "Cylinder radius must be > 0.");
        }

        if (length <= 0.0)
        {
            throw std::invalid_argument(
                "Cylinder length must be > 0.");
        }

        if (angularResolution < 8)
        {
            throw std::invalid_argument(
                "Cylinder angular resolution must be >= 8.");
        }

        // ----------------------------------------------------
        // Generate VTK cylinder
        // ----------------------------------------------------

        auto cylinder =
            vtkSmartPointer<vtkCylinderSource>::New();

        cylinder->SetRadius(radius);
        cylinder->SetHeight(length);
        cylinder->SetResolution(angularResolution);

        cylinder->SetCenter(0.0, 0.0, 0.0);

        // We need a CLOSED solid.
        cylinder->CappingOn();

        // ----------------------------------------------------
        // Rotate VTK cylinder from y-axis to z-axis
        // ----------------------------------------------------

        auto transform =
            vtkSmartPointer<vtkTransform>::New();

        transform->RotateX(90.0);

        auto transformFilter =
            vtkSmartPointer<vtkTransformPolyDataFilter>::New();

        transformFilter->SetTransform(transform);
        transformFilter->SetInputConnection(
            cylinder->GetOutputPort());

        // ----------------------------------------------------
        // Convert all faces to triangles
        // ----------------------------------------------------

        auto triangleFilter =
            vtkSmartPointer<vtkTriangleFilter>::New();

        triangleFilter->SetInputConnection(
            transformFilter->GetOutputPort());

        triangleFilter->PassLinesOff();
        triangleFilter->PassVertsOff();

        // ----------------------------------------------------
        // Merge coincident points
        // ----------------------------------------------------

        auto clean =
            vtkSmartPointer<vtkCleanPolyData>::New();

        clean->SetInputConnection(
            triangleFilter->GetOutputPort());

        clean->PointMergingOn();
        clean->Update();

        auto surface =
            vtkSmartPointer<vtkPolyData>::New();

        surface->DeepCopy(clean->GetOutput());

        // ----------------------------------------------------
        // Ensure destination directory exists
        // ----------------------------------------------------

        const std::filesystem::path path(outputFile);

        if (path.has_parent_path())
        {
            std::filesystem::create_directories(
                path.parent_path());
        }

        // ----------------------------------------------------
        // Write STL
        // ----------------------------------------------------

        auto writer =
            vtkSmartPointer<vtkSTLWriter>::New();

        writer->SetFileName(outputFile.c_str());
        writer->SetInputData(surface);

        // Binary is much smaller than ASCII STL.
        writer->SetFileTypeToBinary();

        if (writer->Write() == 0)
        {
            throw std::runtime_error(
                "Failed to write STL: " + outputFile);
        }

        std::cout
            << "\nGenerated cylinder STL:\n"
            << "  file   = " << outputFile << '\n'
            << "  radius = " << radius << " m\n"
            << "  length = " << length << " m\n";

        return surface;
    }


    // ========================================================
    // Load STL
    // ========================================================
    //
    // STL itself contains no physical unit metadata.
    //
    // scaleToMeters:
    //
    //     STL already in meters:
    //
    //         scaleToMeters = 1.0
    //
    //     STL coordinates in millimeters:
    //
    //         scaleToMeters = 1e-3
    //
    // ========================================================

    static vtkSmartPointer<vtkPolyData>
    loadSTL(
        const std::string& filename,
        double scaleToMeters = 1.0)
    {
        if (!std::filesystem::exists(filename))
        {
            throw std::runtime_error(
                "STL file does not exist: " + filename);
        }

        if (scaleToMeters <= 0.0)
        {
            throw std::invalid_argument(
                "scaleToMeters must be > 0.");
        }

        auto reader =
            vtkSmartPointer<vtkSTLReader>::New();

        reader->SetFileName(filename.c_str());

        // Merge duplicate STL vertices where possible.
        reader->MergingOn();

        reader->Update();

        if (reader->GetOutput()->GetNumberOfPoints() == 0)
        {
            throw std::runtime_error(
                "STL contains no points: " + filename);
        }

        // ----------------------------------------------------
        // Convert incoming coordinates to SI meters
        // ----------------------------------------------------

        auto transform =
            vtkSmartPointer<vtkTransform>::New();

        transform->Scale(
            scaleToMeters,
            scaleToMeters,
            scaleToMeters);

        auto transformFilter =
            vtkSmartPointer<vtkTransformPolyDataFilter>::New();

        transformFilter->SetTransform(transform);
        transformFilter->SetInputConnection(
            reader->GetOutputPort());

        // ----------------------------------------------------
        // Ensure triangle-only representation
        // ----------------------------------------------------

        auto triangleFilter =
            vtkSmartPointer<vtkTriangleFilter>::New();

        triangleFilter->SetInputConnection(
            transformFilter->GetOutputPort());

        triangleFilter->PassLinesOff();
        triangleFilter->PassVertsOff();

        // ----------------------------------------------------
        // Clean duplicate points
        // ----------------------------------------------------

        auto clean =
            vtkSmartPointer<vtkCleanPolyData>::New();

        clean->SetInputConnection(
            triangleFilter->GetOutputPort());

        clean->PointMergingOn();
        clean->Update();

        auto surface =
            vtkSmartPointer<vtkPolyData>::New();

        surface->DeepCopy(clean->GetOutput());

        std::cout
            << "\nLoaded STL:\n"
            << "  file = " << filename << '\n'
            << "  scale-to-meters = "
            << scaleToMeters << '\n';

        return surface;
    }


    // ========================================================
    // Extract information from STL
    // ========================================================

    static CavityInfo
    gatherCavityInfo(vtkPolyData* surface)
    {
        if (!surface)
        {
            throw std::invalid_argument(
                "gatherCavityInfo(): surface is null.");
        }

        if (surface->GetNumberOfPoints() == 0)
        {
            throw std::runtime_error(
                "gatherCavityInfo(): empty surface.");
        }

        CavityInfo info;

        // ----------------------------------------------------
        // Bounds
        // ----------------------------------------------------

        double b[6];
        surface->GetBounds(b);

        for (int q = 0; q < 6; ++q)
        {
            info.bounds[q] = b[q];
        }

        info.dimensions =
        {
            b[1] - b[0],
            b[3] - b[2],
            b[5] - b[4]
        };

        info.center =
        {
            0.5 * (b[0] + b[1]),
            0.5 * (b[2] + b[3]),
            0.5 * (b[4] + b[5])
        };

        info.radiusX =
            0.5 * info.dimensions[0];

        info.radiusY =
            0.5 * info.dimensions[1];

        info.approximateRadius =
            0.5 * (info.radiusX + info.radiusY);

        info.axialLengthZ =
            info.dimensions[2];

        // ----------------------------------------------------
        // Topology
        // ----------------------------------------------------

        info.numberOfPoints =
            surface->GetNumberOfPoints();

        info.numberOfTriangles =
            surface->GetNumberOfPolys();

        info.closedAndManifold =
            vtkSelectEnclosedPoints::
                IsSurfaceClosed(surface) > 0;

        // ----------------------------------------------------
        // Explicitly inspect bad edges
        // ----------------------------------------------------

        auto edges =
            vtkSmartPointer<vtkFeatureEdges>::New();

        edges->SetInputData(surface);

        edges->BoundaryEdgesOn();
        edges->NonManifoldEdgesOn();

        edges->FeatureEdgesOff();
        edges->ManifoldEdgesOff();

        edges->Update();

        info.problematicEdgeCount =
            edges->GetOutput()->GetNumberOfCells();

        // ----------------------------------------------------
        // Surface area + enclosed volume
        // ----------------------------------------------------

        auto mass =
            vtkSmartPointer<vtkMassProperties>::New();

        mass->SetInputData(surface);
        mass->Update();

        info.volume =
            mass->GetVolume();

        info.surfaceArea =
            mass->GetSurfaceArea();

        return info;
    }


    // ========================================================
    // Voxelize STL interior
    // ========================================================
    //
    // This is intentionally a simple Cartesian center-sampling
    // prototype.
    //
    // If the center of a voxel lies inside the STL, that voxel
    // becomes part of the quantum interaction domain.
    //
    // Later we can improve boundary voxels using:
    //
    //   - volume fractions,
    //   - signed-distance functions,
    //   - adaptive grids,
    //   - tetrahedral meshes.
    //
    // ========================================================

    static VoxelDomain
    voxelizeInterior(
        vtkPolyData* surface,
        double dx,
        double dy,
        double dz)
    {
        if (!surface)
        {
            throw std::invalid_argument(
                "voxelizeInterior(): surface is null.");
        }

        if (dx <= 0.0 ||
            dy <= 0.0 ||
            dz <= 0.0)
        {
            throw std::invalid_argument(
                "Voxel spacing must be > 0.");
        }

        if (vtkSelectEnclosedPoints::
                IsSurfaceClosed(surface) <= 0)
        {
            throw std::runtime_error(
                "Cannot voxelize STL: "
                "surface is not closed/manifold.");
        }

        double b[6];
        surface->GetBounds(b);

        VoxelDomain domain;

        domain.dx = dx;
        domain.dy = dy;
        domain.dz = dz;

        domain.nx =
            static_cast<std::size_t>(
                std::ceil((b[1] - b[0]) / dx));

        domain.ny =
            static_cast<std::size_t>(
                std::ceil((b[3] - b[2]) / dy));

        domain.nz =
            static_cast<std::size_t>(
                std::ceil((b[5] - b[4]) / dz));

        // ----------------------------------------------------
        // Point-inside-surface accelerator
        // ----------------------------------------------------

        auto inside =
            vtkSmartPointer<vtkSelectEnclosedPoints>::New();

        // Fraction of bounding-box diagonal.
        inside->SetTolerance(1.0e-7);

        // Builds the internal search structures once.
        inside->Initialize(surface);

        const double dV = dx * dy * dz;

        std::uint64_t id = 0;

        // ----------------------------------------------------
        // Sample voxel centers
        // ----------------------------------------------------

        for (std::size_t k = 0;
             k < domain.nz;
             ++k)
        {
            const double z =
                b[4] +
                (static_cast<double>(k) + 0.5) * dz;

            if (z > b[5])
            {
                continue;
            }

            for (std::size_t j = 0;
                 j < domain.ny;
                 ++j)
            {
                const double y =
                    b[2] +
                    (static_cast<double>(j) + 0.5) * dy;

                if (y > b[3])
                {
                    continue;
                }

                for (std::size_t i = 0;
                     i < domain.nx;
                     ++i)
                {
                    const double x =
                        b[0] +
                        (static_cast<double>(i) + 0.5) * dx;

                    if (x > b[1])
                    {
                        continue;
                    }

                    if (!inside->IsInsideSurface(x, y, z))
                    {
                        continue;
                    }

                    DomainCell cell;

                    cell.id = id++;

                    cell.i = i;
                    cell.j = j;
                    cell.k = k;

                    cell.x = x;
                    cell.y = y;
                    cell.z = z;

                    cell.volume = dV;

                    // ----------------------------------------
                    // Prototype quantum initialization
                    // ----------------------------------------
                    //
                    // rho_11 = 1
                    //
                    // Later this initial condition will be
                    // determined by Li2 state populations.
                    //
                    cell.rho.setInitialState1();

                    domain.cells.push_back(cell);
                }
            }
        }

        inside->Complete();

        return domain;
    }


    // ========================================================
    // Pretty-print cavity information
    // ========================================================

    static void printCavityInfo(
        const CavityInfo& info)
    {
        std::cout
            << "\n"
            << "========================================\n"
            << "       CAVITY / DOMAIN INFORMATION\n"
            << "========================================\n";

        std::cout
            << std::scientific
            << std::setprecision(8);

        std::cout
            << "Bounds [m]\n"
            << "  x : "
            << info.bounds[0]
            << "  ->  "
            << info.bounds[1]
            << '\n'

            << "  y : "
            << info.bounds[2]
            << "  ->  "
            << info.bounds[3]
            << '\n'

            << "  z : "
            << info.bounds[4]
            << "  ->  "
            << info.bounds[5]
            << '\n';

        std::cout
            << "\nCenter [m]\n"
            << "  ("
            << info.center[0] << ", "
            << info.center[1] << ", "
            << info.center[2] << ")\n";

        std::cout
            << "\nDimensions [m]\n"
            << "  dx = " << info.dimensions[0] << '\n'
            << "  dy = " << info.dimensions[1] << '\n'
            << "  dz = " << info.dimensions[2] << '\n';

        std::cout
            << "\nCylinder-like interpretation\n"
            << "  Rx       = "
            << info.radiusX << " m\n"

            << "  Ry       = "
            << info.radiusY << " m\n"

            << "  R approx = "
            << info.approximateRadius << " m\n"

            << "  Lz       = "
            << info.axialLengthZ << " m\n";

        std::cout
            << "\nMesh information\n"
            << "  Points      = "
            << info.numberOfPoints << '\n'

            << "  Triangles   = "
            << info.numberOfTriangles << '\n'

            << "  Bad edges   = "
            << info.problematicEdgeCount << '\n';

        std::cout
            << "\nTopology\n"
            << "  Closed/manifold = "
            << (info.closedAndManifold
                    ? "YES"
                    : "NO")
            << '\n';

        std::cout
            << "\nGeometric properties\n"
            << "  Volume       = "
            << info.volume << " m^3\n"

            << "  Surface area = "
            << info.surfaceArea << " m^2\n";

        std::cout
            << "========================================\n";
    }


    // ========================================================
    // Print voxelization information
    // ========================================================

    static void printVoxelInfo(
        const VoxelDomain& domain,
        double stlVolume)
    {
        const double voxelVol =
            domain.totalVoxelizedVolume();

        const double relativeError =
            (stlVolume > 0.0)
                ? std::abs(voxelVol - stlVolume)
                    / stlVolume
                : 0.0;

        std::cout
            << "\n"
            << "========================================\n"
            << "          VOXEL DOMAIN INFO\n"
            << "========================================\n";

        std::cout
            << std::scientific
            << std::setprecision(8);

        std::cout
            << "Spacing\n"
            << "  dx = " << domain.dx << " m\n"
            << "  dy = " << domain.dy << " m\n"
            << "  dz = " << domain.dz << " m\n";

        std::cout
            << "\nBounding-grid dimensions\n"
            << "  nx = " << domain.nx << '\n'
            << "  ny = " << domain.ny << '\n'
            << "  nz = " << domain.nz << '\n';

        std::cout
            << "\nInterior cells\n"
            << "  N = "
            << domain.cells.size()
            << '\n';

        std::cout
            << "\nVolumes\n"
            << "  STL volume   = "
            << stlVolume << " m^3\n"

            << "  voxel volume = "
            << voxelVol << " m^3\n"

            << "  relative difference = "
            << 100.0 * relativeError
            << " %\n";

        std::cout
            << "========================================\n";
    }
};

} // namespace lmh



// ============================================================
// First standalone domain test
// ============================================================
//
// For the prototype we intentionally keep main() here.
//
// Later:
//
//      domain_model.cc
//      domain_model.h
//      quantum_solver.cc
//      laser_model.cc
//      main.cc
//
// ============================================================

int main()
{
    try
    {
        // ----------------------------------------------------
        // Prototype geometry
        // ----------------------------------------------------
        //
        // R = 15 mm
        // L = 40 mm
        //
        // Cylinder extends:
        //
        //     z = -20 mm ... +20 mm
        //
        // ----------------------------------------------------

        constexpr double radius =
            15.0e-3;

        constexpr double length =
            40.0e-3;

        // 1 mm initial quantum-domain resolution.
        constexpr double spacing =
            1.0e-3;

        const std::string stlFile =
            "geometry/prototype_cylinder.stl";

        // ----------------------------------------------------
        // STEP 1
        //
        // Generate STL.
        // ----------------------------------------------------

        lmh::DomainModel::generateCylinderSTL(
            radius,
            length,
            stlFile,
            256);

        // ----------------------------------------------------
        // STEP 2
        //
        // Load STL again.
        //
        // This deliberately tests our STL loader rather than
        // simply reusing the in-memory generated surface.
        // ----------------------------------------------------

        auto cavity =
            lmh::DomainModel::loadSTL(
                stlFile,
                1.0); // already meters

        // ----------------------------------------------------
        // STEP 3
        //
        // Gather geometry.
        // ----------------------------------------------------

        const auto info =
            lmh::DomainModel::gatherCavityInfo(
                cavity);

        lmh::DomainModel::printCavityInfo(
            info);

        // ----------------------------------------------------
        // Analytic check for our prototype cylinder
        // ----------------------------------------------------

        const double analyticVolume =
            std::numbers::pi
            * radius
            * radius
            * length;

        std::cout
            << "\nAnalytic cylinder volume = "
            << std::scientific
            << analyticVolume
            << " m^3\n";

        std::cout
            << "STL volume error = "
            << 100.0
                * std::abs(
                    info.volume
                    - analyticVolume)
                / analyticVolume
            << " %\n";

        // ----------------------------------------------------
        // STEP 4
        //
        // Fill STL interior with voxels.
        // ----------------------------------------------------

        auto domain =
            lmh::DomainModel::voxelizeInterior(
                cavity,
                spacing,
                spacing,
                spacing);

        lmh::DomainModel::printVoxelInfo(
            domain,
            info.volume);

        // ----------------------------------------------------
        // STEP 5
        //
        // Inspect a few cells.
        // ----------------------------------------------------

        std::cout
            << "\nFirst few interaction cells:\n";

        const std::size_t nPrint =
            std::min<std::size_t>(
                5,
                domain.cells.size());

        for (std::size_t q = 0;
             q < nPrint;
             ++q)
        {
            const auto& c =
                domain.cells[q];

            std::cout
                << "cell "
                << c.id
                << "  r = ("
                << c.x << ", "
                << c.y << ", "
                << c.z << ") m"
                << "  Tr(rho) = "
                << c.rho.trace()
                << '\n';
        }

        std::cout
            << "\nDomain prototype completed successfully.\n";

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr
            << "\nERROR: "
            << e.what()
            << '\n';

        return 1;
    }
}