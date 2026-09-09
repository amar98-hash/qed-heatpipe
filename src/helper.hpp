#pragma once

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace helper
{

// ============================================================
// Simple 2D numerical array
// ============================================================

template <typename T>
struct Matrix2D
{
    std::vector<std::vector<T>> data;

    std::size_t rows = 0;
    std::size_t cols = 0;

    T& operator()(std::size_t row,
                  std::size_t col)
    {
        return data.at(row).at(col);
    }

    const T& operator()(std::size_t row,
                        std::size_t col) const
    {
        return data.at(row).at(col);
    }
};


// ============================================================
// Read a REAL-valued 2D matrix from text file
// ============================================================

inline Matrix2D<double> readMatrix2D(
    const std::string& filename)
{
    std::ifstream file(filename);

    if (!file)
    {
        throw std::runtime_error(
            "Could not open matrix file: "
            + filename
        );
    }

    Matrix2D<double> matrix;

    std::string line;

    std::size_t expectedCols = 0;

    while (std::getline(file, line))
    {
        // Remove comments
        const auto commentPosition =
            line.find('#');

        if (commentPosition != std::string::npos)
        {
            line =
                line.substr(
                    0,
                    commentPosition
                );
        }

        // Allow comma-separated values
        std::replace(
            line.begin(),
            line.end(),
            ',',
            ' '
        );

        std::stringstream stream(line);

        std::vector<double> row;

        double value;

        while (stream >> value)
        {
            row.push_back(value);
        }

        if (row.empty())
            continue;

        if (expectedCols == 0)
        {
            expectedCols =
                row.size();
        }
        else if (row.size() != expectedCols)
        {
            throw std::runtime_error(
                "Matrix file is not rectangular: "
                + filename
            );
        }

        matrix.data.push_back(
            std::move(row)
        );
    }

    if (matrix.data.empty())
    {
        throw std::runtime_error(
            "No numerical matrix data found in: "
            + filename
        );
    }

    matrix.rows =
        matrix.data.size();

    matrix.cols =
        matrix.data.front().size();

    return matrix;
}


// ============================================================
// Generic matrix writer
//
// Works for types that support operator<<
//
// Examples:
//     double
//     float
//     int
//     std::complex<double>
//
// ============================================================

template <typename T>
inline void writeMatrix2D(
    const std::string& filename,
    const Matrix2D<T>& matrix,
    int precision = 12)
{
    std::ofstream file(filename);

    if (!file)
    {
        throw std::runtime_error(
            "Could not open output file: "
            + filename
        );
    }

    file
        << std::scientific
        << std::setprecision(precision);

    for (std::size_t i = 0;
         i < matrix.rows;
         ++i)
    {
        for (std::size_t j = 0;
             j < matrix.cols;
             ++j)
        {
            file
                << matrix(i, j);

            if (j + 1 < matrix.cols)
            {
                file << ' ';
            }
        }

        file << '\n';
    }
}

// ============================================================
// Generic matrix printer
//
// Prints matrix to stdout.
//
// Works for:
//     double
//     float
//     int
//     std::complex<double>
// ============================================================

template <typename T>
inline void printMatrix2D(
    const Matrix2D<T>& matrix,
    int precision = 6)
{
    std::cout
        << std::scientific
        << std::setprecision(precision);

    std::cout
        << "Matrix: "
        << matrix.rows
        << " x "
        << matrix.cols
        << '\n';

    for (std::size_t i = 0;
         i < matrix.rows;
         ++i)
    {
        for (std::size_t j = 0;
             j < matrix.cols;
             ++j)
        {
            std::cout
                << matrix(i, j);

            if (j + 1 < matrix.cols)
            {
                std::cout << ' ';
            }
        }

        std::cout << '\n';
    }
}


template <typename T>
inline void writeVectorField2D(
    const std::string& filename,
    const Matrix2D<T>& Fx,
    const Matrix2D<T>& Fy,
    const Matrix2D<T>& Fz,
    int precision = 12)
{
    if (Fx.rows != Fy.rows ||
        Fx.rows != Fz.rows ||
        Fx.cols != Fy.cols ||
        Fx.cols != Fz.cols)
    {
        throw std::runtime_error(
            "Vector field matrices must have identical dimensions."
        );
    }

    std::ofstream file(filename);

    if (!file)
    {
        throw std::runtime_error(
            "Could not open output file: " + filename
        );
    }

    file
        << std::scientific
        << std::setprecision(precision);

    file << "# row col Fx Fy Fz\n";

    for (std::size_t i = 0; i < Fx.rows; ++i)
    {
        for (std::size_t j = 0; j < Fx.cols; ++j)
        {
            file
                << i << ' '
                << j << ' '
                << Fx(i, j) << ' '
                << Fy(i, j) << ' '
                << Fz(i, j)
                << '\n';
        }
    }
}

} // namespace helper