#ifndef LINEAR_SOLVERS_HPP
#define LINEAR_SOLVERS_HPP

#include <vector>

std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>>
lu_factorization_partial_pivot(std::vector<std::vector<double>> A);

std::vector<double> solve_lu_system_by_substitution(const std::vector<std::vector<double>> &L,
                                                    const std::vector<std::vector<double>> &U,
                                                    const std::vector<std::vector<double>> &P,
                                                    const std::vector<double> &b);

// Solve AX = b using LU factorization
std::vector<double> solve_system_with_LU(const std::vector<std::vector<double>> &A, const std::vector<double> &b);

#endif // LINEAR_SOLVERS_HPP