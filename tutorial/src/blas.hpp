#ifndef BLAS_HPP
#define BLAS_HPP

#include <vector>

void print_matrix_simple(const std::vector<std::vector<double>> &matrix);
void print_matrix(const std::vector<std::vector<double>> &matrix);
double frobenius_norm(const std::vector<std::vector<double>> &matrix);
double infinity_norm(const std::vector<std::vector<double>> &matrix);
double l2_norm(const std::vector<double> &v);
std::vector<double> multiply_matrix_vector(const std::vector<std::vector<double>> &A, const std::vector<double> &x);
std::vector<double> subtract_vectors(const std::vector<double> &x, const std::vector<double> &y);
std::vector<std::vector<double>> subtract_matrices(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);
std::vector<std::vector<double>> matrix_multiply(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);
std::vector<std::vector<double>> transpose_matrix(const std::vector<std::vector<double>> &matrix);

#endif // BLAS_HPP