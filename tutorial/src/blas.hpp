#ifndef BLAS_HPP
#define BLAS_HPP

#include <vector>

using namespace std;

void print_matrix(const vector<vector<double>> &matrix);
double frobenius_norm(const vector<vector<double>> &matrix);
double infinity_norm(const vector<vector<double>> &matrix);
double l2_norm(const vector<double> &v);
vector<double> multiply_matrix_vector(const vector<vector<double>> &A, const vector<double> &x);
vector<double> subtract_vectors(const vector<double> &x, const vector<double> &y);
vector<vector<double>> subtract_matrices(const vector<vector<double>> &A, const vector<vector<double>> &B);
vector<vector<double>> matrix_multiply(const vector<vector<double>> &A, const vector<vector<double>> &B);

#endif // BLAS_HPP