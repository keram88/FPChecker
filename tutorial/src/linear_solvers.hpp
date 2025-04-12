#ifndef LINEAR_SOLVERS_HPP
#define LINEAR_SOLVERS_HPP

#include <vector>

using namespace std;

tuple<vector<vector<double>>, vector<vector<double>>, vector<vector<double>>>
lu_factorization_partial_pivot(vector<vector<double>> A);
vector<double> solve_lu_system(const vector<vector<double>> &L,
                               const vector<vector<double>> &U,
                               const vector<vector<double>> &P,
                               const vector<double> &b);

#endif // LINEAR_SOLVERS_HPP