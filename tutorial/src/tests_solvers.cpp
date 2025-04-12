#include <iostream>
#include <vector>

#include "linear_solvers.hpp"
#include "blas.hpp"

using namespace std;

int main()
{
    // Example matrix 1
    vector<vector<double>> A1 =
        {{1, 2, 3},
         {4, 5, 6},
         {7, 8, 9}};

    vector<vector<double>> A2 =
        {{0.4, 3, 0, 0.2},
         {0.1, 17, 1, 0.1},
         {21, 0.3, 0.5, 0.1},
         {1, 1, 1, 0.1}};

    vector<vector<double>> A3 =
        {{2, -1, 0, 0},
         {-1, 2, -1, 0},
         {0, -1, 2, -1},
         {1, 0, -1, 2}};

    vector<vector<vector<double>>> matrix_list = {A1, A2, A3};

    for (auto A : matrix_list)
    {
        auto [L, U, P] = lu_factorization_partial_pivot(A);

        cout << "Matrix A:" << endl;
        print_matrix(A);
        cout << "Matrix L:" << endl;
        print_matrix(L);
        cout << "Matrix U:" << endl;
        print_matrix(U);
        cout << "Matrix P:" << endl;
        print_matrix(P);

        vector<vector<double>> LU = matrix_multiply(L, U);
        cout << "Matrix LU:" << endl;
        print_matrix(LU);

        vector<vector<double>> PA = matrix_multiply(P, A);
        cout << "Matrix PA:" << endl;
        print_matrix(PA);

        vector<vector<double>> difference = subtract_matrices(PA, LU);
        double frobenius_norm_diff = frobenius_norm(difference);
        cout << ">>> Frobenius norm of the difference PA-LU: " << frobenius_norm_diff << endl;

        // Solve inear system Ax = 1
        vector<double> b(A.size(), 1.0);
        vector<double> x = solve_lu_system(L, U, P, b);
        for (size_t i = 0; i < x.size(); ++i)
            cout << "x[" << i << "]: " << x[i] << endl;

        // Residual
        auto a = multiply_matrix_vector(A, x);
        auto residual = subtract_vectors(a, b);
        auto norm = l2_norm(residual);
        cout << "Residual norm: " << norm << endl;
    }

    return 0;
}