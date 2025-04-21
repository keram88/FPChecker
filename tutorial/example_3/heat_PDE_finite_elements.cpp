#include <vector>
#include <cmath>
#include <map>
#include <array>
#include <iostream>
#include <iomanip>   // For printing matrix with precision
#include <algorithm> // For std::find
#include <fstream>

#include "../common/linear_solvers.hpp"

// Structure to hold the mesh data
struct MeshData
{
    std::vector<std::array<double, 2>> points;
    std::vector<std::array<int, 3>> triangles;
};

// Helper function to generate linearly spaced values (equivalent to numpy.linspace)
std::vector<double> linspace(double start, double end, int num)
{
    std::vector<double> linspaced;
    if (num == 0)
    {
        return linspaced;
    }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }
    double delta = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    return linspaced;
}

MeshData generate_L_mesh(int N)
{
    MeshData mesh;

    // Generate grid points
    std::vector<double> x = linspace(0, 1, N);
    std::vector<double> y = linspace(0, 1, N);

    std::map<std::pair<int, int>, int> point_indices;
    int index = 0;

    // Create points only in the "L" shape (exclude upper-right region)
    for (int j = 0; j < N; ++j)
    {
        for (int i = 0; i < N; ++i)
        {
            // Exclude upper-right region
            if (!(x[i] > 0.5 && y[j] > 0.5))
            {
                mesh.points.push_back({x[i], y[j]});
                point_indices[{i, j}] = index;
                index++;
            }
        }
    }

    // Create triangles
    for (int j = 0; j < N - 1; ++j)
    {
        for (int i = 0; i < N - 1; ++i)
        {
            // Skip cells in the upper-right corner or those adjacent to it in the excluded region
            if ((x[i] > 0.5 && y[j] > 0.5) || (x[i + 1] > 0.5 && y[j] > 0.5) ||
                (x[i] > 0.5 && y[j + 1] > 0.5) || (x[i + 1] > 0.5 && y[j + 1] > 0.5))
            {
                continue;
            }

            // Get the indices of the four points in the current cell
            auto it1 = point_indices.find({i, j});
            auto it2 = point_indices.find({i + 1, j});
            auto it3 = point_indices.find({i, j + 1});
            auto it4 = point_indices.find({i + 1, j + 1});

            // Check if all four points of the quadrilateral exist in the L-shape
            if (it1 != point_indices.end() && it2 != point_indices.end() &&
                it3 != point_indices.end() && it4 != point_indices.end())
            {

                int p1 = it1->second; // Bottom-left
                int p2 = it2->second; // Bottom-right
                int p3 = it3->second; // Top-left
                int p4 = it4->second; // Top-right

                // Create lower-left triangle (p1, p2, p3) - counter-clockwise
                mesh.triangles.push_back({p1, p2, p3});
                // Create upper-right triangle (p2, p4, p3) - counter-clockwise
                mesh.triangles.push_back({p2, p4, p3});
            }
        }
    }

    return mesh;
}

// Checks areas are positive
double signed_area(double x_1, double y_1, double x_2, double y_2, double x_3, double y_3)
{
    double A = 0.5 * (x_1 * (y_2 - y_3) + x_2 * (y_3 - y_1) + x_3 * (y_1 - y_2));
    return A;
}

// Check if the resulting matrix is symmetric
bool is_symmetric(const std::vector<std::vector<double>> &matrix, double tolerance = 1e-9)
{
    // Check if the matrix is a square matrix
    if (matrix.empty() || matrix.size() != matrix[0].size())
    {
        return false;
    }
    // Check if the matrix is equal to its transpose within a tolerance
    size_t n = matrix.size();
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            if (std::abs(matrix[i][j] - matrix[j][i]) > tolerance)
            {
                return false;
            }
        }
    }
    return true;
}

// Function to save a double vector to a file
// Each element is written on a new line.
bool save_temperature_vector(const std::vector<double> &T, const std::string &filename)
{
    std::ofstream outfile(filename);

    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    // Set precision for writing double values
    outfile << std::fixed << std::setprecision(10);

    for (double value : T)
    {
        outfile << value << std::endl;
    }

    outfile.close();
    return true;
}

// Function to save coordinates (vector of arrays of 2 doubles) to a file
bool save_coords(const std::vector<std::array<double, 2>> &coords, const std::string &filename)
{
    std::ofstream outfile(filename);

    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    // Optional: Write the number of points and dimension on the first line
    // outfile << coords.size() << " " << 2 << std::endl;

    // Set precision for writing double values
    outfile << std::fixed << std::setprecision(10);

    for (const auto &point : coords)
    {
        outfile << point[0] << " " << point[1] << std::endl;
    }

    outfile.close();
    return true;
}

// Function to save triangles (vector of arrays of 3 ints) to a file
bool save_triangles(const std::vector<std::array<int, 3>> &triangles, const std::string &filename)
{
    std::ofstream outfile(filename);

    if (!outfile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    // Optional: Write the number of triangles and vertices per triangle on the first line
    // outfile << triangles.size() << " " << 3 << std::endl;

    for (const auto &triangle : triangles)
    {
        outfile << triangle[0] << " " << triangle[1] << " " << triangle[2] << std::endl;
    }

    outfile.close();
    return true;
}

int main()
{
    int N = 20; // Example value for N
    MeshData mesh = generate_L_mesh(N);

    // Print points
    std::cout << "Points:" << std::endl;
    for (const auto &point : mesh.points)
    {
        std::cout << "(" << point[0] << ", " << point[1] << ")" << std::endl;
    }

    // Print triangles
    std::cout << "\nTriangles:" << std::endl;
    for (const auto &triangle : mesh.triangles)
    {
        std::cout << "[" << triangle[0] << ", " << triangle[1] << ", " << triangle[2] << "]" << std::endl;
    }

    const auto &coords = mesh.points;
    const auto &triangles = mesh.triangles;

    // For large-scale FEM, a sparse matrix library (e.g., Eigen, SuiteSparse) would be more efficient.
    size_t num_nodes = coords.size();
    std::vector<std::vector<double>> K_global(num_nodes, std::vector<double>(num_nodes, 0.0));
    std::vector<double> f_global(num_nodes, 0.0); // f_global is initialized to zeros

    int nodes_per_elem = 3; // For triangles

    for (const auto &tria_indices : triangles)
    {
        int n1 = tria_indices[0];
        int n2 = tria_indices[1];
        int n3 = tria_indices[2];

        const auto &coord_1 = coords[n1];
        const auto &coord_2 = coords[n2];
        const auto &coord_3 = coords[n3];

        double x_1 = coord_1[0];
        double y_1 = coord_1[1];
        double x_2 = coord_2[0];
        double y_2 = coord_2[1];
        double x_3 = coord_3[0];
        double y_3 = coord_3[1];

        double A = signed_area(x_1, y_1, x_2, y_2, x_3, y_3);
        if (A <= 0)
        {
            // Handle error: non-positive area.
            // This might indicate incorrect node ordering or mesh issues.
            std::cerr << "Error: Non-positive signed area for a triangle." << std::endl;
            exit(1);
        }

        // Compute [B^e]
        std::vector<std::vector<double>> B(2, std::vector<double>(3, 0.0));
        double c = 1.0 / (2.0 * A);
        B[0][0] = y_2 - y_3;
        B[0][1] = y_3 - y_1;
        B[0][2] = y_1 - y_2;
        B[1][0] = x_3 - x_2;
        B[1][1] = x_1 - x_3;
        B[1][2] = x_2 - x_1;

        for (int row = 0; row < 2; ++row)
        {
            for (int col = 0; col < 3; ++col)
            {
                B[row][col] *= c;
            }
        }

        // Compute [L^e]
        // L_T is a matrix that maps local element degrees of freedom to global degrees of freedom.
        std::vector<std::vector<double>>
            L_T(num_nodes, std::vector<double>(nodes_per_elem, 0.0));
        L_T[n1][0] = 1.0;
        L_T[n2][1] = 1.0;
        L_T[n3][2] = 1.0;

        // L is the transpose of L_T
        std::vector<std::vector<double>> L(nodes_per_elem, std::vector<double>(num_nodes, 0.0));
        for (size_t i = 0; i < num_nodes; ++i)
        {
            for (int j = 0; j < nodes_per_elem; ++j)
            {
                L[j][i] = L_T[i][j];
            }
        }

        // Compute [k^e]
        // k = B.T * B * A
        std::vector<std::vector<double>> B_T(3, std::vector<double>(2, 0.0));
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                B_T[j][i] = B[i][j];
            }
        }

        std::vector<std::vector<double>> B_T_B(3, std::vector<double>(3, 0.0));
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int l = 0; l < 2; ++l)
                {
                    B_T_B[i][j] += B_T[i][l] * B[l][j];
                }
            }
        }

        std::vector<std::vector<double>> k(3, std::vector<double>(3, 0.0));
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                k[i][j] = B_T_B[i][j] * A;
            }
        }

        // Compute [L^e]^T [k^e] [L^e] and add to K_global
        // This is the assembly process. The element stiffness matrix k (3x3)
        // contributes to the global stiffness matrix K_global at the indices
        // corresponding to the nodes of the current triangle (n1, n2, n3).
        std::array<int, 3> global_indices = {n1, n2, n3};
        for (int i = 0; i < nodes_per_elem; ++i)
        {
            for (int j = 0; j < nodes_per_elem; ++j)
            {
                K_global[global_indices[i]][global_indices[j]] += k[i][j];
            }
        }

        // In this specific code, f_global is all zeros, so no force assembly is shown.
        // If there were body forces or boundary conditions contributing to f_global,
        // the corresponding terms would be calculated here and added to f_global.
    }

    // Optional: Check if K_global is symmetric (for verification)
    if (!is_symmetric(K_global))
    {
        std::cerr << "Warning: K_global matrix is not symmetric." << std::endl;
    }

    // Example: Print a portion of K_global and f_global
    std::cout << std::fixed << std::setprecision(6); // For formatted output

    // std::cout << K_global.size() << " x " << K_global[0].size() << std::endl;

    std::cout << "\nK_global:" << std::endl;
    for (size_t i = 0; i < K_global.size(); ++i)
    {
        for (size_t j = 0; j < K_global[i].size(); ++j)
        {
            std::cout << K_global[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "\nf_global:" << std::endl;
    for (size_t i = 0; i < f_global.size(); ++i)
    {
        std::cout << f_global[i] << std::endl;
    }

    std::map<int, double> T_known;
    // Boundary Conditions:
    double BC_temp = 100.0; // boundary condition temperature

    // Apply boundary conditions on the left side (nodes where original grid i % N == 0)
    // Need to iterate through original grid indices to identify these nodes
    for (const auto &triangle : mesh.triangles)
    {
        for (size_t i = 0; i < triangle.size(); ++i)
        {
            if (triangle[i] % N == 0) // on the side
            {
                T_known[triangle[i]] = BC_temp;
            }
        }
    }

    // Change the boundary to be ascending
    // Constant to multiply the keys by
    double constant = 3.0;
    // Sorting the keys and updating values in T_known
    std::vector<int> sorted_keys;
    for (const auto &pair : T_known)
    {
        sorted_keys.push_back(pair.first);
    }
    std::sort(sorted_keys.begin(), sorted_keys.end());

    for (int key : sorted_keys)
    {
        T_known[key] = static_cast<double>(key) * constant;
    }

    // Add some boundary conditions on the right side (where x == 1.0)
    double BC_right_side = 1.0;
    for (size_t node_index = 0; node_index < coords.size(); ++node_index)
    {
        if (coords[node_index][0] == 1.0)
        { // Right side
            T_known[node_index] = BC_right_side;
        }
    }

    // Applying Boundary Conditions (Direct Substitution)
    // Identify unknown nodes
    std::vector<int> unknown_nodes;
    for (size_t i = 0; i < num_nodes; ++i)
    {
        if (T_known.find(i) == T_known.end())
        { // If node i is NOT in T_known
            unknown_nodes.push_back(i);
        }
    }

    // Reduced System
    size_t num_unknown_nodes = unknown_nodes.size();
    std::vector<std::vector<double>> K_reduced(num_unknown_nodes, std::vector<double>(num_unknown_nodes, 0.0));
    std::vector<double> f_reduced(num_unknown_nodes, 0.0);

    // Populate K_reduced and initial f_reduced
    for (size_t i = 0; i < num_unknown_nodes; ++i)
    {
        f_reduced[i] = f_global[unknown_nodes[i]];
        for (size_t j = 0; j < num_unknown_nodes; ++j)
        {
            K_reduced[i][j] = K_global[unknown_nodes[i]][unknown_nodes[j]];
        }
    }

    // Incorporate known temperatures into f_reduced
    for (const auto &pair : T_known)
    {
        int node = pair.first;
        double temp = pair.second;

        for (size_t i = 0; i < num_unknown_nodes; ++i)
        {
            int unknown_node = unknown_nodes[i];
            f_reduced[i] -= K_global[unknown_node][node] * temp;
        }
    }

    // Print K_reduced and f_reduced (optional)
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "K_reduced:" << std::endl;
    for (size_t i = 0; i < num_unknown_nodes; ++i)
    {
        for (size_t j = 0; j < num_unknown_nodes; ++j)
        {
            std::cout << K_reduced[i][j] << (j == num_unknown_nodes - 1 ? "" : " ");
        }
        std::cout << std::endl;
    }

    std::cout << "\nf_reduced:" << std::endl;
    for (size_t i = 0; i < num_unknown_nodes; ++i)
    {
        std::cout << f_reduced[i] << std::endl;
    }

    // Solve linear system
    auto T_unknown = solve_system_with_LU(K_reduced, f_reduced);

    // Print solution
    std::cout << "Solution T_unknown:" << std::endl;
    for (size_t i = 0; i < T_unknown.size(); ++i)
        std::cout << "T_unknown[" << i << "]: " << T_unknown[i] << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;

    // Reconstruct Full Temperature Vector
    std::vector<double> T(num_nodes, 0.0);

    // Assign known temperatures
    for (const auto &pair : T_known)
    {
        T[pair.first] = pair.second;
    }

    // Assign unknown temperatures
    for (size_t i = 0; i < num_unknown_nodes; ++i)
    {
        T[unknown_nodes[i]] = T_unknown[i];
    }

    // Print the full temperature vector T
    std::cout << "\nFull Temperature Vector T:" << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < num_nodes; ++i)
    {
        std::cout << "T[" << i << "] = " << std::setw(10) << T[i] << std::endl;
    }

    if (save_temperature_vector(T, "temperature_output.txt"))
    {
        std::cout << "Temperature vector saved successfully to temperature_output.txt" << std::endl;
    }
    else
    {
        std::cerr << "Failed to save temperature vector." << std::endl;
    }

    if (save_coords(mesh.points, "coords.txt"))
    {
        std::cout << "Coordinates saved successfully to coords.txt" << std::endl;
    }
    else
    {
        std::cerr << "Failed to save coordinates." << std::endl;
    }

    if (save_triangles(mesh.triangles, "triangles.txt"))
    {
        std::cout << "Triangles saved successfully to triangles.txt" << std::endl;
    }
    else
    {
        std::cerr << "Failed to save triangles." << std::endl;
    }

    return 0;
}