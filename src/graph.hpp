#pragma once
#include "types.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace propane
{
/** Contains the fields and coefficients as well as useful methods for construction and error checking.
 * Contains an adjacency matrix of the problem where the nonzeros are the value of the coefficient J_ij between S_i and S_j.
 * The field h_i is that felt by S_i.
 */
class Graph {
    Eigen::SparseMatrix<EdgeType> adjacent_;
    Eigen::Matrix<EdgeType, Eigen::Dynamic, 1> fields_;
    bool field_nonzero_;
public:
/**  Resizes adjacency matrix and field vector.
 *  Zeros out field vector and reserves memory for adjacency matrix.
 */
    void Resize(int no_vertices, int no_couplers);
/** Sets the field on a particular vertex.
 */
    void SetField(int vertex, EdgeType field);
/** Sets a non-zero coefficient between S_from and S_to.
 * This is unidirecational and must be set in both directions in most cases.
 */
    void AddEdge(int from, int to, EdgeType coupler);
/** Checks validity of the adjacency matrix.
 * Checks for zero rows and asymmetry.
 */
    bool IsConsistent() const;
/** Gets number of vertices.
 */
    int size() const;
/** Number of edges with bidirectional asumption.
 */
    int edges() const;
/** Returns true if a field has been set.
 */
    bool has_field() const;
/** Reference to adjacency matrix.
 * This will be const in the future.
 */
    Eigen::SparseMatrix<EdgeType>& Adjacent();
/** Reference to field vector.
 * This will be const in the future.
 */
     Eigen::Matrix<EdgeType, Eigen::Dynamic, 1>& Fields();
};
}