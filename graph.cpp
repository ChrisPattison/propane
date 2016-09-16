#include "graph.hpp"

void Graph::Resize(int no_vertices, int no_couplers) {
    fields_.resize(no_vertices);
    fields_.setZero();
    adjacent_.resize(no_vertices, no_vertices);
    adjacent_.reserve(no_couplers);
}

void Graph::SetField(int vertex, EdgeType field) {
    fields_(vertex) = field;
}

void Graph::AddEdge(int from, int to, EdgeType coupler) {
    adjacent_.insert(from, to) = coupler;
}

bool Graph::IsConsistent() const {
    bool not_consistent = false;
    for(std::size_t k = 0; k < adjacent_.outerSize(); ++k) {
        bool zero = true;
        for(Eigen::SparseMatrix<EdgeType>::InnerIterator it(adjacent_, k); it; ++it) {
            if(std::abs(it.value()) > kEpsilon) {
                zero = false;
                break;
            }
        }
        if(zero) {
            not_consistent = true;
        }
    }
    return not_consistent | adjacent_.isApprox(adjacent_.transpose());
}

Eigen::SparseMatrix<EdgeType>& Graph::Adjacent() {
    return adjacent_;
}

Eigen::Matrix<VertexType, Eigen::Dynamic, 1>& Graph::Fields() {
    return fields_;
}

int Graph::size() const {
    return fields_.size();
}