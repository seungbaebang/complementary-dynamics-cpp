// #include <Eigen/Core>
// #include <Eigen/Sparse>
#include "lumped_mass_matrix.h"


// Sparse Version LBS
void lumped_mass_matrix(const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & T,
        Eigen::SparseMatrix<double>& M) {

        // lumped mass matrix
        Eigen::SparseMatrix<double> Ms;
        igl::massmatrix(V, T, igl::MASSMATRIX_TYPE_DEFAULT, Ms);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(V.rows()*V.cols());
        M.resize(V.rows()*V.cols(),V.rows()*V.cols());
        for (int k=0; k < Ms.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Ms,k); it; ++it) {
                Eigen::Triplet<double> trplt1(3 * it.row() + 0, 3 * it.col() + 0, it.value());
                Eigen::Triplet<double> trplt2(3 * it.row() + 1, 3 * it.col() + 1, it.value());
                Eigen::Triplet<double> trplt3(3 * it.row() + 2, 3 * it.col() + 2, it.value());

                triplets.emplace_back(trplt1);
                triplets.emplace_back(trplt2);
                triplets.emplace_back(trplt3);
            }
        }
        M.setFromTriplets(triplets.begin(), triplets.end());

    }