#ifndef CREATE_MASK_MATRIX_H
#define CREATE_MASK_MATRIX_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/unique.h>
#include <igl/diag.h>
#include <igl/massmatrix.h>

#include <Eigen/Core>
#include <Eigen/Sparse>

void create_mask_matrix(const int &mask, const Eigen::MatrixXd& V,
		const Eigen::MatrixXi &T, const Eigen::MatrixXi& F,
		Eigen::SparseMatrix<double>& phi);

void create_binary_mask_matrix(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& T, const Eigen::MatrixXi& F,
    Eigen::SparseMatrix<double>& phi);

void create_surface_mask_matrix(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::SparseMatrix<double>& phi);

void create_poisson_mask_matrix(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& T,
    Eigen::SparseMatrix<double>& phi);

#endif