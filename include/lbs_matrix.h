#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/lbs_matrix.h>

// Sparse Version LBS
void lbs_matrix( const Eigen::MatrixXd & V,
		 const Eigen::MatrixXd & W,
		 Eigen::SparseMatrix<double>& M);


// Seungbae's hand written LBS column
// output vertex order [x1 y1 z1 x2 y2 z2 x3 y3 z3]' instead of [x1 x2 x3 y1 y2 y3 z1 z2 z3]'
void lbs_matrix_column( const Eigen::MatrixXd & V,
		 const Eigen::MatrixXd & W,
		 Eigen::SparseMatrix<double>& M);