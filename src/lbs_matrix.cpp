#include "lbs_matrix.h"

void lbs_matrix(const Eigen::MatrixXd & V,
		 const Eigen::MatrixXd & W,
		 Eigen::SparseMatrix<double>& M)
{
	Eigen::MatrixXd Md;
	igl::lbs_matrix(V, W, Md);
	M = Md.sparseView();
}


void lbs_matrix_column(const Eigen::MatrixXd & V,
		 const Eigen::MatrixXd & W,
		 Eigen::SparseMatrix<double>& M)
{
	std::vector<Eigen::Triplet<double>> coefficients;
	 // number of mesh vertices
	 int n = V.rows();
	 assert(n == W.rows());
	 // dimension of mesh
	 int dim = V.cols();
	 // number of handles
	 int m = W.cols();
	 M.resize(n*dim,m*dim*(dim+1));
	 // loop over coordinates of mesh vertices
	 for(int x = 0; x < dim; x++)
	 {
	  // loop over mesh vertices
	  for(int j = 0; j < n; j++)
	  {
	   // loop over handles
	   for(int i = 0; i < m; i++)
	   {
	    // loop over cols of affine transformations
	    for(int c = 0; c < (dim+1); c++)
	    {
	     double value = W(j,i);
	     if(c<dim)
	     {
	      value *= V(j,c);
	     }
	     //coefficients.push_back(Eigen::Triplet<double>(x*n + j,x*m + c*m*dim + i,value));
	     coefficients.push_back(Eigen::Triplet<double>(dim*j+x, x*m + c*m*dim + i,value));
	     //M.insert(x*n + j,x*m + c*m*dim + i) = value;
	    }
	   }
	  }
	 }
	 M.setFromTriplets(coefficients.begin(), coefficients.end());
	 //M.makeCompressed();
}