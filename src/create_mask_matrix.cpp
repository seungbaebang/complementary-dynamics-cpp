
#include "create_mask_matrix.h"
#include <iostream>

enum MaskType
{
    IDENTITY = 1, 
    BINARY = 2, 
    SURFACE = 3,
    POISSON = 4
};

void create_mask_matrix(const int &mask, const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& T, const Eigen::MatrixXi& F,
    Eigen::SparseMatrix<double>& phi) {

        // tet mass
        Eigen::SparseMatrix<double> Mt;
        igl::massmatrix(V,T,igl::MASSMATRIX_TYPE_DEFAULT,Mt);
        Eigen::VectorXd Mtdiag;
        igl::diag(Mt, Mtdiag);

        // surface mass
        Eigen::SparseMatrix<double> Ms;
        igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT, Ms);
        Eigen::VectorXd Msdiag;
        igl::diag(Ms, Msdiag);

        // surface vertices
        Eigen::VectorXi face_vertices;
        igl::unique(F, face_vertices);


        phi.resize(3*V.rows(),3*V.rows());
        std::vector<Eigen::Triplet<double>> coeff;


        switch (mask) {
            case IDENTITY:
                std::cout << "Identity\n";
                for (int i = 0; i < V.rows(); ++i) {
                    coeff.push_back(Eigen::Triplet<double>(3*i+0, 3*i+0, 1));
                    coeff.push_back(Eigen::Triplet<double>(3*i+1, 3*i+1, 1));
                    coeff.push_back(Eigen::Triplet<double>(3*i+2, 3*i+2, 1));
                }
                phi.setFromTriplets(coeff.begin(), coeff.end());
                break;
            case BINARY:
                std::cout << "Binary\n";
                for (int i = 0; i < face_vertices.size(); ++i) {
                    int idx = face_vertices(i);
                    coeff.push_back(Eigen::Triplet<double>(3*idx+0, 3*idx+0, Mtdiag(idx)));
                    coeff.push_back(Eigen::Triplet<double>(3*idx+1, 3*idx+1, Mtdiag(idx)));
                    coeff.push_back(Eigen::Triplet<double>(3*idx+2, 3*idx+2, Mtdiag(idx)));
                }
                phi.setFromTriplets(coeff.begin(), coeff.end());
                break;
            case SURFACE:
                std::cout << "Surface\n";
                for (int i = 0; i < face_vertices.size(); ++i) {
                    int idx = face_vertices(i);
                    coeff.push_back(Eigen::Triplet<double>(3*idx+0, 3*idx+0, Msdiag(idx)));
                    coeff.push_back(Eigen::Triplet<double>(3*idx+1, 3*idx+1, Msdiag(idx)));
                    coeff.push_back(Eigen::Triplet<double>(3*idx+2, 3*idx+2, Msdiag(idx)));
                }
                phi.setFromTriplets(coeff.begin(), coeff.end());
                break;
            case POISSON:
                std::cout << "Invalid Selection\n";
                break;
            default:
                std::cout << "Invalid Selection\n";
                break;
        }


    }




// binary mask matrix
void create_binary_mask_matrix(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& T, const Eigen::MatrixXi& F,
    Eigen::SparseMatrix<double>& phi) {

    Eigen::SparseMatrix<double> Mm;
    igl::massmatrix(V,T,igl::MASSMATRIX_TYPE_DEFAULT,Mm);
    Eigen::VectorXd Mdiag;
    igl::diag(Mm, Mdiag);

    // construct W matrix
    phi.resize(3*V.rows(),3*V.rows());
    Eigen::VectorXi face_vertices;
    igl::unique(F, face_vertices);
    std::vector<Eigen::Triplet<double>> coeff;
    for (int i = 0; i < face_vertices.size(); ++i) {
        int idx = face_vertices(i);
        coeff.push_back(Eigen::Triplet<double>(3*idx+0, 3*idx+0, Mdiag(idx)));
        coeff.push_back(Eigen::Triplet<double>(3*idx+1, 3*idx+1, Mdiag(idx)));
        coeff.push_back(Eigen::Triplet<double>(3*idx+2, 3*idx+2, Mdiag(idx)));
    }
    phi.setFromTriplets(coeff.begin(), coeff.end());

}


// Ms matrix
void create_surface_mask_matrix(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::SparseMatrix<double>& phi) {

    Eigen::SparseMatrix<double> Mm;
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_DEFAULT, Mm);
    Eigen::VectorXd Mdiag;
    igl::diag(Mm, Mdiag);

    // construct W matrix
    phi.resize(3*V.rows(),3*V.rows());
    Eigen::VectorXi face_vertices;
    igl::unique(F, face_vertices);
    std::vector<Eigen::Triplet<double>> coeff;
    for (int i = 0; i < face_vertices.size(); ++i) {
        int idx = face_vertices(i);
        coeff.push_back(Eigen::Triplet<double>(3*idx+0, 3*idx+0, Mdiag(idx)));
        coeff.push_back(Eigen::Triplet<double>(3*idx+1, 3*idx+1, Mdiag(idx)));
        coeff.push_back(Eigen::Triplet<double>(3*idx+2, 3*idx+2, Mdiag(idx)));
    }
    phi.setFromTriplets(coeff.begin(), coeff.end());

}

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/boundary_facets.h>
#include <igl/diag.h>
#include <igl/min_quad_with_fixed.h>

// possion solve
void create_poisson_mask_matrix(const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& T,
    Eigen::SparseMatrix<double>& phi) {
	Eigen::SparseMatrix<double> L,M;
	Eigen::VectorXd ones = Eigen::VectorXd::Constant(V.rows(),1,1);
	igl::cotmatrix(V,T,L);
	igl::massmatrix(V,T,igl::MASSMATRIX_TYPE_DEFAULT,M);
    Eigen::MatrixXi F;
    igl::boundary_facets(T,F);
    Eigen::VectorXi b;
    igl::unique(F, b);
    Eigen::VectorXd bc = Eigen::VectorXd::Zero(b.rows(),1);
    Eigen::VectorXd Z;
    Eigen::SparseMatrix<double> Q = -L;
    Eigen::VectorXd l = M*ones;
	igl::min_quad_with_fixed(Q,l,b,bc,
			Eigen::SparseMatrix<double>(),
			Eigen::VectorXd(),
			false,
			Z);
    Z = (M*Z).eval();
    Eigen::VectorXd ZZ(Z.rows()*3,1);	
    for(int i =0;i<Z.rows();i++)
    {
        ZZ(i*3+0) = Z(i);
        ZZ(i*3+1) = Z(i);
        ZZ(i*3+2) = Z(i);
    }
    igl::diag(ZZ,phi);

}
