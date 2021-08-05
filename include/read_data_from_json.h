#ifndef READ_MODEL_AND_RIG_FROM_JSON_H
#define READ_MODEL_AND_RIG_FROM_JSON_H

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Sparse>
//#include "Skeleton.h"
#include <vector>
#include <utility> // std::pair


// Implementation

#include "json.hpp"

#include <igl/forward_kinematics.h>
#include <igl/directed_edge_parents.h>
#include <igl/readOBJ.h>
#include <igl/readMESH.h>
#include <igl/readDMAT.h>
#include <igl/dirname.h>
#include <igl/png/readPNG.h>
#include <igl/readTGF.h>
#include <igl/mat_max.h>
#include <igl/mat_min.h>

#include <igl/PI.h>
#include <Eigen/Geometry>
#include <fstream>
#include <iostream>
#include <cassert>
#include <dirent.h>


#if defined(WIN32) || defined(_WIN32)
	#define PATH_SEPARATOR std::string("\\")
	#else
#define PATH_SEPARATOR std::string("/")
#endif


void convert_vector_mask_to_sparsematrix(int dim, Eigen::VectorXd& V, Eigen::SparseMatrix<double>& SM)
{
	int num_size = V.size();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.resize(num_size*dim);
    SM.resize(num_size*dim, num_size*dim);
    for (int k=0; k < num_size; ++k) {

        triplets[3*k]   = Eigen::Triplet<double>(3*k, 3*k,V(k));
        triplets[3*k+1] = Eigen::Triplet<double>(3*k+1, 3*k+1,V(k));
        triplets[3*k+2] = Eigen::Triplet<double>(3*k+2, 3*k+2,V(k));
    }
    SM.setFromTriplets(triplets.begin(), triplets.end());
}


void load_sparse_matrix(std::string file_path, Eigen::SparseMatrix<double>& SM)
{
	std::vector<Eigen::Triplet<double>> coefficients;
	FILE* file;
	file = fopen(file_path.c_str(), "rb");
	int num_row, num_col, num_coeff;
	fscanf(file, "%d%d%d\n", &num_row, &num_col, &num_coeff);
	coefficients.resize(num_coeff);
	for(int i=0; i<num_coeff; i++){
		int ri, ci;
		double val;
		fscanf(file, "%d%d%lf", &ri,&ci,&val);
		coefficients[i] = Eigen::Triplet<double>(ri,ci, val);
	}
	fclose(file);

	SM.resize(num_row, num_col);
	SM.setFromTriplets(coefficients.begin(), coefficients.end());
}

void normalize_obj(Eigen::MatrixXd& V, const Eigen::RowVectorXd& center, const double& scale){
    V = V.rowwise()-center;
    V = V/scale;
}

void read_mesh_files_from_directory(const std::string& directory, std::vector<std::string>& mesh_list){

	using namespace std;

	std::vector<std::string> file_list;

	DIR *d;
	struct dirent *dir;
	int i=0;
	d = opendir(directory.c_str());
	if (d)
	{
		while ((dir = readdir(d)) != NULL)
		{
			i++;
			file_list.push_back(dir->d_name);
		}
		closedir(d);
	}
	for(int i=0; i<file_list.size(); i++){
		std::string file = file_list[i];
		size_t last_index = file.find_last_of(".");
		if(file.substr(last_index+1) == "mesh"){
			std::string name = file.substr(0,last_index);
			mesh_list.push_back(name);
		}
	}
}

void read_blendshape_list(const std::string& directory, const std::string& list_file,
		std::vector<std::string>& name_list,
		std::vector<Eigen::MatrixXd>& V_list)
{
	FILE* file;
	file= fopen(list_file.c_str(), "rb");
	int num_list;
	fscanf(file, "%d\n", &num_list);


	std::vector<int> reordered_id(num_list);
	for(int i=0; i<num_list; i++)
	{
		//std::string bs_name;
		char bs_char[256];
		fscanf(file, "%s\n", bs_char);
		std::string bs_name(bs_char);
		for(int j=0; j<name_list.size(); j++){
			if(bs_name==name_list[j]){
				reordered_id[i]=j;
			}
		}
	}
	std::vector<std::string> reorder_name_list(num_list);
	for(int i=0; i<reordered_id.size(); i++){
		reorder_name_list[i] = name_list[reordered_id[i]];
	}
	name_list = reorder_name_list;
	fclose(file);

	V_list.resize(name_list.size());
	for(int j=0; j<name_list.size(); j++){

		std::string mesh_path = directory+"/"+name_list[j]+".mesh";
		Eigen::MatrixXd BV;
		Eigen::MatrixXi BT, BF;
		igl::readMESH(mesh_path, BV,BT,BF);
		V_list[j] = BV;
	}
}




void read_blendshape_anim(const std::string& anim_file, std::vector<Eigen::VectorXd>& w_list)
{
	FILE* file;
	file= fopen(anim_file.c_str(), "rb");
	int num_bs, num_frame;
	fscanf(file, "%d %d\n", &num_bs, &num_frame);

	w_list.resize(num_frame);
	for(int i=0; i<num_frame; i++){
		Eigen::VectorXd w(num_bs);
		for(int j=0; j<num_bs; j++){
			double val;
			fscanf(file, "%lf", &val);
			w(j) = val;
		}
		w_list[i] = w;
	}
	fclose(file);
}


void construct_blendshape_mat(const Eigen::MatrixXd& V,
		const std::vector<Eigen::MatrixXd>& BV_list,
		Eigen::RowVectorXd Vcenter, double Vbound,
		Eigen::SparseMatrix<double>& BS)
{
	Eigen::MatrixXd VT = V.transpose();
	Eigen::VectorXd VCol(Eigen::Map<Eigen::VectorXd>(VT.data(), V.cols()*V.rows()));
	std::vector<Eigen::Triplet<double> > coefficients;
	for(int j=0; j<BV_list.size(); j++){
		Eigen::MatrixXd BV = BV_list[j];

		normalize_obj(BV, Vcenter, Vbound);

		Eigen::MatrixXd BVT = BV.transpose();

		Eigen::VectorXd bVcol = Eigen::Map<Eigen::VectorXd>(BVT.data(), BV.cols()*BV.rows());

		Eigen::VectorXd deltaV = bVcol-VCol;
		for(int i=0; i<deltaV.size(); i++){
			if(deltaV(i)>1e-12){
				coefficients.push_back(Eigen::Triplet<double>(i,j, deltaV(i)));
			}
		}
	}
	BS.resize(VCol.size(), BV_list.size());
	BS.setFromTriplets(coefficients.begin(), coefficients.end());
}

void euler_to_quat(const Eigen::Vector3d& euler, Eigen::Quaterniond& q) {

	q = Eigen::AngleAxisd(euler(2), Eigen::Vector3d::UnitZ())
			* Eigen::AngleAxisd(euler(1), Eigen::Vector3d::UnitY())
			* Eigen::AngleAxisd(euler(0), Eigen::Vector3d::UnitX());
}


void euler_to_quat(const Eigen::Vector3d& euler, const Eigen::Affine3d& a, Eigen::Quaterniond& q) {
	q = Eigen::AngleAxisd(euler(2), a.rotation().col(2))
			* Eigen::AngleAxisd(euler(1), a.rotation().col(1))
			* Eigen::AngleAxisd(euler(0), a.rotation().col(0));
}



void read_bone_anim(std::string anim_file,
		const Eigen::MatrixXd& C,
		const Eigen::MatrixXi& BE,
		const Eigen::VectorXi& P,
		std::vector<Eigen::MatrixXd>& T_list)
{
	FILE* file;
	file= fopen(anim_file.c_str(), "rb");

	double degree2radian = igl::PI/180.0;

	int num_bone, num_frame;
	double val = 0;

	typedef std::vector<Eigen::Quaterniond,
			Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


	fscanf(file, "%d %d\n", &num_bone, &num_frame);
	T_list.resize(num_frame);
	RotationList rot_list(num_bone);
	std::vector<Eigen::Vector3d> tran_list(num_bone);
	std::vector<Eigen::Affine3d> rest_list(num_bone);

	Eigen::Vector3d root_rest_tran;
	Eigen::Affine3d root_affine;
	for (int j = 0; j < 3; j++) {
		fscanf(file, "%lf", &val);
		root_rest_tran(j) = val;
	}
	Eigen::Vector3d euler;
	for (int j = 0; j < 3; j++) {
		fscanf(file, "%lf", &val);
		euler(j) = val;
	}
	euler = euler.array() * degree2radian;
	Eigen::Quaterniond q;
	euler_to_quat(euler, q);
	root_affine = Eigen::Affine3d::Identity();
	root_affine.rotate(q);

	for (int i = 0; i < num_bone; i++) {
		Eigen::Vector3d euler;
		for (int j = 0; j < 3; j++) {
			fscanf(file, "%lf", &val);
			euler(j) = val;
		}
		euler = euler.array() * degree2radian;
		Eigen::Quaterniond q;
		euler_to_quat(euler, q);
		Eigen::Affine3d a = Eigen::Affine3d::Identity();
		a.rotate(q);
		rest_list[i] = a;
	}
	//Eigen::Affine3d root_affine = rest_list[0];
	Eigen::Vector3d root_tran, root_rot;
	for (int k = 0; k < num_frame; k++) {
		for (int j = 0; j < 3; j++) {
			fscanf(file, "%lf", &val);
			root_rot(j) = val;
		}
		for (int j = 0; j < 3; j++) {
			fscanf(file, "%lf", &val);
			root_tran(j) = val;
		}
		root_rot = root_rot.array() * degree2radian;
		Eigen::Quaterniond root_q;
		euler_to_quat(root_rot, root_affine, root_q);

		for (int i = 0; i < num_bone; i++) {
			Eigen::Vector3d euler;
			for (int j = 0; j < 3; j++) {
				fscanf(file, "%lf", &val);
				euler(j) = val;
			}
			euler = euler.array() * degree2radian;
			Eigen::Quaterniond q;
			euler_to_quat(euler, rest_list[i], q);

			if(P(i)==-1){
				int root_cnt = 0;
				for(int ii=0; ii<num_bone; ii++)
					if(P(ii)==-1)
						root_cnt++;
				if(root_cnt>1)
					q = q*root_q;
				tran_list[i] = root_tran - root_rest_tran;
			}
			else
				tran_list[i] = Eigen::Vector3d::Zero();
			rot_list[i] = q;
		}
		RotationList vQ;
		std::vector<Eigen::Vector3d> vT;
		igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

		Eigen::MatrixXd T(num_bone * 4, 3);
		for (int i = 0; i < num_bone; i++) {
			Eigen::Affine3d a = Eigen::Affine3d::Identity();
			a.translate(vT[i]);
			a.rotate(vQ[i]);
			T.block(i * 4, 0, 4, 3) = a.matrix().transpose().block(0, 0, 4, 3);
		}
		T_list[k] = T;
	}
	fclose(file);
}


void read_bone_anim(std::string anim_file,
		const Eigen::MatrixXd& C,
		const Eigen::MatrixXi& BE,
		const Eigen::VectorXi& P,
		const Eigen::VectorXd& center,
		const double& scale,
		std::vector<Eigen::MatrixXd>& T_list)
{
	FILE* file;
	file= fopen(anim_file.c_str(), "rb");

	double degree2radian = igl::PI/180.0;

	int num_bone, num_frame;
	double val = 0;

	typedef std::vector<Eigen::Quaterniond,
			Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


	fscanf(file, "%d %d\n", &num_bone, &num_frame);
	T_list.resize(num_frame);
	RotationList rot_list(num_bone);
	std::vector<Eigen::Vector3d> tran_list(num_bone);
	std::vector<Eigen::Affine3d> rest_list(num_bone);

	Eigen::Vector3d root_rest_tran;
	Eigen::Affine3d root_affine;
	for (int j = 0; j < 3; j++) {
		fscanf(file, "%lf", &val);
		root_rest_tran(j) = val;
	}
	root_rest_tran = (root_rest_tran-center)/scale;

	Eigen::Vector3d euler;
	for (int j = 0; j < 3; j++) {
		fscanf(file, "%lf", &val);
		euler(j) = val;
	}
	euler = euler.array() * degree2radian;
	Eigen::Quaterniond q;
	euler_to_quat(euler, q);
	root_affine = Eigen::Affine3d::Identity();
	root_affine.rotate(q);

	for (int i = 0; i < num_bone; i++) {
		Eigen::Vector3d euler;
		for (int j = 0; j < 3; j++) {
			fscanf(file, "%lf", &val);
			euler(j) = val;
		}
		euler = euler.array() * degree2radian;
		Eigen::Quaterniond q;
		euler_to_quat(euler, q);
		Eigen::Affine3d a = Eigen::Affine3d::Identity();
		a.rotate(q);
		rest_list[i] = a;
	}
	//Eigen::Affine3d root_affine = rest_list[0];
	Eigen::Vector3d root_tran, root_rot;
	for (int k = 0; k < num_frame; k++) {
		for (int j = 0; j < 3; j++) {
			fscanf(file, "%lf", &val);
			root_rot(j) = val;
		}
		for (int j = 0; j < 3; j++) {
			fscanf(file, "%lf", &val);
			root_tran(j) = val;
		}
		root_tran = (root_tran - center)/scale;
		root_rot = root_rot.array() * degree2radian;
		Eigen::Quaterniond root_q;
		euler_to_quat(root_rot, root_affine, root_q);

		for (int i = 0; i < num_bone; i++) {
			Eigen::Vector3d euler;
			for (int j = 0; j < 3; j++) {
				fscanf(file, "%lf", &val);
				euler(j) = val;
			}
			euler = euler.array() * degree2radian;
			Eigen::Quaterniond q;
			euler_to_quat(euler, rest_list[i], q);

			if(P(i)==-1){
				int root_cnt = 0;
				for(int ii=0; ii<num_bone; ii++)
					if(P(ii)==-1)
						root_cnt++;
				if(root_cnt>1)
					q = q*root_q;
				tran_list[i] = root_tran - root_rest_tran;
			}
			else
				tran_list[i] = Eigen::Vector3d::Zero();
			rot_list[i] = q;
		}
		RotationList vQ;
		std::vector<Eigen::Vector3d> vT;
		igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

		Eigen::MatrixXd T(num_bone * 4, 3);
		for (int i = 0; i < num_bone; i++) {
			Eigen::Affine3d a = Eigen::Affine3d::Identity();
			a.translate(vT[i]);
			a.rotate(vQ[i]);
			T.block(i * 4, 0, 4, 3) = a.matrix().transpose().block(0, 0, 4, 3);
		}
		T_list[k] = T;
	}
	fclose(file);
}



void read_pnt_anim(std::string anim_file,
		//const Eigen::MatrixXd& C,
		//int dim,
		const Eigen::VectorXd& center,
		const double& scale,
		std::vector<Eigen::MatrixXd>& T_list)
{
	FILE* file;
	file= fopen(anim_file.c_str(), "rb");

	int dim = 3;//= C.cols();

	double degree2radian = igl::PI/180.0;

	int num_bone, num_frame;
	double val = 0;

	typedef std::vector<Eigen::Quaterniond,
			Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


	fscanf(file, "%d %d\n", &num_bone, &num_frame);
	T_list.resize(num_frame);
	//std::cout<<"num_bone: "<<num_bone<<", num_frame: "<<num_frame<<std::endl;

	std::vector<Eigen::Vector3d> rest_tran_list(num_bone);
	std::vector<Eigen::Affine3d> rest_affine_list(num_bone);

	for(int i=0; i<num_bone ; i++){
		Eigen::Vector3d rest_tran;
		for (int j = 0; j < dim; j++) {
			fscanf(file, "%lf", &val);
			rest_tran(j) = val;
		}
		rest_tran_list[i] = (rest_tran-center)/scale;

		Eigen::Vector3d euler;
		for (int j = 0; j < dim; j++) {
			fscanf(file, "%lf", &val);
			euler(j) = val;
		}
		euler = euler.array() * degree2radian;
		Eigen::Quaterniond q;
		euler_to_quat(euler, q);
		Eigen::Affine3d root_affine = Eigen::Affine3d::Identity();
		root_affine.rotate(q);
		rest_affine_list[i] = root_affine;
	}
	RotationList rot_list(num_bone);
	std::vector<Eigen::Vector3d> tran_list(num_bone);

	for (int k = 0; k < num_frame; k++) {
		for (int i = 0; i < num_bone; i++) {
			Eigen::Vector3d tran;
			for (int j = 0; j < dim; j++) {
				fscanf(file, "%lf", &val);
				tran(j) = val;
			}
			tran = (tran - center)/scale;
			Eigen::Vector3d euler;
			for (int j = 0; j < dim; j++) {
				fscanf(file, "%lf", &val);
				euler(j) = val;
			}
			euler = euler.array() * degree2radian;
			Eigen::Quaterniond q;
			euler_to_quat(euler, rest_affine_list[i], q);

			tran_list[i] = tran - rest_tran_list[i];
			//tran_list[i] = (tran - rest_tran_list[i] - center)/scale;
			rot_list[i] = q;
		}
		RotationList vQ = rot_list;
		std::vector<Eigen::Vector3d> vT= tran_list;
		//igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

		Eigen::MatrixXd T(num_bone * (dim+1), dim);
		for (int i = 0; i < num_bone; i++) {
			Eigen::Affine3d a = Eigen::Affine3d::Identity();
			a.translate(vT[i]);
			a.rotate(vQ[i]);
			T.block(i * (dim+1), 0, dim+1, dim) = a.matrix().transpose().block(0, 0, dim+1, dim);
		}
		T_list[k] = T;
	}
	fclose(file);

}


void read_pnt_anim(std::string anim_file,
    const Eigen::MatrixXd& C,
    const Eigen::VectorXd& center,
    const double& scale,
    std::vector<Eigen::MatrixXd>& T_list)
{
  FILE* file;
  file= fopen(anim_file.c_str(), "rb");

  int dim = 3;//= C.cols();

  double degree2radian = igl::PI/180.0;

  int num_bone, num_frame;
  double val = 0;

  typedef std::vector<Eigen::Quaterniond,
      Eigen::aligned_allocator<Eigen::Quaterniond> > RotationList;


  fscanf(file, "%d %d\n", &num_bone, &num_frame);
  T_list.resize(num_frame);
  //std::cout<<"num_bone: "<<num_bone<<", num_frame: "<<num_frame<<std::endl;

  std::vector<Eigen::Vector3d> rest_tran_list(num_bone);
  std::vector<Eigen::Affine3d> rest_affine_list(num_bone);

  for(int i=0; i<num_bone ; i++){
    Eigen::Vector3d rest_tran;
    for (int j = 0; j < dim; j++) {
      fscanf(file, "%lf", &val);
      rest_tran(j) = val;
    }
    rest_tran_list[i] = (rest_tran-center)/scale;

    Eigen::Vector3d euler;
    for (int j = 0; j < dim; j++) {
      fscanf(file, "%lf", &val);
      euler(j) = val;
    }
    euler = euler.array() * degree2radian;
    Eigen::Quaterniond q;
    euler_to_quat(euler, q);
    Eigen::Affine3d root_affine = Eigen::Affine3d::Identity();
    root_affine.rotate(q);
    rest_affine_list[i] = root_affine;
    //std::cout<<"i: "<<i<<", mat: "<<root_affine.matrix()<<std::endl;
  }
  RotationList rot_list(num_bone);
  std::vector<Eigen::Vector3d> tran_list(num_bone);

  for (int k = 0; k < num_frame; k++) {
    for (int i = 0; i < num_bone; i++) {
      Eigen::Vector3d tran;
      for (int j = 0; j < dim; j++) {
        fscanf(file, "%lf", &val);
        tran(j) = val;
      }
      tran = (tran - center)/scale;
      Eigen::Vector3d euler;
      for (int j = 0; j < dim; j++) {
        fscanf(file, "%lf", &val);
        euler(j) = val;
      }
      euler = euler.array() * degree2radian;
      Eigen::Quaterniond q;
      euler_to_quat(euler, rest_affine_list[i], q);
      //std::cout<<"1 i: "<<i<<", mat: "<<rest_affine_list[i].rotation()<<std::endl;

      q = rest_affine_list[i].rotation().transpose()*q;

      //euler_to_quat(euler, q);

      //std::cout<<"q w: "<<q.w()<<"< vec: "<<q.vec()<<std::endl;

      const Eigen::Vector3d cn = C.row(i).transpose();

      tran_list[i] = tran - rest_tran_list[i] + cn-q*cn;
      //tran_list[i] = tran - rest_tran_list[i];
      rot_list[i] = q;
    }
    RotationList vQ = rot_list;
    std::vector<Eigen::Vector3d> vT= tran_list;
    //igl::forward_kinematics(C, BE, P, rot_list, tran_list, vQ, vT);

    Eigen::MatrixXd T(num_bone * (dim+1), dim);
    for (int i = 0; i < num_bone; i++) {
      Eigen::Affine3d a = Eigen::Affine3d::Identity();
      a.translate(vT[i]);
      a.rotate(vQ[i]);
      T.block(i * (dim+1), 0, dim+1, dim) = a.matrix().transpose().block(0, 0, dim+1, dim);
    }
    T_list[k] = T;
  }
  fclose(file);

}

void read_cage_anim(std::string anim_file, 
		const Eigen::VectorXd& center,
		const double& scale,
		std::vector<Eigen::MatrixXd>& T_list)
{
	FILE* file;
	file= fopen(anim_file.c_str(), "rb");

	int dim = 3;//= C.cols();

	int num_bone, num_frame;
	double val = 0;

	fscanf(file, "%d %d\n", &num_bone, &num_frame);
	T_list.resize(num_frame);

	Eigen::MatrixXd C(num_bone, dim);
	for (int k = 0; k < num_frame; k++) {
		for (int i = 0; i < num_bone; i++) {
			Eigen::Vector3d tran;
			for (int j = 0; j < dim; j++) {
				fscanf(file, "%lf", &val);
				tran(j) = val;
			}
			tran = (tran - center)/scale;
			C.row(i) = tran.transpose();
		}
		T_list[k] = C;
	}
	fclose(file);
}





// bool read_data_from_json(const std::string &filename,
// 		std::string& model_name,
// 		Eigen::MatrixXd &TV,
// 		Eigen::MatrixXi &TT,
// 		Eigen::MatrixXd &SV,
// 		Eigen::MatrixXi &SF,
// 		Eigen::MatrixXd &TC,
// 		Eigen::MatrixXi &FTC,
// 		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& CR,
// 		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& CG,
// 		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& CB,
// 		Eigen::SparseMatrix<double>& T2S,
// 		Eigen::MatrixXd &CV,
// 		Eigen::VectorXi &PI,
// 		Eigen::MatrixXi &BE,
// 		Eigen::MatrixXd &W,
// 		Eigen::SparseMatrix<double>& SMsk,
// 		std::vector<Eigen::MatrixXd>& T_list,
// 		double& density,
// 		double& YM,
// 		double& mu,
// 		double& fps,
// 		Eigen::RowVectorXd& Vcenter,
// 		double& Vbound
// 		) {
// 	using json = nlohmann::json;

// 	std::ifstream infile(filename);
// 	if (!infile)
// 		return false;
// 	json j;
// 	infile >> j;


// 	const std::string dir = igl::dirname(filename) + PATH_SEPARATOR;

// 	if(j.count("model"))
// 		model_name = j["model"];
// 	else{
// 		std::cerr<<"no model specified"<<std::endl;
// 		return false;
// 	}
	
// 	if(j.count("mesh"))
// 	{
// 		std::string mesh_filename = j["mesh"];
// 		Eigen::MatrixXi TF;
// 		igl::readMESH(dir + mesh_filename, TV, TT, TF);
// 	}
// 	else{
// 		std::cerr<<"no mesh specified"<<std::endl;
// 		return false;
// 	}

// 	//igl::readOBJ(dir + obj_filename, SV, SF);
// 	if(j.count("obj"))
// 	{
// 		std::string obj_filename = j["obj"];
// 		Eigen::MatrixXd N;
// 		Eigen::MatrixXi FN;
// 		igl::readOBJ(dir + obj_filename, SV, TC, N, SF, FTC, FN);
// 	}

// 	if(j.count("texture"))
// 	{
// 		std::string tex_filename = j["texture"];
// 		Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> CA;
// 		igl::png::readPNG(dir + tex_filename,CR,CG,CB,CA);
// 	}

// 	if(j.count("tet2surf"))
// 	{
// 		std::string tet2surf_filename = j["tet2surf"];
// 		load_sparse_matrix(dir + tet2surf_filename,T2S);
// 	}
	
// 	if(j.count("handle"))
// 	{
// 		std::string handle_filename = j["handle"];
// 		Eigen::MatrixXi E, CE, PE;
// 		igl::readTGF(dir + handle_filename, CV, E, PI, BE, CE, PE);
// 	}
// 	else{
// 		std::cerr<<"no handle specified"<<std::endl;
// 		return false;
// 	}

// 	if(j.count("weights"))
// 	{
// 		std::string weights_filename = j["weights"];
// 		igl::readDMAT(dir + weights_filename, W);
// 	}
// 	else{
// 		std::cerr<<"no weights specified"<<std::endl;
// 		return false;
// 	}

// 	if(j.count("mask"))
// 	{
// 		std::string mask_filename = j["mask"];
// 		Eigen::VectorXd Msk;
// 		igl::readDMAT(dir + mask_filename, Msk);
// 	    convert_vector_mask_to_sparsematrix(TV.cols(), Msk, SMsk);
// 	}


// 	{

// 		Vcenter = TV.colwise().minCoeff();
// 		Vbound = TV.maxCoeff();
        
// 	}
//     normalize_obj(TV, Vcenter, Vbound);
//     normalize_obj(SV, Vcenter, Vbound);
//     normalize_obj(CV, Vcenter, Vbound);


// 	if(j.count("animation"))
// 	{
// 		std::string anim_filename = j["animation"];
// 		if(PI.size()>0){
// 			if(j.count("cage")){
// 				std::cout<<"read cage handle"<<std::endl;
// 				read_cage_anim(dir+anim_filename, Vcenter.transpose(), Vbound, T_list);
// 			}
// 			else{
// 				std::cout<<"read pnt animation"<<std::endl;
// 				read_pnt_anim(dir + anim_filename, Vcenter.transpose(), Vbound, T_list);
// 			}
// 		}
// 		else if(BE.rows()>0){
// 			std::cout<<"read bone animation"<<std::endl;
// 			Eigen::VectorXi P;
// 			igl::directed_edge_parents(BE, P);
// 			read_bone_anim(dir + anim_filename, CV, BE, P, Vcenter.transpose(), Vbound, T_list);
// 		}
		
// 		else{
// 			std::cerr<<"doesn't support animation"<<std::endl;
// 			return false;
// 		}
// 	}
// 	else{
// 		std::cerr<<"no animation specified"<<std::endl;
// 		return false;
// 	}
// 	if(j.count("density"))
// 		density = j["density"];
// 	if(j.count("YM"))
// 		YM = j["YM"];
// 	if(j.count("mu"))
// 		mu = j["mu"];
// 	if(j.count("fps"))
// 		fps = j["fps"];

// 	return true;
// }

// bool read_data_from_json(const std::string &filename,
// 		std::string& model_name,
// 		Eigen::MatrixXd &TV,
// 		Eigen::MatrixXi &TT,
// 		Eigen::MatrixXd &SV,
// 		Eigen::MatrixXi &SF,
// 		Eigen::MatrixXd &TC,
// 		Eigen::MatrixXi &FTC,
// 		Eigen::SparseMatrix<double>& T2S,
// 		Eigen::MatrixXd &CV,
// 		Eigen::VectorXi &PI,
// 		Eigen::MatrixXi &BE,
// 		Eigen::MatrixXd &W,
// 		Eigen::SparseMatrix<double>& SMsk,
// 		std::vector<Eigen::MatrixXd>& T_list,
// 		double& density,
// 		double& YM,
// 		double& mu,
// 		double& fps,
// 		Eigen::RowVectorXd& Vcenter,
// 		double& Vbound
// 		) 
// {
// 	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> CR, CG, CB;
// 	return read_data_from_json(filename,model_name,TV,TT,SV,SF,TC,FTC,CR,CG,CB,
// 			T2S,CV,PI,BE,W,SMsk,T_list,density,YM,mu,fps, Vcenter, Vbound);
// }



// bool read_data_from_json(const std::string &filename,
// 		std::string& model_name,
// 		Eigen::MatrixXd &TV,
// 		Eigen::MatrixXi &TT,
// 		Eigen::MatrixXd &SV,
// 		Eigen::MatrixXi &SF, 
// 		Eigen::SparseMatrix<double>& T2S,
// 		Eigen::MatrixXd &CV,
// 		Eigen::VectorXi &PI,
// 		Eigen::MatrixXi &BE,
// 		Eigen::MatrixXd &W,
// 		Eigen::SparseMatrix<double>& SMsk,
// 		std::vector<Eigen::MatrixXd>& T_list,
// 		double& density,
// 		double& YM,
// 		double& mu,
// 		double& fps,
// 		Eigen::RowVectorXd& Vcenter,
// 		double& Vbound
// 		) 
// {

// 	Eigen::MatrixXd TC;
// 	Eigen::MatrixXi FTC;
// 	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> CR, CG, CB;

// 	return read_data_from_json(filename,model_name,TV,TT,SV,SF,TC,FTC,CR,CG,CB,
// 			T2S,CV,PI,BE,W,SMsk,T_list,density,YM,mu,fps, Vcenter,Vbound);

// }

// bool read_data_from_json(const std::string &filename,
// 		std::string& model_name,
// 		Eigen::MatrixXd &TV,
// 		Eigen::MatrixXi &TT,
// 		Eigen::MatrixXd &SV,
// 		Eigen::MatrixXi &SF,
// 		Eigen::MatrixXd &CV,
// 		Eigen::VectorXi &PI,
// 		Eigen::MatrixXi &BE,
// 		Eigen::MatrixXd &W,
// 		Eigen::SparseMatrix<double>& SMsk,
// 		std::vector<Eigen::MatrixXd>& T_list,
// 		double& density,
// 		double& YM,
// 		double& mu,
// 		double& fps,
// 		Eigen::RowVectorXd& Vcenter,
// 		double& Vbound
// 		)
// {


// 	Eigen::MatrixXd TC;
// 	Eigen::MatrixXi FTC;
// 	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> CR, CG, CB;
// 	Eigen::SparseMatrix<double> T2S;

// 	return read_data_from_json(filename,model_name,TV,TT,SV,SF,TC,FTC,CR,CG,CB,
// 			T2S,CV,PI,BE,W,SMsk,T_list,density,YM,mu,fps, Vcenter, Vbound);
// }


// bool read_data_from_json(const std::string &filename,
// 		std::string& model_name,
// 		Eigen::MatrixXd &TV,
// 		Eigen::MatrixXi &TT,
// 		Eigen::MatrixXd &CV,
// 		Eigen::VectorXi &PI,
// 		Eigen::MatrixXi &BE,
// 		Eigen::MatrixXd &W,
// 		std::vector<Eigen::MatrixXd>& T_list,
// 		double& density,
// 		double& YM,
// 		double& mu,
// 		double& fps,
// 		Eigen::RowVectorXd& Vcenter,
// 		double& Vbound
// 		)
// {

// 	Eigen::MatrixXd SV;
// 	Eigen::MatrixXi SF;
// 	Eigen::MatrixXd TC;
// 	Eigen::MatrixXi FTC;
// 	Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> CR, CG, CB;
// 	Eigen::SparseMatrix<double> T2S;
// 	Eigen::SparseMatrix<double> SMsk;

// 	return read_data_from_json(filename,model_name,TV,TT,SV,SF,TC,FTC,CR,CG,CB,
// 			T2S,CV,PI,BE,W,SMsk,T_list,density,YM,mu,fps, Vcenter, Vbound);

// }


// bool read_blendshape_data_from_json(const std::string &filename,
// 		std::string& model_name,
// 		Eigen::MatrixXd &TV,
// 		Eigen::MatrixXi &TT,
// 		Eigen::MatrixXd &SV,
// 		Eigen::MatrixXi &SF,
// 		Eigen::MatrixXd &TC,
// 		Eigen::MatrixXi &FTC,
// 		Eigen::SparseMatrix<double>& BS,
// 		Eigen::SparseMatrix<double>& T2S,
// 		Eigen::SparseMatrix<double>& SMsk,
// 		Eigen::SparseMatrix<double>& UeSMsk,
// 		std::vector<Eigen::VectorXd>& T_list,
// 		double& density,
// 		double& YM,
// 		double& mu,
// 		double& fps,
// 		Eigen::RowVectorXd& Vcenter,
// 		double& Vbound
// 		) {
// 	using json = nlohmann::json;

// 	std::ifstream infile(filename);
// 	if (!infile)
// 		return false;
// 	json j;
// 	infile >> j;


// 	const std::string dir = igl::dirname(filename) + PATH_SEPARATOR;

// 	if(j.count("model"))
// 		model_name = j["model"];
// 	else{
// 		std::cerr<<"no model specified"<<std::endl;
// 		return false;
// 	}
	
// 	if(j.count("mesh"))
// 	{
// 		std::string mesh_filename = j["mesh"];
// 		Eigen::MatrixXi TF;
// 		igl::readMESH(dir + mesh_filename, TV, TT, TF);
// 	}
// 	else{
// 		std::cerr<<"no mesh specified"<<std::endl;
// 		return false;
// 	}

// 	//igl::readOBJ(dir + obj_filename, SV, SF);
// 	if(j.count("obj"))
// 	{
// 		std::string obj_filename = j["obj"];
// 		Eigen::MatrixXd N;
// 		Eigen::MatrixXi FN;
// 		igl::readOBJ(dir + obj_filename, SV, TC, N, SF, FTC, FN);
// 	}

// 	if(j.count("tet2surf"))
// 	{
// 		std::string tet2surf_filename = j["tet2surf"];
// 		load_sparse_matrix(dir + tet2surf_filename,T2S);
// 	}

// 	if(j.count("mask"))
// 	{
// 		std::string mask_filename = j["mask"];
// 		Eigen::VectorXd Msk;
// 		igl::readDMAT(dir + mask_filename, Msk);
// 	    convert_vector_mask_to_sparsematrix(TV.cols(), Msk, SMsk);
// 	}
// 	if(j.count("Ue_mask"))
// 	{
// 		std::string mask_filename = j["Ue_mask"];
// 		Eigen::VectorXd Msk;
// 		igl::readDMAT(dir + mask_filename, Msk);
// 	    convert_vector_mask_to_sparsematrix(TV.cols(), Msk, UeSMsk);
// 	}


// 	{
// 		// move the whole mesh to the center of the coordinate and normalize it
// 		Vcenter = TV.colwise().minCoeff();
// 		Vbound = TV.maxCoeff();
// 	}
//     normalize_obj(TV, Vcenter, Vbound);
//     normalize_obj(SV, Vcenter, Vbound);

// 	if(j.count("bs_list"))
// 	{
// 		std::string bs_list_file = j["bs_list"];
// 		std::string bs_folder = j["bs_folder"];
// 		std::string animfile = j["animation"];

// 		std::vector<std::string> mesh_name_list;
// 		read_mesh_files_from_directory(dir+bs_folder, mesh_name_list);

// 		std::vector<Eigen::MatrixXd> BV_list;
// 		read_blendshape_list(dir+bs_folder, dir+bs_list_file, mesh_name_list, BV_list);
// 		read_blendshape_anim(dir+animfile, T_list);

// 		//Eigen::SparseMatrix<double> BS;
// 		construct_blendshape_mat(TV, BV_list, Vcenter, Vbound, BS);
// 	}
// 	std::cout<<"2"<<std::endl;

	
// 	// if(j.count("animation"))
// 	// {
// 	// 	std::string anim_filename = j["animation"];
// 	// 	if(PI.size()>0){
// 	// 		std::cout<<"read pnt animation"<<std::endl;
// 	// 		read_pnt_anim(dir + anim_filename, Vcenter.transpose(), Vbound, T_list);
// 	// 	}
// 	// 	else if(BE.rows()>0){
// 	// 		std::cout<<"read bone animation"<<std::endl;
// 	// 		Eigen::VectorXi P;
// 	// 		igl::directed_edge_parents(BE, P);
// 	// 		read_bone_anim(dir + anim_filename, CV, BE, P, Vcenter.transpose(), Vbound, T_list);
// 	// 	}
// 	// 	else{
// 	// 		std::cerr<<"doesn't support animation"<<std::endl;
// 	// 		return false;
// 	// 	}
// 	// }
// 	// else{
// 	// 	std::cerr<<"no animation specified"<<std::endl;
// 	// 	return false;
// 	// }
// 	if(j.count("density"))
// 		density = j["density"];
// 	if(j.count("YM"))
// 		YM = j["YM"];
// 	if(j.count("mu"))
// 		mu = j["mu"];
// 	if(j.count("fps"))
// 		fps = j["fps"];
// 	return true;
// }


#endif