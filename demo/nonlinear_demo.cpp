//igl
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readMESH.h>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/boundary_facets.h>
#include <igl/lbs_matrix.h>
#include <igl/deform_skeleton.h>
#include <igl/mat_max.h>
#include <igl/mat_min.h>
#include <igl/PI.h>
#include <igl/volume.h>
#include <igl/cat.h>

// #include <igl/matlab/matlabinterface.h>


//internal
#include <read_data_from_json.h>
#include <lumped_mass_matrix.h>
#include <lbs_matrix.h>
#include <create_mask_matrix.h>
#include <line_search.h>

//Bartels
#include <linear_tetmesh_dphi_dX.h>
#include <linear_tetmesh_arap_dq.h>
#include <linear_tetmesh_arap_dq2.h>
#include <linear_tetmesh_arap_q.h>

#include <Eigen/Sparse>
#include <Eigen/Core>

#include <iostream>

std::string folder_path = "/Users/seungbae/Sources/GitHub/ExteriorRigSpace_cpp/anim_data";

//#include "include/read_data_from_json.h"
const Eigen::Vector3d red(255./255.,0./255.,0./255.);

int frame = 0;

void emu_to_lame(const double& E, const double& p, double& lambda, double& mu)
{
  lambda = (E*p)/((1.0+p)*(1.0-2.0*p));
  mu = E/(2.0*(1.0+p));
}

bool read_json_data(const std::string& filename,
Eigen::MatrixXd& V, Eigen::MatrixXi& T, Eigen::MatrixXi& F,
Eigen::MatrixXd& C, Eigen::VectorXi& PI, Eigen::MatrixXi& BE,
Eigen::MatrixXd& W, std::vector<Eigen::MatrixXd>& T_list,
double& dt, int& k, double& YM, double& pr, double& scale)
{
  using json = nlohmann::json;

	std::ifstream infile(filename);
	if (!infile)
		return false;
	json j;
	infile >> j;


	const std::string dir = igl::dirname(filename) + PATH_SEPARATOR;
  std::string model_name;
	if(j.count("model"))
		model_name = j["model"];
	else{
		std::cerr<<"no model specified"<<std::endl;
		return false;
	}
  	if(j.count("mesh_file"))
	{
		std::string mesh_filename = j["mesh_file"];
		Eigen::MatrixXi TF;
		igl::readMESH(dir+mesh_filename, V, T, F);
    std::cout<<"dir+mesh_filename: "<<dir+mesh_filename<<std::endl;
	}
  	if(j.count("handle_file"))
	{
		std::string handle_filename = j["handle_file"];
		Eigen::MatrixXi E, CE, PE;
		igl::readTGF(dir+handle_filename, C, E, PI, BE, CE, PE);
	}
  	if(j.count("weight_file"))
	{
		std::string weights_filename = j["weight_file"];
		igl::readDMAT(dir+weights_filename, W);
	}

  if(j.count("dt"))
		dt = j["dt"];
	if(j.count("k"))
		k = j["k"];
	if(j.count("YM"))
		YM = j["YM"];
	if(j.count("pr"))
		pr = j["pr"];
	if(j.count("scale"))
		scale = j["scale"];

  Eigen::RowVectorXd Vcenter;
  double Vbound;
  Eigen::VectorXd Vmax, Vmin;
  Eigen::VectorXi I;
  igl::mat_max(V, 1, Vmax, I);
  igl::mat_min(V, 1, Vmin, I);
  Vcenter = (Vmax+Vmin)/2.0;
  //Vcenter = V.colwise().minCoeff();
  V = V.rowwise()-Vcenter;
  Vbound = V.maxCoeff();
  V = V/Vbound;

  C = C.rowwise()-Vcenter;
  C = C/Vbound;


  if(j.count("anim_file"))
	{
		std::string anim_filename = j["anim_file"];
    read_pnt_anim(dir+anim_filename, C, Vcenter.transpose(), Vbound, T_list);
	}

}

Eigen::VectorXd vectorize(const Eigen::MatrixXd& M)
{
  Eigen::MatrixXd MT = M.transpose();
  return Eigen::Map<Eigen::VectorXd>(MT.data(),M.rows()*M.cols());
}

Eigen::MatrixXd matrixize(Eigen::VectorXd V)
{
  Eigen::MatrixXd MT = Eigen::Map<Eigen::MatrixXd>(V.data(), 3, V.size()/3);
  return MT.transpose();
}


int main(int argc, char *argv[])
{
  Eigen::MatrixXd V, C, W, VM;
  Eigen::MatrixXi T,F, E, BE;
  Eigen::VectorXi PI;

  std::vector<Eigen::MatrixXd> T_list;
  double dt, YM, pr, scale;
  int k;

  std::string model = argc>1?std::string(argv[1]):"sphere";
  std::string json_path = folder_path+"/"+model+"/"+model+".json";

  //std::string filename = "/Users/seungbae/Sources/GitHub/ExteriorRigSpace_cpp/anim_data/sphere/sphere.json";

  Eigen::MatrixXd SV,TC; Eigen::MatrixXi SF,FTC,FN;
  {
    std::string obj_path = folder_path+"/"+model+"/"+model+".obj";
    Eigen::MatrixXd N; 
    igl::readOBJ(obj_path,SV,TC,N,SF,FTC,FN);
    std::cout<<"N row: "<<N.rows()<<", col: "<<N.cols()<<std::endl;
  }

  read_json_data(json_path,V,T,F,C,PI,BE,W,T_list,dt,k,YM,pr,scale);

  YM = argc>2?std::stod(argv[2]):YM;
  pr = argc>3?std::stod(argv[3]):pr;

  double lambda, mu;
  emu_to_lame(YM,pr,lambda,mu);


  Eigen::MatrixXd params(T.rows(),2);
  params.col(0) = 0.5*lambda*Eigen::VectorXd::Ones(T.rows());
  params.col(1) = mu*Eigen::VectorXd::Ones(T.rows());



  igl::boundary_facets(T, F);
  F = F.rowwise().reverse().eval();

	igl::lbs_matrix(V,W,VM);


  Eigen::VectorXd vol;

  igl::volume(V,T,vol);

  Eigen::MatrixXd dX;
  sim::linear_tetmesh_dphi_dX(dX,V,T);

  Eigen::SparseMatrix<double> M;

  lumped_mass_matrix(V,T,M);
  M = M*1000;

  Eigen::SparseMatrix<double> A;
  lbs_matrix_column(V,W,A);




  Eigen::SparseMatrix<double> phi;
  create_poisson_mask_matrix(V,T,phi);

  Eigen::SparseMatrix<double> Aeq = A.transpose()*M*phi;
  Eigen::VectorXd Beq = Eigen::VectorXd::Zero(A.cols());

  Eigen::MatrixXd T_new = T_list[0];

  Eigen::MatrixXd Vr = VM*T_new;
  Eigen::MatrixXd U = Vr-V;

  int nc = V.rows()*V.cols();

  Eigen::VectorXd UCol = vectorize(U);

  Eigen::VectorXd VCol = vectorize(V);

  Eigen::VectorXd UdCol = Eigen::VectorXd::Zero(nc);
  Eigen::VectorXd UcCol = Eigen::VectorXd::Zero(nc);

  // Engine* engine;
  // igl::matlab::mlinit(&engine);
  // igl::matlab::mlsetmatrix(&engine,"Vrc",Vr);
  // igl::matlab::mleval(&engine,"save('Vrc.mat','Vrc')");
  std::vector<Eigen::MatrixXd> Vn_list(T_list.size());
  int max_iter = 20;
  for(int ai=0; ai<T_list.size(); ai++){
    std::cout<<"ai: "<<ai<<std::endl;

    Eigen::MatrixXd TM = T_list[ai];

    Eigen::VectorXd UCol0 = UCol;
    Eigen::VectorXd UdCol0 = UdCol;
    Eigen::VectorXd UcCol0 = UcCol;

    Vr = VM*TM;
    Eigen::MatrixXd Ur = Vr-V;
    Eigen::VectorXd UrCol = vectorize(Ur);

    for(int i=0; i<max_iter; i++)
    {

      Eigen::VectorXd q = VCol+UrCol+UcCol;

      Eigen::VectorXd G;
      Eigen::SparseMatrix<double> K;


      sim::linear_tetmesh_arap_dq(G,V,T,q,dX,vol,params);

      sim::linear_tetmesh_arap_dq2(K,V,T,q,dX,vol,params);

      double e = sim::linear_tetmesh_arap_q(V, T, q, dX, vol, params);
      //std::cout<<"e: "<<e<<std::endl;
      
      Eigen::VectorXd tmp_g = M/(dt*dt) * (UrCol+UcCol) - M*(UCol0/(dt*dt)+UdCol0/dt) + G;
      Eigen::SparseMatrix<double> tmp_H = M/(dt*dt) + K;
      tmp_H = 0.5 * (tmp_H+Eigen::SparseMatrix<double>(tmp_H.transpose()));

      Eigen::SparseMatrix<double> AR1, AR2, AA;
      igl::cat(2,tmp_H,Eigen::SparseMatrix<double>(Aeq.transpose()),AR1);
      Eigen::SparseMatrix<double> SI(Aeq.rows(), Aeq.rows());
      SI.setZero();
      igl::cat(2,Aeq,SI,AR2);
      igl::cat(1, AR1,AR2, AA);

      Eigen::VectorXd b(tmp_g.size()+Beq.size());
      b.head(tmp_g.size())=-tmp_g;
      b.tail(Beq.size())=Beq;

      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double > > ldlt(AA);
      Eigen::VectorXd dUc = ldlt.solve(b).head(tmp_g.size());

      if(tmp_g.transpose()*dUc > -1e-6)
        break;

      std::function<double(const Eigen::VectorXd &, const Eigen::VectorXd &)> f;
      f = [&](const Eigen::VectorXd &UrColi, const Eigen::VectorXd &UcColi){
        return
        0.5*(UrColi+UcColi-UCol0-dt*UdCol0).transpose()*M/(dt*dt)*(UrColi+UcColi-UCol0-dt*UdCol0)
        + sim::linear_tetmesh_arap_q(V, T, VCol+UrColi+UcColi, dX, vol, params);
      };


      double alpha = line_search(f,tmp_g,dUc,UrCol,UcCol);
      UcCol = UcCol + alpha * dUc;
      //std::cout<<"i: "<<i<<", alpha: "<<alpha<<std::endl;
    }
    UCol = UrCol + UcCol;
    UdCol = (UCol-UCol0)/dt;

    Eigen::MatrixXd Vn = V + matrixize(UCol);
    Vn_list[ai] = Vn;

    std::string numstr = std::to_string(ai);
    int strlen = numstr.length();

    std::string seqstr;
    for(int k=0; k<(4-strlen); k++){
        seqstr.append("0");
    }
    seqstr.append(numstr);

    std::string save_path = folder_path+"/"+model+"/output/"+model+seqstr+".obj";
    //std::string save_path = "/Users/seungbae/Sources/GitHub/ExteriorRigSpace_cpp/data/sphere/sphere"+seqstr+".obj";
    //igl::writeOBJ(save_path, Vn, F);

    Eigen::MatrixXd SVn = Vn.block(0,0,SV.rows(),3);
    Eigen::MatrixXd SN;
    //igl::per_vertex_normals(SVn, SF,SN);
    igl::per_corner_normals(SVn, SF, 20, SN);
    igl::writeOBJ(save_path,SVn,SF,Eigen::MatrixXd(),Eigen::MatrixXi(),TC,FTC);
    std::cout<<"SN row: "<<SN.rows()<<", col: "<<SN.cols()<<std::endl;


  }

// //   Eigen::MatrixXd U = V;
// 	igl::opengl::glfw::Viewer viewer;
// 	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool
// 	{
// 		if(viewer.core().is_animating)
// 		{
// 			Eigen::MatrixXd T = T_list[frame];
//       Eigen::MatrixXd Vn = Vn_list[frame];
// 			//U = VM*T;
// 			frame++;

//       viewer.data().set_vertices(Vn);
//       viewer.data().compute_normals();

//       if(frame == T_list.size()){
//         frame = 0;
//         //viewer.core().is_animating = false;
//       }
// 		}
// 		return false;
// 	};
// viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
// {
//   if(key==' '){
//     viewer.core().is_animating = !viewer.core().is_animating;
//   }
//   return false;
// };
// 	viewer.data().set_mesh(V, F);
// 	viewer.core().is_animating = false;
// 	viewer.core().animation_max_fps = 24.;
// 	viewer.data().show_lines = true;
// 	viewer.data().show_overlay_depth = false;
//   viewer.launch();
}
