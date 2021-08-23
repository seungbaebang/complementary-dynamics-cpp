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
#include <igl/png/writePNG.h>

//internal
#include <json.hpp>
#include <read_data_from_json.h>
#include <lumped_mass_matrix.h>
#include <lbs_matrix.h>
#include <create_mask_matrix.h>
#include <line_search.h>
#include <util.h>

//Bartels
#include <linear_tetmesh_dphi_dX.h>
#include <linear_tetmesh_arap_dq.h>
#include <linear_tetmesh_arap_dq2.h>
#include <linear_tetmesh_arap_q.h>
#include <linear_tetmesh_neohookean_dq.h>
#include <linear_tetmesh_neohookean_dq2.h>
#include <linear_tetmesh_neohookean_q.h>
#include <linear_tetmesh_stvk_dq.h>
#include <linear_tetmesh_stvk_dq2.h>
#include <linear_tetmesh_stvk_q.h>
#include <linear_tetmesh_corotational_dq.h>
#include <linear_tetmesh_corotational_dq2.h>
#include <linear_tetmesh_corotational_q.h>
#include <simple_psd_fix.h>
#include <Eigen/Sparse>
#include <Eigen/Core>

#include <iostream>
#include <filesystem>


const Eigen::Vector3d red(255./255.,0./255.,0./255.);
int frame = 0;

void export_png_seq(igl::opengl::glfw::Viewer& viewer, const std::vector<Eigen::MatrixXd>& Vn_list, 
                    const std::string& json_path, const std::string& physic_model)
{
  frame = 0;
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool
	{
		if(viewer.core().is_animating)
		{
      Eigen::MatrixXd Vn = Vn_list[frame];
      // viewer.data().set_mesh(Vn,F);
      viewer.data().set_vertices(Vn);
      viewer.data().compute_normals();

      const std::string dir = igl::dirname(json_path) + PATH_SEPARATOR;;
      size_t dot_found = json_path.find_last_of(".");
      std::string model = std::string(json_path.begin()+dir.size(),json_path.begin()+dot_found);

      std::string output_folder = "../showcases";
      std::string output_path = output_folder+"/"+model+"_"+physic_model;
      if(!std::filesystem::exists(output_path)){
        std::filesystem::create_directory(output_path);
      }

      std::string numstr = std::to_string(frame);
      int strlen = numstr.length();

      std::string seqstr;
      for(int k=0; k<(4-strlen); k++){
          seqstr.append("0");
      }
      seqstr.append(numstr);
      std::string save_path = output_path+"/"+model+"_"+physic_model+"_"+seqstr+".png";

      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R(1280,800);
      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G(1280,800);
      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B(1280,800);
      Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A(1280,800);

      // Draw the scene in the buffers
      viewer.core().draw_buffer(viewer.data(),true,R,G,B,A);
      igl::png::writePNG(R,G,B,A,save_path);

			frame++;
      if(frame == Vn_list.size()){
        viewer.core().is_animating = !viewer.core().is_animating;
        frame = 0;
      }
		}
		return false;
	};
}

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V,C,W,VM; //V: vertices of tet-mesh, C: joint positions, W: skinning weight, VM: lbs matrix
  Eigen::MatrixXi T,F,BE; //T: tet indices  of tet-mesh, F: face indices of tet-mesh, BE: edge indices of bone handle
  Eigen::VectorXi PI; //PI: point indices of point handle

  std::vector<Eigen::MatrixXd> TF_list; //list of transformation
  double dt, YM, pr, scale; //dt: time step, YM: young's modulus, pr: poisson ratio, scale: custom scale of mesh

  std::string json_path = argc>1?std::string(argv[1]):"../examples/sphere/sphere.json";
  std::string physic_model;
  read_json_data(json_path,V,T,F,C,PI,BE,W,TF_list,dt,YM,pr,scale,physic_model);

  std::string str = argc>2?std::string(argv[2]):"";
  bool export_objs = false;
  if(str=="export")
    export_objs=true;

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

  Eigen::MatrixXd TF = TF_list[0];

  Eigen::MatrixXd Vr = VM*TF;
  Eigen::MatrixXd U = Vr-V;

  int nc = V.rows()*V.cols();

  Eigen::VectorXd UCol = vectorize(U);

  Eigen::VectorXd VCol = vectorize(V);

  Eigen::VectorXd UdCol = Eigen::VectorXd::Zero(nc);
  Eigen::VectorXd UcCol = Eigen::VectorXd::Zero(nc);

  std::vector<Eigen::MatrixXd> Vn_list(TF_list.size());
  int max_iter = 20;
  for(int ai=0; ai<TF_list.size(); ai++){
    std::cout<<"frame: "<<ai<<std::endl;

    Eigen::MatrixXd TF = TF_list[ai];

    Eigen::VectorXd UCol0 = UCol;
    Eigen::VectorXd UdCol0 = UdCol;
    Eigen::VectorXd UcCol0 = UcCol;

    Vr = VM*TF;
    Eigen::MatrixXd Ur = Vr-V;
    Eigen::VectorXd UrCol = vectorize(Ur);

    for(int i=0; i<max_iter; i++)
    {
      Eigen::VectorXd q = VCol+UrCol+UcCol;
      Eigen::VectorXd G;
      Eigen::SparseMatrix<double> K;

      if(physic_model=="arap"){
        sim::linear_tetmesh_arap_dq(G,V,T,q,dX,vol,params);
        sim::linear_tetmesh_arap_dq2(K,V,T,q,dX,vol,params,[](auto &a) {sim::simple_psd_fix(a, 1e-3);});
      }
      else if(physic_model=="neohookean"){
        sim::linear_tetmesh_neohookean_dq(G,V,T,q,dX,vol,params);
        sim::linear_tetmesh_neohookean_dq2(K,V,T,q,dX,vol,params,[](auto &a) {sim::simple_psd_fix(a, 1e-3);});
      }
      else if(physic_model=="stvk"){
        sim::linear_tetmesh_stvk_dq(G,V,T,q,dX,vol,params);
        sim::linear_tetmesh_stvk_dq2(K,V,T,q,dX,vol,params,[](auto &a) {sim::simple_psd_fix(a, 1e-3);});
      }
      else if(physic_model=="corotational"){
        sim::linear_tetmesh_corotational_dq(G,V,T,q,dX,vol,params);
        sim::linear_tetmesh_corotational_dq2(K,V,T,q,dX,vol,params,[](auto &a) {sim::simple_psd_fix(a, 1e-3);});
      }
      Eigen::VectorXd tmp_g = M/(dt*dt) * (UrCol+UcCol) - M*(UCol0/(dt*dt)+UdCol0/dt) + G;
      Eigen::SparseMatrix<double> tmp_H = M/(dt*dt) + K;
      tmp_H = 0.5 * (tmp_H+Eigen::SparseMatrix<double>(tmp_H.transpose()));

      Eigen::SparseMatrix<double> AR1, AR2, AA;
      igl::cat(2,tmp_H,Eigen::SparseMatrix<double>(Aeq.transpose()),AR1);
      Eigen::SparseMatrix<double> SI(Aeq.rows(), Aeq.rows());
      SI.setZero();
      igl::cat(2,Aeq,SI,AR2);
      igl::cat(1,AR1,AR2,AA);

      Eigen::VectorXd b(tmp_g.size()+Beq.size());
      b.head(tmp_g.size())=-tmp_g;
      b.tail(Beq.size())=Beq;

      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double > > ldlt(AA);
      Eigen::VectorXd dUc = ldlt.solve(b).head(tmp_g.size());

      if(tmp_g.transpose()*dUc > -1e-6)
        break;

      std::function<double(const Eigen::VectorXd &, const Eigen::VectorXd &)> f;
      f = [&](const Eigen::VectorXd &UrColi, const Eigen::VectorXd &UcColi){
        double e = 0.5*(UrColi+UcColi-UCol0-dt*UdCol0).transpose()*M/(dt*dt)*(UrColi+UcColi-UCol0-dt*UdCol0);
        if(physic_model=="arap"){
          return e + sim::linear_tetmesh_arap_q(V, T, VCol+UrColi+UcColi, dX, vol, params);
        }
        else if(physic_model=="neohookean"){
          return e + sim::linear_tetmesh_neohookean_q(V, T, VCol+UrColi+UcColi, dX, vol, params);
        }
        else if(physic_model=="stvk"){
          return e + sim::linear_tetmesh_stvk_q(V, T, VCol+UrColi+UcColi, dX, vol, params);
        }
        else if(physic_model=="corotational"){
          return e + sim::linear_tetmesh_corotational_q(V, T, VCol+UrColi+UcColi, dX, vol, params);
        }
      };
      double alpha = line_search(f,tmp_g,dUc,UrCol,UcCol);
      UcCol = UcCol + alpha * dUc;
    }
    UCol = UrCol + UcCol;
    UdCol = (UCol-UCol0)/dt;

    Eigen::MatrixXd Vn = V + matrixize(UCol);
    Vn_list[ai] = Vn;

  }
  if(export_objs)
  {
    const std::string dir = igl::dirname(json_path) + PATH_SEPARATOR;;
    size_t dot_found = json_path.find_last_of(".");
    std::string model = std::string(json_path.begin()+dir.size(),json_path.begin()+dot_found);

    std::string output_folder = "../output";
    if(!std::filesystem::exists(output_folder)){
      std::filesystem::create_directory(output_folder);
    }
    std::string output_path = output_folder+"/"+model+"_"+physic_model;
    if(!std::filesystem::exists(output_path)){
      std::filesystem::create_directory(output_path);
    }
    for(int ai=0; ai<Vn_list.size(); ai++)
    {
      Eigen::MatrixXd Vn = Vn_list[ai];
      std::string numstr = std::to_string(ai);
      int strlen = numstr.length();

      std::string seqstr;
      for(int k=0; k<(4-strlen); k++){
          seqstr.append("0");
      }
      seqstr.append(numstr);
      std::string save_path = output_path+"/"+model+"_"+physic_model+"_"+seqstr+".obj";
      igl::writeOBJ(save_path,Vn,F);
    }
  }


	igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.core().align_camera_center(V);
  viewer.core().is_animating = true;
  viewer.core().animation_max_fps = 24.;
  viewer.data().show_lines = true;
  viewer.data().show_overlay_depth = false;
  viewer.core().background_color=Eigen::Vector4f::Ones();
  viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer &) -> bool
	{
		if(viewer.core().is_animating)
		{
      Eigen::MatrixXd Vn = Vn_list[frame];
      viewer.data().set_vertices(Vn);
      viewer.data().compute_normals();
			frame++;
      if(frame == Vn_list.size()){
        //viewer.core().is_animating = !viewer.core().is_animating;
        frame = 0;
      }
		}
		return false;
	};


  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if(key==' '){
      viewer.core().is_animating = !viewer.core().is_animating;
    }
    if(key=='e' || key=='E'){
      std::cout<<"export png seq"<<std::endl;
      viewer.core().is_animating = true;
      export_png_seq(viewer,Vn_list,json_path,physic_model);
    }
    return false;
  };
  viewer.launch();
}
