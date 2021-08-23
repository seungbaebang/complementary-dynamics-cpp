#ifndef READ_DATA_FROM_JSON_H
#define READ_DATA_FROM_JSON_H

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

bool read_json_data(const std::string& filename, 
Eigen::MatrixXd& V, Eigen::MatrixXi& T, Eigen::MatrixXi& F,
Eigen::MatrixXd& C, Eigen::VectorXi& PI, Eigen::MatrixXi& BE,
Eigen::MatrixXd& W, std::vector<Eigen::MatrixXd>& T_list,
double& dt, double& YM, double& pr, double& scale, std::string& physic_model);


#endif