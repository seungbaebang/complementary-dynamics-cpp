#include "line_search.h"
#include <iostream>

// #include <igl/matlab/matlabinterface.h>

// double line_search(
//     const std::function<double(const Eigen::VectorXd &)> & f,
//     const std::function<void(Eigen::VectorXd &)> & proj_z,
//     const Eigen::VectorXd & z,
//     const Eigen::VectorXd & dz,
//     const double max_step) {
//     /////////////////////////////////////////////////////////////////////////////
//     double optimal_step_sigma = max_step;
//     double Energy = f(z);
//     Eigen::VectorXd next_step = z+optimal_step_sigma*dz;
//     proj_z(next_step); 
    
//     while (f(next_step) > Energy && optimal_step_sigma != 0){
//         // we decrease σ by a constant factor
//         optimal_step_sigma = optimal_step_sigma*0.5;
//         next_step = z-optimal_step_sigma*dz;
//         proj_z(next_step);
//     }
//     return optimal_step_sigma;
//     /////////////////////////////////////////////////////////////////////////////
// }

double line_search(
    const std::function<double(const Eigen::VectorXd &, const Eigen::VectorXd &)> & f,
    const Eigen::VectorXd & tmp_g,
    const Eigen::VectorXd & dUc,
    const Eigen::VectorXd & Ur,
    const Eigen::VectorXd & Uc) 
    {
    double alpha = 1;
    double p = 0.5;
    double c = 1e-8;

    double f0 = f(Ur,Uc);
    double s = f0 + c * tmp_g.transpose() * dUc; 


    while (alpha > c){
      Eigen::VectorXd Uc_tmp = Uc + alpha * dUc;

      if (f(Ur,Uc_tmp) <= s){
        break;
      }
      alpha = alpha * p;
    }
    return alpha;

}