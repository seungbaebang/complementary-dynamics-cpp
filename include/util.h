#ifndef UTIL_H
#define UTIL_H



void emu_to_lame(const double& E, const double& p, double& lambda, double& mu)
{
  lambda = (E*p)/((1.0+p)*(1.0-2.0*p));
  mu = E/(2.0*(1.0+p));
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

#endif 