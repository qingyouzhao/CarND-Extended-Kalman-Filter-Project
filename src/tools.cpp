#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;


Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   * done
   */
  using namespace std;
  // initialize return value
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  const bool b_valid_data = estimations.size() != 0 && estimations.size() == ground_truth.size();
  // 1st check for zero
  if (!b_valid_data)
  {
    cout << "Calculate RMSE require the estimation and ground truth to be the same non-zero size" << endl;
    assert(b_valid_data);
  }

  // Accumulate square error
  for (unsigned int i = 0; i < estimations.size(); i++)
  {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // mean
  rmse /= estimations.size();
  // root
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   * done
   */
  MatrixXd Hj(3, 4);
  
  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // check division by zero
  const bool b_devide_by_zero = (px * px + py * py) <= 1e-8;
  if (b_devide_by_zero)
  {
    cout << "Division by zero, returning initialized matrix" << endl;
    assert(!b_devide_by_zero);
  }
  // compute the Jacobian matrix
  const double pxpy2 = px * px + py * py;
  const double pxpy3_2 = pow(pxpy2, 1.5f);
  Hj << px / sqrt(pxpy2),             py / sqrt(pxpy2),               0,                0,
        -py / pxpy2,                  px / pxpy2,                     0,                0,
        py*(vx*py - vy*px) / pxpy3_2, px*(vy*px - vx * py) / pxpy3_2, px / sqrt(pxpy2), py / sqrt(pxpy2);

  return Hj;
}
