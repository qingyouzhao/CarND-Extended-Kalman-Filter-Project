#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd F_t = F_.transpose();
  P_ = F_ * P_ * F_t;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */

  VectorXd y = z - H_ * x_;
  MatrixXd H_t = H_.transpose();
  MatrixXd S = H_ * P_ * H_t + R_;
  MatrixXd S_inverse = S.inverse();
  MatrixXd K = P_ * H_t * S_inverse;
  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(K.rows(), K.cols()) - K * H_)*P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  MatrixXd H_j = Tools::CalculateJacobian(x_);
  VectorXd h_prime = LinearStateToPolar(x_);
  VectorXd y = z - h_prime;

  MatrixXd H_j_t = H_j.transpose();
  MatrixXd S = H_j * P_ * H_j_t + R_;
  MatrixXd S_inverse = S.inverse();
  MatrixXd K = P_ * H_j_t * S_inverse;
  
  x_ = x_ + K * y;
  P_ = (MatrixXd::Identity(K.rows(), K.cols()) - K * H_)*P_;
}

VectorXd KalmanFilter::LinearStateToPolar(const VectorXd &x)
{
  // Extract x with variable name for readability
  const float & px = x(0);
  const float & py = x(1);
  const float & vx = x(2);
  const float & vy = x(3);

  VectorXd h(3);
  h(0) = sqrt(px *px + py*py);
  h(1) = atan(py/px);
  h(2) = (px*vx+py*vy)/h(0);

  return h;
}
