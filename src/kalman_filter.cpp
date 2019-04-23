#include "kalman_filter.h"
#include <iostream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

#define PRINT_INTERMEDIATE_STEPS 0

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
   * TODO: predict the state, done
   */
  x_ = F_ * x_;
  MatrixXd F_t = F_.transpose();
  P_ = F_ * P_ * F_t + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations, done
   */
  cout << "Update: " << endl;
  VectorXd y = z - H_ * x_;
#if PRINT_INTERMEDIATE_STEPS
  cout << "Update: y \n" << y << endl;
#endif
  MatrixXd H_t = H_.transpose();
#if PRINT_INTERMEDIATE_STEPS
  cout << "Update: Ht \n" << H_t << endl;
  cout << "Update: P_ \n" << P_ << endl;
  cout << "Update: R_ \n" << R_ << endl;
#endif
  MatrixXd S = H_ * P_ * H_t + R_;

#if PRINT_INTERMEDIATE_STEPS
  cout << "Update: S \n" << S << endl;
#endif
  MatrixXd S_inverse = S.inverse();
  MatrixXd K = P_ * H_t * S_inverse;
#if PRINT_INTERMEDIATE_STEPS
  cout << "Update: K \n" << K << endl;
#endif

  x_ = x_ + K * y;
#if PRINT_INTERMEDIATE_STEPS
  cout << "Update: x_ \n" << x_ << endl;
#endif
  MatrixXd I = MatrixXd::Identity(P_.rows(), P_.cols());

#if PRINT_INTERMEDIATE_STEPS
  cout << "Update: I \n" << I << endl;
  cout << "Update: K \n" << K << endl;
  cout << "Update: H_ \n" << H_ << endl;
#endif

  P_ = (I - K * H_)*P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations, done
   */
  cout << "UpdateEKF: " << endl;

  VectorXd h_prime = LinearStateToPolar(x_);
  VectorXd y = z - h_prime;
  // normalize y 
  while (y[1] < -M_PI)
  {
    y[1] += 2.0*M_PI;
  }
  while (y[1] >= M_PI)
  {
    y[1] -= 2.0*M_PI;
  }

  MatrixXd H_t = H_.transpose();
  MatrixXd S = H_ * P_ * H_t + R_;
  MatrixXd S_inverse = S.inverse();
  MatrixXd K = P_ * H_t * S_inverse;
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(P_.rows(), P_.cols());
  P_ = (I - K * H_)*P_;
}

VectorXd KalmanFilter::LinearStateToPolar(const VectorXd &x)
{
  // Extract x with variable name for readability
  const double & px = x(0);
  const double & py = x(1);
  const double & vx = x(2);
  const double & vy = x(3);

  VectorXd h(3);
  h(0) = sqrt(px *px + py*py);
  h(1) = atan2(py,px);
  h(2) = (px*vx+py*vy)/h(0);

  return h;
}
