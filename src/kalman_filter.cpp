#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;
// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in,
                        MatrixXd &P_in,
                        MatrixXd &F_in,
                        MatrixXd &HL_in,
                        MatrixXd &RL_in,
                        MatrixXd &RR_in,
                        MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    HL_ = HL_in;
    RL_ = RL_in;
    RR_ = RR_in;
    Q_ = Q_in;
    I_ = MatrixXd::Identity(4, 4);
}

void KalmanFilter::Predict() {
    /**
      predict the state
    */

    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
    /**
      * update the state by using Kalman Filter equations
    */
    // Update the distribution
    MatrixXd HLt_ = HL_.transpose();

    VectorXd y_ = z - (HL_ * x_);
    MatrixXd S = (HL_ * P_ * HLt_) + RL_;
    MatrixXd K_ = P_ * HLt_ * S.inverse();

    // Update the state
    x_ = x_ + (K_ * y_);
    P_ = (I_ - K_ * HL_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
      * update the state by using Extended Kalman Filter equations
    */
    if (x_[0] == 0 and x_[1] == 0) return;

    VectorXd h_x = tools.CalculateCoordinates(x_);
    MatrixXd Hj_ = tools.CalculateJacobian(x_);
    MatrixXd Hjt_ = Hj_.transpose();

    // Compute the error in state vector
    VectorXd y_ = z - h_x;

    // Normalizing y
    tools.NormalizePhi(y_);


    // Update the distribution
    MatrixXd S = Hj_ * P_ * Hjt_ + RR_;
    MatrixXd K_ = P_ * Hjt_ * S.inverse();

    // Update the state
    x_ = x_ + (K_ * y_);
    P_ = (I_ - (K_ * Hj_)) * P_;
}
