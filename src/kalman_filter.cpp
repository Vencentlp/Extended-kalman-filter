#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
    * predict the state
  */
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	/*
	  * update the state by using Kalman Filter equations
	*/
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd S_inv = S.inverse();
	MatrixXd K = P_ * Ht * S_inv;
	x_ = x_ + (K * y);
	long size = x_.size();
	MatrixXd I_ = MatrixXd::Identity(size, size);
	P_ = (I_ - K * H_)*P_;
}
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
	
	float px = x_(0);
	float py = x_(1);
	float vx = x_(2);
	float vy = x_(3);
	//Pre-computer a set of terms to avoid repeat calculatio	
	
	float c1 = px * px + py * py;
	if (c1 < 0.0001)
	{	
	px = 0.0001;
	py = 0.0001;
	};	
	if (px < 0.0001){
	px = 0.0001;
	};	
	float c2 = sqrt(c1);
	float c3 = atan2(py, px);
	
	VectorXd hc_(3);
	hc_ << c2, c3, (px*vy + py * vx) / c2;
	//float r = z(0);
	//float phi = z(1);
	//float r_r = z(2);
	//VectorXd z_(3);
	//z << r, phi, r_r;
	VectorXd y = z - hc_;
	while ((y(1)<-M_PI) || (y(1)>M_PI))
	{
     		if (y(1) < -M_PI)
		{
			y(1) = y(1) + M_PI;
		}
		else
		{
			y(1) = y(1) - M_PI;
		}
	}
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd S_inv = S.inverse();
	MatrixXd K = P_ * Ht * S_inv;
	x_ = x_ + (K * y);
	long size = x_.size();
	MatrixXd I_ = MatrixXd::Identity(size, size);
	P_ = (I_ - K * H_)*P_;
}
