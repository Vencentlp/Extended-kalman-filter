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
	
	double px = x_(0);
	double py = x_(1);
	double vx = x_(2);
	double vy = x_(3);
	//Pre-computer a set of terms to avoid repeat calculatio	
	double phi;
	double rho;
	double rho_dot;
	
	double rho = sqrt(px*px + py * py);

	// avoid division by zero
	if (fabs(px) < 0.0001) {
		cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
	}
	else {
		phi = atan2(py, px);
	}

	// avoid division by zero
	if (rho < 0.0001) {
		cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
	}
	else {
		rho_dot = (px*vx + py * vy) / rho;
	}

	VectorXd h = VectorXd(3);
	h << rho, phi, rho_dot; // For radar H * x becomes h(x)

	VectorXd y = z - h; // Using h instead of Jacobian Hj_ here!

	
	//double r = z(0);
	//double phi = z(1);
	//double r_r = z(2);
	//VectorXd z_(3);
	//z << r, phi, r_r;
	//VectorXd y = z - hc_;
	while ((y(1)<-M_PI) || (y(1)>M_PI))
	{
     		if (y(1) < -M_PI)
		{
			y(1) = y(1) + 2*M_PI;
		}
		else
		{
			y(1) = y(1) - 2*M_PI;
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
