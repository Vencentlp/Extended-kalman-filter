#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
	VectorXd RMSE(4);
	RMSE << 0, 0, 0, 0;
	if (estimations.size() != ground_truth.size()) 
	{
		cout << "Error - Calculate RMSE () - The groud truth and estimation vector must have same size";
		return RMSE;
	}
	for (uint i = 0; i < estimations.size(); i++)
	{
		VectorXd Diff = estimations[i] - ground_truth[i];
		VectorXd residual = Diff.array() * Diff.array();
		RMSE = RMSE + residual;
	}
	return RMSE;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
	MatrixXd Hj(3, 4);
	//VectorXd h_c(3);
  //Recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
	//Pre-computer a set of terms to avoid repeat calculation
	float c1 = px * px + py * py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);
	
	if (c1 < 0.0001) 
	{
		return Hj;
	}
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy * px) / c3, px*(px*vy - py * vx) / c3, px / c2, py / c2;
	return Hj;


}
