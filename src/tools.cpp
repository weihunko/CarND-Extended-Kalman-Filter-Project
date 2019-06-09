#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	if(estimations.size() != ground_truth.size() || estimations.size() == 0){
		std::cout << "invalid estimations or ground truth data" << std::endl;
		return rmse;
	}

	for(int i=0; i<estimations.size(); i++){
		VectorXd residual = estimations[i] - ground_truth[i];
		// std::cout << "estimations x " << estimations[i](0) << std::endl;
		// std::cout << "ground truth x " <<  ground_truth[i](0) << std::endl;
		residual = residual.array()*residual.array();
		rmse = rmse + residual;
	}

	rmse = rmse/estimations.size();
	rmse = rmse.array().sqrt();
	std::cout << "x rmse " << rmse(0) << std::endl;
	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
	MatrixXd Hj(3,4);
 	// recover state parameters
	float px = x_state(0);
	float py = x_state(1);
 	float vx = x_state(2);
	float vy = x_state(3);

  // check division by zero
  float px2py2 = px*px + py*py;
  if(fabs(px2py2) < 0.0001){
  	std::cout << "invalid division by zero, return zero Hj" << std::endl;
   	return Hj;
  }
  // compute the Jacobian matrix
	float px2py2_square = sqrt(px2py2);
	float px2py2_sq1_5  = pow(px2py2, 1.5);
	Hj(0, 0) = px/px2py2_square;
	Hj(0, 1) = py/px2py2_square;
	Hj(0, 2) = 0;
	Hj(0, 3) = 0;
	Hj(1, 0) = -py/px2py2;
	Hj(1, 1) = px/px2py2;
	Hj(1, 2) = 0;
	Hj(1, 3) = 0;
	Hj(2, 0) = py*(vx*py - vy*px)/px2py2_sq1_5;
	Hj(2, 1) = px*(vy*px - vx*py)/px2py2_sq1_5;
	Hj(2, 2) = px/px2py2_square;
	Hj(2, 3) = py/px2py2_square;

	return Hj;
}
