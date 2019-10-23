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
   */
  VectorXd RMSE(4);
  RMSE<<0,0,0,0;
  	// check the validity of the following inputs:
	// the estimation vector size should not be zero
	// the estimation vector size should equal ground truth vector size
    if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
    cout << "CalculateRMSE() - Error - Invalid inputs." << endl;
    return RMSE;
    }
  	// Accumulate squared residuals
  	for(int i = 0; i < estimations.size(); i++){
    	VectorXd residual = estimations[i] - ground_truth[i];
    	residual = residual.array()*residual.array(); //coefficient-wise multiplication
    	RMSE += residual;
    }
  	// Compute the mean
  	RMSE = RMSE/estimations.size();

  	// Compute the squared root
  	RMSE = RMSE.array().sqrt();

  return RMSE;  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  float d1 = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  float el00 = x_(0)/d1;
  float el01 = x_(1)/d1;
  float el02 = 0.0;
  float el03 = 0.0;
  float el10 = -x_(1)/(d1*d1);
  float el11 = x_(0)/(d1*d1);
  float el12 = 0.0;
  float el13 = 0.0;
  float el20 = (x_(1)*(x_(2)*x_(1) - x_(3)*x_(0))) / d1*d1*d1;
  float el21 = -x_(0)*el20/x_(1);
  float el22 = el00;
  float el23 = el01;
  MatrixXd Hj(3,4);
  Hj<<el00, el01, el02, el03,
  	  el10, el11, el12, el13,
  	  el20, el21, el22, el23;
  return Hj;
}
