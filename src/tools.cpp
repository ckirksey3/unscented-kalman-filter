#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() < 1) {
    throw std::invalid_argument( "Error: Estimation vector size should not be zero");
  }
  else if(estimations.size() != ground_truth.size()) {
    throw std::invalid_argument( "Error: Estimation vector size should equal ground truth size");
  } else {
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
      VectorXd current_estimate = estimations[i];
      VectorXd cur_ground_truth = ground_truth[i];
      VectorXd residual = current_estimate - cur_ground_truth;
      for(int j=0; j<residual.size(); ++j) {
        residual[j] = residual[j] * residual[j];
      }
      rmse += residual;
    }

    //calculate the mean
    rmse /= estimations.size();

    //calculate the squared root
    for(int i=0; i<rmse.size(); ++i) {
      rmse[i] = sqrt(rmse[i]);
    }

    //return the result
    return rmse;
  }
}
