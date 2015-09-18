//
// Created by Liangsheng Zhang on 9/18/15.
//

#ifndef MBL_V1_DISORDER_TRANSITION_FUNC_H
#define MBL_V1_DISORDER_TRANSITION_FUNC_H

/**
 ** It implements various functions useful for disorder_transition
 **/

#include <vector>
#include <utility>
#include <Eigen/Dense>

void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::VectorXd& vector);
void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::MatrixXd& matrix);
void Middle_Half_Extract(const std::vector< std::pair<double,int> >& pair_index,  Eigen::MatrixXcd& matrix);

#endif //MBL_V1_DISORDER_TRANSITION_FUNC_H
