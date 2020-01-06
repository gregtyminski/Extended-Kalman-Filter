#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
        // cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    // accumulate squared residuals
    for (unsigned int i = 0; i < estimations.size(); ++i) {
        VectorXd residual = estimations[i] - ground_truth[i];
        residual = residual.array() * residual.array();
        rmse += residual;
    }

    // calculate the mean
    rmse = rmse / estimations.size();
    // calculate the squared root
    rmse = rmse.array().sqrt();
    // return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
    MatrixXd Hj(3, 4);
    // recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    float px_py_2 = px * px + py * py;
    float sqrt_px_py_2 = sqrt(px_py_2);
    float vxpy_vypx = vx * py - vy * px;
    float vypx_vxpy = -vxpy_vypx;
    float px_py_2_32 = pow(px_py_2, 3) / sqrt_px_py_2;

    // check division by zero
    if (fabs(px_py_2) < 0.00001) {
        return Hj;
    }

    // compute the Jacobian matrix
    Hj <<   (px / sqrt_px_py_2),            (py / sqrt_px_py_2),            0, 0,
            (-py / px_py_2),                (px / px_py_2),                 0, 0,
            (py * vxpy_vypx / px_py_2_32),  (px * vypx_vxpy / px_py_2_32),  (px / sqrt_px_py_2), (py / sqrt_px_py_2);

    return Hj;
}

VectorXd Tools::Polar2Cartesian(const MeasurementPackage &measurement_pack) {
    VectorXd result = VectorXd(4);
    double rho = measurement_pack.raw_measurements_[0]; // range
    double phi = measurement_pack.raw_measurements_[1]; // bearing
    double rho_dot = measurement_pack.raw_measurements_[2]; // velocity

    double x = verifyZero(rho * cos(phi));
    double y = verifyZero(rho * sin(phi));

    double vx = rho_dot * cos(phi);
    double vy = rho_dot * sin(phi);

    result << x, y, vx, vy;
    return result;
}

double Tools::verifyZero(double value) {
    double treshold = 0.001;
    if (value < treshold){
        return treshold;
    }
    return value;
}
