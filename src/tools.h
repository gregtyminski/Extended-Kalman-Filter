#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

class Tools {
public:
    /**
     * Constructor.
     */
    Tools();

    /**
     * Destructor.
     */
    virtual ~Tools();

    /**
     * A helper method to calculate RMSE.
     */
    Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                  const std::vector<Eigen::VectorXd> &ground_truth);

    /**
     * A helper method to calculate Jacobians.
     */
    Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd &x_state);

    /**
     * Method converting Polar coordinates to Cartesian ones.
     * @param measurement_pack
     * @return
     */
    Eigen::VectorXd Polar2Cartesian(const MeasurementPackage &measurement_pack);

    /**
     * Verifies if value is smaller than Tools::epsilon
     * @param value
     * @return
     */
    double verifyZero(double value);
};

#endif  // TOOLS_H_
