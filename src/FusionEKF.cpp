#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    /**
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
     */
    ekf_ = KalmanFilter();
    noise_ax = 9;
    noise_ay = 9;

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    tools = Tools();
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    /**
     * Initialization
     */
    if (!is_initialized_) {
        /**
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * You'll need to convert radar from polar to cartesian coordinates.
         */

        // first measurement
        cout << "EKF: " << endl;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            // Convert radar from polar to cartesian coordinates and initialize state.
            cout << "EKF init RADAR" << endl;
            ekf_.x_ = tools.Polar2Cartesian(measurement_pack);
            ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
            ekf_.R_ = R_radar_;
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            // Initialize state.
            cout << "EKF init LASER" << endl;
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
            ekf_.H_ = H_laser_;
            ekf_.R_ = R_laser_;
        }

        // done initializing, no need to predict or update
        previous_timestamp_ = measurement_pack.timestamp_;
        is_initialized_ = true;
        return;
    }

    /**
     * Prediction
     */

    /**
     * Update the state transition matrix F according to the new elapsed time.
     * Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    // find time elapsed in seconds
    double dt =
            (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // divide by 1M to change from [ms] to [s]
    previous_timestamp_ = measurement_pack.timestamp_;

    // init F
    // 1. Modify the F matrix so that the time is integrated
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, dt, 0,
            0, 1, 0, dt,
            0, 0, 1, 0,
            0, 0, 0, 1;

    // 2. Set the process covariance matrix Q
    float dt_4, dt_3, dt_2;
    dt_2 = pow(dt, 2);
    dt_3 = dt_2 * dt;
    dt_4 = dt_3 * dt;
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << (noise_ax * dt_4 / 4), 0, (noise_ax * dt_3 / 2), 0,
            0, (noise_ay * dt_4 / 4), 0, (noise_ay * dt_3 / 2),
            (noise_ax * dt_3 / 2), 0, (noise_ax * dt_2), 0,
            0, (noise_ay * dt_3 / 2), 0, (noise_ay * dt_2);

    ekf_.Predict();

    /**
     * Update
     */

    /**
     * - Use the sensor type to perform the update step.
     * - Update the state and covariance matrices.
     */
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
    return;
}
