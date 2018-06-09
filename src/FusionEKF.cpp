#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matric
    x_ = VectorXd(4);
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_      = MatrixXd(3, 4);
    P_ = MatrixXd(4, 4);
    F_ = MatrixXd(4, 4);
    Q_ = MatrixXd(4, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    /**
    TODO:
      * Finish initializing the FusionEKF.
      * Set the process and measurement noises
    */

    // Initialize measurement matrix for laser measurements
    H_laser_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    // Initialize transition matrix
    F_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;

    // Initializing emtpy matrix Q_

    //Set the process and measurement noises
    noise_ax = 9.0;
    noise_ay = 9.0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

    cout << "Got Measurement: " << measurement_pack.raw_measurements_ << endl;

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            double rho = measurement_pack.raw_measurements_[0];
            double phi = measurement_pack.raw_measurements_[1];

            x_ << cos(phi) * rho, sin(phi) * rho, 0.0, 0.0;

        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x_ <<   measurement_pack.raw_measurements_[0],
                    measurement_pack.raw_measurements_[1],
                    0.0,
                    0.0;
        }

        // Intializing EKF object
        ekf_.Init(x_, /*x_in*/
                  P_, /*P_in*/
                  F_, /*F_in*/
                  H_laser_, /*H_in*/
                  R_laser_, /*R_in*/
                  R_radar_, /*R_in*/
                  Q_); /*Q_in*/

        // Recording the timestamp
        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    // Computing delta in time
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    //
    if (dt > 0.00) {

        // Updating the state transition matrix F
        ekf_.F_(0, 2) = dt;
        ekf_.F_(1, 3) = dt;

        // Update the process noise covariance matrix.
        float dt2 = dt * dt;
        float dt3 = dt2 * dt;
        float dt4 = dt3 * dt;
        float dt4over4 = dt4 / 4.;
        float dt3over2 = dt3 / 2.;
        ekf_.Q_ << dt4over4 * noise_ax, 0, dt3over2 * noise_ax, 0,
                0, dt4over4 * noise_ay, 0, dt3over2 * noise_ay,
                dt3over2 * noise_ax, 0, dt2 * noise_ax, 0,
                0, dt3over2 * noise_ay, 0, dt2 * noise_ay;

        ekf_.Predict();
    }

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

        // Radar updates
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    } else {

        // Laser updates
        ekf_.Update(measurement_pack.raw_measurements_);

    }

    // print the output
    //cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
}
