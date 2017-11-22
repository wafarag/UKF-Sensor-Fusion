#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // Set flag to false
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Initialize measurement dimension, radar can measure r, phi, and r_dot
  n_z_radar_ = 3;

  // Initialize measurement dimension, laser can measure positions Px and Py
  n_z_laser_ = 2;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // measurement covariance matrix S for Radar
  S_radar_ = MatrixXd(n_z_radar_, n_z_radar_);

  // measurement covariance matrix S for Laser
  S_laser_ = MatrixXd(n_z_laser_, n_z_laser_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_      = 1.0;        // Tuned using Trial and Error

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_  = 0.6;        // Tuned using Trial and Error

  // Laser measurement noise standard deviation position1 in m
  std_laspx_  = 0.15;       // Kept as is

  // Laser measurement noise standard deviation position2 in m
  std_laspy_  = 0.15;       // Kept as is

  // Radar measurement noise standard deviation radius in m
  std_radr_   = 0.3;        // Kept as is

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;      // Kept as is

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_  = 0.3;       // Kept as is

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Initialize sigma point matrix
  Xsig_ = MatrixXd(n_x_, (2 * n_x_ + 1));
  Xsig_.fill(0.0);

  // Initialize augmented sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, (2 * n_aug_ + 1));
  Xsig_aug_.fill(0.0);

  // Initialize matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, (2 * n_aug_ + 1));
  Xsig_pred_.fill(0.0);

  // create matrix for sigma points in measurement space
  Zsig_radar_ = MatrixXd(n_z_radar_, (2 * n_aug_ + 1));
  Zsig_radar_.fill(0.0);

  // create matrix for sigma points in measurement space
  Zsig_laser_ = MatrixXd(n_z_laser_, (2 * n_aug_ + 1));
  Zsig_laser_.fill(0.0);

  // Initialize mean predicted measurement vector
  z_pred_radar_ = VectorXd(n_z_radar_);
  z_pred_radar_.fill(0.0);

  // Initialize mean predicted measurement vector for Laser
  z_pred_laser_ = VectorXd(n_z_laser_);
  z_pred_laser_.fill(0.0);

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<(2*n_aug_+1); i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // Radar measurement noise covariance matrix
  R_radar_ = MatrixXd(n_z_radar_,n_z_radar_);
  //measurement covariance matrix - radar
  // Tuned by trial and error
  R_radar_ << std_radr_*std_radr_, 0,                       0,
              0,                   std_radphi_*std_radphi_, 0,
              0,                   0,                       std_radrd_*std_radrd_;

  // measurement covariance matrix - laser
  R_laser_   = MatrixXd(n_z_laser_, n_z_laser_);
  // Tuned by trial and error
  R_laser_ << std_laspx_*std_laspx_  ,0,
              0                      ,std_laspy_*std_laspy_;

   // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  // Initialize the State Covariance Matrix
  // Tuned by Trial and Error
  float p11 = 0.1;
  float p22 = 0.005;
  float p33 = 0.05;
  float p44 = 0.01;
  float p55 = 0.05;

  P_ <<    p11,  0,    0,    0,    0,
           0,    p22,  0,    0,    0,
           0,    0,    p33,  0,    0,
           0,    0,    0,    p44,  0,
           0,    0,    0,    0,    p55;

   // Initialize NIS for radar
  NIS_radar_ = 0.0;

  // Initialize NIS for laser
  NIS_laser_ = 0.0;

  // Initialize the summation of radar NIS records
  NIS_radar_sum = 0.0;

  // Initialize the counter of radar readings
  Radar_records_counter = 0;

  // Initialize the summation of laser NIS records
  NIS_laser_sum = 0.0;

  // Initialize the counter of laser readings
  Laser_records_counter = 0;

  // Initialize the counter of laser readings above DOF threshold
  Laser_DOF_threshold_counter = 0;

  // Initialize the counter of radar readings above DOF threshold
  Radar_DOF_threshold_counter = 0;

  // Initialize sampling time for measurements
  delta_t = 0.1;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /// ****************************************************************
  /// ***************** Initialization Step **************************
  /// ****************************************************************
  if (!is_initialized_) {

  /// Initialize the states for both RADAR and LiDAR
  if ((measurement_pack.sensor_type_ == MeasurementPackage::RADAR) and use_radar_) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

     /** px is adjacent to theta and is calculated as rho * cos(theta)
      py is opposed to theta and is calculated as rho * sin(theta)   */

      float roa = measurement_pack.raw_measurements_[0];
      float theta = measurement_pack.raw_measurements_[1];
      float roa_dot = measurement_pack.raw_measurements_[2];

      float Px = roa * cos(theta); //
      float Py = roa * sin(theta); //

      //cout<< "theta = " << theta;

      /** vx = ro_dot * cos(phi); and vy = ro_dot * sin(phi)    */

      float Vx = roa_dot * cos(theta); // Approximate Initialization
      float Vy = roa_dot * sin(theta); // Approximate Initialization
      float Vtotal = sqrt(Vx*Vx+Vy*Vy);
      float yaw_angle = 0.0;            // in rad
      float yaw_rate  = 0.0;

      x_ << Px, Py, Vtotal, yaw_angle, yaw_rate;

      Radar_records_counter++;

      // Update time-stamp to be used for the next iteration
      previous_timestamp_ = measurement_pack.timestamp_;
      // done initializing, no need to predict or update
      is_initialized_ = true;
    }
    else if ((measurement_pack.sensor_type_ == MeasurementPackage::LASER) and use_laser_) {
      /**
      Initialize state.
      */
      float Px = measurement_pack.raw_measurements_[0];
      float Py = measurement_pack.raw_measurements_[1];

      float Vx =  5.2;  // Tuning by trial and error
      float Vy =  0.05; // Tuning by trial and error
      float Vtotal = sqrt(Vx*Vx+Vy*Vy);
      float yaw_angle = 0.0;
      float yaw_rate  = 0.0;

      x_ << Px, Py, Vtotal, yaw_angle, yaw_rate;

      Laser_records_counter++;

      // Update time-stamp to be used for the next iteration
      previous_timestamp_ = measurement_pack.timestamp_;
      // done initializing, no need to predict or update
      is_initialized_ = true;

    }

    return;
  }

  /// ****************************************************************
  /// ***************** Prediction Step ******************************
  /// ****************************************************************

  // Calculate Delta Time (dt) in Seconds
  //delta_t = (measurement_pack.timestamp_ - previous_timestamp_)/1e6;
  // Update time-stamp to be used for the next iteration
  //previous_timestamp_ = measurement_pack.timestamp_;

  if ((measurement_pack.sensor_type_ == MeasurementPackage::RADAR) and use_radar_) {

     // Calculate Delta Time (dt) in Seconds
     delta_t = (measurement_pack.timestamp_ - previous_timestamp_)/1e6;
     // Update time-stamp to be used for the next iteration
     previous_timestamp_ = measurement_pack.timestamp_;

     Prediction(measurement_pack);
   }
   else if ((measurement_pack.sensor_type_ == MeasurementPackage::LASER) and use_laser_) {

     // Calculate Delta Time (dt) in Seconds
     delta_t = (measurement_pack.timestamp_ - previous_timestamp_)/1e6;
     // Update time-stamp to be used for the next iteration
     previous_timestamp_ = measurement_pack.timestamp_;

     Prediction(measurement_pack);
   }

  /// ****************************************************************
  /// ***************** Update Step **********************************
  /// ****************************************************************
  if ((measurement_pack.sensor_type_ == MeasurementPackage::RADAR) and use_radar_)
  {
    // Radar updates
    UpdateRadar(measurement_pack);

  } else if ((measurement_pack.sensor_type_ == MeasurementPackage::LASER) and use_laser_)
  {
    // Laser updates
    UpdateLidar(measurement_pack);
  }

}
///*/////////////////////////////////////////////////////////////////////////////////
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(MeasurementPackage measurement_pack) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //GenerateSigmaPoints();
  AugmentedSigmaPoints();
  SigmaPointPrediction();
  PredictMeanAndCovariance();

}
///*/////////////////////////////////////////////////////////////////////////////////
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Laser updates
  PredictLaserMeasurement();
  //
  UpdateState_Laser(meas_package);
}
///*/////////////////////////////////////////////////////////////////////////////////
/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Radar updates
  PredictRadarMeasurement();
  //
  UpdateState_Radar(meas_package);
}
///*/////////////////////////////////////////////////////////////////////////////////
void UKF::GenerateSigmaPoints() {

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig_.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig_.col(i+1)      = x_ + sqrt(lambda_ + n_x_) * A.col(i);
    Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
  }

  //print result
  //std::cout << "Xsig = " << std::endl << Xsig << std::endl;

}
///*/////////////////////////////////////////////////////////////////////////////////
void UKF::AugmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug_.head(n_x_) = x_;
  x_aug_(5)         = 0;
  x_aug_(6)         = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //print result
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
}
///*/////////////////////////////////////////////////////////////////////////////////
void UKF::SigmaPointPrediction() {

  //predict sigma points
  for (int i = 0; i< (2*n_aug_+1); i++)
  {
    //extract values for better readability
    double p_x      = Xsig_aug_(0,i);
    double p_y      = Xsig_aug_(1,i);
    double v        = Xsig_aug_(2,i);
    double yaw      = Xsig_aug_(3,i);
    double yawd     = Xsig_aug_(4,i);
    double nu_a     = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + (v/yawd) * ( sin(yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + (v/yawd) * ( cos(yaw) - cos(yaw + yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p    = v;
    double yaw_p  = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p  = v_p  + nu_a*delta_t;

    yaw_p  = yaw_p  + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

}

///*/////////////////////////////////////////////////////////////////////////////////
void UKF::PredictMeanAndCovariance() {

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.0*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.0*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  //std::cout << "Predicted state" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
  //std::cout << P_ << std::endl;

}

///*/////////////////////////////////////////////////////////////////////////////////
void UKF::PredictRadarMeasurement() {

  //transform sigma points into measurement space
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points

    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double vx = cos(yaw)*v;
    double vy = sin(yaw)*v;

    // measurement model
    Zsig_radar_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_radar_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_radar_(2,i) = (p_x*vx + p_y*vy) / Zsig_radar_(0,i);           //r_dot
    //Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  z_pred_radar_.fill(0.0);
  for (int i=0; i < (2*n_aug_+1); i++) {
      z_pred_radar_ = z_pred_radar_ + weights_(i) * Zsig_radar_.col(i);
  }

  S_radar_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.0*M_PI;

    S_radar_ = S_radar_ + weights_(i) * z_diff * z_diff.transpose();
  }

  S_radar_ = S_radar_ + R_radar_;

  //print result
  //std::cout << "z_pred: " << std::endl << z_pred_radar_ << std::endl;
  //std::cout << "S: " << std::endl << S_radar_ << std::endl;

}

///*/////////////////////////////////////////////////////////////////////////////////
void UKF::PredictLaserMeasurement() {

  //transform sigma points into measurement space
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points

    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig_laser_(0,i) = p_x;         //Px
    Zsig_laser_(1,i) = p_y;         //Py
  }

  z_pred_laser_.fill(0.0);
  for (int i=0; i < (2*n_aug_+1); i++) {
      z_pred_laser_ = z_pred_laser_ + weights_(i) * Zsig_laser_.col(i);
  }

  S_laser_.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_laser_.col(i) - z_pred_laser_;

    //angle normalization
    //while (z_diff(1)> M_PI) z_diff(1) -= 2.0*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1) += 2.0*M_PI;

    S_laser_ = S_laser_ + weights_(i) * z_diff * z_diff.transpose();
  }

  S_laser_ = S_laser_ + R_laser_;

  //print result
  //std::cout << "z_pred: " << std::endl << z_pred_radar_ << std::endl;
  //std::cout << "S: " << std::endl << S_radar_ << std::endl;

}

///*/////////////////////////////////////////////////////////////////////////////////
void UKF::UpdateState_Radar(const MeasurementPackage &measurement_pack) {

  //create example vector for incoming radar measurement
  VectorXd z = measurement_pack.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.0*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.0 * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.0 * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_radar_.inverse();

  //residual
  VectorXd z_diff = z - z_pred_radar_;

  //angle normalization
  while (z_diff(1)>  M_PI) z_diff(1) -= 2.0*M_PI;
  while (z_diff(1)< -M_PI) z_diff(1) += 2.0*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_radar_* K.transpose();

  ///*   NIS Calculations and some statistics  /////////////////
  #define Radar_3DoF_95percent_Threshold 7.8

  // Calculate NIS
  NIS_radar_ = z_diff.transpose() * S_radar_.inverse() * z_diff;

  NIS_radar_sum = NIS_radar_sum + NIS_radar_;

  Radar_records_counter++;

  if (NIS_radar_ > Radar_3DoF_95percent_Threshold)
      Radar_DOF_threshold_counter++;

  //print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}
///*/////////////////////////////////////////////////////////////////////////////////
void UKF::UpdateState_Laser(const MeasurementPackage &measurement_pack) {

  //create example vector for incoming radar measurement
  VectorXd z = measurement_pack.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_laser_);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < (2 * n_aug_ + 1); i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_laser_.col(i) - z_pred_laser_;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.0 * M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.0 * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_laser_.inverse();

  //residual
  VectorXd z_diff = z - z_pred_laser_;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_laser_ * K.transpose();


  ///*   NIS Calculations and some statistics  /////////////////
  #define Laser_2DoF_95percent_Threshold 5.99

  // Calculate NIS
  NIS_laser_ = z_diff.transpose() * S_laser_.inverse() * z_diff;

  NIS_laser_sum = NIS_laser_sum + NIS_radar_;

  Laser_records_counter++;

  if (NIS_laser_ > Laser_2DoF_95percent_Threshold)
      Laser_DOF_threshold_counter++;

  //print result
  //std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;

}
