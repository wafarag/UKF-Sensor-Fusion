#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* saving "k-1" sample time
  long long previous_timestamp_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* measurement covariance matrix S for Radar
  MatrixXd S_radar_;

  ///* measurement covariance matrix S for Laser
  MatrixXd S_laser_;

  ///* generated sigma points matrix
  MatrixXd Xsig_;

  ///* generated augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix (as columns)
  MatrixXd Xsig_pred_;

  ///* matrix for sigma points in measurement space - Radar
  MatrixXd Zsig_radar_;

  ///* matrix for sigma points in measurement space - Laser
  MatrixXd Zsig_laser_;

  ///* Initialize mean predicted measurement vector for Radar
  VectorXd z_pred_radar_;

  ///* Initialize mean predicted measurement vector for Laser
  VectorXd z_pred_laser_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* time difference between measurements in sec
  double delta_t; //time diff in sec

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* define measurement dimension, radar can measure r, phi, and r_dot
  int n_z_radar_;

  ///* define measurement dimension, lidar (laser) can measure positions Px and Py
  int n_z_laser_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the summation of radar NIS records
  double NIS_radar_sum;

  ///* the counter of radar readings
  int Radar_records_counter;

  ///* the counter of radar readings above DOF threshold
  int Radar_DOF_threshold_counter;

  ///* the current NIS for laser
  double NIS_laser_;

  ///* the summation of laser NIS records
  double NIS_laser_sum;

  ///* the counter of laser readings
  int Laser_records_counter;

  ///* the counter of laser readings above DOF threshold
  int Laser_DOF_threshold_counter;

  ///* Radar measurement noise covariance matrix
  MatrixXd R_radar_;

  ///* Lidar measurement noise covariance matrix
  MatrixXd R_laser_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void GenerateSigmaPoints();

  void AugmentedSigmaPoints();

  void SigmaPointPrediction();

  void PredictMeanAndCovariance();

  void PredictRadarMeasurement();

  void PredictLaserMeasurement();

  void UpdateState_Radar(const MeasurementPackage &measurement_pack);

  void UpdateState_Laser(const MeasurementPackage &measurement_pack);
};

#endif /* UKF_H */
