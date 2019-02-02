#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include <iostream>
#include <string>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

class UKF {
  private:
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);
  void PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t);
  void PredictMeanAndCovariance();
  void NIS(MatrixXd S, VectorXd z_pred, VectorXd z_meas, bool type);


  public:
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
  void Prediction(double delta_t);

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


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // if this is true, initialisation will be done based on two measurements
  bool use_second_init_;

  // if this is true, NIS values will be calculated and saved to files
  bool use_NIS_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;

  //variables for second initialisation loop - for calculating teta and speed: x_(2) and x_(3)
  bool is_second_initialized_;
  int iter_;
  float x_init_1_;
  float y_init_1_;
  float x_init_2_;
  float y_init_2_;

  float NIS_radar_;
  float NIS_laser_;

  const double NIS_RADAR_050 = 7.815;
  const double NIS_LASER_050 = 5.991;

  int NIS_radar_total_;
  int NIS_laser_total_;
  int NIS_radar_overlimit_;
  int NIS_laser_overlimit_;

};

#endif  // UKF_H
