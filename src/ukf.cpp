#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;

  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  // create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);

  // set weights
  weights.fill(1/(2*(lambda + n_aug)));
  weights(0) = lambda/(lambda + n_aug); 
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
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
    cout << "UKF Init: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = measurement_pack.raw_measurements_(0);
      float teta = measurement_pack.raw_measurements_(1);
      ekf_.x_(0) = ro * cos(teta);
      ekf_.x_(1) = ro * sin(teta);
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
      ekf_.x_(4) = 0;
      }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0;
      ekf_.x_(3) = 0;
      ekf_.x_(4) = 0;
      }
      
    time_us_ = measurement_pack.timestamp_;

    is_initialized_ = true;


    }else{

      double delta_t = (measurement_pack.timestamp_ - time_us_)/1000000.0;
      time_us_ = measurement_pack.timestamp_;

      Prediction(delta_t);

      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
        UpdateRadar(measurement_pack);
      else
        UpdateLidar(measurement_pack);
    } 
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  /**  
  * Augmentation 
  */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(&Xsig_aug);
  /**  
  * Sigma point prediction 
  */
  //MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  PredictSigmaPoints(&Xsig_pred);
  /**  
  * the state x and the state covariance matrix P prediction 
  */
  PredictMeanAndCovariance();

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2;
  VectorXd z_out = VectorXd(2);
  MatrixXd S_out = MatrixXd(2, 2);
  PredictRadarMeasurement(z_out, S_out);
  UpdateStateLaser();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = 3;
  VectorXd z_out = VectorXd(3);
  MatrixXd S_out = MatrixXd(3, 3);
  PredictRadarMeasurement(z_out, S_out);
  UpdateStateRadar();
}



void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){
    // create augmented mean vector
  VectorXd x_aug = VectorXd(7);

    x_aug.head(5) = x;
    x_aug(5) = 0;
    x_aug(6) = 0;

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x, n_x) = P;
    MatrixXd Q = MatrixXd(2, 2);
    Q << std_a*std_a, 0,
        0, std_yawdd*std_yawdd;
    P_aug.bottomRightCorner(2, 2) = Q;

  // create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
    for (int i = 0; i<n_aug; i++){
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda+n_aug)*A.col(i);
        Xsig_aug.col(n_aug+i+1) = x_aug - sqrt(lambda+n_aug)*A.col(i);
    }

  *Xsig_out = Xsig_aug;

}


void UKF::PredictSigmaPoints(MatrixXd* Xsig_out){
// predict sigma points

  for (int i = 0; i < Xsig_aug.cols(); i++){
      double v_k = Xsig_aug(2, i);
      double yaw_k = Xsig_aug(3, i);
      double yawd_k = Xsig_aug(4, i);
      double std_a = Xsig_aug(5, i);
      double std_yawd = Xsig_aug(6, i);
      
      if (yawd_k < 0.001){
          Xsig_pred(0, i) = Xsig_aug(0, i) + v_k * cos(yaw_k) * delta_t;
          Xsig_pred(1, i) = Xsig_aug(1, i) + v_k * sin(yaw_k) * delta_t;
      }else{
          Xsig_pred(0, i) = Xsig_aug(0, i) + (v_k/yawd_k)*(sin(yaw_k + yawd_k * delta_t)-sin(yaw_k));
          Xsig_pred(1, i) = Xsig_aug(1, i) + (v_k/yawd_k)*(-cos(yaw_k + yawd_k * delta_t)+cos(yaw_k));
      }
      Xsig_pred(2, i) = Xsig_aug(2, i);
      Xsig_pred(3, i) = Xsig_aug(3, i) + yawd_k * delta_t;
      Xsig_pred(4, i) = Xsig_aug(4, i);
      
      //add uncertainty
      Xsig_pred(0, i) += (1./2.)*delta_t*delta_t*cos(yaw_k)*std_a;
      Xsig_pred(1, i) += (1./2.)*delta_t*delta_t*sin(yaw_k)*std_a;
      Xsig_pred(2, i) += delta_t * std_a;
      Xsig_pred(3, i) += (1./2.)*delta_t*delta_t*std_yawd;
      Xsig_pred(4, i) += delta_t * std_yawd;
  }
}

void UKF::PredictMeanAndCovariance(){
  // predict state mean
    x = Xsig_pred * weights;
    
    P.fill(0);
    // predict state covariance matrix
    for (int i = 0; i<Xsig_pred.cols(); i++){
        Xsig_pred.col(i) -= x;
        while (Xsig_pred.col(i)(3)> M_PI) Xsig_pred.col(i)(3)-=2.*M_PI;
        while (Xsig_pred.col(i)(3)<-M_PI) Xsig_pred.col(i)(3)+=2.*M_PI;
        P += weights(i) * Xsig_pred.col(i) * Xsig_pred.col(i).transpose();
    }
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out){
      // create matrix for sigma points in measurement space
    int n_z = 3;

    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    
    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);



    // transform sigma points into measurement space
    for (int i = 0; i < Xsig_pred.cols(); i++){
        double p_x = Xsig_pred(0, i);
        double p_y = Xsig_pred(1, i);
        double v = Xsig_pred(2, i);
        double yaw = Xsig_pred(3, i);
        Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
        Zsig(1, i) = atan2(p_y, p_x);
        Zsig(2, i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v)/Zsig(0, i);
    }
    // calculate mean predicted measurement
    z_pred = Zsig * weights;
    
    // calculate innovation covariance matrix S
    MatrixXd R = MatrixXd(n_z,n_z);
    R << std_radr*std_radr, 0, 0,
    0, std_radphi*std_radphi, 0,
    0, 0, std_radrd*std_radrd;
    S = R;
    for (int i = 0; i < Zsig.cols(); i++){
        Zsig.col(i) -= z_pred;
        while (Zsig.col(i)(1)> M_PI) Zsig.col(i)(1)-=2.*M_PI;
        while (Zsig.col(i)(1)<-M_PI) Zsig.col(i)(1)+=2.*M_PI;
        S += weights(i) * Zsig.col(i) * Zsig.col(i).transpose();
    }
}

void UKF::PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out){
    int n_z = 2;
      // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
    
    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);

    // transform sigma points into measurement space
    for (int i = 0; i < Xsig_pred.cols(); i++){
        Zsig(0, i) = Xsig_pred(0, i);
        Zsig(1, i) = Xsig_pred(1, i);
    }
    // calculate mean predicted measurement
    z_pred = Zsig * weights;
    
    // calculate innovation covariance matrix S
    MatrixXd R = MatrixXd(n_z, n_z);
    R << std_laspx_*std_laspx_, 0,
    0, std_laspy_*std_laspy_;

    S = R;
    for (int i = 0; i < Zsig.cols(); i++){
        Zsig.col(i) -= z_pred;
        S += weights(i) * Zsig.col(i) * Zsig.col(i).transpose();
    }
}

void UKF::UpdateStateRadar(){
  // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    /**
     * Student part begin
     */
    
    // calculate cross correlation matrix
    Tc.fill(0);
    for (int i = 0; i < Zsig.cols(); i++){
        Xsig_pred.col(i) -= x;
        while (Xsig_pred.col(i)(3)> M_PI) Xsig_pred.col(i)(3)-=2.*M_PI;
        while (Xsig_pred.col(i)(3)<-M_PI) Xsig_pred.col(i)(3)+=2.*M_PI;
        Zsig.col(i) -= z_pred;
        while (Zsig.col(i)(1)> M_PI) Zsig.col(i)(1)-=2.*M_PI;
        while (Zsig.col(i)(1)<-M_PI) Zsig.col(i)(1)+=2.*M_PI;
        
        Tc += weights(i) * Xsig_pred.col(i) * Zsig.col(i).transpose();
    }
    
    // calculate Kalman gain K;
    MatrixXd K = MatrixXd(n_x, n_z);
    K = Tc * S.inverse();
    // update state mean and covariance matrix
    VectorXd y = z - z_pred;
    while (y(1)> M_PI) y(1) -=2.*M_PI;
    while (y(1)<-M_PI) y(1) +=2.*M_PI;
    x += K * y;
    
    P -= K * S * K.transpose();
    

}

void UKF::UpdateStateLaser(){
  // create matrix for cross correlation Tc

    MatrixXd Tc = MatrixXd(n_x, n_z);
    
    /**
     * Student part begin
     */
    
    // calculate cross correlation matrix
    Tc.fill(0);
    for (int i = 0; i < Zsig.cols(); i++){
        Xsig_pred.col(i) -= x;
        while (Xsig_pred.col(i)(3)> M_PI) Xsig_pred.col(i)(3)-=2.*M_PI;
        while (Xsig_pred.col(i)(3)<-M_PI) Xsig_pred.col(i)(3)+=2.*M_PI;
        Zsig.col(i) -= z_pred;
        
        Tc += weights(i) * Xsig_pred.col(i) * Zsig.col(i).transpose();
    }
    
    // calculate Kalman gain K;
    MatrixXd K = MatrixXd(n_x, n_z);
    K = Tc * S.inverse();
    // update state mean and covariance matrix
    VectorXd y = z - z_pred;
    x += K * y;
    
    P -= K * S * K.transpose();
    
}
