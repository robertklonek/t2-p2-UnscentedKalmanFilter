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

  // if this is true, initialisation will be done based on two measurements
  use_second_init_ = false;

  // if this is true, NIS values will be calculated and saved to files
  use_NIS_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.34;
  
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
  is_second_initialized_ = false;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  P_ << 0.1, 0,   0, 0, 0,
        0,   0.1, 0, 0, 0,
        0,   0,   1, 0, 0,
        0,   0,   0, 1, 0,
        0,   0,   0, 0, 1;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // create vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  // set weights
  weights_.fill(1/(2*(lambda_ + n_aug_)));
  weights_(0) = lambda_/(lambda_ + n_aug_); 

  iter_ = 0;

  NIS_radar_total_ = 0;
  NIS_laser_total_ = 0;
  NIS_radar_overlimit_ = 0;
  NIS_laser_overlimit_ = 0;
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


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_(0);
      float teta = meas_package.raw_measurements_(1);
      x_(0) = ro * cos(teta);
      x_(1) = ro * sin(teta);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
      }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
      }

    //second init
    if (use_second_init_){
    x_init_1_ = x_(0);
    y_init_1_ = x_(1);
    }

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    }else if (!is_second_initialized_ and use_second_init_){
    /**
      Second initialization based on two measurements data for calculation init speed and teta
    */
      cout << "UKF Second : " << endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_(0);
      float teta = meas_package.raw_measurements_(1);
      x_(0) = ro * cos(teta);
      x_(1) = ro * sin(teta);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
      }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      x_(2) = 0;
      x_(3) = 0;
      x_(4) = 0;
      }

      x_init_2_ = x_(0);
      y_init_2_ = x_(1);

       cout <<"x_init_1: "<<x_init_1_<< endl;
       cout <<"x_init_2: "<<x_init_2_<< endl;
       cout <<"y_init_1: "<<y_init_1_<< endl;
       cout <<"y_init_2: "<<y_init_2_<< endl;

      double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
      time_us_ = meas_package.timestamp_;

      float x_diff = x_init_2_ - x_init_1_;
      float y_diff = y_init_2_ - y_init_1_;
      //teta
      x_(3) = atan2(y_diff, x_diff);
      while (x_(3)> M_PI) x_(3)-=2.*M_PI;
      while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
      //v_a
      x_(2) = sqrt(x_diff*x_diff + y_diff*y_diff)/delta_t;

      cout <<"v speed: "<<x_(2)<< endl;
      cout <<"teta: "<<x_(3)<< endl;

    is_second_initialized_ = true;

    }else{
      iter_++;
      cout << "UKF loop: " <<iter_<< endl;



      double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
      time_us_ = meas_package.timestamp_;

      Prediction(delta_t);
      cout << "UKF prediction done: " << endl;

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
        UpdateRadar(meas_package);
        cout << "Radar Update done: " << endl;}
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
        UpdateLidar(meas_package);
        cout << "Lidar Update done: " << endl;}
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
  PredictSigmaPoints(Xsig_aug, delta_t);
  /**  
  * the state x and the state covariance matrix P prediction 
  */
  PredictMeanAndCovariance();
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out){


  // create augmented mean vector
  VectorXd x_aug = VectorXd(7);

    x_aug.head(5) = x_;
    x_aug(5) = 0;
    x_aug(6) = 0;

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    MatrixXd Q = MatrixXd(2, 2);
    Q <<  std_a_*std_a_,  0,
          0,              std_yawdd_*std_yawdd_;

    P_aug.bottomRightCorner(2, 2) = Q;

  // create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
    for (int i = 0; i<n_aug_; i++){
        Xsig_aug.col(i+1) = x_aug + sqrt(lambda_+n_aug_)*A.col(i);
        Xsig_aug.col(n_aug_+i+1) = x_aug - sqrt(lambda_+n_aug_)*A.col(i);
    }

  *Xsig_out = Xsig_aug;

}


void UKF::PredictSigmaPoints(MatrixXd Xsig_aug, double delta_t){
  // predict sigma points

  for (int i = 0; i < Xsig_aug.cols(); i++){
      double v_k = Xsig_aug(2, i);
      double yaw_k = Xsig_aug(3, i);
      double yawd_k = Xsig_aug(4, i);
      double std_a = Xsig_aug(5, i);
      double std_yawd = Xsig_aug(6, i);
      
      if (yawd_k < 0.001){
          Xsig_pred_(0, i) = Xsig_aug(0, i) + v_k * cos(yaw_k) * delta_t;
          Xsig_pred_(1, i) = Xsig_aug(1, i) + v_k * sin(yaw_k) * delta_t;
      }else{
          Xsig_pred_(0, i) = Xsig_aug(0, i) + (v_k/yawd_k)*(sin(yaw_k + yawd_k * delta_t)-sin(yaw_k));
          Xsig_pred_(1, i) = Xsig_aug(1, i) + (v_k/yawd_k)*(-cos(yaw_k + yawd_k * delta_t)+cos(yaw_k));
      }
      Xsig_pred_(2, i) = Xsig_aug(2, i);
      Xsig_pred_(3, i) = Xsig_aug(3, i) + yawd_k * delta_t;
      Xsig_pred_(4, i) = Xsig_aug(4, i);
      
      //add uncertainty
      Xsig_pred_(0, i) += (1./2.)*delta_t*delta_t*cos(yaw_k)*std_a;
      Xsig_pred_(1, i) += (1./2.)*delta_t*delta_t*sin(yaw_k)*std_a;
      Xsig_pred_(2, i) += delta_t * std_a;
      Xsig_pred_(3, i) += (1./2.)*delta_t*delta_t*std_yawd;
      Xsig_pred_(4, i) += delta_t * std_yawd;
  }

  
}

void UKF::PredictMeanAndCovariance(){
  // predict state mean
  //cout <<"Xsig_pred_: "<< Xsig_pred_<< endl;

    x_ = Xsig_pred_ * weights_;
    P_.fill(0);
    // predict state covariance matrix
    for (int i = 0; i<Xsig_pred_.cols(); i++){
        VectorXd diff_col = Xsig_pred_.col(i) - x_;
        while (diff_col(3)> M_PI) diff_col(3)-=2.*M_PI;
        while (diff_col(3)<-M_PI) diff_col(3)+=2.*M_PI;
        P_ += weights_[i] * diff_col * diff_col.transpose();
    }
}


void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //measurement z
  VectorXd z = VectorXd(n_z);
  z <<  meas_package.raw_measurements_(0),
        meas_package.raw_measurements_(1);

  /**
  * Predict Laser measurements
  */

  // transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); i++){
      Zsig(0, i) = Xsig_pred_(0, i);
      Zsig(1, i) = Xsig_pred_(1, i);
  }
  // calculate mean predicted measurement

  z_pred = Zsig * weights_;
  
  // calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<  std_laspx_*std_laspx_,  0,
        0,                      std_laspy_*std_laspy_;

  S = R;
  for (int i = 0; i < Zsig.cols(); i++){
      VectorXd diff_col = Zsig.col(i) - z_pred;
      S += weights_(i) * diff_col * diff_col.transpose();
  }

  /**
   *Update Laser
   */
  // create matrix for cross correlation Tc

  NIS(S, z_pred, z, false);


  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  // calculate cross correlation matrix
  Tc.fill(0);
  for (int i = 0; i < Zsig.cols(); i++){
      VectorXd diff_col = Xsig_pred_.col(i) - x_;
      while (diff_col(3)> M_PI) diff_col(3)-=2.*M_PI;
      while (diff_col(3)<-M_PI) diff_col(3)+=2.*M_PI;
      Zsig.col(i) -= z_pred;
      
      Tc += weights_(i) * diff_col * Zsig.col(i).transpose();
  }
  
  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  // update state mean and covariance matrix
  VectorXd y = z - z_pred;
  x_ += K * y;
  
  P_ -= K * S * K.transpose();

  cout <<"x_ updated: "<< x_<< endl;
  cout <<"P_ updated: "<< P_<< endl;
}

void UKF::NIS(MatrixXd S, VectorXd z_pred, VectorXd z_meas, bool type){
  
  VectorXd z_diff = z_meas - z_pred;

  cout<<z_diff.transpose() * S.inverse() * z_diff<<endl;


  VectorXd NIS = z_diff.transpose() * S.inverse() * z_diff;
  cout<<"NIS: "<<NIS(0);

  FILE * pFile;
  
  if (type){//radar  
    pFile = fopen ("Radar_NIS.txt","a+");
    if (pFile!=NULL)
  {
    //fputs (NIS[0], pFile);
    fprintf (pFile, "%f\n",NIS(0));
    fclose (pFile);
  }

  NIS_radar_total_++;
  if (NIS(0) > NIS_RADAR_050) NIS_radar_overlimit_++;
  cout<<"Radar NIS: "<<float(NIS_radar_overlimit_)/float(NIS_radar_total_);


  }else{

      pFile = fopen ("Laser_NIS.txt","a+");
        if (pFile!=NULL)
      {
        //fputs (NIS[0], pFile);
        fprintf (pFile, "%f\n",NIS(0));
        fclose (pFile);
      }

    NIS_laser_total_++;
    if (NIS(0) > NIS_LASER_050) NIS_laser_overlimit_++;

    cout<<"Laser NIS: "<<float(NIS_laser_overlimit_)/float(NIS_laser_total_);
  }

  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
      // create matrix for sigma points in measurement space
  int n_z = 3;

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

    //measurement z
  VectorXd z = VectorXd(n_z);
  z <<  meas_package.raw_measurements_(0), 
        meas_package.raw_measurements_(1),
        meas_package.raw_measurements_(2);

  cout<<"z measurement: "<<z<<endl;

  /**
   *Predict Radar measurements
   */

  // transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); i++){
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3, i);
      Zsig(0, i) = sqrt(p_x*p_x + p_y*p_y);
      Zsig(1, i) = atan2(p_y, p_x);
      Zsig(2, i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v)/Zsig(0, i);
  }
  // calculate mean predicted measurement
  z_pred = Zsig * weights_;
  cout<<"z_pred: "<<z_pred<<endl;
  
  // calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R <<  std_radr_*std_radr_,  0,                        0,
        0,                    std_radphi_*std_radphi_,  0,
        0,                    0,                        std_radrd_*std_radrd_;
  S = R;
  for (int i = 0; i < Zsig.cols(); i++){
      
    VectorXd diff_col = Zsig.col(i) - z_pred;
    while (diff_col(1)> M_PI) diff_col(1)-=2.*M_PI;
    while (diff_col(1)<-M_PI) diff_col(1)+=2.*M_PI;
    S += weights_(i) * diff_col * diff_col.transpose();
  }

  /**
   *Update Radar
   */
  NIS(S, z_pred, z, true);

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  
  /**
   * Student part begin
   */
  
  // calculate cross correlation matrix
  Tc.fill(0);
  for (int i = 0; i < Zsig.cols(); i++){
      Xsig_pred_.col(i) -= x_;
      while (Xsig_pred_.col(i)(3)> M_PI) Xsig_pred_.col(i)(3)-=2.*M_PI;
      while (Xsig_pred_.col(i)(3)<-M_PI) Xsig_pred_.col(i)(3)+=2.*M_PI;
      Zsig.col(i) -= z_pred;
      while (Zsig.col(i)(1)> M_PI) Zsig.col(i)(1)-=2.*M_PI;
      while (Zsig.col(i)(1)<-M_PI) Zsig.col(i)(1)+=2.*M_PI;
      
      Tc += weights_(i) * Xsig_pred_.col(i) * Zsig.col(i).transpose();
  }


  
  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();
  // update state mean and covariance matrix
  VectorXd y = z - z_pred;
  while (y(1)> M_PI) y(1) -=2.*M_PI;
  while (y(1)<-M_PI) y(1) +=2.*M_PI;
  x_ += K * y;
  
  P_ -= K * S * K.transpose();

  cout <<"x_ updated: "<< x_<< endl;
  cout <<"P_ updated: "<< P_<< endl;

}


