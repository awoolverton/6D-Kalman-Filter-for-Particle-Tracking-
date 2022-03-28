#include <map>
#include <vector>
#include <algorithm>
#include <TVector3.h>

//=====================================================
//Kalman Filter for 3D particle tracking
class KalmanFilter{
public:
  ///constructor
  KalmanFilter(const TMatrixD &Fin,
	       const TMatrixD &Hin,
	       const TMatrixD &Qin,
	       const TMatrixD &Rin,
	       double Delta_t);

  KalmanFilter(const TMatrixD &Fin,
	       const TMatrixD &Hin,
	       const TMatrixD &Qin,
	       const TMatrixD &Rin,
	       double Delta_t, int max);
  
  ///destructor
  ~KalmanFilter(){};

  // Initialize the filter with a guess for initial states.
  void init(double t0, const TVectorD &x_in, const TMatrixD &P_in);
  void SetConfidance(double confidance); 

  void RunFilter(map<int, TVectorD> &measurments);
  void Reset();
  
  vector<TVectorD> GetStateVectors();
  vector<TVectorD> GetResiduals();
  vector<int> GetOutliers();
 
  double GetChi2();
  
private:

  /*
  * Create a Kalman filter with the specified matrices.
  *   F - System dynamics matrix
  *   H - Observation model
  *   Q - Process noise covariance
  *   R - Measurement noise covariance
  *   P - Estimate error covariance
  */

  const TMatrixD F;
  const TMatrixD H;
  const TMatrixD Q;
  const TMatrixD R;

  const TMatrixD P0;
  const TVectorD x0;
  Double_t Chi2;
  
  vector<TVectorD> states;
  vector<TMatrixD> covariance;
  vector<TVectorD> residuals;
  
  vector<int> outliers;
  bool initialize;

  double t_last;
  double t, dt;    // estimated time and time step 
  int n,m;         // dimention of the filter

  bool check_outlier = true;
  double alpha = 0.05; // 95% confidance
  
  bool Check_Residual(const TVectorD &y, const TMatrixD &C);
  pair<TVectorD, TMatrixD> Predict(const TVectorD &xin, const TMatrixD &Pin);
  tuple<bool, TVectorD, TMatrixD> Update(const TVectorD &z, const TVectorD &x, const TMatrixD &P);
};
